#include "native_reference.hpp"
#include "legacy_provider_reference.hpp"
#include <numeric>

namespace native_reference {
namespace io=snapshot_io;
namespace gp=gather_push;
namespace {
struct CaptureClock {double origin,relative;};
CaptureClock capture_clock(const std::filesystem::path& fixture) {
    // Exact original thermal origin, read from its frozen binary header rather
    // than inferred by subtracting rounded clocks. See original Thermal schema.
    std::ifstream thermal("/home/treefor/sim-workspace/powerflow/outputs/heating-history-full.bin",std::ios::binary);
    std::uint64_t magic=0,n=0;thermal.read(reinterpret_cast<char*>(&magic),8);thermal.read(reinterpret_cast<char*>(&n),8);
    io::require(magic==0x544845524d303031ULL && n>=2 && n<=10000,"Invalid frozen thermal origin header");
    thermal.seekg(48);double origin;thermal.read(reinterpret_cast<char*>(&origin),8);
    io::require(bool(thermal) && std::isfinite(origin),"Missing exact thermal origin");
    const auto diagnostic=io::read<GapPicDiagnostics>(fixture/"diagnostics_before.bin");
    io::require(diagnostic.size()==1,"Expected actual captured diagnostic clock");
    return {origin,diagnostic[0].time_s};
}
}
RectField::RectField(const pic_batch::HostSnapshot& s,const std::vector<double>& data)
    :nx(s.parameters.self_nx),ny(s.parameters.self_ny),width(s.parameters.width),height(s.parameters.height),
     mesh(mfem::Mesh::MakeCartesian2D(s.parameters.self_nx,s.parameters.self_ny,mfem::Element::QUADRILATERAL,true,s.parameters.width,s.parameters.height,false)),
     collection(1,2),space(&mesh,&collection),field(&space),gather(field) {
    io::require(data.size()==std::size_t(field.Size()),"Actual field nodal shape changed");
    std::copy(data.begin(),data.end(),field.HostWrite());
    std::vector<int> elements(mesh.GetNE());std::iota(elements.begin(),elements.end(),0);
    const auto cells=io::export_grid(field,elements);
    io::require(cells.size()==s.self_cells.size() && std::memcmp(cells.data(),s.self_cells.data(),cells.size()*sizeof(gp::Cell))==0,"Recreated MFEM field layout differs from actual capture");
}
PicVector RectField::electric(double x,double y) {
    // Dimensions come from the exact regular field space. Retain original
    // GapPic locate expression/order and point domain checks.
    io::require(std::isfinite(x) && std::isfinite(y) && x>=0 && x<=width && y>=0 && y<=height,"Point outside native reference domain");
    const double sx=x/width*nx,sy=y/height*ny;
    const int ix=std::min(static_cast<int>(sx),nx-1),iy=std::min(static_cast<int>(sy),ny-1);
    mfem::IntegrationPoint ip;ip.Set2(sx-ix,sy-iy);
    const auto gradient=gather.gradient(iy*nx+ix,ip);return {-gradient[0],-gradient[1],0};
}
double RectField::value(double x,double y) {
    const double sx=x/width*nx,sy=y/height*ny;
    const int ix=std::clamp(static_cast<int>(sx),0,nx-1),iy=std::clamp(static_cast<int>(sy),0,ny-1);
    mfem::IntegrationPoint ip;ip.Set2(sx-ix,sy-iy);return field.GetValue(iy*nx+ix,ip);
}
NativeReference::NativeReference(const pic_batch::HostSnapshot& s,const std::filesystem::path& fixture,const std::filesystem::path& voltage_csv)
    :p_(s.parameters),self_(s,s.phi),unit_(s,s.psi),
     native_(std::make_unique<Fields>("/home/treefor/projects/Power-Flow/runs/restore_field_mpi8_domain40mm")),voltage_(std::make_unique<Voltage>(voltage_csv)) {
    io::require(p_.self_nx==360 && p_.self_ny==360,"Frozen native timing requires actual360x360 fixture");
    const auto clock=capture_clock(fixture);origin_=clock.origin;relative_midpoint_=clock.relative+.5*p_.dt;absolute_midpoint_=origin_+relative_midpoint_;
    native_->set(absolute_midpoint_);const auto V=voltage_->at(absolute_midpoint_);
    io::require(native_->fraction==p_.fraction && V==p_.voltage,"Actual selected fraction or voltage differs");
    io::require(io::values(*native_->lower->az)==s.az0 && io::values(*native_->upper->az)==s.az1 && io::values(*native_->lower->rate)==s.rate0 && io::values(*native_->upper->rate)==s.rate1 && native_->lower->xs==s.xs && native_->lower->ys==s.ys,"Frozen native fields differ from actual capture");
    electric_=[this](double x,double y,double t){return native_->value(x,y,origin_+t,false);};
    magnetic_=[this](double x,double y,double t){return native_->value(x,y,origin_+t,true);};
    axial_=[this](double x,double y,double z,double t) {
        double V=voltage_->at(origin_+t),g=(z-p_.short_z)/p_.depth;
        auto E=unit_.electric(x,y);for(int k=0;k<2;++k)E[k]*=V*(2*g-1);
        E[2]=-2*V/p_.depth*(unit_.value(x,y)-.5);return E;
    };
}
NativeReference::~NativeReference()=default;
std::pair<PicVector,PicVector> NativeReference::fields(const GapParticle& p) {
    auto E=self_.electric(p.position_m[0],p.position_m[1]);PicVector B{};
    if(electric_) {auto extra=electric_(p.position_m[0],p.position_m[1],relative_midpoint_);for(int i=0;i<3;++i)E[i]+=extra[i];}
    if(axial_) {auto extra=axial_(p.position_m[0],p.position_m[1],p.position_m[2],relative_midpoint_);for(int i=0;i<3;++i)E[i]+=extra[i];}
    if(magnetic_)B=magnetic_(p.position_m[0],p.position_m[1],relative_midpoint_);
    return {E,B};
}
void NativeReference::force_update(std::span<GapParticle> particles) {
    for(auto& p:particles) {const auto [E,B]=fields(p);p.proper_velocity_m_s=pic_boris(p.proper_velocity_m_s,p.charge_C/p.mass_kg,E,B,p_.dt);}
}
void NativeReference::generate(std::span<const GapParticle> particles,std::span<gp::Outputs> output) {
    io::require(particles.size()==output.size(),"Native reference output size mismatch");
    for(std::size_t i=0;i<particles.size();++i) {
        const auto& p=particles[i];const auto [E,B]=fields(p);
        output[i]={io::vec(E),io::vec(B),io::vec(pic_boris(p.proper_velocity_m_s,p.charge_C/p.mass_kg,E,B,p_.dt)),0};
    }
}
}
