#include "gather_push_core.hpp"
#include "snapshot_io.hpp"
#include "batch_compare.hpp"
#include "native_heating_fields.hpp"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <set>
#include <sstream>

namespace gp=gather_push;
using Clock=std::chrono::steady_clock;
static_assert(sizeof(gp::Cell)==80 && sizeof(gp::Inputs)==56 && sizeof(gp::Outputs)==80 && sizeof(gp::Tables)==152);
namespace {
constexpr double E_abs=1e-6,B_abs=1e-13,u_abs=1e-5,field_rel=5e-12,u_rel=1e-11,weighted_limit=1e-12;
using snapshot_io::require;using snapshot_io::read;using snapshot_io::write;
using snapshot_io::vec;using snapshot_io::values;using snapshot_io::export_grid;
struct Saved {
    GapPicOptions options;GapPicDiagnostics ledger;std::vector<GapParticle> particles;
};
Saved saved(const std::filesystem::path& path) {
    std::ifstream f(path,std::ios::binary);Saved s;std::uint64_t magic=0,id=0,n=0;
    const auto get=[&](auto& x){f.read(reinterpret_cast<char*>(&x),sizeof(x));};
    get(magic);get(s.options);get(s.ledger);get(id);get(n);
    require(bool(f) && magic==0x5241494b50494332ULL && n>0 && n<=2000000 && id>0,"Invalid frozen-ABI saved state");
    std::set<std::uint64_t> ids;
    for(std::uint64_t i=0;i<n;++i) {
        GapParticle p;get(p);
        require(bool(f) && p.id<id && ids.insert(p.id).second && std::isfinite(p.weight) && p.weight>0 &&
                std::isfinite(p.mass_kg) && p.mass_kg>0 && std::isfinite(p.charge_C) &&
                std::isfinite(p.charge_C*p.weight) && std::isfinite(p.charge_C/p.mass_kg) &&
                gp::finite(vec(p.position_m)) && gp::finite(vec(p.proper_velocity_m_s)) &&
                gp::finite(gp::norm2(vec(p.proper_velocity_m_s))),"Invalid saved particle");
        require(p.position_m[0]>=0 && p.position_m[0]<=s.options.width_m && p.position_m[1]>=0 &&
                p.position_m[1]<=s.options.gap_m,"Saved particle outside domain");
        s.particles.push_back(p);
    }
    require(bool(f) && f.peek()==std::char_traits<char>::eof(),"Truncated/trailing saved state");return s;
}
std::pair<double,double> gradient(const mfem::GridFunction& f,int e,const mfem::IntegrationPoint& ip) {
    auto* tr=f.FESpace()->GetMesh()->GetElementTransformation(e);tr->SetIntPoint(&ip);
    mfem::Vector result(2);f.GetGradient(*tr,result);return {result[0],result[1]};
}
void boris_tests() {
    for(int i=0;i<4096;++i) {
        const PicVector u={i%2?0.:double(i)*1e6,std::sin(i)*1e9,std::cos(i)*3e8};
        const PicVector E={std::sin(i*.3)*1e8,0,std::cos(i*.7)*1e7};
        const PicVector B={i%3?0.:2.,std::cos(i)*.1,std::sin(i)*.1};
        const double qm=i%2 ? -pic_e/pic_me : pic_e/(1836*pic_me),dt=1e-15*(1+i%10);
        const auto ref=pic_boris(u,qm,E,B,dt);gp::Vec actual{};
        require(gp::boris(vec(u),qm,vec(E),vec(B),dt,actual)==0,"Shared Boris failed finite input");
        require(std::memcmp(ref.data(),actual.v,3*sizeof(double))==0,"Shared host Boris differs from original result bits");
    }
    gp::Vec out{};const double nan=std::numeric_limits<double>::quiet_NaN(),inf=std::numeric_limits<double>::infinity();
    require(gp::boris({},1,{}, {},0,out)==gp::invalid_input,"Nonpositive dt accepted");
    require(gp::boris({{nan,0,0}},1,{}, {},1,out)==gp::invalid_input,"NaN velocity accepted");
    require(gp::boris({},inf,{}, {},1,out)==gp::invalid_input,"Infinite q/m accepted");
    require(gp::boris({},1,{{1e308,0,0}}, {},1e308,out)==gp::boris_overflow,"Overflow accepted");
}
void geometry_tests() {
    auto mesh=mfem::Mesh::MakeCartesian2D(3,4,mfem::Element::QUADRILATERAL,true,.009,.009,false);
    for(int i=0;i<mesh.GetNV();++i) {auto* p=mesh.GetVertex(i);p[0]+=.2*p[1];p[1]+=p[0]*p[0];}
    mfem::H1_FECollection collection(1,2);mfem::FiniteElementSpace space(&mesh,&collection);mfem::GridFunction f(&space);
    for(int i=0;i<f.Size();++i) f[i]=std::sin(i)*100;
    std::vector<int> elements(mesh.GetNE());std::iota(elements.begin(),elements.end(),0);
    auto cells=export_grid(f,elements);auto data=values(f);gp::Grid grid{cells.data(),3,4};
    for(int e:elements)for(double x:{0.,.001,.25,.5,.999,1.})for(double y:{0.,.001,.37,.5,.999,1.}) {
        mfem::IntegrationPoint ip;ip.Set2(x,y);auto [dx,dy]=gradient(f,e,ip);
        const auto actual=gp::sample(grid,data.data(),e,x,y);
        require(actual.status==0 && std::abs(actual.value-f.GetValue(e,ip))<1e-10 &&
                std::abs(actual.dx-dx)<1e-7 && std::abs(actual.dy-dy)<1e-7,"Bilinear geometry comparison failed");
    }
    double axes[]={-1,0,2};require(gp::bin(axes,2,0)==1 && gp::bin(axes,2,2)==1 && gp::bin(axes,2,-1)==0 &&
                                 gp::bin(axes,2,-2)==-1 && gp::bin(axes,2,std::numeric_limits<double>::quiet_NaN())==-1,
                                 "Native bin edge semantics changed");
    require(gp::sample(grid,data.data(),-1,.5,.5).status==gp::outside_domain,"Invalid cell accepted");
    mfem::H1_FECollection quadratic(2,2);mfem::FiniteElementSpace high_space(&mesh,&quadratic);mfem::GridFunction high(&high_space);high=0;
    bool rejected=false;try{export_grid(high,elements);}catch(const std::exception&){rejected=true;}
    require(rejected,"Higher order field accepted");
}
using batch_compare::Comparison;using batch_compare::compare;using batch_compare::print_comparison;
void validate_result_gate() {
    std::vector<gp::Outputs> ref(1),actual(1);std::vector<GapParticle> p(1);
    require(compare(ref,actual,p).failures==0,"Identity comparison failed");
    actual[0].E.v[0]=std::numeric_limits<double>::quiet_NaN();require(compare(ref,actual,p).failures>0,"NaN bypassed comparison");
    actual=ref;actual[0].u.v[0]=1;require(compare(ref,actual,p).failures>0,"Wrong component bypassed comparison");
    actual=ref;actual[0].status=gp::invalid_input;require(compare(ref,actual,p).failures>0,"Device error status bypassed comparison");
    actual=ref;ref[0].E.v[1]=std::numeric_limits<double>::quiet_NaN();require(compare(ref,actual,p).failures>0,"Invalid reference bypassed comparison");
}
} // namespace
int main(int argc,char** argv) {
 try {
    require(argc>=3,"Usage: gather_push_fixture make DIR [CHECKPOINT] | check DIR DEVICE_OUTPUT");
    const std::filesystem::path out=argv[2];
    if(std::string(argv[1])=="check") {
        require(argc==4,"Expected GPU output filename");
        const auto c=compare(read<gp::Outputs>(out/"reference.bin"),read<gp::Outputs>(argv[3]),read<GapParticle>(out/"particles.bin"));
        print_comparison(c);return c.failures?1:0;
    }
    require(std::string(argv[1])=="make" && argc<=4,"Unknown fixture command");
    boris_tests();geometry_tests();validate_result_gate();
    const auto state=saved(argc==4?argv[3]:"<workspace>/outputs/pic-profile/seed-current/checkpoint_403642.bin");
    const auto& o=state.options;
    require(o.nx>0 && o.nx<=4096 && o.ny>0 && o.ny<=4096 && std::isfinite(o.width_m) && o.width_m>0 &&
            std::isfinite(o.gap_m) && o.gap_m>0 && std::isfinite(o.depth_m) && o.depth_m>0 && std::isfinite(o.short_z_m),"Invalid saved geometry");
    auto mesh=mfem::Mesh::MakeCartesian2D(o.nx,o.ny,mfem::Element::QUADRILATERAL,true,o.width_m,o.gap_m,false);
    mfem::H1_FECollection collection(1,2);mfem::FiniteElementSpace space(&mesh,&collection);mfem::GridFunction phi(&space),psi(&space);
    // Checkpoints do NOT save potential. These two fields are explicitly
    // synthetic; actual native Az/rate and all particle metadata are saved data.
    mfem::FunctionCoefficient voltage([&](const mfem::Vector& x){return -8e4*x[1]/o.gap_m+1000*std::sin(7*x[0]/o.width_m)*std::cos(3*x[1]/o.gap_m);});
    mfem::FunctionCoefficient unit([&](const mfem::Vector& x){return x[1]/o.gap_m;});
    phi.ProjectCoefficient(voltage);psi.ProjectCoefficient(unit);
    const std::filesystem::path native_dir="/home/treefor/projects/Power-Flow/runs/restore_field_mpi8_domain40mm";
    std::ifstream index(native_dir/"em_field_snapshots.csv");std::string line;std::getline(index,line);std::vector<double> times;
    while(std::getline(index,line))times.push_back(std::stod(split(line).at(0)));
    constexpr double dt=5e-15;const double t=7.49788055492e-10+state.ledger.time_s+.5*dt;
    require(times.size()>=2 && std::isfinite(t) && t>=times.front() && t<=times.back(),"Native time outside coverage");
    const double lookup=t+32*std::numeric_limits<double>::epsilon()*std::max(std::abs(t),1e-12);
    const int k=std::clamp(int(std::upper_bound(times.begin(),times.end(),lookup)-times.begin())-1,0,int(times.size())-2);
    const double fraction=std::clamp((t-times[k])/(times[k+1]-times[k]),0.,1.);
    NativeField lower(native_dir,times[k]),upper(native_dir,times[k+1]);
    std::vector<int> elements(mesh.GetNE());std::iota(elements.begin(),elements.end(),0);
    auto self_cells=export_grid(phi,elements),native_cells=export_grid(*lower.az,lower.elements);
    const auto upper_cells=export_grid(*upper.az,upper.elements);
    require(native_cells.size()==upper_cells.size() && std::memcmp(native_cells.data(),upper_cells.data(),native_cells.size()*sizeof(gp::Cell))==0,"Native snapshot layout changed");
    for(const auto* field:{lower.rate.get(),upper.rate.get()}) {
        const auto cells=export_grid(*field,lower.elements);
        require(cells.size()==native_cells.size() && std::memcmp(cells.data(),native_cells.data(),cells.size()*sizeof(gp::Cell))==0,"Rate snapshot layout changed");
    }
    auto ph=values(phi),ps=values(psi),az0=values(*lower.az),az1=values(*upper.az),rate0=values(*lower.rate),rate1=values(*upper.rate);
    gp::Tables tables{{self_cells.data(),o.nx,o.ny},{native_cells.data(),int(lower.xs.size())-1,int(lower.ys.size())-1},
                      ph.data(),ps.data(),az0.data(),az1.data(),rate0.data(),rate1.data(),lower.xs.data(),lower.ys.data(),
                      o.width_m,o.gap_m,o.depth_m,o.short_z_m,fraction,o.upper_potential_V,dt};
    Q1ScalarGather phi_gather(phi),unit_gather(psi);
    std::vector<gp::Inputs> inputs;std::vector<gp::Outputs> reference,actual;
    for(const auto& p:state.particles)inputs.push_back({p.position_m[0],p.position_m[1],p.position_m[2],vec(p.proper_velocity_m_s),p.charge_C/p.mass_kg});
    const auto begin=Clock::now();
    for(const auto& p:inputs) {
        const double sx=p.x/o.width_m*o.nx,sy=p.y/o.gap_m*o.ny;
        const int ix=std::min(int(sx),o.nx-1),iy=std::min(int(sy),o.ny-1),e=iy*o.nx+ix;
        mfem::IntegrationPoint ip;ip.Set2(sx-ix,sy-iy);
        auto [ex,ey]=phi_gather.gradient(e,ip);PicVector E={-ex,-ey,0};
        auto a=lower.electric(p.x,p.y),b=upper.electric(p.x,p.y);
        for(int j=0;j<3;++j)E[j]+=(1-fraction)*a[j]+fraction*b[j];
        auto [ux,uy]=unit_gather.gradient(e,ip);PicVector extra={-ux,-uy,0};
        const double g=(p.z-o.short_z_m)/o.depth_m;
        for(int j=0;j<2;++j)extra[j]*=tables.voltage*(2*g-1);
        extra[2]=-2*tables.voltage/o.depth_m*(psi.GetValue(e,ip)-.5);
        for(int j=0;j<3;++j)E[j]+=extra[j];
        a=lower.magnetic(p.x,p.y);b=upper.magnetic(p.x,p.y);
        for(int j=0;j<3;++j)a[j]=(1-fraction)*a[j]+fraction*b[j];
        reference.push_back({vec(E),vec(a),vec(pic_boris(vec(p.u),p.qm,E,a,dt)),0});
    }
    const double reference_ms=std::chrono::duration<double,std::milli>(Clock::now()-begin).count();
    const auto shared_begin=Clock::now();for(const auto& p:inputs)actual.push_back(gp::evaluate(tables,p));
    const double shared_ms=std::chrono::duration<double,std::milli>(Clock::now()-shared_begin).count();
    const auto comparison=compare(reference,actual,state.particles);print_comparison(comparison);
    require(comparison.failures==0,"CPU fixture comparison failed");
    // Validate all domain/cell boundary guards independently of the saved sample.
    for(double x:{0.,o.width_m/o.nx,o.width_m})for(double y:{0.,o.gap_m/o.ny,o.gap_m}) {
        auto p=inputs.front();p.x=x;p.y=y;require(gp::evaluate(tables,p).status==0,"Domain endpoint rejected");
    }
    auto invalid=inputs.front();invalid.x=-1;require(gp::evaluate(tables,invalid).status==gp::outside_domain,"Outside input accepted");
    invalid=inputs.front();invalid.z=std::numeric_limits<double>::quiet_NaN();require(gp::evaluate(tables,invalid).status==gp::invalid_input,"NaN position accepted");
    std::filesystem::create_directories(out);
    std::vector<double> pack_ms,commit_ms;auto committed=state.particles;
    for(int trial=0;trial<9;++trial) {
        const auto pack_begin=Clock::now();
        for(std::size_t i=0;i<inputs.size();++i) {
            const auto& p=state.particles[i];
            inputs[i]={p.position_m[0],p.position_m[1],p.position_m[2],vec(p.proper_velocity_m_s),p.charge_C/p.mass_kg};
        }
        std::atomic_signal_fence(std::memory_order_seq_cst);
        const auto pack_end=Clock::now();
        // All-or-nothing host validation before any velocity is committed.
        for(const auto& r:actual)require(r.status==0 && gp::finite(r.E) && gp::finite(r.B) && gp::finite(r.u) && gp::finite(gp::norm2(r.u)),"Invalid result before commit");
        for(std::size_t i=0;i<actual.size();++i)committed[i].proper_velocity_m_s=vec(actual[i].u);
        std::atomic_signal_fence(std::memory_order_seq_cst);
        const auto commit_end=Clock::now();
        if(trial>=2) {
            pack_ms.push_back(std::chrono::duration<double,std::milli>(pack_end-pack_begin).count());
            commit_ms.push_back(std::chrono::duration<double,std::milli>(commit_end-pack_end).count());
        }
    }
    // Saved identity, position, charge, mass and macro weight remain host-owned.
    for(std::size_t i=0;i<committed.size();++i) {
        auto original=state.particles[i],pushed=committed[i];pushed.proper_velocity_m_s=original.proper_velocity_m_s;
        require(pushed.id==original.id && pushed.position_m==original.position_m && pushed.charge_C==original.charge_C &&
                pushed.mass_kg==original.mass_kg && pushed.weight==original.weight,"Host commit changed particle ownership or weight");
    }
    std::ofstream host_timing(out/"host-staging.json");host_timing<<std::setprecision(17)<<"{\"logical_cpu\":6,\"warmups\":2,\"pack_ms\":[";
    for(std::size_t i=0;i<pack_ms.size();++i)host_timing<<(i?",":"")<<pack_ms[i];
    host_timing<<"],\"finite_validation_and_velocity_commit_ms\":[";
    for(std::size_t i=0;i<commit_ms.size();++i)host_timing<<(i?",":"")<<commit_ms[i];
    host_timing<<"],\"packing\":\"saved GapParticle AoS to56B Inputs\",\"commit\":\"validate E/B/u/status for entire batch, then overwrite only proper velocity\"}\n";
    require(bool(host_timing),"Host staging timing write failed");
    write(out/"self_cells.bin",self_cells);write(out/"native_cells.bin",native_cells);
    write(out/"phi.bin",ph);write(out/"psi.bin",ps);write(out/"az0.bin",az0);write(out/"az1.bin",az1);
    write(out/"rate0.bin",rate0);write(out/"rate1.bin",rate1);write(out/"xs.bin",lower.xs);write(out/"ys.bin",lower.ys);
    write(out/"inputs.bin",inputs);write(out/"reference.bin",reference);write(out/"shared_cpu.bin",actual);write(out/"particles.bin",state.particles);
    std::ofstream manifest(out/"manifest.json");manifest<<std::setprecision(17)
        <<"{\n\"schema\":1,\"particles\":"<<inputs.size()<<",\"self_nx\":"<<o.nx<<",\"self_ny\":"<<o.ny
        <<",\"native_nx\":"<<tables.native.nx<<",\"native_ny\":"<<tables.native.ny
        <<",\"width\":"<<o.width_m<<",\"height\":"<<o.gap_m<<",\"depth\":"<<o.depth_m<<",\"short_z\":"<<o.short_z_m
        <<",\"fraction\":"<<fraction<<",\"voltage\":"<<tables.voltage<<",\"dt\":"<<dt
        <<",\"absolute_native_time_s\":"<<t<<",\"native_lower_time_s\":"<<times[k]<<",\"native_upper_time_s\":"<<times[k+1]
        <<",\"self_phi_and_unit_psi\":\"synthetic; checkpoint contains no potential; NOT a coupled replay qualification\""
        <<",\"particle_positions\":\"saved full-time positions; no midpoint drift performed\""
        <<",\"native_fields\":\"unchanged actual MFEM Az and dAz_dt snapshots\""
        <<",\"reference_host_ms_single_sample\":"<<reference_ms<<",\"shared_host_ms_single_sample\":"<<shared_ms
        <<",\"host_boris_bitwise_cases\":4096,\"geometry_sample_cases\":432,\"gpu_executed\":false"
        <<",\"gates\":{\"E_atol\":"<<E_abs<<",\"B_atol\":"<<B_abs<<",\"u_atol\":"<<u_abs
        <<",\"field_rtol\":"<<field_rel<<",\"u_rtol\":"<<u_rel<<",\"weighted_relative_limit\":"<<weighted_limit<<"}\n}\n";
    require(bool(manifest),"Manifest write failed");
    std::cout<<"Fixture written: "<<inputs.size()<<" saved particles; original host Boris bitwise cases4096; CPU only\n";
 }catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 2;}
}
