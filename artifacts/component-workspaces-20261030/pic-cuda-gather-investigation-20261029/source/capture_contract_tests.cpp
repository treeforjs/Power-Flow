#include "capture_hooks.hpp"
#include "snapshot_io.hpp"
#include <cstring>
#include <iostream>

using namespace raikou::mfem_backend;
struct Counts {int E=0,axial=0,B=0;};
void fields(GapPic& pic,Counts& count) {
    pic.set_external_fields([&](double,double,double){++count.E;return PicVector{1,2,3};},
                            [&](double,double,double){++count.B;return PicVector{.1,.2,.3};});
    pic.set_external_electric_3d([&](double,double,double,double){++count.axial;return PicVector{4,5,6};});
}
void native_writer(const void* context,const std::filesystem::path&,double,std::ostream& manifest) {
    const auto& count=*static_cast<const Counts*>(context);
    snapshot_io::require(count.E==8 && count.axial==8 && count.B==8,"Diagnostic made extra field callbacks");
    manifest<<",\"test_only_mock_native_export\":true";
}
int main(int argc,char** argv) {
    try {
        snapshot_io::require(argc==2,"Expected fresh test output directory");
        GapPicOptions options;options.nx=4;options.ny=5;options.upper_potential_V=1;
        GapPic reference(options),captured(options),unit(options);
        for(int i=0;i<8;++i) {
            GapParticle p;p.position_m={options.width_m*(.2+.08*i),.4*options.gap_m,0};
            p.proper_velocity_m_s={1e6,2e6,3e6};p.weight=1e-9*(i+1);
            reference.add_particle(p);captured.add_particle(p);
        }
        Counts original,observed;fields(reference,original);fields(captured,observed);
        reference.advance(5e-15);unit.solve_fields();
        pic_capture::configure(argv[1],unit.potential(),0,&observed,native_writer);
        captured.advance(5e-15);
        snapshot_io::require(original.E==observed.E && original.axial==observed.axial && original.B==observed.B,"Callback counts changed");
        snapshot_io::require(reference.particles().size()==captured.particles().size() &&
            std::memcmp(reference.particles().data(),captured.particles().data(),reference.particles().size()*sizeof(GapParticle))==0,"Capture changed particle bits");
        const auto a=reference.diagnostics(),b=captured.diagnostics();
        snapshot_io::require(std::memcmp(&a,&b,sizeof(a))==0,"Capture changed diagnostic bits");
        snapshot_io::require(std::memcmp(reference.potential().HostRead(),captured.potential().HostRead(),reference.potential().Size()*sizeof(double))==0,"Capture changed potential bits");
        snapshot_io::require(reference.horizontal_transfers_C()==captured.horizontal_transfers_C() && reference.vertical_transfers_C()==captured.vertical_transfers_C(),"Capture changed continuity transfers");
        std::cout<<"{\"passed\":true,\"GPU_executed\":false,\"capture_on_off_full_step_bitwise_equal\":true,\"extra_field_callbacks\":0}\n";
    }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
