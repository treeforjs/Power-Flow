#include "evolving_force.hpp"
#include "cuda_driver_mock.hpp"
#include "bench_profile.hpp"
// Original provider code copied verbatim in already accepted timing evidence.
#include "/home/treefor/projects/runtimes/pic-cuda-gather-investigation-20261029/backend-v2/host-timing/legacy_provider_reference.hpp"
#include <iostream>
#include <limits>
using namespace raikou::mfem_backend;
using namespace pic_batch;
namespace io=snapshot_io;namespace gp=gather_push;namespace mc=mock_cuda;
namespace {
constexpr double seed_time=2.6478441131125727e-08;
struct Scalar {
    mfem::Mesh mesh;mfem::H1_FECollection fec;mfem::FiniteElementSpace space;mfem::GridFunction field;Q1ScalarGather gather;
    Scalar(const HostSnapshot& s,const std::vector<double>& data)
        :mesh(mfem::Mesh::MakeCartesian2D(s.parameters.self_nx,s.parameters.self_ny,mfem::Element::QUADRILATERAL,true,s.parameters.width,s.parameters.height,false)),
         fec(1,2),space(&mesh,&fec),field(&space),gather(field){io::require(std::size_t(field.Size())==data.size(),"Scalar shape");std::copy(data.begin(),data.end(),field.HostWrite());}
    std::pair<PicVector,double> at(const GapParticle& p,const Parameters& params) {
        const double sx=p.position_m[0]/params.width*params.self_nx,sy=p.position_m[1]/params.height*params.self_ny;
        const int ix=std::min(int(sx),params.self_nx-1),iy=std::min(int(sy),params.self_ny-1);
        mfem::IntegrationPoint ip;ip.Set2(sx-ix,sy-iy);const int e=iy*params.self_nx+ix;
        const auto grad=gather.gradient(e,ip);return {{-grad[0],-grad[1],0},field.GetValue(e,ip)};
    }
};
pic_evolving::NativeSelection select(native_reference::Fields& f,native_reference::Voltage& v,Scalar& unit,double t) {
    f.set(t);const double V=v.at(t);
    return {&unit.field,f.lower->az.get(),f.upper->az.get(),f.lower->rate.get(),f.upper->rate.get(),&f.lower->elements,&f.lower->xs,&f.lower->ys,
        1,1,std::uint64_t(f.current),std::uint64_t(f.current+1),0,t,f.times[f.current],f.times[f.current+1],f.fraction,V};
}
}
int main(int argc,char** argv) {
    try {
        io::require(argc==4,"Expected CASE ACTUAL_FIXTURE OUTPUT");const std::string mode=argv[1];const std::filesystem::path out=argv[3];std::filesystem::create_directories(out);
        mc::reset();
        if(mode=="empty") {
            GapPicOptions options;options.nx=options.ny=2;GapPic pic(options);int providers=0,callbacks=0;
            pic.set_external_fields([&](double,double,double){++callbacks;return PicVector{};},[&](double,double,double){++callbacks;return PicVector{};});
            pic.set_external_electric_3d([&](double,double,double,double){++callbacks;return PicVector{};});
            pic_evolving::bound_restart(0,seed_time,seed_time+1e-12,437411);
            pic_evolving::configure(pic,pic_evolving::Mode::OptimizedCandidate,"not-used.ptx",out,[&](double){++providers;return pic_evolving::NativeSelection{};});
            pic.advance(1e-15);
            io::require(pic.particles().empty() && pic.diagnostics().time_s==1e-15 && providers==0 && callbacks==0 && mc::calls["init"]==0 && mc::calls["launch"]==0,
                        "Empty midpoint advanced provider/callback/CUDA or changed clock rule");
            io::require(pic_bench::entries["solve_fields_total"].calls==2,"Empty advance skipped either guarded field solve");
            pic_evolving::observe_step(pic,437412,1e-15);
            std::array<double,4> spans{};std::array<std::span<const double>,4> empty{};pic_evolving::observe_sources(spans,empty,empty);
            pic_evolving::shutdown();
        } else {
            auto s=HostSnapshot::load(argv[2],true);Scalar self(s,s.phi),unit(s,s.psi);
            native_reference::Fields fields("/home/treefor/projects/Power-Flow/runs/restore_field_mpi8_domain40mm");
            native_reference::Voltage voltage("/home/treefor/sim-workspace/powerflow/outputs/gap-pulse-midplane110kV.csv");
            constexpr double midpoint=2.6478443519271536e-08;
            fields.set(midpoint);
            if(mode=="frame-boundary") {
                const int before=fields.current;const double boundary=fields.times.at(before+1);
                const double fudge=64*std::numeric_limits<double>::epsilon()*std::max(std::abs(boundary),1e-12);
                const std::array<double,5> times={boundary-fudge,boundary,std::nextafter(boundary,0.),std::nextafter(boundary,std::numeric_limits<double>::infinity()),boundary+fudge};
                for(std::size_t t=0;t<times.size();++t) {
                    const auto n=select(fields,voltage,unit,times[t]);
                    io::require(fields.current==(t==0?before:before+1),"Original epsilon/upper_bound frame boundary rule changed");
                    s.native_cells=io::export_grid(*n.az0,*n.elements);s.az0=io::values(*n.az0);s.az1=io::values(*n.az1);s.rate0=io::values(*n.rate0);s.rate1=io::values(*n.rate1);
                    s.xs=*n.xs;s.ys=*n.ys;s.parameters.fraction=n.fraction;s.parameters.voltage=n.voltage;
                    std::vector<gp::Outputs> original,shared;std::vector<GapParticle> particles;
                    for(std::size_t j=0;j<16;++j) {
                        const auto index=j*(s.particles.size()-1)/15;const auto& p=s.particles[index];
                        auto E=self.at(p,s.parameters).first;auto extra=fields.value(p.position_m[0],p.position_m[1],times[t],false);
                        for(int k=0;k<3;++k)E[k]+=extra[k];
                        auto [axial,psi]=unit.at(p,s.parameters);const double g=(p.position_m[2]-s.parameters.short_z)/s.parameters.depth;
                        for(int k=0;k<2;++k)axial[k]*=n.voltage*(2*g-1);
                        axial[2]=-2*n.voltage/s.parameters.depth*(psi-.5);for(int k=0;k<3;++k)E[k]+=axial[k];
                        const auto B=fields.value(p.position_m[0],p.position_m[1],times[t],true);
                        original.push_back({io::vec(E),io::vec(B),io::vec(pic_boris(p.proper_velocity_m_s,p.charge_C/p.mass_kg,E,B,s.parameters.dt)),0});
                        shared.push_back(gp::evaluate(s.tables(),s.inputs[index]));particles.push_back(p);
                    }
                    const auto c=batch_compare::compare(original,shared,particles);io::require(c.failures==0,"CPU frame-boundary original native/axial/Boris parity failed");
                    io::require(fields.last_time==times[t] && voltage.last_time==times[t],"Provider time advanced unexpectedly");
                }
                io::require(mc::calls["init"]==0 && mc::calls["launch"]==0,"Frame CPU parity initialized real/mock CUDA");
            } else if(mode=="native-epoch-refresh") {
                GapPicOptions options;options.nx=options.ny=360;options.width_m=options.gap_m=.009;options.depth_m=.025;options.short_z_m=-.0125;
                GapPic pic(options);pic_evolving::bound_restart(0,seed_time,seed_time+1e-12,437411);bool bad_unit_epoch=false;
                pic_evolving::configure(pic,pic_evolving::Mode::OptimizedCandidate,"/home/treefor/projects/runtimes/pic-cuda-gather-investigation-20261029/build/gather_push.ptx",out,
                    [&](double t){auto n=select(fields,voltage,unit,t);if(bad_unit_epoch)n.unit_epoch=2;return n;});
                const double boundary=fields.times.at(fields.current+1);const std::array<double,2> times={midpoint,boundary+1e-21};
                for(double t:times) {
                    auto live=s.particles;pic_evolving::force(&pic,live,self.field,options,t,s.parameters.dt,{});
                    const auto n=select(fields,voltage,unit,t);
                    s.native_cells=io::export_grid(*n.az0,*n.elements);s.az0=io::values(*n.az0);s.az1=io::values(*n.az1);s.rate0=io::values(*n.rate0);s.rate1=io::values(*n.rate1);
                    s.parameters.fraction=n.fraction;s.parameters.voltage=n.voltage;
                    for(std::size_t i=0;i<live.size();++i) {
                        const auto expected=gp::evaluate(s.tables(),s.inputs[i]);
                        for(int k=0;k<3;++k)io::require(std::abs(live[i].proper_velocity_m_s[k]-expected.u.v[k])<=batch_compare::u_abs+batch_compare::u_rel*std::abs(expected.u.v[k]),"Native epoch refresh used stale uploaded fields");
                    }
                }
                io::require(mc::calls["context_sync"]==2 && mc::calls["launch"]==2,"Frame epoch did not explicitly refresh static data exactly once");
                auto live=s.particles;const auto before=live;bad_unit_epoch=true;bool rejected=false;
                try{pic_evolving::force(&pic,live,self.field,options,times.back(),s.parameters.dt,{});}catch(const std::exception&){rejected=true;}
                io::require(rejected && mc::calls["launch"]==2 && std::memcmp(live.data(),before.data(),live.size()*sizeof(GapParticle))==0,"Unit epoch mismatch dispatched or partially committed");
                pic_evolving::shutdown("expected CPU driver-mock unit epoch rejection");
                io::require(mc::live_allocations()==0 && mc::unsafe_releases()==0 && mc::primary_refs()==1,"Epoch rejection leaked or disturbed borrowed primary");
            } else if(mode=="failed-shadow-transaction") {
                GapPicOptions options;options.nx=options.ny=360;options.width_m=options.gap_m=.009;options.depth_m=.025;options.short_z_m=-.0125;
                GapPic pic(options);pic_evolving::bound_restart(0,seed_time,seed_time+1e-12,437411);
                pic_evolving::configure(pic,pic_evolving::Mode::Shadow,"/home/treefor/projects/runtimes/pic-cuda-gather-investigation-20261029/build/gather_push.ptx",out,
                    [&](double t){return select(fields,voltage,unit,t);});
                auto live=s.particles;const auto before=live;auto wrong=s.reference;wrong.back().u.v[0]+=1e8;bool rejected=false;
                try{pic_evolving::force(&pic,live,self.field,options,midpoint,s.parameters.dt,wrong);}catch(const std::exception& error){rejected=std::string(error.what()).find("shadow independent")!=std::string::npos;}
                io::require(rejected && std::memcmp(live.data(),before.data(),live.size()*sizeof(GapParticle))==0,"Failed shadow gate partially committed particles");
                bool retry_rejected=false;try{pic_evolving::begin_force(&pic,midpoint,s.parameters.dt);}catch(const std::exception&){retry_rejected=true;}
                io::require(retry_rejected,"Failed kick permitted retry/native fallback");
                pic_evolving::shutdown("expected CPU driver-mock shadow rejection");
                io::require(mc::live_allocations()==0 && mc::unsafe_releases()==0 && mc::primary_refs()==1 && !std::filesystem::exists(out/"restart.meta") && std::filesystem::file_size(out/"evolving-state.bin")==0,
                            "Failed shadow leaked resources or saved partial step");
            } else throw std::invalid_argument("Unknown CPU parity case");
        }
        std::cout<<"{\"passed\":true,\"case\":"<<std::quoted(mode)<<",\"real_CUDA_initialized_or_executed\":false,\"driver_mock_only\":true}\n";
    }catch(const std::exception& error){try{pic_evolving::shutdown(error.what());}catch(...){}std::cerr<<error.what()<<'\n';return 1;}
}
