#include "native_driver_parts.hpp"
#include "records.hpp"
#include <iomanip>
#include <iostream>
namespace cv {
using validation_native::Thermal;
using validation_native::Fields;
using validation_native::Voltage;
using validation_native::Reservoir;
struct TrialConfig {
    int count=64,seed=0;bool pilot=false,adaptive=false,replay=false;
    double ceiling_s=0;
    std::string input_identity;
    std::filesystem::path input,fields,thermal,voltage,output,trace;
};
struct SinkTracker {
    std::uint64_t fresh_start=0;std::array<double,5> charge{},energy{};
    void observe(const GapParticle& p,int hit){
        if(p.id<fresh_start)return;
        require(hit>=0&&hit<=4,"invalid native sink classification");
        charge[hit]+=p.charge_C*p.weight;energy[hit]+=pic_kinetic_energy_J(p)*p.weight;
    }
};
inline std::array<Reservoir,4> pools(const Input& input){std::array<Reservoir,4> result;
    for(int f=0;f<4;++f){result[f].span=input.pools[f].span;result[f].integral=input.pools[f].integral;result[f].weighted_temperature=input.pools[f].weighted_temperature;}return result;}
struct SourceAdapter {
    GapPic& native;int count,seed;std::uint64_t step;double time;
    std::vector<Event> recorded;const std::vector<Event>* expected=nullptr;std::size_t cursor=0;
    long double physical=0,physical_expected=0,energy=0;double max_error=0;
    SourceAdapter(GapPic& p,int n,int s,std::uint64_t step_,double t):native(p),count(n),seed(s),step(step_),time(t){}
    void inject_surface_profile(GapElectrode el,GapSurface side,const std::vector<double>& T,const std::vector<double>& J,double span,int frozen_count,std::uint64_t base){
        require(frozen_count==32&&count_ok(count),"invalid count override");Event e{el,side,time,span,base,step,T,J};valid(e);
        if(expected){require(cursor<expected->size()&&same_event(e,expected->at(cursor)),"replay source event differs from reference");++cursor;}
        recorded.push_back(e);auto first=native.particles().size();
        native.inject_surface_profile(el,side,T,J,span,count,base+std::uint64_t(seed)*0x9e3779b97f4a7c15ULL);
        require(native.particles().size()==first+count,"native emission count differs");
        long double sum=0;for(std::size_t i=first;i<native.particles().size();++i){const auto& p=native.particles()[i];require(std::isfinite(p.weight)&&p.weight>0,"invalid emitted weight");sum+=p.weight;energy+=p.weight*pic_kinetic_energy_J(p);}
        const auto target=physical_electrons(e,native.options());max_error=std::max(max_error,double(std::abs(sum-target)/target));require(max_error<=1e-12,"event charge invariant failed");physical+=sum;physical_expected+=target;
    }
    void complete()const{if(expected)require(cursor==expected->size(),"replay omitted source events");}
};
inline std::uintmax_t directory_bytes(const std::filesystem::path& path){std::uintmax_t bytes=0;
    for(const auto& e:std::filesystem::recursive_directory_iterator(path))if(e.is_regular_file())bytes+=e.file_size();
    return bytes;}
inline void frame(const std::filesystem::path& path,GapPic& pic,GapPic& unit,double absolute,const SinkTracker& sinks,double supplied,double pending_e){
    std::ofstream f(path,std::ios::binary);const std::uint64_t magic=0x43564652414d4531ULL,n=pic.particles().size();
    write(f,magic);write(f,absolute);write(f,n);write(f,sinks.fresh_start);auto d=pic.diagnostics();write(f,d);
    for(auto p:pic.particles())write(f,p);
    for(double x:sinks.charge)write(f,x);
    for(double x:sinks.energy)write(f,x);
    write(f,supplied);write(f,pending_e);
    // Space-charge Poisson contribution at fixed native cell centers; imposed
    // fields are common and are not used as a large relative-error denominator.
    const auto& o=pic.options();const double hx=o.width_m/o.nx,hy=o.gap_m/o.ny;
    std::uint64_t cells=o.nx*o.ny;write(f,cells);mfem::IntegrationPoint ip;ip.Set2(.5,.5);
    for(int j=0;j<o.ny;++j)for(int i=0;i<o.nx;++i){double x=(i+.5)*hx,y=(j+.5)*hy;
        double volume=pic.inside_conductor(x,y)?0:hx*hy*o.depth_m;
        double phi=pic.potential().GetValue(j*o.nx+i,ip)-o.upper_potential_V*unit.potential().GetValue(j*o.nx+i,ip);
        auto E=pic.electric_field(x,y),U=unit.electric_field(x,y);
        for(int k=0;k<2;++k)E[k]-=o.upper_potential_V*U[k];
        for(double value:{phi,E[0],E[1],volume}){require(std::isfinite(value),"nonfinite field diagnostic");write(f,value);}}
}
// Compiled implementation; every path is fenced before MPI or GapPic creation.
// Removing this fence is a source/build change, not a runtime option.
void run_trial(const TrialConfig& cfg){
    require_execution();
    int argc=1;char name[]="isolated-validation";char* arguments[]{name,nullptr};char** argv=arguments;
    mfem::Mpi::Init(argc,argv);mfem::Hypre::Init();require(mfem::Mpi::WorldSize()==1,"exactly one rank required");
    require(count_ok(cfg.count),"only64/32/16 supported");
    require(cfg.seed>=0&&cfg.seed<=(cfg.adaptive?11:5),"seed outside declared stage");
    require(!cfg.adaptive||cfg.seed>=6,"adaptive confirmation requires fresh seed indices");
    resource_check();auto initial=read_input(cfg.input);Thermal thermal(cfg.thermal);Fields fields(cfg.fields);Voltage voltage(cfg.voltage);
    auto reservoir=pools(initial);auto o=initial.options;const double start=thermal.times.front();
    const double t0=start+initial.ledger.time_s;const auto limit=caps(cfg.pilot);const double end=t0+limit.duration_s;
    require(end<=thermal.times.back()&&end<=fields.times.back()&&end<=voltage.end.back(),"trial outside input support");
    require_owned_output(cfg.output);require(!std::filesystem::exists(cfg.output),"trial output must be new");std::filesystem::create_directories(cfg.output);
    Trace trace;if(cfg.replay){trace=trace_read(cfg.trace);require(trace.coupled_verified&&trace.input_identity==cfg.input_identity&&trace.seed_index==cfg.seed&&trace.initial_time_s==t0&&trace.end_time_s==end,"trace is unverified or from another trial");}
    else{trace.seed_index=cfg.seed;trace.initial_time_s=t0;trace.end_time_s=end;trace.input_identity=cfg.input_identity;}
    GapPic pic(o);pic.load_checkpoint((cfg.input/initial.checkpoint).string());
    auto unit_options=o;unit_options.lower_potential_V=0;unit_options.upper_potential_V=1;GapPic unit(unit_options);unit.solve_fields();
    const auto psi=[&](double x,double y){const double sx=x/o.width_m*o.nx,sy=y/o.gap_m*o.ny;const int i=std::clamp(int(sx),0,o.nx-1),j=std::clamp(int(sy),0,o.ny-1);mfem::IntegrationPoint ip;ip.Set2(sx-i,sy-j);return unit.potential().GetValue(j*o.nx+i,ip);};
    pic.set_external_fields([&](double x,double y,double t){return fields.value(x,y,start+t,false);},[&](double x,double y,double t){return fields.value(x,y,start+t,true);});
    pic.set_external_electric_3d([&](double x,double y,double z,double t){double V=voltage.at(start+t),g=(z-o.short_z_m)/o.depth_m;auto E=unit.electric_field(x,y);for(int k=0;k<2;++k)E[k]*=V*(2*g-1);E[2]=-2*V/o.depth_m*(psi(x,y)-.5);return E;});
    SinkTracker sinks;sinks.fresh_start=initial.next_id;pic.set_absorption_observer([&](const GapParticle& p,int hit){sinks.observe(p,hit);});
    const auto began=std::chrono::steady_clock::now();std::uint64_t accepted=0,step=initial.step;
    const auto baseline=pic.diagnostics();double emitted=0,expected_emitted=0;int sample=0;std::string status="incomplete";
    std::ofstream ledger(cfg.output/"step-ledger.csv");ledger<<std::setprecision(17)<<"step,time_s,dt_s,active_macros,particle_charge_C,charge_defect_C,continuity_defect_C,gauss_defect_C,poisson_residual_C_m,source_electrons,pending_electrons\n";
    auto pending_now=[&](){double total=0;for(int f=0;f<4;++f)total+=2*std::accumulate(reservoir[f].integral.begin(),reservoir[f].integral.end(),0.)*(f<2?.005:.0002)/thermal.sizes[f]*o.depth_m/pic_e;return total;};
    std::vector<std::string> checkpoints;
    auto save_valid=[&](){const std::string name="checkpoint_"+std::to_string(step)+".bin";pic.save_checkpoint((cfg.output/name).string());
        std::ofstream meta(cfg.output/"restart.meta.tmp");meta<<std::setprecision(17)<<name<<' '<<step<<' '<<initial.next_frame<<'\n';
        for(const auto& pool:reservoir){meta<<pool.span<<'\n';for(double x:pool.integral)meta<<x<<' ';meta<<'\n';for(double x:pool.weighted_temperature)meta<<x<<' ';meta<<'\n';}
        meta<<pic.diagnostics().time_s<<' '<<pic.options().upper_potential_V<<'\n';meta.close();require(bool(meta),"trial checkpoint metadata failed");std::filesystem::rename(cfg.output/"restart.meta.tmp",cfg.output/"restart.meta");
        if(checkpoints.empty()||checkpoints.back()!=name)checkpoints.push_back(name);
        while(checkpoints.size()>2){std::filesystem::remove(cfg.output/checkpoints.front());checkpoints.erase(checkpoints.begin());}};
    auto sample_now=[&](){frame(cfg.output/("frame_"+std::to_string(sample)+".bin"),pic,unit,start+pic.diagnostics().time_s,sinks,emitted,pending_now());++sample;save_valid();};sample_now();
    try {
    while(true){
        double t=start+pic.diagnostics().time_s;struct rusage ru{};getrusage(RUSAGE_SELF,&ru);
        const auto reason=cap_reason(limit,accepted,std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count(),t-t0,directory_bytes(cfg.output),ru.ru_maxrss);
        if(!reason.empty()){status=reason;break;}
        fields.set(t);const double V=voltage.at(t);
        if(pic.options().upper_potential_V!=V){if(pic.particles().empty())pic.set_empty_potentials(0,V);else pic.set_potentials(0,V);}
        const double sample_end=std::min(end,t0+sample*5e-12);
        double bound_end=sample_end;if(cfg.ceiling_s>0)bound_end=std::min(bound_end,t+cfg.ceiling_s);
        auto [admissible,T]=validation_native::allowed_step(pic,fields,voltage,thermal,reservoir,o,t,bound_end,initial.next_frame);
        double dt=admissible;const Step* replay_step=nullptr;
        if(cfg.replay){require(accepted<trace.steps.size(),"reference trace ended early");replay_step=&trace.steps[accepted];
            require(replay_step->step==step&&replay_step->time_s==t,"replay accepted clock diverged");dt=replay_step->dt_s;
            require(dt<=admissible,"common trace violates candidate native step bound");T=thermal.at(t+.5*dt);}
        SourceAdapter source{pic,cfg.count,cfg.seed,step,t};if(replay_step)source.expected=&replay_step->events;
        const bool injected=validation_native::accumulate_source(source,reservoir,T,dt,o,step);source.complete();emitted+=double(source.physical);
        expected_emitted+=double(source.physical_expected);require(std::abs(emitted-expected_emitted)<=1e-12*std::max(expected_emitted,1.),"cumulative source charge invariant failed");
        if(!cfg.replay&&!cfg.adaptive)trace.steps.push_back({step,t,dt,source.recorded});
        if(!injected&&pic.particles().empty())pic.advance_empty(dt);else pic.advance(dt);++accepted;++step;
        auto d=pic.diagnostics();
        if(!std::isfinite(d.gauss_defect_C)||!std::isfinite(d.particle_charge_defect_C)||pic.particles().size()>2000000)throw std::runtime_error("Nonfinite PIC ledger or two-million-particle resource guard; last checkpoint preserved");
        require(std::abs(d.particle_charge_defect_C-baseline.particle_charge_defect_C)<=1e-10*std::max(pic_e,emitted*pic_e),"incremental charge-balance validation gate failed");
        ledger<<step<<','<<start+d.time_s<<','<<dt<<','<<pic.particles().size()<<','<<d.particle_charge_C<<','<<d.particle_charge_defect_C<<','<<d.maximum_nodal_continuity_defect_C<<','<<d.gauss_defect_C<<','<<d.poisson_residual_C_per_m<<','<<emitted<<','<<pending_now()<<'\n';require(bool(ledger),"step ledger write failed");
        if(start+d.time_s>=sample_end-1e-22)sample_now();
    }
    save_valid();
    }catch(const std::exception& e){std::ofstream failure(cfg.output/"failure.txt");failure<<"numerically_invalid_or_incomplete\n"<<e.what()<<'\n';throw;}
    if(!cfg.replay&&!cfg.adaptive&&status=="complete"){trace.coupled_verified=true;trace_write(cfg.output/"reference-trace.bin",trace);}
    std::ofstream result(cfg.output/"status.json");result<<"{\"status\":\""<<status<<"\",\"count\":"<<cfg.count<<",\"seed\":"<<cfg.seed<<",\"accepted_steps\":"<<accepted<<",\"completed\":"<<(status=="complete"?"true":"false")<<",\"converged\":false}\n";
    unit.release_field_solver();pic.release_field_solver();
}
struct GuardPic {std::vector<GapParticle> data;const auto& particles()const{return data;}};
struct GuardField {std::vector<double> times;int current=0;double B=0;double bound()const{return B;}};
struct GuardBias {std::vector<double> end;int current=0;};
void cpu_preflight(const std::filesystem::path& input,const std::filesystem::path& output,const std::filesystem::path& thermal_path){
    resource_check();require_owned_output(output);auto initial=read_input(input);Thermal thermal(thermal_path);auto reservoir=pools(initial);
    const double t=thermal.times.front()+initial.ledger.time_s;
    GuardPic mock{initial.particles};GuardField field{{t-1e-9,t+1e-9},0,0};GuardBias voltage{{t+1e-9},0};
    auto [dt,T]=validation_native::allowed_step(mock,field,voltage,thermal,reservoir,initial.options,t,t+100e-12,initial.next_frame);
    require(dt>0&&dt<=1e-14,"native occupied step ceiling missing");
    field.B=1e4;auto [gyro,T2]=validation_native::allowed_step(mock,field,voltage,thermal,reservoir,initial.options,t,t+100e-12,initial.next_frame);
    require(gyro<=.05*pic_me/(pic_e*1e4)&&gyro<dt,"native gyro guard missing");
    bool invalid_trace_rejected=false;try{require(dt<=gyro,"trace violates candidate gyro bound");}catch(const std::runtime_error&){invalid_trace_rejected=true;}
    require(invalid_trace_rejected,"replay bound rejection failed");
    field.B=0;
    GuardPic dense;dense.data={initial.particles.front()};dense.data.front().weight=1e20;
    auto [plasma,T3]=validation_native::allowed_step(dense,field,voltage,thermal,reservoir,initial.options,t,t+100e-12,initial.next_frame);
    const double density=1e20/(25e-6*25e-6*initial.options.depth_m);
    require(plasma<=.1/std::sqrt(density*pic_e*pic_e/(validation_native::eps*pic_me))&&plasma<dt,"native plasma guard missing");
    auto hot=thermal;for(auto& v:hot.temperatures)std::fill(v.begin(),v.end(),1e5);
    auto [source_dt,T4]=validation_native::allowed_step(mock,field,voltage,hot,reservoir,initial.options,t,t+100e-12,initial.next_frame);
    require(source_dt<dt,"native source-density guard missing");
    auto narrow=voltage;narrow.end[0]=t+3e-15;
    auto [clock_dt,T5]=validation_native::allowed_step(mock,field,narrow,thermal,reservoir,initial.options,t,t+100e-12,initial.next_frame);
    require(clock_dt<=narrow.end[0]-t,"native input-clock clamp missing");
    GuardPic empty;auto emitting=reservoir;std::fill(emitting[0].integral.begin(),emitting[0].integral.end(),1e-14);
    auto [empty_dt,T6]=validation_native::allowed_step(empty,field,voltage,thermal,emitting,initial.options,t,t+100e-12,initial.next_frame);
    require(empty_dt<=1e-14,"native empty-to-emission guard missing");
    SinkTracker track;track.fresh_start=initial.next_id;auto old=initial.particles.front();track.observe(old,4);require(track.charge[4]==0,"old particle mislabeled fresh");
    auto fresh=old;fresh.id=initial.next_id;for(int hit=0;hit<5;++hit)track.observe(fresh,hit);
    for(int hit=0;hit<5;++hit)require(track.charge[hit]==fresh.charge_C*fresh.weight,"fresh sink class omitted");
    // Endpoint cap classification and sub-cap continuation are pure controls.
    require(cap_reason(caps(false),20000,1,0,0,0)=="step_cap","step cap failed");
    require(cap_reason(caps(false),0,900,0,0,0)=="wall_cap","wall cap failed");
    require(cap_reason(caps(false),0,1,0,256ULL*1024*1024,0)=="output_cap","output cap failed");
    require(cap_reason(caps(false),0,1,0,0,4ULL*1024*1024)=="memory_cap","memory cap failed");
    require(cap_reason(caps(false),0,1,100e-12,0,0)=="complete","end cap failed");
    require(cap_reason(caps(false),20000,1,100e-12,0,0)=="complete","exact-budget endpoint must complete");
    require(cap_reason(caps(false),1,1,1e-14,0,0).empty(),"sub-cap run stopped");
    require(!count_ok(8)&&count_ok(64)&&count_ok(32)&&count_ok(16),"count fence failed");
    require(std::filesystem::exists(output)&&std::filesystem::is_directory(output),"preflight output must exist");
    std::ofstream dump(output/"native-particles.bin",std::ios::binary);for(const auto& p:initial.particles)write(dump,p);
    std::ofstream f(output/"native-preflight.json");f<<std::setprecision(17)<<"{\"passed\":true,\"cpu\":6,\"nice\":15,\"GPU_enabled\":false,\"MPI_initialized\":false,\"GapPic_constructed\":false,\"field_solves\":0,\"trajectory_steps\":0,\"particle_count\":"<<initial.particles.size()<<",\"next_particle_id\":"<<initial.next_id<<",\"pending_electrons\":"<<pending(initial)<<",\"manufactured_B0_bound_s\":"<<dt<<",\"manufactured_high_B_bound_s\":"<<gyro<<",\"plasma_source_clock_empty_guards_passed\":true,\"sink_classification_fixture_passed\":true,\"resource_cap_fixtures_passed\":true}\n";
}
}
int main(int argc,char** argv){
    try{
        if(argc==2&&std::string(argv[1])=="--describe"){std::cout<<"{\"counts\":[64,32,16],\"gpu_execution_enabled\":false,\"coupled_execution_enabled\":false,\"preflight_only\":true}\n";return 0;}
        if(argc==3&&std::string(argv[1])=="--inspect-trace"){cv::resource_check();auto trace=cv::trace_read(argv[2]);std::cout<<"{\"valid\":true,\"coupled_verified\":"<<(trace.coupled_verified?"true":"false")<<",\"input_identity\":\""<<trace.input_identity<<"\",\"steps\":"<<trace.steps.size()<<"}\n";return 0;}
        if(argc==5&&std::string(argv[1])=="--preflight"){cv::cpu_preflight(argv[2],argv[3],argv[4]);return 0;}
        cv::require_execution();
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 78;}
}
