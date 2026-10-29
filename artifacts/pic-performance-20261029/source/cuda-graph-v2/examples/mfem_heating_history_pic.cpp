#include "native_heating_fields.hpp"
#include <array>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <limits>
#include <numeric>
#include "bench_profile.hpp"
#include "cudss_graph_record.hpp"
namespace {
constexpr double eps=8.8541878188e-12;
struct Voltage {
    std::vector<double> start,end,value;
    int current=0;
    explicit Voltage(const std::filesystem::path& path) {
        std::ifstream in(path);std::string line;std::getline(in,line);
        while(std::getline(in,line)){auto r=split(line);start.push_back(std::stod(r.at(0)));end.push_back(std::stod(r.at(1)));value.push_back(std::stod(r.at(2)));}
        if(start.empty())throw std::runtime_error("Missing pulsed voltage intervals");
        for(std::size_t i=0;i<start.size();++i)if(!std::isfinite(start[i]) || !std::isfinite(end[i]) || !std::isfinite(value[i]) || end[i]<=start[i] || value[i]>=0 || value[i] < -1e6 || (i && std::abs(start[i]-end[i-1])>1e-21))throw std::runtime_error("Invalid pulsed voltage intervals");
    }
    double at(double t) {
        if(t<start.front()-1e-20 || t>end.back()+1e-20)throw std::runtime_error("Voltage extrapolation prohibited");
        const double lookup=t+32*std::numeric_limits<double>::epsilon()*std::max(std::abs(t),1e-12);
        current=std::clamp(int(std::upper_bound(start.begin(),start.end(),lookup)-start.begin())-1,0,int(start.size())-1);
        return value[current];
    }
};
struct Thermal {
    std::vector<double> times;
    std::array<std::vector<double>,4> temperatures;
    std::array<std::uint64_t,4> sizes{};
    explicit Thermal(const std::filesystem::path& path) {
        std::ifstream in(path,std::ios::binary);std::uint64_t magic=0,n=0;
        in.read(reinterpret_cast<char*>(&magic),8);in.read(reinterpret_cast<char*>(&n),8);
        in.read(reinterpret_cast<char*>(sizes.data()),32);
        if(!in || magic!=0x544845524d303031ULL || n<2 || n>10000)throw std::runtime_error("Invalid thermal history header");
        times.resize(n);in.read(reinterpret_cast<char*>(times.data()),8*n);
        for(int f=0;f<4;++f){if(sizes[f]<1 || sizes[f]>4096)throw std::runtime_error("Invalid thermal patch count");temperatures[f].resize(n*sizes[f]);in.read(reinterpret_cast<char*>(temperatures[f].data()),8*n*sizes[f]);}
        if(!in || in.peek()!=std::char_traits<char>::eof())throw std::runtime_error("Truncated thermal history");
        for(std::size_t i=0;i<times.size();++i)if(!std::isfinite(times[i]) || (i && times[i]<=times[i-1]))throw std::runtime_error("Invalid thermal clock");
        for(const auto& v:temperatures)for(double T:v)if(!std::isfinite(T) || T<=0 || T>1e6)throw std::runtime_error("Invalid absolute temperature");
    }
    std::array<std::vector<double>,4> at(double t) const {
        if(t<times.front()-1e-20 || t>times.back()+1e-20)throw std::runtime_error("Temperature extrapolation prohibited");
        auto k=std::clamp(int(std::upper_bound(times.begin(),times.end(),t)-times.begin())-1,0,int(times.size())-2);
        double q=std::clamp((t-times[k])/(times[k+1]-times[k]),0.,1.);std::array<std::vector<double>,4> result;
        for(int f=0;f<4;++f){result[f].resize(sizes[f]);for(std::size_t i=0;i<sizes[f];++i)result[f][i]=(1-q)*temperatures[f][k*sizes[f]+i]+q*temperatures[f][(k+1)*sizes[f]+i];}
        return result;
    }
};
struct Fields {
    std::filesystem::path directory;
    std::vector<double> times;
    std::vector<double> bounds;
    std::unique_ptr<NativeField> lower,upper;
    int current=-1;double fraction=0;
    explicit Fields(const std::filesystem::path& path):directory(path) {
        std::ifstream in(path/"em_field_snapshots.csv");std::string line;std::getline(in,line);
        while(std::getline(in,line))times.push_back(std::stod(split(line).at(0)));
        if(times.size()<2)throw std::runtime_error("Missing field history");bounds.resize(times.size(),-1);
    }
    static double maximum(NativeField& f) {
        double value=0;mfem::Vector g(2);
        // Q1 gradient norm attains its maximum at an element corner.
        for(int e=0;e<f.mesh.GetNE();++e)for(double x:{0.,1.})for(double y:{0.,1.}) {
            mfem::IntegrationPoint ip;ip.Set2(x,y);auto* tr=f.mesh.GetElementTransformation(e);tr->SetIntPoint(&ip);f.az->GetGradient(*tr,g);value=std::max(value,g.Norml2());
        }
        return value;
    }
    void set(double t) {
        if(t<times.front()-1e-20 || t>times.back()+1e-20)throw std::runtime_error("Field extrapolation prohibited");
        const double lookup=t+32*std::numeric_limits<double>::epsilon()*std::max(std::abs(t),1e-12);
        int k=std::clamp(int(std::upper_bound(times.begin(),times.end(),lookup)-times.begin())-1,0,int(times.size())-2);
        if(k!=current) {
            if(k==current+1 && upper)lower=std::move(upper);else lower=std::make_unique<NativeField>(directory,times[k]);
            upper=std::make_unique<NativeField>(directory,times[k+1]);current=k;
            if(bounds[k]<0)bounds[k]=maximum(*lower);if(bounds[k+1]<0)bounds[k+1]=maximum(*upper);
        }
        fraction=std::clamp((t-times[k])/(times[k+1]-times[k]),0.,1.);
    }
    double bound() const {return std::max(bounds[current],bounds[current+1]);}
    PicVector value(double x,double y,double t,bool magnetic) {
        set(t);auto a=magnetic?lower->magnetic(x,y):lower->electric(x,y);auto b=magnetic?upper->magnetic(x,y):upper->electric(x,y);
        for(int i=0;i<3;++i)a[i]=(1-fraction)*a[i]+fraction*b[i];return a;
    }
};
struct Reservoir {std::vector<double> integral,weighted_temperature;double span=0;};
void snapshot(GapPic& pic,const std::filesystem::path& out,int frame,double time) {
    const auto& o=pic.options();const double h=o.width_m/o.nx,volume=h*h*o.depth_m;std::vector<double> density(o.nx*o.ny,0.),midplane(o.nx*o.ny,0.);
    double represented=0;
    for(const auto& p:pic.particles()){int i=std::clamp(int(p.position_m[0]/h),0,o.nx-1),j=std::clamp(int(p.position_m[1]/h),0,o.ny-1);density[j*o.nx+i]+=p.weight;if(std::abs(p.position_m[2])<.0005)midplane[j*o.nx+i]+=p.weight;represented+=p.weight;}
    double mapped=std::accumulate(density.begin(),density.end(),0.);
    if(std::abs(mapped-represented)>1e-11*std::max(represented,1.))throw std::runtime_error("Snapshot lost represented charge");
    std::string stem="frame_"+std::to_string(frame);std::ofstream f(out/(stem+".csv.tmp"));f<<std::setprecision(17)<<"time_s,x_m,y_m,potential_V,Ex_V_m,Ey_V_m,electron_density_cm3,height_average_electron_density_cm3,solid\n";
    mfem::IntegrationPoint ip;ip.Set2(.5,.5);
    for(int j=0;j<o.ny;++j)for(int i=0;i<o.nx;++i){double x=(i+.5)*h,y=(j+.5)*h;auto E=pic.electric_field(x,y);f<<time<<','<<x<<','<<y<<','<<pic.potential().GetValue(j*o.nx+i,ip)<<','<<E[0]<<','<<E[1]<<','<<midplane[j*o.nx+i]/(h*h*.001)/1e6<<','<<density[j*o.nx+i]/volume/1e6<<','<<pic.inside_conductor(x,y)<<'\n';}
    f.close();if(!f)throw std::runtime_error("Frame write failed");std::filesystem::rename(out/(stem+".csv.tmp"),out/(stem+".csv"));
}
}
int main(int argc,char** argv) {
    mfem::Mpi::Init(argc,argv);mfem::Hypre::Init();
    try {
        if(mfem::Mpi::WorldRank()==0)pic_graph::RunRecord::instance().begin("heating_history_pic_gpu");
        if(argc!=6)throw std::invalid_argument("Usage: heating-history-pic FIELD_DIR THERMAL_BINARY OUTPUT END_NS VOLTAGE_INTERVAL_CSV");
        if(mfem::Mpi::WorldRank()!=0){pic_parallel_workers();pic_parallel_record("");return 0;}
        Thermal thermal(argv[2]);Fields fields(argv[1]);Voltage voltage(argv[5]);const std::filesystem::path out(argv[3]);
        const double start=thermal.times.front(),end=std::stod(argv[4])*1e-9;
        if(!std::isfinite(end) || end<=start || end>thermal.times.back()+1e-20 || start<fields.times.front()-1e-20 || end>fields.times.back()+1e-20)throw std::runtime_error("Full run outside supported input clocks");
        const bool resume=std::filesystem::exists(out/"restart.meta");
        if(std::filesystem::exists(out) && !resume)throw std::runtime_error("Output exists without a restart record");
        std::filesystem::create_directories(out);
        GapPicOptions o;o.nx=o.ny=360;o.width_m=o.gap_m=.009;o.depth_m=.025;o.finite_foils=true;o.distributed_fields=true;o.finite_axial_extent=true;
        o.lower_foil={.002,.007,.0018,.002};o.upper_foil={.002,.007,.007,.0072};o.upper_potential_V=voltage.at(start);
        GapPic pic(o);std::array<Reservoir,4> reservoir;
        for(int f=0;f<4;++f){reservoir[f].integral.resize(thermal.sizes[f],0.);reservoir[f].weighted_temperature.resize(thermal.sizes[f],0.);}
        std::uint64_t step=0;int next_frame=1;std::string checkpoint;
        if(resume) {
            std::ifstream meta(out/"restart.meta");meta>>checkpoint>>step>>next_frame;
            for(auto& pool:reservoir){meta>>pool.span;for(auto& v:pool.integral)meta>>v;for(auto& v:pool.weighted_temperature)meta>>v;}
            if(!meta)throw std::runtime_error("Invalid reservoir restart");double elapsed=0,stored_voltage=0;meta>>elapsed>>stored_voltage;if(!meta)throw std::runtime_error("Missing checkpoint clock");
            pic.set_potentials(0,stored_voltage);pic.load_checkpoint((out/checkpoint).string());
        }else pic.solve_fields();
        // Midplane of a25mm line: feed at+12.5mm, ideal short at-12.5mm.
        // Local gap bias is half the feed bias. Symmetric conductor potentials
        // give Ez_scalar=-V_feed/H*(psi_upper-.5) at the observation plane.
        auto unit_options=o;unit_options.lower_potential_V=0;unit_options.upper_potential_V=1;
        GapPic unit_bias(unit_options);unit_bias.solve_fields();
        const auto psi=[&](double x,double y) {
            const double sx=x/o.width_m*o.nx,sy=y/o.gap_m*o.ny;
            const int i=std::clamp(int(sx),0,o.nx-1),j=std::clamp(int(sy),0,o.ny-1);
            mfem::IntegrationPoint ip;ip.Set2(sx-i,sy-j);
            return unit_bias.potential().GetValue(j*o.nx+i,ip);
        };
        pic.set_external_fields([&](double x,double y,double t){return fields.value(x,y,start+t,false);},[&](double x,double y,double t){return fields.value(x,y,start+t,true);});
        pic.set_external_electric_3d([&](double x,double y,double z,double t){
            double V=voltage.at(start+t),g=(z-o.short_z_m)/o.depth_m;
            auto E=unit_bias.electric_field(x,y);for(int k=0;k<2;++k)E[k]*=V*(2*g-1);
            E[2]=-2*V/o.depth_m*(psi(x,y)-.5);return E;
        });
        std::ofstream scope(out/"scope.json");scope<<std::setprecision(17)<<"{\"status\":\"continuous_heating_history_pic_gpu_fields\",\"physical_prediction\":false,\"converged\":false,\"start_time_s\":"<<start<<",\"requested_end_time_s\":"<<end<<",\"nx\":360,\"ny\":360,\"cell_size_m\":0.000025,\"observation_plane_z_m\":0,\"target_height_m\":0.025,\"short_z_m\":-0.0125,\"feed_z_m\":0.0125,\"short_scope\":\"out-of-plane5x5mm conducting short absorbs particles; feed and exterior escapes recorded; full3D field fringing not resolved\",\"depth_m\":0.025,\"finite_axial_extent\":true,\"density_sampling_depth_m\":0.001,\"space_charge_model\":\"height-averaged2D Poisson with3D particle transport and imposed linear axial voltage profile\",\"mpi_ranks\":1,\"lower_potential_V\":0,\"voltage_model\":\"midpoint pulse half of220kV feed pulse; ideal linear axial drop to0atshort; symmetric scalar longitudinal gradient plus native inductive Ez\",\"initial_upper_potential_V\":"<<o.upper_potential_V<<",\"initial_electrons\":\"empty at first supported heating timestamp; prior electron population unknown\",\"thermal_source\":\"295 high resolution upper boundary states; mirrored lower reconstruction verified on small case\",\"fields\":\"linear temporal interpolation of coarse regenerated native Az curl and inductive Ez; conditional AK voltage\",\"source_law\":\"patch Richardson phi4.3eV ideal A, no Schottky or tunneling\",\"source_reservoir\":\"sub-electron face supplies accumulated before emission; less than one physical electron per face/electrode can be delayed\",\"step_rule\":\"gyroangle<=0.05, dt<=10fs with particles, plasma-frequency and particle-displacement limits; empty vacuum<=1ps\",\"limitations\":\"constant-property heating beyond material validity; no phase change, collisions, gas feedback or sheath/grid convergence\",\"density_units\":\"cm^-3\",\"checkpoint_format\":\"same frozen executable ABI only; input SHA checked by launcher\",\"parallelism\":\"root particles/assembly/output; one CPU rank; optional CUDA cuDSS field solve\"}\n";scope.close();
        std::ofstream history(out/"history.csv",resume?std::ios::app:std::ios::out);history<<std::setprecision(17);if(!resume)history<<"time_s,step,macro_particles,gauss_defect_C,particle_charge_defect_C,nodal_continuity_defect_C,field_energy_J,kinetic_energy_J,source_reservoir_electrons,short_absorbed_charge_C,short_absorbed_energy_J,axial_escaped_charge_C\n";
        auto wall=std::chrono::steady_clock::now(),saved=wall,reported=wall;
        const auto clock=[]{return std::chrono::steady_clock::now();};
        const auto save=[&] {
            std::string name="checkpoint_"+std::to_string(step)+".bin";pic.save_checkpoint((out/name).string());
            std::ofstream meta(out/"restart.meta.tmp");meta<<std::setprecision(17)<<name<<' '<<step<<' '<<next_frame<<'\n';
            for(const auto& pool:reservoir){meta<<pool.span<<'\n';for(double v:pool.integral)meta<<v<<' ';meta<<'\n';for(double v:pool.weighted_temperature)meta<<v<<' ';meta<<'\n';}
            meta<<pic.diagnostics().time_s<<' '<<pic.options().upper_potential_V<<'\n';meta.close();if(!meta)throw std::runtime_error("Reservoir checkpoint failed");std::filesystem::rename(out/"restart.meta.tmp",out/"restart.meta");
            if(!checkpoint.empty() && checkpoint!=name)std::filesystem::remove(out/checkpoint);checkpoint=name;saved=clock();
        };
        const auto progress=[&](bool done) {
            auto d=pic.diagnostics();double pending=0;
            for(int f=0;f<4;++f)pending+=2*std::accumulate(reservoir[f].integral.begin(),reservoir[f].integral.end(),0.)*(f<2?.005:.0002)/thermal.sizes[f]*o.depth_m/pic_e;
            history<<start+d.time_s<<','<<step<<','<<pic.particles().size()<<','<<d.gauss_defect_C<<','<<d.particle_charge_defect_C<<','<<d.maximum_nodal_continuity_defect_C<<','<<d.field_energy_J<<','<<d.particle_kinetic_energy_J<<','<<pending<<','<<d.short_absorbed_charge_C<<','<<d.short_absorbed_energy_J<<','<<d.axial_escaped_charge_C<<'\n';history.flush();
            std::ofstream p(out/"progress.json.tmp");p<<std::setprecision(17)<<"{\"phase\":\""<<(done?"completed":"running")<<"\",\"full_history_run\":"<<(end>=60e-9-1e-20?"true":"false")<<",\"checkpoint_benchmark\":"<<(end>=60e-9-1e-20?"false":"true")<<",\"accepted_time_s\":"<<start+d.time_s<<",\"start_time_s\":"<<start<<",\"requested_end_time_s\":"<<end<<",\"steps\":"<<step<<",\"macro_particles\":"<<pic.particles().size()<<",\"mpi_ranks\":1,\"upper_potential_V\":"<<pic.options().upper_potential_V<<",\"wall_time_since_restart_s\":"<<std::chrono::duration<double>(clock()-wall).count()<<",\"source_reservoir_electrons\":"<<pending<<"}\n";p.close();std::filesystem::rename(out/"progress.json.tmp",out/"progress.json");
            std::cout<<std::setprecision(12)<<"Accepted heating clock "<<(start+d.time_s)*1e9<<" ns / "<<end*1e9<<" ns; steps "<<step<<" particles "<<pic.particles().size()<<'\n'<<std::flush;reported=clock();
        };
        if(!resume){snapshot(pic,out,0,start);save();}progress(false);
        while(start+pic.diagnostics().time_s<end-1e-22) {
            const double t=start+pic.diagnostics().time_s;fields.set(t);const double boundary_V=voltage.at(t);
            if(pic.options().upper_potential_V!=boundary_V){if(pic.particles().empty())pic.set_empty_potentials(0,boundary_V);else pic.set_potentials(0,boundary_V);}
            double dt=pic.particles().empty()?1e-12:std::min(1e-14,.05*pic_me/(pic_e*std::max(fields.bound(),1e-20)));
            if(!pic.particles().empty()) {
                std::vector<double> count(o.nx*o.ny,0.);double vmax=0;
                for(const auto& p:pic.particles()){int i=std::clamp(int(p.position_m[0]/25e-6),0,o.nx-1),j=std::clamp(int(p.position_m[1]/25e-6),0,o.ny-1);count[j*o.nx+i]+=p.weight;double u=std::hypot(p.proper_velocity_m_s[0],p.proper_velocity_m_s[1],p.proper_velocity_m_s[2]);vmax=std::max(vmax,u/std::sqrt(1+u*u/(pic_c*pic_c)));}
                double density=*std::max_element(count.begin(),count.end())/(25e-6*25e-6*o.depth_m);
                dt=std::min(dt,.1/std::sqrt(std::max(density,1e-100)*pic_e*pic_e/(eps*pic_me)));
                if(vmax>0)dt=std::min(dt,.2*25e-6/vmax);
            }
            dt=std::min({dt,end-t,fields.times[fields.current+1]-t,voltage.end[voltage.current]-t,std::max(1e-22,next_frame*1e-9-t)});
            if(dt<=1e-23 || !std::isfinite(dt))throw std::runtime_error("Nonprogressing time step");
            auto T=thermal.at(t+.5*dt);double jmax=0;
            for(const auto& v:T)for(double vT:v)jmax=std::max(jmax,pic_richardson_supply_A_m2(vT,4.3,1201732.289489034));
            // Source-induced density also bounds the new explicit plasma step.
            double source_bound=std::cbrt(.01*eps*pic_me*25e-6/(std::max(jmax,1e-100)*pic_e));
            if(dt>source_bound){dt=source_bound;T=thermal.at(t+.5*dt);}
            if(pic.particles().empty()) {
                bool will_emit=false;
                for(int f=0;f<4;++f) {
                    double sum=std::accumulate(reservoir[f].integral.begin(),reservoir[f].integral.end(),0.);
                    for(double vT:T[f])sum+=pic_richardson_supply_A_m2(vT,4.3,1201732.289489034)*dt;
                    if(sum*(f<2?.005:.0002)/T[f].size()*o.depth_m/pic_e>=1.)will_emit=true;
                }
                if(will_emit) {
                    dt=std::min(dt,std::min(1e-14,.05*pic_me/(pic_e*std::max(fields.bound(),1e-20))));
                    T=thermal.at(t+.5*dt);
                }
            }
            bool injected=false;
            for(int f=0;f<4;++f) {
                auto& pool=reservoir[f];pool.span+=dt;
                for(std::size_t i=0;i<T[f].size();++i){double amount=pic_richardson_supply_A_m2(T[f][i],4.3,1201732.289489034)*dt;pool.integral[i]+=amount;pool.weighted_temperature[i]+=amount*T[f][i];}
                double ne=std::accumulate(pool.integral.begin(),pool.integral.end(),0.)*(f<2?.005:.0002)/T[f].size()*o.depth_m/pic_e;
                if(ne<1.)continue;
                auto temp=T[f],J=pool.integral;
                for(std::size_t i=0;i<J.size();++i){if(J[i]>0)temp[i]=pool.weighted_temperature[i]/J[i];J[i]/=pool.span;}
                const GapSurface side=f==0?GapSurface::Back:f==1?GapSurface::Gap:f==2?GapSurface::Left:GapSurface::Right;
                pic.inject_surface_profile(GapElectrode::Upper,side,temp,J,pool.span,32,1000000ULL*step+f*100);
                if(f>=2){std::reverse(temp.begin(),temp.end());std::reverse(J.begin(),J.end());}
                pic.inject_surface_profile(GapElectrode::Lower,side,temp,J,pool.span,32,1000000ULL*step+f*100+50);
                std::fill(pool.integral.begin(),pool.integral.end(),0.);std::fill(pool.weighted_temperature.begin(),pool.weighted_temperature.end(),0.);pool.span=0;injected=true;
            }
            if(!injected && pic.particles().empty())pic.advance_empty(dt);else pic.advance(dt);++step;
            auto d=pic.diagnostics();
            if(!std::isfinite(d.gauss_defect_C) || !std::isfinite(d.particle_charge_defect_C) || pic.particles().size()>2000000)throw std::runtime_error("Nonfinite PIC ledger or two-million-particle resource guard; last checkpoint preserved");
            if(start+d.time_s>=next_frame*1e-9-1e-21){snapshot(pic,out,next_frame,start+d.time_s);++next_frame;save();}
            if(std::chrono::duration<double>(clock()-reported).count()>10)progress(false);
            if(std::chrono::duration<double>(clock()-saved).count()>60)save();
        }
        snapshot(pic,out,next_frame,end);save();progress(true);pic_bench::write((out/"benchmark_profile.json").string());pic_bench::resources((out/"process_resources.json").string());unit_bias.release_field_solver();pic.release_field_solver();pic_parallel_finish();pic_parallel_record(out.string());
        pic_graph::RunRecord::instance().complete();
    }catch(const std::exception& e){if(mfem::Mpi::WorldRank()==0)pic_graph::RunRecord::instance().fail(std::size_t(-1),e.what());std::cerr<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);return 1;}
}
