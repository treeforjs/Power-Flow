#include "native_slice.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <sys/resource.h>
#include <sched.h>
#include <unistd.h>

namespace {
constexpr int seeds=32;
constexpr std::array<int,4> counts={64,32,16,8};
constexpr double window=100e-12, coarse_dt=1e-14, slab_fraction=.04;
struct Event {
    GapElectrode electrode; GapSurface surface;
    double time,span; std::uint64_t seed,step;
    std::vector<double> T,J;
};
struct Capture {
    std::vector<Event> events;
    double time=0; std::uint64_t step=0;
    void inject_surface_profile(GapElectrode electrode,GapSurface surface,
        const std::vector<double>& T,const std::vector<double>& J,double span,int count,std::uint64_t seed) {
        if(count!=32)throw std::runtime_error("Frozen source count changed");
        events.push_back({electrode,surface,time,span,seed,step,T,J});
    }
};
struct Initial {
    GapPicOptions options; std::array<Reservoir,4> pools;
    double time=0; std::uint64_t step=0;
};
Initial read_initial(const std::filesystem::path& dir,const Thermal& thermal) {
    Initial initial;
    std::ifstream meta(dir/"restart.meta");std::string name;int next_frame=0;
    meta>>name>>initial.step>>next_frame;
    for(int f=0;f<4;++f) {
        auto& pool=initial.pools[f];pool.integral.resize(thermal.sizes[f]);pool.weighted_temperature.resize(thermal.sizes[f]);
        meta>>pool.span;for(auto& x:pool.integral)meta>>x;for(auto& x:pool.weighted_temperature)meta>>x;
    }
    double elapsed=0,voltage=0;meta>>elapsed>>voltage;
    if(!meta || std::filesystem::path(name).filename()!=name)throw std::runtime_error("Bad copied reservoir metadata");
    std::ifstream checkpoint(dir/name,std::ios::binary);std::uint64_t magic=0;GapPicDiagnostics ledger;
    checkpoint.read(reinterpret_cast<char*>(&magic),sizeof(magic));
    checkpoint.read(reinterpret_cast<char*>(&initial.options),sizeof(initial.options));
    checkpoint.read(reinterpret_cast<char*>(&ledger),sizeof(ledger));
    if(!checkpoint || magic!=0x5241494b50494332ULL || ledger.time_s!=elapsed || initial.options.upper_potential_V!=voltage ||
       sizeof(initial.options)!=168 || sizeof(ledger)!=200 || sizeof(GapParticle)!=80)
        throw std::runtime_error("Frozen checkpoint header/ABI mismatch");
    initial.time=thermal.times.front()+elapsed;return initial;
}
int face(GapSurface surface) {
    return surface==GapSurface::Back?0:surface==GapSurface::Gap?1:surface==GapSurface::Left?2:3;
}
int group(const Event& e){return face(e.surface)*2+(e.electrode==GapElectrode::Lower?1:0);}
double length(const Event& e,const GapPicOptions& o) {
    const auto f=e.electrode==GapElectrode::Lower?o.lower_foil:o.upper_foil;
    return face(e.surface)<2?f.right_m-f.left_m:f.top_m-f.bottom_m;
}
double pending(const std::array<Reservoir,4>& pools,const Thermal& thermal,const GapPicOptions& o) {
    double result=0;for(int f=0;f<4;++f)result+=2*std::accumulate(pools[f].integral.begin(),pools[f].integral.end(),0.)*(f<2?.005:.0002)/thermal.sizes[f]*o.depth_m/pic_e;
    return result;
}
Capture generate(const Initial& initial,Thermal& thermal,double dtbase,double& remaining) {
    Capture capture;auto reservoir=initial.pools;
    const auto& o=initial.options;double t=initial.time;std::uint64_t step=initial.step;
    const int steps=static_cast<int>(std::llround(window/dtbase));
    for(int k=0;k<steps;++k) {
        double next=initial.time+(k+1)*dtbase,dt=next-t;
        auto T=thermal.at(t+.5*dt);capture.time=next;capture.step=step;
        accumulate_native(capture,reservoir,T,dt,o,thermal,step);
        t=next;++step;
    }
    remaining=pending(reservoir,thermal,o);return capture;
}
struct Model {
    double N=0,total=0,mean_energy=0,var_energy=0,trace_patch=0,miss_mass=0;
    std::vector<double> p;
};
Model model(const Event& e,const GapPicOptions& o,int M) {
    Model m;m.total=std::accumulate(e.J.begin(),e.J.end(),0.);
    m.N=m.total*length(e,o)/e.J.size()*o.depth_m*e.span/pic_e;
    double second=0,p2=0;
    for(std::size_t i=0;i<e.J.size();++i) {
        double p=e.J[i]/m.total,kt=1.380649e-23*e.T[i]/pic_e;
        m.p.push_back(p);m.mean_energy+=p*2*kt;second+=p*6*kt*kt;p2+=p*p;
        m.miss_mass+=p*std::pow(1-p,M);
    }
    m.var_energy=second-m.mean_energy*m.mean_energy;m.trace_patch=m.N*m.N*(1-p2)/M;return m;
}
int patch_index(const Event& e,const GapParticle& p,const GapPicOptions& o) {
    const auto f=e.electrode==GapElectrode::Lower?o.lower_foil:o.upper_foil;
    const double fraction=face(e.surface)<2?(p.position_m[0]-f.left_m)/(f.right_m-f.left_m):(p.position_m[1]-f.bottom_m)/(f.top_m-f.bottom_m);
    const int i=static_cast<int>(fraction*e.J.size());
    if(i<0 || i>=static_cast<int>(e.J.size()))throw std::runtime_error("Launch outside face");
    return i;
}
void check_direction(const Event& e,const GapParticle& p,const GapPicOptions& o) {
    const int axis=face(e.surface)<2?1:0;
    double sign=e.electrode==GapElectrode::Lower?1:-1;
    if(e.surface==GapSurface::Back)sign*=-1;
    if(e.surface==GapSurface::Left)sign=-1;
    if(e.surface==GapSurface::Right)sign=1;
    if(p.proper_velocity_m_s[axis]*sign<=0 || p.position_m[2]<o.short_z_m || p.position_m[2]>o.feed_z_m)
        throw std::runtime_error("Invalid inward velocity/axial launch");
}
struct Stats {
    long double N=0,energy=0,nr_energy=0,slab=0,weight2=0,patch_l2=0,miss_mass=0;
    std::uint64_t samples=0,events=0;
};
struct PatchStats {
    long double expected=0,predicted_var=0,expected_hits=0,variance_hits=0;
    std::uint64_t hits=0;
    std::array<long double,seeds> realized{};
};
struct Expected {long double N=0,energy=0,var_energy=0,var_slab=0,trace_patch=0,miss_mass=0;std::uint64_t events=0;};
void write_events(const Capture& capture,const std::filesystem::path& out,const GapPicOptions& o) {
    std::ofstream bin(out/"events.bin",std::ios::binary),csv(out/"events.csv");
    auto put=[&](const auto& x){bin.write(reinterpret_cast<const char*>(&x),sizeof(x));};
    const std::uint64_t magic=0x5352434556543031ULL,n=capture.events.size();put(magic);put(n);
    csv<<std::setprecision(17)<<"event,absolute_time_s,step,face,electrode,span_s,patches,native_seed,physical_electrons,expected_mean_energy_eV,effective_patches\n";
    int index=0;for(const auto& e:capture.events) {
        std::uint64_t size=e.J.size();int f=face(e.surface),lower=e.electrode==GapElectrode::Lower;
        put(e.time);put(e.span);put(e.seed);put(e.step);put(f);put(lower);put(size);
        bin.write(reinterpret_cast<const char*>(e.T.data()),size*sizeof(double));bin.write(reinterpret_cast<const char*>(e.J.data()),size*sizeof(double));
        auto m=model(e,o,32);double p2=0;for(double p:m.p)p2+=p*p;
        csv<<index++<<','<<e.time<<','<<e.step<<','<<f<<','<<lower<<','<<e.span<<','<<size<<','<<e.seed<<','<<m.N<<','<<m.mean_energy<<','<<1/p2<<'\n';
    }
    if(!bin || !csv)throw std::runtime_error("Event write failed");
}
}
int main(int argc,char** argv) {
 try {
    const auto started=std::chrono::steady_clock::now();
    if(argc!=4)throw std::runtime_error("Usage: validate THERMAL_BINARY COPIED_CHECKPOINT_DIR OUTPUT");
    cpu_set_t mask;CPU_ZERO(&mask);sched_getaffinity(0,sizeof(mask),&mask);
    if(CPU_COUNT(&mask)!=1 || !CPU_ISSET(6,&mask) || getpriority(PRIO_PROCESS,0)!=15)
        throw std::runtime_error("Must run on CPU6 at nice15");
    const std::filesystem::path out=argv[3];Thermal thermal(argv[1]);auto initial=read_initial(argv[2],thermal);
    double remaining=0,fine_remaining=0;auto capture=generate(initial,thermal,coarse_dt,remaining);
    auto fine=generate(initial,thermal,coarse_dt/2,fine_remaining);
    double coarse_supply=remaining,fine_supply=fine_remaining;
    for(const auto& e:capture.events)coarse_supply+=model(e,initial.options,32).N;
    for(const auto& e:fine.events)fine_supply+=model(e,initial.options,32).N;
    write_events(capture,out,initial.options);
    std::ofstream rows(out/"seed-metrics.csv"),expect(out/"expected-metrics.csv"),patches(out/"patch-metrics.csv");
    rows<<std::setprecision(17)<<"M,seed,group,events,samples,physical_electrons,source_energy_eV,nonrel_source_energy_eV,slab_physical_electrons,weight_effective_samples,patch_squared_error,missed_probability_mass\n";
    expect<<std::setprecision(17)<<"M,group,events,physical_electrons,source_energy_eV,variance_source_energy_eV2,slab_physical_electrons,variance_slab_electrons2,expected_patch_squared_error,expected_missed_probability_mass\n";
    patches<<std::setprecision(17)<<"M,group,patch,expected_charge_electrons,predicted_variance,ensemble_mean,ensemble_variance,z_mean,nonzero_seeds,ensemble_expected_hits,ensemble_variance_hits,ensemble_hits\n";
    double max_event_charge_error=0,max_prefix_weight_error=0;std::uint64_t prefix_checks=0,zero_patch_draws=0;
    for(int M:counts) {
        std::array<Expected,9> expected{};std::array<std::vector<PatchStats>,8> patch_stats;
        for(int g=0;g<8;++g)patch_stats[g].resize(thermal.sizes[g/2]);
        for(const auto& e:capture.events) {
            const auto m=model(e,initial.options,M);const int g=group(e);
            for(int h:{g,8}) {auto& x=expected[h];x.N+=m.N;x.energy+=m.N*m.mean_energy;x.var_energy+=m.N*m.N*m.var_energy/M;
                x.var_slab+=m.N*m.N*slab_fraction*(1-slab_fraction)/M;x.trace_patch+=m.trace_patch;x.miss_mass+=m.miss_mass;++x.events;}
            for(std::size_t i=0;i<m.p.size();++i){auto& p=patch_stats[g][i];p.expected+=m.N*m.p[i];p.predicted_var+=m.N*m.N*m.p[i]*(1-m.p[i])/M;
                p.expected_hits+=seeds*M*m.p[i];p.variance_hits+=seeds*M*m.p[i]*(1-m.p[i]);}
        }
        for(int seed_index=0;seed_index<seeds;++seed_index) {
            std::array<Stats,9> stats{};
            for(const auto& e:capture.events) {
                if(std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()>50)throw std::runtime_error("50 second source-only wall cap reached; incomplete");
                const auto m=model(e,initial.options,M);const int g=group(e);
                const std::uint64_t seed=e.seed+std::uint64_t(seed_index)*0x9e3779b97f4a7c15ULL;
                GapPic pic(initial.options);pic.inject_surface_profile(e.electrode,e.surface,e.T,e.J,e.span,M,seed);
                GapPic ref(initial.options);ref.inject_surface_profile(e.electrode,e.surface,e.T,e.J,e.span,64,seed);
                std::vector<int> hits(e.J.size());Stats event;
                if(pic.particles_.size()!=static_cast<std::size_t>(M))throw std::runtime_error("Bad emitted macro count");
                for(int n=0;n<M;++n) {
                    const auto& p=pic.particles_[n];const auto& r=ref.particles_[n];
                    if(p.position_m!=r.position_m || p.proper_velocity_m_s!=r.proper_velocity_m_s || p.charge_C!=r.charge_C || p.mass_kg!=r.mass_kg)
                        throw std::runtime_error("Sample prefix differs across counts");
                    max_prefix_weight_error=std::max(max_prefix_weight_error,std::abs(p.weight*M-r.weight*64)/(r.weight*64));++prefix_checks;
                    int i=patch_index(e,p,initial.options);++hits[i];zero_patch_draws+=e.J[i]==0;check_direction(e,p,initial.options);
                    double u2=norm2(p.proper_velocity_m_s),v2=u2/(1+u2/(pic_c*pic_c));
                    event.N+=p.weight;event.weight2+=p.weight*p.weight;event.energy+=p.weight*pic_kinetic_energy_J(p)/pic_e;
                    event.nr_energy+=p.weight*.5*pic_me*v2/pic_e;
                    if(std::abs(p.position_m[2])<.0005)event.slab+=p.weight;
                    patch_stats[g][i].realized[seed_index]+=p.weight;
                    ++patch_stats[g][i].hits;
                }
                max_event_charge_error=std::max(max_event_charge_error,double(std::abs(event.N-m.N)/m.N));
                for(std::size_t i=0;i<hits.size();++i) {
                    double delta=m.N*(double(hits[i])/M-m.p[i]);event.patch_l2+=delta*delta;
                    if(hits[i]==0)event.miss_mass+=m.p[i];
                }
                for(int h:{g,8}) {auto& s=stats[h];s.N+=event.N;s.energy+=event.energy;s.nr_energy+=event.nr_energy;s.slab+=event.slab;
                    s.weight2+=event.weight2;s.patch_l2+=event.patch_l2;s.miss_mass+=event.miss_mass;s.samples+=M;++s.events;}
            }
            for(int g=0;g<9;++g) {const auto& s=stats[g];rows<<M<<','<<seed_index<<','<<g<<','<<s.events<<','<<s.samples<<','<<s.N<<','<<s.energy<<','<<s.nr_energy<<','<<s.slab<<','<<(s.weight2>0?s.N*s.N/s.weight2:0)<<','<<s.patch_l2<<','<<(s.events?s.miss_mass/s.events:0)<<'\n';}
        }
        for(int g=0;g<9;++g){const auto& x=expected[g];expect<<M<<','<<g<<','<<x.events<<','<<x.N<<','<<x.energy<<','<<x.var_energy<<','<<slab_fraction*x.N<<','<<x.var_slab<<','<<x.trace_patch<<','<<(x.events?x.miss_mass/x.events:0)<<'\n';}
        for(int g=0;g<8;++g)for(std::size_t i=0;i<patch_stats[g].size();++i){const auto& x=patch_stats[g][i];long double mean=0,var=0;int nonzero=0;
            for(auto v:x.realized){mean+=v/seeds;nonzero+=v>0;}for(auto v:x.realized)var+=(v-mean)*(v-mean)/(seeds-1);
            double z=x.predicted_var>0?double((mean-x.expected)/std::sqrt(x.predicted_var/seeds)):0;
            patches<<M<<','<<g<<','<<i<<','<<x.expected<<','<<x.predicted_var<<','<<mean<<','<<var<<','<<z<<','<<nonzero<<','<<x.expected_hits<<','<<x.variance_hits<<','<<x.hits<<'\n';}
        std::cout<<"Completed M="<<M<<" with "<<seeds<<" seeds\n"<<std::flush;
    }
    // Explicit zeros, orientations, and the no-emission path supplement real thermal events.
    std::uint64_t synthetic_draws=0;
    for(int M:counts)for(int s=0;s<seeds;++s)for(auto side:{GapSurface::Back,GapSurface::Gap,GapSurface::Left,GapSurface::Right})for(auto electrode:{GapElectrode::Upper,GapElectrode::Lower}) {
        Event e{electrode,side,0,1e-12,std::uint64_t(s),0,{1000,2000,3000,2500},{0,1,3,6}};
        GapPic pic(initial.options);pic.inject_surface_profile(e.electrode,e.surface,e.T,e.J,e.span,M,e.seed);
        for(const auto& p:pic.particles_) {zero_patch_draws+=patch_index(e,p,initial.options)==0;check_direction(e,p,initial.options);++synthetic_draws;}
        pic.particles_.clear();pic.inject_surface_profile(e.electrode,e.surface,e.T,std::vector<double>(4,0),e.span,M,e.seed);
        if(!pic.particles_.empty())throw std::runtime_error("Zero source emitted particles");
    }
    if(max_event_charge_error>1e-12 || max_prefix_weight_error>1e-12 || zero_patch_draws)throw std::runtime_error("Source validation guard failed");
    std::ofstream result(out/"run-summary.json");result<<std::setprecision(17)
        <<"{\n\"complete\":true,\n\"source_only\":true,\n\"coupled_convergence\":false,\n\"seeds\":"<<seeds
        <<",\n\"cpu\":6,\n\"nice\":15,\n\"threads\":1,\n\"absolute_start_s\":"<<initial.time<<",\n\"duration_s\":"<<window
        <<",\n\"source_quadrature_dt_s\":"<<coarse_dt<<",\n\"events\":"<<capture.events.size()<<",\n\"fine_quadrature_events\":"<<fine.events.size()
        <<",\n\"initial_pending_electrons\":"<<pending(initial.pools,thermal,initial.options)<<",\n\"terminal_pending_electrons\":"<<remaining
        <<",\n\"fine_terminal_pending_electrons\":"<<fine_remaining
        <<",\n\"coarse_emitted_electrons\":"<<coarse_supply-remaining<<",\n\"fine_emitted_electrons\":"<<fine_supply-fine_remaining
        <<",\n\"coarse_emitted_plus_pending_electrons\":"<<coarse_supply<<",\n\"fine_emitted_plus_pending_electrons\":"<<fine_supply
        <<",\n\"quadrature_relative_difference\":"<<std::abs(coarse_supply-fine_supply)/fine_supply
        <<",\n\"max_event_charge_relative_error\":"<<max_event_charge_error<<",\n\"prefix_checks\":"<<prefix_checks
        <<",\n\"max_prefix_scaled_weight_relative_error\":"<<max_prefix_weight_error<<",\n\"zero_current_patch_draws\":"<<zero_patch_draws
        <<",\n\"synthetic_orientation_draws\":"<<synthetic_draws<<",\n\"wall_s\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<"\n}\n";
    if(!rows || !expect || !patches || !result)throw std::runtime_error("Metrics write failed");
    return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 2;}
}
