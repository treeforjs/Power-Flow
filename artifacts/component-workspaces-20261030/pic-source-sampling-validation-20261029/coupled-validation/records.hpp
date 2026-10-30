#pragma once
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include <sys/resource.h>
#include <sched.h>
namespace cv {
using namespace raikou::mfem_backend;
inline constexpr bool gpu_execution_enabled=false;
inline constexpr bool coupled_execution_enabled=false;
inline constexpr std::array<int,3> counts{64,32,16};
inline void require(bool ok,const std::string& s){if(!ok)throw std::runtime_error(s);}
inline bool count_ok(int n){return std::find(counts.begin(),counts.end(),n)!=counts.end();}
inline void resource_check(){
    cpu_set_t mask;CPU_ZERO(&mask);require(sched_getaffinity(0,sizeof(mask),&mask)==0,"affinity unavailable");
    require(CPU_COUNT(&mask)==1&&CPU_ISSET(6,&mask),"CPU6 only");
    require(getpriority(PRIO_PROCESS,0)==15,"nice15 required");
    std::ifstream status("/proc/self/status");std::string line;bool found=false;
    while(std::getline(status,line))if(line.rfind("Threads:",0)==0){std::istringstream row(line.substr(8));int threads=0;row>>threads;require(threads==1,"one preflight thread required");found=true;}
    require(found,"thread-count evidence unavailable");
}
template<class T> void read(std::istream& f,T& x){f.read(reinterpret_cast<char*>(&x),sizeof(x));require(bool(f),"truncated record");}
template<class T> void write(std::ostream& f,const T& x){f.write(reinterpret_cast<const char*>(&x),sizeof(x));require(bool(f),"record write failed");}
inline void vector_read(std::istream& f,std::vector<double>& x,std::uint64_t n){
    require(n>0&&n<=4096,"invalid patch count");x.resize(n);f.read(reinterpret_cast<char*>(x.data()),n*sizeof(double));require(bool(f),"truncated patches");
}
inline void vector_write(std::ostream& f,const std::vector<double>& x){f.write(reinterpret_cast<const char*>(x.data()),x.size()*sizeof(double));require(bool(f),"patch write failed");}
struct Pool {std::vector<double> integral,weighted_temperature;double span=0;};
struct Input {
    GapPicOptions options;GapPicDiagnostics ledger;std::uint64_t next_id=0,step=0;int next_frame=0;
    std::array<Pool,4> pools;std::vector<GapParticle> particles;std::string checkpoint;
};
inline Input read_input(const std::filesystem::path& directory){
    static_assert(sizeof(GapPicOptions)==168&&sizeof(GapPicDiagnostics)==200&&sizeof(GapParticle)==80);
    Input r;std::ifstream meta(directory/"restart.meta");meta>>r.checkpoint>>r.step>>r.next_frame;
    require(bool(meta)&&std::filesystem::path(r.checkpoint).filename()==r.checkpoint,"unsafe checkpoint name");
    std::ifstream f(directory/r.checkpoint,std::ios::binary);std::uint64_t magic=0,n=0;
    read(f,magic);read(f,r.options);read(f,r.ledger);read(f,r.next_id);read(f,n);
    require(magic==0x5241494b50494332ULL&&n<=2000000&&r.next_id>0,"checkpoint ABI/header invalid");
    require(std::filesystem::file_size(directory/r.checkpoint)==392+80*n,"checkpoint byte length invalid");
    const auto& o=r.options;
    std::array<double,25> ledger_values{};std::memcpy(ledger_values.data(),&r.ledger,sizeof(r.ledger));
    for(double x:ledger_values)require(std::isfinite(x),"nonfinite checkpoint ledger");
    require(o.nx==360&&o.ny==360&&o.width_m==.009&&o.gap_m==.009&&o.depth_m==.025,
            "frozen grid/depth changed");
    require(o.finite_foils&&o.finite_axial_extent&&o.distributed_fields&&o.short_z_m==-.0125&&o.feed_z_m==.0125,"frozen geometry/backend flags changed");
    std::set<std::uint64_t> ids;r.particles.resize(n);
    for(auto& p:r.particles){
        read(f,p);require(p.id<r.next_id&&ids.insert(p.id).second,"invalid/duplicate particle id");
        require(p.charge_C==-pic_e&&p.mass_kg==pic_me&&std::isfinite(p.weight)&&p.weight>0,"invalid species/weight");
        for(double x:p.position_m)require(std::isfinite(x),"nonfinite position");
        for(double u:p.proper_velocity_m_s)require(std::isfinite(u),"nonfinite velocity");
        const auto& u=p.proper_velocity_m_s;require(std::isfinite(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]),"velocity norm overflow");
        require(p.position_m[0]>=0&&p.position_m[0]<=o.width_m&&p.position_m[1]>=0&&p.position_m[1]<=o.gap_m&&p.position_m[2]>=o.short_z_m&&p.position_m[2]<=o.feed_z_m,"particle outside domain");
        for(const auto& foil:{o.lower_foil,o.upper_foil})require(!(p.position_m[0]>foil.left_m&&p.position_m[0]<foil.right_m&&p.position_m[1]>foil.bottom_m&&p.position_m[1]<foil.top_m),"particle inside foil");
    }
    require(f.peek()==std::char_traits<char>::eof(),"trailing checkpoint bytes");
    const std::array<std::size_t,4> sizes{2048,2048,256,256};
    for(int k=0;k<4;++k){auto& p=r.pools[k];p.integral.resize(sizes[k]);p.weighted_temperature.resize(sizes[k]);meta>>p.span;
        require(std::isfinite(p.span)&&p.span>=0,"invalid reservoir span");
        for(auto& x:p.integral){meta>>x;require(std::isfinite(x)&&x>=0,"invalid reservoir integral");}
        for(auto& x:p.weighted_temperature){meta>>x;require(std::isfinite(x)&&x>=0,"invalid weighted temperature");}
        for(std::size_t i=0;i<sizes[k];++i)require(p.integral[i]==0?p.weighted_temperature[i]==0:(p.span>0&&p.weighted_temperature[i]/p.integral[i]>0&&p.weighted_temperature[i]/p.integral[i]<=1e6),"reservoir temperature/span inconsistent");
    }
    double t=0,v=0;meta>>t>>v;require(bool(meta)&&t==r.ledger.time_s&&v==o.upper_potential_V,"restart clock/voltage mismatch");
    std::string trailing;require(!(meta>>trailing),"trailing restart tokens");
    return r;
}
inline double pending(const Input& r){double n=0;for(int f=0;f<4;++f)n+=2*std::accumulate(r.pools[f].integral.begin(),r.pools[f].integral.end(),0.)*(f<2?.005:.0002)/r.pools[f].integral.size()*r.options.depth_m/pic_e;return n;}
struct Event {
    GapElectrode electrode=GapElectrode::Upper;GapSurface surface=GapSurface::Gap;
    double injection_time_s=0,span_s=0;std::uint64_t native_seed=0,step=0;
    std::vector<double> T,J;
};
inline void valid(const Event& e){
    require(std::isfinite(e.injection_time_s)&&std::isfinite(e.span_s)&&e.span_s>0,"invalid source time/span");
    require(e.electrode==GapElectrode::Lower||e.electrode==GapElectrode::Upper,"unknown electrode");
    require(int(e.surface)>=0&&int(e.surface)<=3,"unknown surface");
    require(!e.T.empty()&&e.T.size()==e.J.size()&&e.T.size()<=4096,"patch shape mismatch");
    double sum=0;for(std::size_t i=0;i<e.T.size();++i){require(std::isfinite(e.T[i])&&e.T[i]>0&&e.T[i]<=1e6&&std::isfinite(e.J[i])&&e.J[i]>=0,"invalid T/J patch");sum+=e.J[i];}
    require(std::isfinite(sum)&&sum>0,"empty/nonfinite source");
}
inline double physical_electrons(const Event& e,const GapPicOptions& o){
    const auto f=e.electrode==GapElectrode::Lower?o.lower_foil:o.upper_foil;
    const double L=(e.surface==GapSurface::Gap||e.surface==GapSurface::Back)?f.right_m-f.left_m:f.top_m-f.bottom_m;
    return std::accumulate(e.J.begin(),e.J.end(),0.)*L/e.J.size()*o.depth_m*e.span_s/pic_e;
}
inline void event_write(std::ostream& out,const Event& e){valid(e);write(out,e.injection_time_s);write(out,e.span_s);write(out,e.native_seed);write(out,e.step);
    auto el=std::uint32_t(e.electrode),side=std::uint32_t(e.surface);write(out,el);write(out,side);auto n=std::uint64_t(e.T.size());write(out,n);vector_write(out,e.T);vector_write(out,e.J);}
inline Event event_read(std::istream& in){Event e;std::uint32_t el=0,side=0;std::uint64_t n=0;
    read(in,e.injection_time_s);read(in,e.span_s);read(in,e.native_seed);read(in,e.step);read(in,el);read(in,side);e.electrode=GapElectrode(el);e.surface=GapSurface(side);read(in,n);vector_read(in,e.T,n);vector_read(in,e.J,n);valid(e);return e;}
inline bool same_event(const Event& a,const Event& b){return a.electrode==b.electrode&&a.surface==b.surface&&a.injection_time_s==b.injection_time_s&&a.span_s==b.span_s&&a.native_seed==b.native_seed&&a.step==b.step&&a.T==b.T&&a.J==b.J;}
struct Step {std::uint64_t step=0;double time_s=0,dt_s=0;std::vector<Event> events;};
struct Trace {
    bool coupled_verified=false;int seed_index=0;double initial_time_s=0,end_time_s=0;
    std::string input_identity;
    std::vector<Step> steps;
};
inline void valid(const Trace& t){
    require(t.seed_index>=0&&t.seed_index<=11&&std::isfinite(t.initial_time_s)&&t.end_time_s>t.initial_time_s,"invalid trace header");
    require(t.input_identity.size()==64&&t.input_identity.find_first_not_of("0123456789abcdef")==std::string::npos,"trace requires bound input identity");
    require(!t.steps.empty()&&t.steps.size()<=20000,"invalid trace step count");
    double clock=t.initial_time_s;std::uint64_t step=t.steps.front().step;
    for(const auto& s:t.steps){const double tol=32*std::numeric_limits<double>::epsilon()*std::max(std::abs(clock),1e-12);
        require(s.step==step++&&std::abs(s.time_s-clock)<=tol&&std::isfinite(s.dt_s)&&s.dt_s>0&&s.dt_s<=1e-12,"trace clock/step invalid");
        require(s.events.size()<=8,"too many face events");for(const auto& e:s.events){valid(e);require(e.step==s.step&&e.injection_time_s==s.time_s,"source injection must use native pre-advance time");}
        clock=s.time_s+s.dt_s;require(clock<=t.end_time_s+1e-22,"trace overshoots end");}
    require(std::abs(clock-t.end_time_s)<=1e-22,"incomplete trace");
}
inline void trace_write(const std::filesystem::path& path,const Trace& t){valid(t);std::ofstream f(path,std::ios::binary);
    const std::uint64_t magic=0x4356545241434531ULL,n=t.steps.size();write(f,magic);std::uint32_t v=t.coupled_verified,seed=t.seed_index;write(f,v);write(f,seed);write(f,t.initial_time_s);write(f,t.end_time_s);write(f,n);f.write(t.input_identity.data(),64);require(bool(f),"trace identity write failed");
    for(const auto& s:t.steps){write(f,s.step);write(f,s.time_s);write(f,s.dt_s);auto count=std::uint64_t(s.events.size());write(f,count);for(const auto& e:s.events)event_write(f,e);}}
inline Trace trace_read(const std::filesystem::path& path){std::ifstream f(path,std::ios::binary);Trace t;std::uint64_t magic=0,n=0;std::uint32_t v=0,seed=0;
    read(f,magic);read(f,v);read(f,seed);read(f,t.initial_time_s);read(f,t.end_time_s);read(f,n);require(magic==0x4356545241434531ULL&&v<=1&&n>0&&n<=20000,"trace header invalid");t.coupled_verified=v;t.seed_index=seed;
    t.input_identity.resize(64);f.read(t.input_identity.data(),64);require(bool(f),"truncated trace identity");
    t.steps.resize(n);for(auto& s:t.steps){read(f,s.step);read(f,s.time_s);read(f,s.dt_s);std::uint64_t count=0;read(f,count);require(count<=8,"bad trace event count");for(std::uint64_t i=0;i<count;++i)s.events.push_back(event_read(f));}
    require(f.peek()==std::char_traits<char>::eof(),"trailing trace data");valid(t);return t;}
struct Caps {
    std::uint64_t steps=20000;double wall_s=900,duration_s=100e-12;
    std::uintmax_t bytes=256ULL*1024*1024,rss_KiB=4ULL*1024*1024;
};
inline Caps caps(bool pilot){return pilot?Caps{2000,90,10e-12}:Caps{};}
inline std::string cap_reason(const Caps& c,std::uint64_t steps,double wall,double elapsed,std::uintmax_t bytes,std::uintmax_t rss){
    if(wall>=c.wall_s)return "wall_cap";
    if(bytes>=c.bytes)return "output_cap";
    if(rss>=c.rss_KiB)return "memory_cap";
    if(elapsed>=c.duration_s-1e-22&&steps<=c.steps)return "complete";
    if(steps>=c.steps)return "step_cap";
    return "";}
inline void require_owned_output(const std::filesystem::path& output){
    const auto root=std::filesystem::canonical("/home/treefor/sim-workspace/powerflow/outputs/pic-profile/source-sampling-validation/coupled-preflight");
    const auto full=std::filesystem::weakly_canonical(output);auto rel=full.lexically_relative(root);
    require(!rel.empty()&&!rel.is_absolute(),"output outside owned directory");
    for(const auto& part:rel)require(part!="..","output escapes owned directory");
}
inline void require_execution(){throw std::runtime_error("Coupled and GPU execution are disabled in this build; no enable flag or environment override exists");}
}
