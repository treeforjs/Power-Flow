#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <optional>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

// Local binary neutral collision scheduling only. No particle mutation, product
// kinematics, field evolution, water data, or claim of physical event admission.
namespace raikou::pic_collisions {
using Vec = std::array<double,3>;
constexpr double c = 299792458.;
constexpr double qe = 1.602176634e-19;
struct Error : std::runtime_error { using std::runtime_error::runtime_error; };
inline void require(bool condition,const char* text){if(!condition)throw Error(text);}
inline bool digest(const std::string& s){return s.size()==64&&s.find_first_not_of("0123456789abcdef")==std::string::npos;}
inline void positive(double x,const char* message){require(std::isfinite(x)&&x>0,message);}
struct Packet {
    std::string id,state;
    std::uint64_t cell=0;
    double weight=0,mass_kg=0;
    Vec proper_velocity_m_s{};
};
struct Cell { std::uint64_t id=0; double volume_m3=0; };
enum class Observable { ExclusiveEvent, InclusiveProductYield, AggregateInelastic };
struct Channel {
    std::string id,projectile_state,target_state,data_sha256,kernel_sha256;
    Observable observable=Observable::AggregateInelastic;
    // Query itself owns threshold and all disjoint support checks. An absent
    // kernel identity cannot be turned into a physical reaction by this clock.
    std::function<double(double)> sigma_m2;
};
struct Clock {
    std::uint64_t step=0,revision=0;
    double begin_s=0,end_s=0;
    std::string state_sha256,mesh_sha256,mechanism_sha256;
};
struct Kinematics { double target_frame_energy_eV=0,moller_flux_m_s=0; };
inline Kinematics relative(const Packet& p,const Packet& n) {
    positive(p.mass_kg,"projectile mass");positive(n.mass_kg,"target mass");
    std::array<long double,3> a{},b{},vp{},vn{},boost{};
    long double aa=0,bb=0,ab=0;
    for(int d=0;d<3;++d){
        require(std::isfinite(p.proper_velocity_m_s[d])&&std::isfinite(n.proper_velocity_m_s[d]),"nonfinite proper velocity");
        a[d]=p.proper_velocity_m_s[d]/static_cast<long double>(c);
        b[d]=n.proper_velocity_m_s[d]/static_cast<long double>(c);
        aa+=a[d]*a[d];bb+=b[d]*b[d];ab+=a[d]*b[d];
    }
    const long double gp=std::sqrt(1+aa),gn=std::sqrt(1+bb);
    long double boosted2=0,difference2=0;
    for(int d=0;d<3;++d){
        boost[d]=a[d]+(ab/(gn+1)-gp)*b[d];boosted2+=boost[d]*boost[d];
        vp[d]=a[d]/gp;vn[d]=b[d]/gn;
        difference2+=(vp[d]-vn[d])*(vp[d]-vn[d]);
    }
    long double cross2=0;
    for(int d=0;d<3;++d){int j=(d+1)%3,k=(d+2)%3;long double v=vp[j]*vn[k]-vp[k]*vn[j];cross2+=v*v;}
    const long double flux2=difference2-cross2;
    require(flux2>=0,"invalid invariant relative flux");
    const long double energy=p.mass_kg*static_cast<long double>(c)*c*boosted2/(std::sqrt(1+boosted2)+1)/qe;
    Kinematics out{static_cast<double>(energy),static_cast<double>(c*std::sqrt(flux2))};
    require(std::isfinite(out.target_frame_energy_eV)&&std::isfinite(out.moller_flux_m_s),"relative kinematics overflow");
    return out;
}
struct Candidate {
    std::string projectile_id,target_id,channel_id,data_sha256,kernel_sha256;
    std::uint64_t cell=0;
    double event_weight=0,projectile_remainder=0,target_remainder=0;
    double rate_s=0,energy_eV=0,flux_m_s=0,sigma_m2=0,volume_m3=0;
};
inline std::vector<Candidate> candidates(const std::vector<Packet>& projectiles,
    const std::vector<Packet>& targets,const std::vector<Cell>& cells,
    const std::vector<Channel>& channels,std::size_t pair_limit=1000000) {
    std::set<std::string> ids,channel_ids;std::set<std::uint64_t> cell_ids;
    for(const auto& cell:cells){require(cell_ids.insert(cell.id).second,"duplicate cell");positive(cell.volume_m3,"physical cell volume");}
    for(const auto* list:{&projectiles,&targets})for(const auto& p:*list){
        require(!p.id.empty()&&!p.state.empty()&&ids.insert(p.id).second,"missing or duplicate particle identity");
        positive(p.weight,"represented weight");positive(p.mass_kg,"particle mass");
        require(cell_ids.count(p.cell),"particle cell absent");
        for(double u:p.proper_velocity_m_s)require(std::isfinite(u),"nonfinite particle momentum");
    }
    for(const auto& ch:channels){
        require(!ch.id.empty()&&channel_ids.insert(ch.id).second,"duplicate channel");
        require(!ch.projectile_state.empty()&&!ch.target_state.empty(),"missing channel reactant states");
        require(ch.observable==Observable::ExclusiveEvent,"nonexclusive observable is not an event cross section");
        require(digest(ch.data_sha256)&&digest(ch.kernel_sha256)&&bool(ch.sigma_m2),"unbound data/product kernel");
    }
    std::vector<Candidate> out;std::size_t pairs=0;
    for(const auto& p:projectiles)for(const auto& n:targets){
        if(p.cell!=n.cell)continue;
        require(++pairs<=pair_limit,"collision pair work limit");
        const auto kin=relative(p,n);
        const auto cell=std::find_if(cells.begin(),cells.end(),[&](const Cell& a){return a.id==p.cell;});
        for(const auto& ch:channels){
            if(ch.projectile_state!=p.state||ch.target_state!=n.state)continue;
            // Query even at zero flux: missing low-energy support is not hidden.
            const double sigma=ch.sigma_m2(kin.target_frame_energy_eV);
            require(std::isfinite(sigma)&&sigma>=0,"negative or nonfinite cross section");
            const double w=std::min(p.weight,n.weight),big=std::max(p.weight,n.weight);
            const long double rate=static_cast<long double>(big)*sigma*kin.moller_flux_m_s/cell->volume_m3;
            require(std::isfinite(rate)&&rate<=std::numeric_limits<double>::max(),"collision propensity overflow");
            const double rounded=static_cast<double>(rate);
            require(rate==0||rounded>0,"collision propensity underflow");
            if(rounded>0)out.push_back({p.id,n.id,ch.id,ch.data_sha256,ch.kernel_sha256,p.cell,w,p.weight-w,n.weight-w,rounded,kin.target_frame_energy_eV,kin.moller_flux_m_s,sigma,cell->volume_m3});
        }
    }
    std::sort(out.begin(),out.end(),[](const Candidate&a,const Candidate&b){
        return std::tie(a.cell,a.projectile_id,a.target_id,a.channel_id)<std::tie(b.cell,b.projectile_id,b.target_id,b.channel_id);});
    return out;
}
class Scheduler;
class Proposal {
    friend class Scheduler;
    const Scheduler* owner_=nullptr;std::uint64_t generation_=0;
    std::mt19937_64 random_after_;Clock clock_;std::optional<Candidate> event_;
    double event_time_s_=0,total_rate_s_=0;
public:
    const std::optional<Candidate>& event()const{return event_;}
    const Clock& clock()const{return clock_;}
    double event_time_s()const{return event_time_s_;}
    double total_rate_s()const{return total_rate_s_;}
};
class Scheduler {
    std::mt19937_64 random_;std::uint64_t generation_=0;
    static double uniform_open(std::mt19937_64& rng){
        // 52-bit midpoint rule ensures strict (0,1) in binary64.
        return (static_cast<double>(rng()>>12)+.5)*0x1p-52;
    }
public:
    explicit Scheduler(std::uint64_t seed):random_(seed){}
    Scheduler(const Scheduler&)=delete;Scheduler& operator=(const Scheduler&)=delete;
    Proposal prepare(const Clock& clock,const std::vector<Candidate>& input)const {
        require(std::isfinite(clock.begin_s)&&clock.begin_s>=0&&std::isfinite(clock.end_s)&&clock.end_s>clock.begin_s,"invalid accepted interval");
        require(digest(clock.state_sha256)&&digest(clock.mesh_sha256)&&digest(clock.mechanism_sha256),"unbound input revision");
        Proposal result;result.owner_=this;result.generation_=generation_;result.random_after_=random_;result.clock_=clock;
        long double total=0;
        std::set<std::tuple<std::uint64_t,std::string,std::string,std::string>> keys;
        auto canonical=input;
        std::sort(canonical.begin(),canonical.end(),[](const Candidate&a,const Candidate&b){
            return std::tie(a.cell,a.projectile_id,a.target_id,a.channel_id)<std::tie(b.cell,b.projectile_id,b.target_id,b.channel_id);});
        for(const auto& a:canonical){
            positive(a.rate_s,"candidate rate");positive(a.event_weight,"event weight");
            require(!a.projectile_id.empty()&&!a.target_id.empty()&&a.projectile_id!=a.target_id&&!a.channel_id.empty(),"candidate identity");
            require(digest(a.data_sha256)&&digest(a.kernel_sha256),"candidate data/product kernel");
            positive(a.flux_m_s,"candidate relative flux");positive(a.sigma_m2,"candidate cross section");positive(a.volume_m3,"candidate physical volume");
            require(std::isfinite(a.energy_eV)&&a.energy_eV>=0,"candidate energy");
            require(std::isfinite(a.projectile_remainder)&&a.projectile_remainder>=0&&std::isfinite(a.target_remainder)&&a.target_remainder>=0,"invalid unreacted weights");
            require(a.projectile_remainder==0||a.target_remainder==0,"binary event must consume the smaller packet");
            const long double expected=(static_cast<long double>(a.event_weight)+std::max(a.projectile_remainder,a.target_remainder))*a.sigma_m2*a.flux_m_s/a.volume_m3;
            require(std::abs(expected-a.rate_s)<=8*std::numeric_limits<double>::epsilon()*expected,"inconsistent candidate propensity");
            require(keys.emplace(a.cell,a.projectile_id,a.target_id,a.channel_id).second,"duplicate candidate");
            total+=a.rate_s;
        }
        require(total<=std::numeric_limits<double>::max(),"total rate overflow");
        result.total_rate_s_=static_cast<double>(total);
        if(total==0)return result;
        const long double wait=-std::log(uniform_open(result.random_after_))/total;
        if(wait>=static_cast<long double>(clock.end_s)-clock.begin_s)return result;
        result.event_time_s_=static_cast<double>(clock.begin_s+wait);
        require(result.event_time_s_>clock.begin_s&&result.event_time_s_<clock.end_s,"collision clock resolution exhausted");
        const long double choice=uniform_open(result.random_after_)*total;
        long double cumulative=0;
        for(const auto& a:canonical){cumulative+=a.rate_s;if(choice<cumulative){result.event_=a;break;}}
        require(result.event_.has_value(),"candidate sum lost selection");
        return result;
    }
    bool can_publish(const Proposal& p)const noexcept{return p.owner_==this&&p.generation_==generation_&&generation_!=std::numeric_limits<std::uint64_t>::max();}
    // RNG publication only. Parent must preflight all physical owners and commit
    // the selected event first within one no-fail common state publication.
    bool try_publish(const Proposal& p)noexcept{
        if(!can_publish(p))return false;
        random_=p.random_after_;++generation_;return true;
    }
    std::string checkpoint()const{
        std::ostringstream o;o<<"weighted-binary-clock-v1 "<<generation_<<' '<<random_;return o.str();
    }
    static Scheduler restore(const std::string& text){
        std::istringstream in(text);std::string schema;std::uint64_t generation;std::mt19937_64 rng;
        require(bool(in>>schema>>generation>>rng)&&schema=="weighted-binary-clock-v1","invalid scheduler checkpoint");
        in>>std::ws;require(in.eof(),"trailing scheduler checkpoint data");return Scheduler(rng,generation);
    }
private:
    Scheduler(std::mt19937_64 r,std::uint64_t g):random_(r),generation_(g){}
};
} // namespace raikou::pic_collisions
