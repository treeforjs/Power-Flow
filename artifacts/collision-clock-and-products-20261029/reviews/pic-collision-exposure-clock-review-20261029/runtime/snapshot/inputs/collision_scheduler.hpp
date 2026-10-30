#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <map>
#include <optional>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
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
    double frozen_physical_stage_s=0;
    double collision_elapsed_s=0;
};
struct Kinematics { double target_frame_energy_eV=0,moller_flux_m_s=0; };
inline Kinematics relative(const Packet& p,const Packet& n) {
    positive(p.mass_kg,"projectile mass");positive(n.mass_kg,"target mass");
    std::array<long double,3> a{},b{},raw_a{},raw_b{},raw_delta{};
    long double aa=0,bb=0,dot=0,radial_numerator=0;
    for(int d=0;d<3;++d){
        require(std::isfinite(p.proper_velocity_m_s[d])&&std::isfinite(n.proper_velocity_m_s[d]),"nonfinite proper velocity");
        raw_a[d]=p.proper_velocity_m_s[d];raw_b[d]=n.proper_velocity_m_s[d];
        raw_delta[d]=raw_a[d]-raw_b[d];
        a[d]=raw_a[d]/static_cast<long double>(c);
        b[d]=raw_b[d]/static_cast<long double>(c);
        aa+=a[d]*a[d];bb+=b[d]*b[d];dot+=a[d]*b[d];
        radial_numerator+=(raw_delta[d]/c)*((raw_a[d]+raw_b[d])/c);
    }
    const long double gp=std::sqrt(1+aa),gn=std::sqrt(1+bb);
    const long double an=std::sqrt(aa),bn=std::sqrt(bb),product=an*bn;
    long double angular=0;
    if(an>0&&bn>0){
        if(dot<=0)angular=product-dot;
        else {
            long double cross2=0;
            for(int d=0;d<3;++d){
                const int j=(d+1)%3,k=(d+2)%3;
                const long double x=raw_a[j]*raw_delta[k],y=raw_a[k]*raw_delta[j];
                // Difference-of-products compensation retains near-parallel
                // cross products without first subtracting normalized vectors.
                const long double cross=((x-y)+(std::fma(raw_a[j],raw_delta[k],-x)-std::fma(raw_a[k],raw_delta[j],-y)))/(static_cast<long double>(c)*c);
                cross2+=cross*cross;
            }
            angular=cross2/(product+dot);
        }
    }
    // gamma_rel-1 = cosh(rapidity1-rapidity2)-1 + |u1||u2|(1-cos theta)/c².
    // Rationalized radial term and squared unit-vector difference avoid subtracting
    // large equal four-vector components for comoving relativistic particles.
    const long double radial=(an+bn)>0?radial_numerator/(an+bn):0;
    const long double dg=radial*radial/(gp*gn+product+1)+angular;
    const long double energy=p.mass_kg*static_cast<long double>(c)*c*dg/qe;
    const long double flux=c*std::sqrt(dg*(dg+2))/(gp*gn);
    Kinematics out{static_cast<double>(energy),static_cast<double>(flux)};
    require(std::isfinite(out.target_frame_energy_eV)&&std::isfinite(out.moller_flux_m_s),"relative kinematics overflow");
    require((energy==0||out.target_frame_energy_eV>0)&&(flux==0||out.moller_flux_m_s>0),"relative kinematics conversion underflow");
    return out;
}
struct Remainder {
    double high=0,low=0;
    double require_single_weight()const{
        require(low==0,"partial packet requires compensated weight storage");return high;
    }
};
inline Remainder difference(double a,double b){
    const double high=a-b,b_virtual=a-high,a_virtual=high+b_virtual;
    return {high,(a-a_virtual)+(b_virtual-b)};
}
struct Candidate {
    std::string projectile_id,target_id,channel_id,data_sha256,kernel_sha256;
    std::uint64_t cell=0;
    double event_weight=0,projectile_before=0,target_before=0;
    Remainder projectile_remainder{},target_remainder{};
    double rate_s=0,energy_eV=0,flux_m_s=0,sigma_m2=0,volume_m3=0;
};
// Caps bound input validation/index construction before internal allocation.
// Candidate string bytes count the copied payload plus five NUL terminators;
// row storage remains separately bounded by candidate_limit. These logical
// budgets are not a claim about allocator overhead or process resident memory.
struct CandidateInputLimits {
    std::size_t projectiles=1000000,targets=1000000,cells=1000000,channels=16384;
    std::size_t input_string_bytes=256*1024*1024;
    std::size_t candidate_string_bytes=256*1024*1024;
};
struct CandidateWork {
    std::size_t input_string_bytes=0,cell_lookups=0,target_index_entries=0;
    std::size_t channel_index_entries=0,projectiles_visited=0,pair_visits=0;
    std::size_t kinematics_calls=0,channel_group_lookups=0,channel_queries=0;
    std::size_t candidate_rows=0,candidate_string_bytes=0;
};
inline std::vector<Candidate> candidates(const std::vector<Packet>& projectiles,
    const std::vector<Packet>& targets,const std::vector<Cell>& cells,
    const std::vector<Channel>& channels,std::size_t pair_limit=1000000,
    std::size_t query_limit=1000000,std::size_t candidate_limit=1000000,
    const CandidateInputLimits& limits=CandidateInputLimits{},CandidateWork* work=nullptr) {
    CandidateWork measured{};
    if(work)*work=measured;
    // All these checks and byte accounting precede map/set/vector allocation.
    require(projectiles.size()<=limits.projectiles,"collision projectile input limit");
    require(targets.size()<=limits.targets,"collision target input limit");
    require(cells.size()<=limits.cells,"collision cell input limit");
    require(channels.size()<=limits.channels,"collision channel input limit");
    std::size_t input_bytes=0;
    auto input_string=[&](const std::string& value){
        require(value.size()<=limits.input_string_bytes-input_bytes,"collision input string byte limit");
        input_bytes+=value.size();
    };
    for(const auto* list:{&projectiles,&targets})for(const auto& p:*list){input_string(p.id);input_string(p.state);}
    for(const auto& ch:channels){input_string(ch.id);input_string(ch.projectile_state);input_string(ch.target_state);input_string(ch.data_sha256);input_string(ch.kernel_sha256);}
    measured.input_string_bytes=input_bytes;if(work)*work=measured;
    struct CellIndex {const Cell* cell;std::vector<const Packet*> targets;};
    std::map<std::uint64_t,CellIndex> cell_index;
    std::set<std::string_view> ids,channel_ids;
    for(const auto& cell:cells){require(cell_index.emplace(cell.id,CellIndex{&cell,{}}).second,"duplicate cell");positive(cell.volume_m3,"physical cell volume");}
    for(const auto* list:{&projectiles,&targets})for(const auto& p:*list){
        require(!p.id.empty()&&!p.state.empty()&&ids.insert(p.id).second,"missing or duplicate particle identity");
        positive(p.weight,"represented weight");positive(p.mass_kg,"particle mass");
        ++measured.cell_lookups;if(work)*work=measured;
        require(cell_index.count(p.cell),"particle cell absent");
        for(double u:p.proper_velocity_m_s)require(std::isfinite(u),"nonfinite particle momentum");
    }
    for(const auto& ch:channels){
        require(!ch.id.empty()&&channel_ids.insert(ch.id).second,"duplicate channel");
        require(!ch.projectile_state.empty()&&!ch.target_state.empty(),"missing channel reactant states");
        require(ch.observable==Observable::ExclusiveEvent,"nonexclusive observable is not an event cross section");
        require(digest(ch.data_sha256)&&digest(ch.kernel_sha256)&&bool(ch.sigma_m2),"unbound data/product kernel");
    }
    // Each local list retains the original target order. Projectiles and matching
    // channel lists also retain input order, preserving the legacy query sequence.
    for(const auto& n:targets){
        ++measured.cell_lookups;
        cell_index.find(n.cell)->second.targets.push_back(&n);
        ++measured.target_index_entries;if(work)*work=measured;
    }
    using States=std::pair<std::string_view,std::string_view>;
    std::map<States,std::vector<const Channel*>> channel_index;
    for(const auto& ch:channels){channel_index[{ch.projectile_state,ch.target_state}].push_back(&ch);++measured.channel_index_entries;}
    if(work)*work=measured;
    std::vector<Candidate> out;std::size_t pairs=0,queries=0,output_bytes=0;
    for(const auto& p:projectiles){
        ++measured.projectiles_visited;++measured.cell_lookups;
        const auto& local=cell_index.find(p.cell)->second;
        if(work)*work=measured;
        for(const Packet* pointer:local.targets){
            const auto& n=*pointer;
            require(pairs<pair_limit,"collision pair work limit");++pairs;
            ++measured.pair_visits;++measured.kinematics_calls;if(work)*work=measured;
            // Preserve even unmatched-state kinematic validation and pair budgets.
            const auto kin=relative(p,n);
            ++measured.channel_group_lookups;if(work)*work=measured;
            const auto group=channel_index.find({p.state,n.state});
            if(group==channel_index.end())continue;
            for(const Channel* ch_pointer:group->second){
                const auto& ch=*ch_pointer;
                require(queries<query_limit,"collision query work limit");++queries;
                ++measured.channel_queries;if(work)*work=measured;
                // Query even at zero flux: missing low-energy support is not hidden.
                const double sigma=ch.sigma_m2(kin.target_frame_energy_eV);
                require(std::isfinite(sigma)&&sigma>=0,"negative or nonfinite cross section");
                const double w=std::min(p.weight,n.weight),big=std::max(p.weight,n.weight);
                const long double rate=static_cast<long double>(big)*sigma*kin.moller_flux_m_s/local.cell->volume_m3;
                require(std::isfinite(rate)&&rate<=std::numeric_limits<double>::max(),"collision propensity overflow");
                const double rounded=static_cast<double>(rate);
                require(rate==0||rounded>0,"collision propensity underflow");
                if(rounded>0){
                    require(out.size()<candidate_limit,"collision candidate allocation limit");
                    std::size_t row_bytes=5; // one terminator for each copied string
                    require(row_bytes<=limits.candidate_string_bytes,"collision candidate string byte limit");
                    for(const std::string* value:{&p.id,&n.id,&ch.id,&ch.data_sha256,&ch.kernel_sha256}){
                        require(value->size()<=limits.candidate_string_bytes-row_bytes,"collision candidate string byte limit");
                        row_bytes+=value->size();
                    }
                    require(row_bytes<=limits.candidate_string_bytes-output_bytes,"collision candidate string byte limit");
                    out.push_back({p.id,n.id,ch.id,ch.data_sha256,ch.kernel_sha256,p.cell,w,p.weight,n.weight,difference(p.weight,w),difference(n.weight,w),rounded,kin.target_frame_energy_eV,kin.moller_flux_m_s,sigma,local.cell->volume_m3});
                    output_bytes+=row_bytes;
                    ++measured.candidate_rows;measured.candidate_string_bytes=output_bytes;if(work)*work=measured;
                }
            }
        }
    }
    std::sort(out.begin(),out.end(),[](const Candidate&a,const Candidate&b){
        return std::tie(a.cell,a.projectile_id,a.target_id,a.channel_id)<std::tie(b.cell,b.projectile_id,b.target_id,b.channel_id);});
    return out;
}

class Scheduler;
class Proposal {
    friend class Scheduler;
    std::shared_ptr<const char> lifetime_;std::uint64_t generation_=0;
    std::mt19937_64 random_after_;Clock clock_;std::optional<Candidate> event_;
    double event_subflow_elapsed_s_=0,total_rate_s_=0;
public:
    const std::optional<Candidate>& event()const{return event_;}
    const Clock& clock()const{return clock_;}
    double event_subflow_elapsed_s()const{return event_subflow_elapsed_s_;}
    double frozen_physical_stage_s()const{return clock_.frozen_physical_stage_s;}
    double total_rate_s()const{return total_rate_s_;}
};
class Scheduler {
    std::mt19937_64 random_;std::uint64_t generation_=0;
    const std::shared_ptr<const char> lifetime_=std::make_shared<const char>(0);
    static double uniform_open(std::mt19937_64& rng){
        // 52-bit midpoint rule ensures strict (0,1) in binary64.
        return (static_cast<double>(rng()>>12)+.5)*0x1p-52;
    }
public:
    explicit Scheduler(std::uint64_t seed):random_(seed){}
    Scheduler(const Scheduler&)=delete;Scheduler& operator=(const Scheduler&)=delete;
    Proposal prepare(const Clock& clock,const std::vector<Candidate>& input,
        std::size_t candidate_limit=1000000)const {
        require(input.size()<=candidate_limit,"scheduler candidate input limit");
        require(std::isfinite(clock.begin_s)&&clock.begin_s>=0&&std::isfinite(clock.end_s)&&clock.end_s>clock.begin_s,"invalid accepted interval");
        require(digest(clock.state_sha256)&&digest(clock.mesh_sha256)&&digest(clock.mechanism_sha256),"unbound input revision");
        require(std::isfinite(clock.frozen_physical_stage_s)&&clock.frozen_physical_stage_s>=clock.begin_s&&clock.frozen_physical_stage_s<=clock.end_s,"frozen physical stage outside accepted interval");
        const double duration=clock.end_s-clock.begin_s;
        require(std::isfinite(clock.collision_elapsed_s)&&clock.collision_elapsed_s>=0&&clock.collision_elapsed_s<duration,"invalid collision-subflow elapsed time");
        Proposal result;result.lifetime_=lifetime_;result.generation_=generation_;result.random_after_=random_;result.clock_=clock;
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
            positive(a.projectile_before,"projectile input weight");positive(a.target_before,"target input weight");
            require(a.event_weight==std::min(a.projectile_before,a.target_before),"binary event must consume the smaller packet");
            const auto ep=difference(a.projectile_before,a.event_weight),en=difference(a.target_before,a.event_weight);
            require(a.projectile_remainder.high==ep.high&&a.projectile_remainder.low==ep.low&&a.target_remainder.high==en.high&&a.target_remainder.low==en.low,"incorrect compensated remainder");
            const long double expected=static_cast<long double>(std::max(a.projectile_before,a.target_before))*a.sigma_m2*a.flux_m_s/a.volume_m3;
            require(std::abs(expected-a.rate_s)<=8*std::numeric_limits<double>::epsilon()*expected,"inconsistent candidate propensity");
            require(keys.emplace(a.cell,a.projectile_id,a.target_id,a.channel_id).second,"duplicate candidate");
            total+=a.rate_s;
        }
        require(total<=std::numeric_limits<double>::max(),"total rate overflow");
        result.total_rate_s_=static_cast<double>(total);
        if(total==0)return result;
        const long double wait=-std::log(uniform_open(result.random_after_))/total;
        if(wait>=static_cast<long double>(duration)-clock.collision_elapsed_s)return result;
        result.event_subflow_elapsed_s_=static_cast<double>(clock.collision_elapsed_s+wait);
        require(result.event_subflow_elapsed_s_>clock.collision_elapsed_s&&result.event_subflow_elapsed_s_<duration,"collision clock resolution exhausted");
        const long double choice=uniform_open(result.random_after_)*total;
        long double cumulative=0;
        for(const auto& a:canonical){cumulative+=a.rate_s;if(choice<cumulative){result.event_=a;break;}}
        require(result.event_.has_value(),"candidate sum lost selection");
        return result;
    }
    bool can_publish(const Proposal& p)const noexcept{return p.lifetime_==lifetime_&&p.generation_==generation_&&generation_!=std::numeric_limits<std::uint64_t>::max();}
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
