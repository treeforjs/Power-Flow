#pragma once
#include "collision_scheduler.hpp"
namespace raikou::streaming_water {
namespace pc=pic_collisions;
using pc::Candidate;using pc::Packet;using pc::Cell;using pc::Channel;
using pc::CandidateInputLimits;using pc::CandidateWork;
using pc::Observable;
using pc::require;using pc::positive;using pc::digest;using pc::relative;using pc::difference;
// Exact frozen indexed builder body plus three nonnumeric evidence hooks.
template<class Observer>
inline std::vector<Candidate> observed_candidates(Observer& observer,const std::vector<Packet>& projectiles,
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
            observer.begin_pair(p,n,*local.cell,kin);
            for(const Channel* ch_pointer:group->second){
                const auto& ch=*ch_pointer;
                require(queries<query_limit,"collision query work limit");++queries;
                ++measured.channel_queries;if(work)*work=measured;
                // Query even at zero flux: missing low-energy support is not hidden.
                const double sigma=ch.sigma_m2(kin.target_frame_energy_eV);
                require(std::isfinite(sigma)&&sigma>=0,"negative or nonfinite cross section");
                observer.channel(ch,kin,sigma);
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
            observer.end_pair();
        }
    }
    std::sort(out.begin(),out.end(),[](const Candidate&a,const Candidate&b){
        return std::tie(a.cell,a.projectile_id,a.target_id,a.channel_id)<std::tie(b.cell,b.projectile_id,b.target_id,b.channel_id);});
    return out;
}
} // namespace raikou::streaming_water
