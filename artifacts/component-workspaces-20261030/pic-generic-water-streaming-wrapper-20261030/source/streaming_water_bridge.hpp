#pragma once
#include "generic_water_bridge.hpp"
#include "streaming_coverage.hpp"
#include "observed_candidates.hpp"
namespace raikou::streaming_water {
namespace g=generic_water;namespace e=chemistry_evidence;using Json=e::Json;
namespace k=g::k;
struct EvidenceRequest {
    Json stream_bindings;
    std::string scope_context_sha256,wrapper_source_sha256,observed_builder_source_sha256;
    e::Limits limits;
};
struct CandidateBatch {
    std::vector<pc::Candidate> rows;Json evidence;pc::CandidateWork work;
    bool complete_mechanism=false;
};
inline void packet_bytes(std::string&bytes,const pc::Packet&p){
    e::word(bytes,p.id);e::word(bytes,p.state);e::u64(bytes,p.cell);
    for(double x:{p.weight,p.mass_kg})e::u64(bytes,std::bit_cast<std::uint64_t>(x));
    for(double x:p.proper_velocity_m_s)e::u64(bytes,std::bit_cast<std::uint64_t>(x));
}
inline std::string pair_id(const pc::Packet&p,const pc::Packet&n,const pc::Cell&cell,const pc::Kinematics&kin){
    std::string bytes="generic-v6-physical-pair-v1";packet_bytes(bytes,p);packet_bytes(bytes,n);e::u64(bytes,cell.id);
    for(double x:{cell.volume_m3,kin.target_frame_energy_eV,kin.moller_flux_m_s})e::u64(bytes,std::bit_cast<std::uint64_t>(x));
    return "physical-pair/"+e::sha256(bytes);
}
class Bridge {
    using Data=g::k::StrictCollisionDataset;
    std::shared_ptr<const Data> data_;
    g::Bridge physical_;
    double M_;
    Json reference_;
    static bool active(const g::Spec&s,double E){
        // Exact frozen-v6 evidence predicate; physical selection stays delegated.
        return E>=s.low_eV&&E<=s.high_eV&&!(s.model=="parent_ground_equal_pair_CM_v1"&&E>=1000);
    }
    void valid_packet(const g::PhysicalPacket&p,bool electron)const {
        pc::require(p.packet.state==(electron?"electron":"H2O")&&p.packet.mass_kg==(electron?g::me:M_),"generic species/catalog mass mismatch");
        pc::require(!p.packet.id.empty()&&p.packet.id.size()<=4096,"generic packet identity cap");
        for(double x:p.position_m)pc::require(std::isfinite(x),"nonfinite collision position");
    }
    struct Memo {
        struct Value {bool known=false;double number=0;std::string reason;};
        const Data& data;double E;
        std::map<std::string,Value> values;
        std::uint64_t diagnostic_queries=0,reuses=0,seeds=0;
        Memo(const Data&d,double energy):data(d),E(energy){}
        void seed(const std::string&id,double value){
            pc::require(values.emplace(id,Value{true,value,{}}).second,"duplicate selected evidence return");++seeds;
        }
        double value(const std::string&id){
            auto it=values.find(id);
            if(it==values.end()){
                pc::require(values.size()<64,"diagnostic lookup memo cap");Value x;
                ++diagnostic_queries;
                try{x.number=data.observable_value_SI(id,E);x.known=true;}
                catch(const g::k::CollisionCoverageError&err){x.reason=err.what();}
                it=values.emplace(id,std::move(x)).first;
            }else ++reuses;
            if(!it->second.known)throw g::k::CollisionCoverageError(it->second.reason);
            return it->second.number;
        }
    };
    Json project(double E,Memo&memo)const {
        pc::require(std::isfinite(E)&&E>=0,"invalid coverage energy");Json selected=Json::array(),holes=Json::array();
        for(const auto&s:physical_.specs()){Json r={{"channel",s.channel},{"model",s.model},{"selected_subset",active(s,E)}};
            try{r["source_sigma_m2"]=memo.value(s.channel);r["source_rate_known"]=true;}
            catch(const k::CollisionCoverageError&x){r["source_sigma_m2"]=nullptr;r["source_rate_known"]=false;r["reason"]=x.what();}
            r["event_kernel_in_subset_support"]=active(s,E);selected.push_back(r);}
        Json unmapped={{"channel","H2O:conditional_outer3:unmapped_2a1"},{"applied_as_event",false}};
        try{unmapped["sigma_m2"]=memo.value("H2O:conditional_outer3:unmapped_2a1");unmapped["status"]="known_unmapped_rate";}
        catch(const k::CollisionCoverageError&x){unmapped["sigma_m2"]=nullptr;unmapped["status"]="unknown_at_energy";unmapped["reason"]=x.what();}
        holes.push_back(unmapped);
        for(const auto& name:{"core_1a1","multiple_ionization_Auger","fragment_species_continuation","unmapped_inelastic_and_angular_channels"})holes.push_back({{"channel",name},{"sigma_m2",nullptr},{"status","unknown_or_unqualified"},{"applied_as_event",false}});
        Json unused=Json::array();
        for(const auto& raw:data_->metadata().at("channels")){
            const auto id=raw.at("id").get<std::string>();bool selected_id=false;for(const auto&s:physical_.specs())if(s.channel==id)selected_id=true;
            if(selected_id)continue;
            Json row={{"channel",id},{"unit_SI",raw.at("unit_SI")},{"applied_as_event",false},{"rate_sum_eligible",false},{"reason","unmapped or diagnostic observable; no exclusive event added"}};
            try{row["source_value_SI"]=memo.value(id);row["source_value_known"]=true;}
            catch(const k::CollisionCoverageError&x){row["source_value_SI"]=nullptr;row["source_value_known"]=false;row["coverage_reason"]=x.what();}
            unused.push_back(row);
        }
        return {{"unapplied_observables",unused},{"energy_eV",E},{"scope","supported_conditional_comparison_subset"},{"complete_mechanism",false},
          {"one_keV_policy","outer3_selected_at_exact_1000_eV"},{"conditional_models",selected},{"missing_rates",holes},
          {"site_policy","retain_scattered_electron_at_primary_and_other_products_at_target; capture_requires_exact_3D_coincidence"}};
    }
    struct Recorder {
        const Bridge& owner;const EvidenceRequest& request;
        Json descriptor;
        e::Stream stream;
        std::optional<Memo> memo;
        e::Query query;
        std::size_t channel_returns=0;
        bool open=false,last_seen=false;
        double last_energy=0;
        Json last_coverage;
        std::uint64_t diagnostic_queries=0,memo_reuses=0,selected_return_seeds=0,last_record_reuses=0;
        static Json input_descriptor(const EvidenceRequest&r){return {{"schema","generic-v6-streaming-query-descriptor-v1"},{"host_query_input_sha256",r.stream_bindings.at("query_input_sha256")},{"scope_context_sha256",r.scope_context_sha256},{"wrapper_source_sha256",r.wrapper_source_sha256},{"observed_builder_source_sha256",r.observed_builder_source_sha256},{"weight_measure","projectile-carrier weight per evaluated electron-H2O pair; repeated pairs are not a unique population census"}};}
        static Json bindings(const EvidenceRequest&r,const Json&descriptor){auto b=r.stream_bindings;b["query_input_sha256"]=e::sha256(descriptor.dump());return b;}
        Recorder(const Bridge&o,const EvidenceRequest&r):owner(o),request(r),descriptor(input_descriptor(r)),stream(o.reference_,bindings(r,descriptor),r.limits){}
        void begin_pair(const pc::Packet&p,const pc::Packet&n,const pc::Cell&cell,const pc::Kinematics&kin){
            pc::require(!open,"overlapping private evidence pair");open=true;channel_returns=0;
            memo.emplace(*owner.data_,kin.target_frame_energy_eV);
            query={stream.queries(),pair_id(p,n,cell,kin),request.scope_context_sha256,kin.target_frame_energy_eV,p.weight};
        }
        void channel(const pc::Channel&ch,const pc::Kinematics&kin,double sigma){
            pc::require(open&&memo&&kin.target_frame_energy_eV==query.energy_eV,"private evidence pair identity");
            const auto&specs=owner.physical_.specs();const auto it=std::find_if(specs.begin(),specs.end(),[&](const auto&s){return s.channel==ch.id;});
            pc::require(it!=specs.end(),"foreign physical evidence channel");
            if(active(*it,query.energy_eV))memo->seed(ch.id,sigma);
            ++channel_returns;
        }
        void end_pair(){
            pc::require(open&&memo&&channel_returns==owner.physical_.specs().size(),"incomplete private evidence pair");
            const bool repeated=last_seen&&std::bit_cast<std::uint64_t>(last_energy)==std::bit_cast<std::uint64_t>(query.energy_eV);
            if(repeated){
                // Cached metadata may not override a changed physical return.
                for(const auto&row:last_coverage.at("conditional_models"))if(row.at("selected_subset").get<bool>()){
                    const auto id=row.at("channel").get<std::string>();const auto it=memo->values.find(id);
                    pc::require(it!=memo->values.end()&&it->second.known&&std::bit_cast<std::uint64_t>(it->second.number)==std::bit_cast<std::uint64_t>(row.at("source_sigma_m2").get<double>()),"cached evidence differs from actual physical return");
                }
                stream.observe(query,last_coverage);++last_record_reuses;
            }else{
                auto coverage=owner.project(query.energy_eV,*memo);stream.observe(query,coverage);last_coverage=std::move(coverage);last_energy=query.energy_eV;last_seen=true;
            }
            diagnostic_queries+=memo->diagnostic_queries;memo_reuses+=memo->reuses;selected_return_seeds+=memo->seeds;
            memo.reset();open=false;
        }
        Json receipt()const {
            pc::require(!open,"unfinished evidence pair");
            return {{"schema","generic-v6-streaming-wrapper-evidence-v1"},{"query_input_descriptor",descriptor},{"stream",stream.snapshot()},
              {"cache_scope","one per-pair raw lookup memo plus one last exact-energy coverage record; no distinct-energy container"},
              {"diagnostic_source_queries",diagnostic_queries},{"diagnostic_lookup_reuses",memo_reuses},{"actual_selected_return_seeds",selected_return_seeds},{"last_energy_record_reuses",last_record_reuses},
              {"counters_exclude_internal_physical_provider_lookups",true},{"accepted_native_completion",false},{"complete_mechanism",false}};
        }
    };
public:
    Bridge(std::shared_ptr<const Data> data,std::string dataset,std::string kernel,double M,double H,g::SitePolicy site)
      :data_(std::move(data)),physical_(data_,std::move(dataset),std::move(kernel),M,H,site),M_(M),reference_(physical_.coverage(0)){}
    const g::Bridge& physical_provider()const{return physical_;}
    g::Outcome sample(const pc::Candidate&a,const g::PhysicalPacket&p,const g::PhysicalPacket&n,std::mt19937_64&rng)const{return physical_.sample(a,p,n,rng);}
    g::Outcome sample(const pc::Candidate&a,const g::PhysicalPacket&p,const g::PhysicalPacket&n,const std::vector<double>&u)const{return physical_.sample(a,p,n,u);}
    CandidateBatch candidates(const std::vector<g::PhysicalPacket>&electrons,const std::vector<g::PhysicalPacket>&water,
       const std::vector<pc::Cell>&cells,g::Scope scope,const EvidenceRequest&request,
       std::size_t pair_cap=1000000,std::size_t query_cap=1000000,std::size_t row_cap=100000)const {
        // Frozen-v6 physical input checks and error ordering; only the old
        // distinct-energy evidence cap is absent at this new entrypoint.
        pc::require(electrons.size()<=100000&&water.size()<=100000,"generic population cap");
        std::size_t bytes=0;constexpr std::size_t byte_cap=64*1024*1024;
        auto validate=[&](const std::vector<g::PhysicalPacket>&population,bool electron){for(const auto&p:population){valid_packet(p,electron);for(const auto*text:{&p.packet.id,&p.packet.state}){pc::require(text->size()<=byte_cap-bytes,"generic input identity byte cap");bytes+=text->size();}}};validate(electrons,true);validate(water,false);
        const auto channels=physical_.comparison_channels(scope);
        pc::require(e::digest(request.scope_context_sha256)&&e::digest(request.wrapper_source_sha256)&&e::digest(request.observed_builder_source_sha256),"unbound streaming wrapper evidence context");
        pc::require(request.stream_bindings.at("dataset_sha256")==physical_.dataset_sha256()&&request.stream_bindings.at("kernel_manifest_sha256")==physical_.kernel_sha256(),"physical/evidence data or kernel mismatch");
        pc::require(request.stream_bindings.at("weight_measure")=="projectile-carrier weight per evaluated electron-H2O pair; repeated pairs are not a unique population census","streaming wrapper evidence weight meaning differs");
        Recorder recorder(*this,request);
        std::vector<pc::Packet> e,n;e.reserve(electrons.size());n.reserve(water.size());
        for(const auto&p:electrons){valid_packet(p,true);e.push_back(p.packet);}for(const auto&p:water){valid_packet(p,false);n.push_back(p.packet);}
        CandidateBatch out;pc::CandidateInputLimits caps;caps.projectiles=caps.targets=100000;
        out.rows=observed_candidates(recorder,e,n,cells,channels,pair_cap,query_cap,row_cap,caps,&out.work);
        out.evidence=recorder.receipt();return out;
    }
};
} // namespace raikou::streaming_water
