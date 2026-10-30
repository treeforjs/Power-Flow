#pragma once
#include "exact_particle_count.h"
#include "raikou/json.h"
#include <openssl/evp.h>
#include <array>
#include <bit>
#include <cstdint>
#include <cfenv>
#ifdef __FAST_MATH__
#error "streaming evidence requires finite/underflow checks; fast math is unsupported"
#endif
#include <limits>
#include <map>
#include <optional>
#include <string>
#include <tuple>

// Evidence only. No cross-section evaluation, channel choice, rate, RNG,
// inventory/field mutation or new accepted owner is implemented here.
namespace raikou::chemistry_evidence {
using Json=raikou::Json;
using Exact=mfem_backend::ExactParticleCount;
inline void require(bool x,const char*s){if(!x)throw std::invalid_argument(s);}
inline void floating_environment(){
    require(std::fegetround()==FE_TONEAREST,"evidence requires round-to-nearest floating environment");
    volatile double tiny=std::numeric_limits<double>::denorm_min(),one=1.;
    const double copied=tiny*one,doubled=tiny+tiny;
    require(std::bit_cast<std::uint64_t>(copied)==1&&std::bit_cast<std::uint64_t>(doubled)==2,"evidence requires preserved binary64 subnormals");
}
inline bool digest(const std::string&s){return s.size()==64&&s.find_first_not_of("0123456789abcdef")==std::string::npos;}
class Digest {
    using Pointer=std::unique_ptr<EVP_MD_CTX,decltype(&EVP_MD_CTX_free)>;
    Pointer p_{EVP_MD_CTX_new(),EVP_MD_CTX_free};
public:
    Digest(){require(bool(p_)&&EVP_DigestInit_ex(p_.get(),EVP_sha256(),nullptr)==1,"SHA initialization");}
    Digest(const Digest&o):Digest(){require(EVP_MD_CTX_copy_ex(p_.get(),o.p_.get())==1,"SHA copy");}
    Digest(Digest&&) noexcept=default;Digest&operator=(Digest&&) noexcept=default;
    void add(const std::string&s){require(EVP_DigestUpdate(p_.get(),s.data(),s.size())==1,"SHA update");}
    std::string hex()const{Digest copy(*this);unsigned char result[EVP_MAX_MD_SIZE];unsigned n=0;require(EVP_DigestFinal_ex(copy.p_.get(),result,&n)==1&&n==32,"SHA final");static constexpr char h[]="0123456789abcdef";std::string s;s.reserve(64);for(unsigned i=0;i<n;++i){s+=h[result[i]>>4];s+=h[result[i]&15];}return s;}
};
inline std::string sha256(const std::string&s){Digest d;d.add(s);return d.hex();}
inline void u64(std::string&s,std::uint64_t n){for(int shift=56;shift>=0;shift-=8)s+=char((n>>shift)&255);}
inline void word(std::string&s,const std::string&x){u64(s,x.size());s+=x;}
inline Json exact_json(const Exact&x){
    struct Free {void operator()(char*p)const noexcept{OPENSSL_free(p);}};
    const auto units=x.copy_units();std::unique_ptr<char,Free> text(BN_bn2hex(units.get()));require(bool(text),"exact evidence encoding");
    return {{"unit_exponent",-1074},{"nonnegative_units_hex",text.get()},{"rounded_binary64",x.rounded_for_rate()}};
}
enum class Category:std::size_t {SelectedPositive,SelectedZero,SelectedUnsupported,UnselectedPositive,UnselectedZero,UnselectedUnsupported,UnqualifiedNoRate,Count};
inline constexpr std::array<const char*,7> category_names={"selected_positive","selected_zero","selected_unsupported","unselected_positive","unselected_zero","unselected_unsupported","unqualified_no_rate"};
struct Key {
    std::string lane,id,model,unit;
    auto tuple()const{return std::tie(lane,id,model,unit);}
    bool operator<(const Key&o)const{return tuple()<o.tuple();}
    bool operator==(const Key&o)const{return tuple()==o.tuple();}
};
struct Observation {Key key;bool selected=false,known=false,unqualified=false;double value=0;};
struct Query {std::uint64_t ordinal=0;std::string id,context_sha256;double energy_eV=0,weight=0;};
struct Limits {
    std::uint64_t queries=1000000,encoded_bytes=16ULL*1024*1024*1024;
    std::size_t record_bytes=262144;
};
struct Stats {
    std::uint64_t count=0,positive_weight_count=0,first=0,last=0;
    Exact weight,energy_product,value_product;
    double minimum_energy=0,maximum_energy=0,minimum_weighted_energy=0,maximum_weighted_energy=0,minimum_value=0,maximum_value=0;
    bool value_seen=false;
    static double product(double a,double b){
        const long double wide=static_cast<long double>(a)*b;const double rounded=a*b;
        require(std::isfinite(rounded)&&(wide==0||rounded>0),"evidence moment product overflow/positive underflow");return rounded;
    }
    static Exact sum(const Exact&old,double x){auto next=old.credited(x);(void)next.rounded_for_rate();return next;}
    void add(const Query&q,const Observation&o){
        require(count!=UINT64_MAX,"evidence category count overflow");
        weight=sum(weight,q.weight);energy_product=sum(energy_product,product(q.weight,q.energy_eV));
        if(o.known){value_product=sum(value_product,product(q.weight,o.value));if(!value_seen){minimum_value=maximum_value=o.value;value_seen=true;}else{minimum_value=std::min(minimum_value,o.value);maximum_value=std::max(maximum_value,o.value);}}
        if(count==0){first=q.ordinal;minimum_energy=maximum_energy=q.energy_eV;}else{minimum_energy=std::min(minimum_energy,q.energy_eV);maximum_energy=std::max(maximum_energy,q.energy_eV);}last=q.ordinal;++count;
        if(q.weight>0){if(positive_weight_count==0)minimum_weighted_energy=maximum_weighted_energy=q.energy_eV;else{minimum_weighted_energy=std::min(minimum_weighted_energy,q.energy_eV);maximum_weighted_energy=std::max(maximum_weighted_energy,q.energy_eV);}++positive_weight_count;}
    }
    Json json()const {
        return {{"count",count},{"positive_weight_count",positive_weight_count},{"weight_exact",exact_json(weight)},
          {"energy_weight_product_sum",exact_json(energy_product)},
          {"source_value_weight_product_sum",value_seen?exact_json(value_product):Json(nullptr)},
          {"energy_eV_min",count?Json(minimum_energy):Json(nullptr)},{"energy_eV_max",count?Json(maximum_energy):Json(nullptr)},
          {"positive_weight_energy_eV_min",positive_weight_count?Json(minimum_weighted_energy):Json(nullptr)},{"positive_weight_energy_eV_max",positive_weight_count?Json(maximum_weighted_energy):Json(nullptr)},
          {"source_value_min",value_seen?Json(minimum_value):Json(nullptr)},{"source_value_max",value_seen?Json(maximum_value):Json(nullptr)},
          {"first_query_ordinal",count?Json(first):Json(nullptr)},{"last_query_ordinal",count?Json(last):Json(nullptr)}};
    }
};
class Stream {
    using Bins=std::array<Stats,static_cast<std::size_t>(Category::Count)>;
    struct State {
        Digest commitment;std::vector<Bins> bins;Exact query_weight;
        std::uint64_t queries=0,observations=0,bytes=0;std::size_t maximum_record_bytes=0;
        explicit State(std::size_t n):bins(n){}
        void swap(State&o) noexcept {std::swap(commitment,o.commitment);bins.swap(o.bins);std::swap(query_weight,o.query_weight);std::swap(queries,o.queries);std::swap(observations,o.observations);std::swap(bytes,o.bytes);std::swap(maximum_record_bytes,o.maximum_record_bytes);}
    };
    std::vector<Key> layout_;Json bindings_,policy_;Limits limits_;State state_;std::string header_;
    static constexpr std::size_t maximum_keys=64,maximum_text_bytes=32768;
    struct Parser {
        std::size_t text_bytes=0;
        const std::string&text(const Json&j,const char*k,std::size_t cap=4096){require(j.contains(k)&&j.at(k).is_string(),"evidence string field");const auto&s=j.at(k).get_ref<const std::string&>();require(s.size()<=cap&&s.size()<=maximum_text_bytes-text_bytes,"evidence text cap");text_bytes+=s.size();return s;}
        bool boolean(const Json&j,const char*k){require(j.contains(k)&&j.at(k).is_boolean(),"evidence boolean field");return j.at(k).get<bool>();}
        static double number(const Json&j,const char*k){require(j.contains(k)&&j.at(k).is_number(),"evidence numeric field");const double x=j.at(k).get<double>();require(std::isfinite(x)&&x>=0,"evidence nonnegative finite operand");return x;}
        void reason(const Json&j,const char*k){(void)text(j,k,8192);}
    };
    static void fields(const Json&j,std::initializer_list<const char*> names){require(j.is_object()&&j.size()==names.size(),"evidence exact record fields");for(auto*n:names)require(j.contains(n),"missing evidence record field");}
    static Category category(const Observation&o){if(o.unqualified)return Category::UnqualifiedNoRate;if(!o.known)return o.selected?Category::SelectedUnsupported:Category::UnselectedUnsupported;if(o.selected)return o.value>0?Category::SelectedPositive:Category::SelectedZero;return o.value>0?Category::UnselectedPositive:Category::UnselectedZero;}
    static Json policy(const Json&j){return {{"scope",j.at("scope")},{"complete_mechanism",j.at("complete_mechanism")},{"one_keV_policy",j.at("one_keV_policy")},{"site_policy",j.at("site_policy")}};}
    static std::vector<Observation> parse(const Json&j){
        fields(j,{"energy_eV","scope","complete_mechanism","one_keV_policy","site_policy","conditional_models","unapplied_observables","missing_rates"});Parser p;(void)Parser::number(j,"energy_eV");
        require(p.text(j,"scope",256)=="supported_conditional_comparison_subset"&&!p.boolean(j,"complete_mechanism"),"evidence cannot certify a complete mechanism");(void)p.text(j,"one_keV_policy",256);(void)p.text(j,"site_policy",512);
        std::vector<Observation> out;out.reserve(maximum_keys);
        for(auto name:{"conditional_models","unapplied_observables","missing_rates"})require(j.at(name).is_array()&&j.at(name).size()<=maximum_keys,"evidence row cap");
        auto add=[&](Observation x){require(out.size()<maximum_keys&&!x.key.id.empty(),"evidence key cap/identity");out.push_back(std::move(x));};
        for(const auto&r:j.at("conditional_models")){
            const bool known=p.boolean(r,"source_rate_known");
            if(known)fields(r,{"channel","model","selected_subset","event_kernel_in_subset_support","source_rate_known","source_sigma_m2"});
            else{fields(r,{"channel","model","selected_subset","event_kernel_in_subset_support","source_rate_known","source_sigma_m2","reason"});p.reason(r,"reason");}
            Observation o;o.key={"conditional_model",p.text(r,"channel"),p.text(r,"model"),"m2"};o.selected=p.boolean(r,"selected_subset");require(o.selected==p.boolean(r,"event_kernel_in_subset_support"),"conditional selection flags disagree");o.known=known;if(known)o.value=Parser::number(r,"source_sigma_m2");else require(r.at("source_sigma_m2").is_null(),"unknown source is not null");add(std::move(o));
        }
        for(const auto&r:j.at("unapplied_observables")){
            const bool known=p.boolean(r,"source_value_known");
            if(known)fields(r,{"channel","unit_SI","applied_as_event","rate_sum_eligible","reason","source_value_known","source_value_SI"});
            else{fields(r,{"channel","unit_SI","applied_as_event","rate_sum_eligible","reason","source_value_known","source_value_SI","coverage_reason"});p.reason(r,"coverage_reason");}
            require(!p.boolean(r,"applied_as_event")&&!p.boolean(r,"rate_sum_eligible"),"diagnostic observable promoted to event");p.reason(r,"reason");
            Observation o;o.key={"unapplied_observable",p.text(r,"channel"),"",p.text(r,"unit_SI",64)};require(o.key.unit=="m2"||o.key.unit=="eV m2","unadmitted evidence unit");o.known=known;if(known)o.value=Parser::number(r,"source_value_SI");else require(r.at("source_value_SI").is_null(),"unknown observable is not null");add(std::move(o));
        }
        for(const auto&r:j.at("missing_rates")){
            const auto status=p.text(r,"status",64);
            if(status=="unknown_at_energy"){fields(r,{"channel","status","sigma_m2","applied_as_event","reason"});p.reason(r,"reason");}
            else fields(r,{"channel","status","sigma_m2","applied_as_event"});
            require(status=="known_unmapped_rate"||status=="unknown_at_energy"||status=="unknown_or_unqualified","unsupported missing-rate status");require(!p.boolean(r,"applied_as_event"),"missing rate promoted to event");
            Observation o;o.key={"missing_rate_contract",p.text(r,"channel"),"","m2"};o.known=status=="known_unmapped_rate";o.unqualified=status=="unknown_or_unqualified";
            if(o.known)o.value=Parser::number(r,"sigma_m2");else require(r.at("sigma_m2").is_null(),"missing rate converted to a value");add(std::move(o));
        }
        std::sort(out.begin(),out.end(),[](const auto&a,const auto&b){return a.key<b.key;});
        for(std::size_t i=1;i<out.size();++i)require(!(out[i].key==out[i-1].key),"duplicate evidence key");
        require(!out.empty(),"empty evidence inventory");return out;
    }
public:
    // Bindings are caller-verified input digests, never owner authentication.
    Stream(const Json&reference,const Json&bindings,Limits limits=Limits{})
      :limits_(limits),state_(0){
        fields(bindings,{"dataset_sha256","provider_source_sha256","kernel_manifest_sha256","kernel_binding_file_sha256","mechanism_contract_sha256","missing_rate_contract_sha256","query_input_sha256","weight_measure","evidence_component_sha256","exact_arithmetic_source_sha256"});
        for(const auto&name:{"dataset_sha256","provider_source_sha256","kernel_manifest_sha256","kernel_binding_file_sha256","mechanism_contract_sha256","missing_rate_contract_sha256","query_input_sha256","evidence_component_sha256","exact_arithmetic_source_sha256"})require(bindings.at(name).is_string()&&digest(bindings.at(name).get_ref<const std::string&>()),"unbound evidence input");
        require(bindings.at("weight_measure").is_string()&&!bindings.at("weight_measure").get_ref<const std::string&>().empty()&&bindings.at("weight_measure").get_ref<const std::string&>().size()<=256,"unbound evidence weight meaning");
        require(limits_.queries>0&&limits_.queries<=100000000&&limits_.encoded_bytes>0&&limits_.record_bytes>0&&limits_.record_bytes<=262144,"evidence admission limits");
        floating_environment();
        const auto observed=parse(reference);bindings_=bindings;policy_=policy(reference);
        for(const auto&o:observed)layout_.push_back(o.key);
        state_=State(layout_.size());
        Json keys=Json::array();for(const auto&k:layout_)keys.push_back({k.lane,k.id,k.model,k.unit});
        header_=Json{{"schema","streaming-generic-v6-evidence-v1"},{"bindings",bindings_},{"policy",policy_},{"layout",keys},{"limits",{{"queries",limits_.queries},{"encoded_bytes",limits_.encoded_bytes},{"record_bytes",limits_.record_bytes}}}}.dump();
        std::string record="CHEM-EVIDENCE-V1";word(record,header_);state_.commitment.add(record);
    }
    Stream(const Stream&)=delete;Stream&operator=(const Stream&)=delete;
    void observe(const Query&q,const Json&coverage){
        floating_environment();
        require(q.ordinal==state_.queries&&state_.queries<limits_.queries,"evidence ordinal/query cap");require(!q.id.empty()&&q.id.size()<=4096&&digest(q.context_sha256),"unbound query identity/context");require(std::isfinite(q.energy_eV)&&q.energy_eV>=0&&std::isfinite(q.weight)&&q.weight>=0,"invalid evidence query operands");
        const auto observations=parse(coverage);require(observations.size()==layout_.size()&&policy(coverage)==policy_,"evidence inventory/policy changed");
        require(std::bit_cast<std::uint64_t>(coverage.at("energy_eV").get<double>())==std::bit_cast<std::uint64_t>(q.energy_eV),"query energy differs from evidence");
        for(std::size_t i=0;i<layout_.size();++i)require(observations[i].key==layout_[i],"evidence key omitted/replaced");
        const auto canonical=coverage.dump();std::string encoded;encoded.reserve(canonical.size()+q.id.size()+128);encoded+='Q';u64(encoded,q.ordinal);word(encoded,q.id);word(encoded,q.context_sha256);u64(encoded,std::bit_cast<std::uint64_t>(q.energy_eV));u64(encoded,std::bit_cast<std::uint64_t>(q.weight));word(encoded,canonical);
        require(encoded.size()<=limits_.record_bytes&&encoded.size()<=limits_.encoded_bytes-state_.bytes,"evidence byte cap");
        require(observations.size()<=UINT64_MAX-state_.observations,"evidence observation count overflow");
        State next=state_;next.query_weight=Stats::sum(next.query_weight,q.weight);
        for(std::size_t i=0;i<observations.size();++i)next.bins[i][static_cast<std::size_t>(category(observations[i]))].add(q,observations[i]);
        next.commitment.add(encoded);++next.queries;next.observations+=observations.size();next.bytes+=encoded.size();next.maximum_record_bytes=std::max(next.maximum_record_bytes,encoded.size());state_.swap(next);
    }
    std::string commitment()const{return state_.commitment.hex();}
    std::uint64_t queries()const noexcept{return state_.queries;}
    std::size_t retained_keys()const noexcept{return layout_.size();}
    Json snapshot()const{
        Json rows=Json::array();for(std::size_t i=0;i<layout_.size();++i){Json bins=Json::object();for(std::size_t c=0;c<category_names.size();++c)bins[category_names[c]]=state_.bins[i][c].json();const auto&k=layout_[i];rows.push_back({{"lane",k.lane},{"id",k.id},{"model",k.model},{"unit",k.unit},{"categories",std::move(bins)}});}
        return {{"canonical_header_json",header_},{"schema","streaming-generic-v6-evidence-receipt-v1"},{"bindings",bindings_},{"policy",policy_},{"queries",state_.queries},{"observations",state_.observations},{"query_weight_exact",exact_json(state_.query_weight)},
          {"content_sha256",commitment()},{"encoded_query_bytes",state_.bytes},{"maximum_encoded_query_bytes",state_.maximum_record_bytes},{"retained_keys",layout_.size()},{"retained_category_slots",layout_.size()*category_names.size()},
          {"query_weight_semantics","per query; repeated particle/pair queries are not a unique population census"},{"moment_semantics","exact sum of per-term binary64 rounded products; positive product underflow or finite-sum overflow refuses atomically"},
          {"lanes_are_additive_event_rates",false},{"complete_mechanism",false},{"rows",std::move(rows)}};
    }
};
} // namespace raikou::chemistry_evidence
