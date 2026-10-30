#pragma once
#include "water_parent_bridge.hpp"
#include "raikou/h2o_dea_kernel.h"
#include "raikou/h2o_two_body_attachment.h"
#include "raikou/h2o_neutral_dissociation.h"
#include "raikou/h2o_common_mass_basis.h"
#include "raikou/h2o_outer_valence.h"
#include "raikou/h2o_elastic_angular.h"
#include <map>
#include <memory>

// Supported conditional comparison subset. No owned particles, field/current
// path, source, physical clock, or accepted RNG publication lives in this port.
namespace raikou::generic_water {
namespace pc=pic_collisions;namespace k=kinetics;namespace lorentz=pic_water;
using Json=raikou::Json;using Vec=pc::Vec;
constexpr double qe=pc::qe,me=lorentz::me,c=pc::c;
enum class Scope { CompareSupportedSubset, RequireCompleteMechanism };
enum class SitePolicy { PreserveReactantSites };
struct State {std::string id;int H=0,O=0,q=0;double mass_kg=0,internal_J=0,reference_J=0;};
struct PhysicalPacket {pc::Packet packet;Vec position_m{};};
struct Spec {
    std::string channel,model;double low_eV=0,high_eV=0;bool captures_primary=false;
    std::size_t variates=0;std::vector<std::vector<std::string>> allowed_product_states;
};
struct Product {State state;Vec proper_velocity_m_s{};bool at_primary=false;double weight=0;};
struct Outcome {
    std::string channel,model,data_sha256,kernel_sha256;double event_weight=0,energy_eV=0,sigma_m2=0;
    pc::Remainder projectile_remainder{},target_remainder{};
    std::vector<Product> products;bool captures_primary=false;int OH_level=-1;
    double weighted_energy_defect_J=0,weighted_energy_bound_J=0,weighted_momentum_bound=0;
    Vec weighted_momentum_defect{};std::size_t random_variates=0;
};
struct CandidateBatch {
    std::vector<pc::Candidate> rows;Json coverage;pc::CandidateWork work;
    bool complete_mechanism=false;
};
class Bridge {
    std::shared_ptr<const k::StrictCollisionDataset> data_;
    std::string dataset_,kernel_;double M_,H_;SitePolicy site_;
    std::vector<State> catalog_;std::map<std::string,State> states_;std::vector<Spec> specs_;
    static bool same(const Vec&a,const Vec&b){return a==b;}
    static double uniform(std::mt19937_64& r){return (static_cast<double>(r()>>12)+.5)*0x1p-52;}
    static Vec isotropic(double u,double v){const double z=2*u-1,a=6.283185307179586476925286766559*v,r=std::sqrt((1-z)*(1+z));return {r*std::cos(a),r*std::sin(a),z};}
    bool active(const Spec&s,double E)const {
        // A single declared policy chooses outer3 at the exact 1keV boundary.
        return E>=s.low_eV&&E<=s.high_eV&&!(s.model=="parent_ground_equal_pair_CM_v1"&&E>=1000);
    }
    double subset_sigma(const Spec&s,double E)const {
        pc::require(std::isfinite(E)&&E>=0,"invalid target-frame energy");
        if(!active(s,E))return 0; // Excluded comparison contribution; see coverage().
        const double sigma=data_->observable_value_SI(s.channel,E);
        if(sigma>0)data_->require_MCC_event_kernel(s.channel,s.model,E);
        return sigma;
    }
    pc::Channel channel(const Spec&s)const {
        // Capture by value; returned callbacks do not depend on Bridge lifetime.
        const auto data=data_;const auto dataset=dataset_,kernel=kernel_;
        return {s.channel,"electron","H2O",dataset,kernel,pc::Observable::ExclusiveEvent,
            [data,s](double E){pc::require(std::isfinite(E)&&E>=0,"invalid target-frame energy");
                if(E<s.low_eV||E>s.high_eV||(s.model=="parent_ground_equal_pair_CM_v1"&&E>=1000))return 0.;
                const double x=data->observable_value_SI(s.channel,E);if(x>0)data->require_MCC_event_kernel(s.channel,s.model,E);return x;}};
    }
    void valid_packet(const PhysicalPacket&p,bool electron)const {
        pc::require(p.packet.state==(electron?"electron":"H2O")&&p.packet.mass_kg==(electron?me:M_),"generic species/catalog mass mismatch");
        pc::require(!p.packet.id.empty()&&p.packet.id.size()<=4096,"generic packet identity cap");
        for(double x:p.position_m)pc::require(std::isfinite(x),"nonfinite collision position");
    }
    const Spec& spec(const std::string&id)const {for(const auto&s:specs_)if(s.channel==id)return s;throw pc::Error("unknown generic channel");}
public:
    // Host must verify the input file and source bytes before supplying hashes;
    // this sampler cannot authenticate a caller's mutable external population.
    Bridge(std::shared_ptr<const k::StrictCollisionDataset> data,std::string dataset,std::string kernel,
           double water_mass_kg,double H_mass_kg,SitePolicy site)
      :data_(std::move(data)),dataset_(std::move(dataset)),kernel_(std::move(kernel)),M_(water_mass_kg),H_(H_mass_kg),site_(site){
        pc::require(data_&&pc::digest(dataset_)&&pc::digest(kernel_),"unbound generic data/kernel");
        pc::require(site_==SitePolicy::PreserveReactantSites,"unsupported collision-site policy");
        const auto b=k::conditional_water_mass_basis(*data_,M_,H_);
        catalog_={{"electron",0,0,-1,me,0,0},{"H2O",2,1,0,M_,0,0},{"H_1s",1,0,0,H_,0,0},
          {"Hm",1,0,-1,b.Hminus_kg,0,-.75*qe},{"OHm_X_v0",1,1,-1,b.OHminus_X_v0_kg,0,3.27*qe},
          {"H2Op_ground",2,1,1,b.H2Oplus_kg,0,12.621*qe}};
        for(int v=0;v<5;++v){const double x=k::hminus_Ram65_vibrational_eV(*data_,v)*qe;catalog_.push_back({"OH_X_v"+std::to_string(v),1,1,0,b.OH_X_v0_kg+x/(c*c),x,5.1*qe});}
        const auto& I=data_->metadata().at("outer3_fragment_ionization_eV");
        const double ih=I.at("H").get<double>()*qe,ioh=I.at("OH").get<double>()*qe;
        catalog_.push_back({"Hp",1,0,1,H_-me+ih/(c*c),0,ih});
        catalog_.push_back({"OHp_X_ground",1,1,1,b.OH_X_v0_kg-me+ioh/(c*c),0,5.1*qe+ioh});
        for(const auto&s:catalog_)pc::require(states_.emplace(s.id,s).second&&std::isfinite(s.mass_kg)&&s.mass_kg>0,"invalid common state catalog");
        specs_={{"H2O:ion_yield:H2O+","parent_ground_equal_pair_CM_v1",13.5,1000,false,2,{{"electron","electron","H2Op_ground"}}},
          {"H2O:attachment_yield:Hminus","Hminus_Ram65_factorized_window_v1",6.25,6.75,true,3,{}},
          {"H2O:attachment_yield:OHminus","OHminus_ground_isotropic_CM_v1",4.5,12.5,true,2,{{"OHm_X_v0","H_1s"}}},
          {"H2O:electronic_A1B1","A1B1_fixed_loss_dipole_ground_v1",9,977.7,false,4,{{"electron","H_1s","OH_X_v0"}}},
          {"H2O:conditional_outer3:parent","outer3_BEB_Tan_parent_v1",1000,10000,false,4,{{"electron","electron","H2Op_ground"}}},
          {"H2O:conditional_outer3:OHp_H","outer3_BEB_Tan_OHp_H_v1",1000,10000,false,4,{{"electron","electron","OHp_X_ground","H_1s"}}},
          {"H2O:conditional_outer3:Hp_OH","outer3_BEB_Tan_Hp_OH_v1",1000,10000,false,4,{{"electron","electron","Hp","OH_X_v0"}}},
          {"H2O:conditional_IAM_NIST64:elastic","IAM_NIST64_piecewise_mu_elastic_v1",1000,10000,false,2,{{"electron","H2O"}}}};
        for(int v=0;v<5;++v)specs_[1].allowed_product_states.push_back({"Hm","OH_X_v"+std::to_string(v)});
        for(const auto&s:specs_){const auto status=data_->kernel_status(s.channel,s.model);pc::require(status.at("status")=="conditional_available"&&status.at("model").at("data_support_eV")==Json::array({s.low_eV,s.high_eV}),"frozen model support differs");}
    }
    const std::vector<State>& catalog()const{return catalog_;}
    const std::vector<Spec>& specs()const{return specs_;}
    const std::string& dataset_sha256()const{return dataset_;}
    const std::string& kernel_sha256()const{return kernel_;}
    Json OH_family_metadata()const {
        return {{"channel",specs_[1].channel},{"model",specs_[1].model},{"dataset_sha256",dataset_},
          {"allowed_product_multisets",specs_[1].allowed_product_states},
          {"OH_vibrational_relative_intensities",data_->kernel_status(specs_[1].channel,specs_[1].model).at("model").at("OH_vibrational_relative_intensities")},
          {"scope","conditional product-state law per one Hminus event; not independent event cross sections"},{"base_rate_multiplicity",1}};
    }
    Json coverage(double E)const {
        pc::require(std::isfinite(E)&&E>=0,"invalid coverage energy");Json selected=Json::array(),holes=Json::array();
        for(const auto&s:specs_){Json r={{"channel",s.channel},{"model",s.model},{"selected_subset",active(s,E)}};
            try{r["source_sigma_m2"]=data_->observable_value_SI(s.channel,E);r["source_rate_known"]=true;}
            catch(const k::CollisionCoverageError&x){r["source_sigma_m2"]=nullptr;r["source_rate_known"]=false;r["reason"]=x.what();}
            r["event_kernel_in_subset_support"]=active(s,E);selected.push_back(r);}
        Json unmapped={{"channel","H2O:conditional_outer3:unmapped_2a1"},{"applied_as_event",false}};
        try{unmapped["sigma_m2"]=data_->observable_value_SI("H2O:conditional_outer3:unmapped_2a1",E);unmapped["status"]="known_unmapped_rate";}
        catch(const k::CollisionCoverageError&x){unmapped["sigma_m2"]=nullptr;unmapped["status"]="unknown_at_energy";unmapped["reason"]=x.what();}
        holes.push_back(unmapped);
        for(const auto& name:{"core_1a1","multiple_ionization_Auger","fragment_species_continuation","unmapped_inelastic_and_angular_channels"})holes.push_back({{"channel",name},{"sigma_m2",nullptr},{"status","unknown_or_unqualified"},{"applied_as_event",false}});
        Json unused=Json::array();
        for(const auto& raw:data_->metadata().at("channels")){
            const auto id=raw.at("id").get<std::string>();bool selected_id=false;for(const auto&s:specs_)if(s.channel==id)selected_id=true;
            if(selected_id)continue;
            Json row={{"channel",id},{"unit_SI",raw.at("unit_SI")},{"applied_as_event",false},{"rate_sum_eligible",false},{"reason","unmapped or diagnostic observable; no exclusive event added"}};
            try{row["source_value_SI"]=data_->observable_value_SI(id,E);row["source_value_known"]=true;}
            catch(const k::CollisionCoverageError&x){row["source_value_SI"]=nullptr;row["source_value_known"]=false;row["coverage_reason"]=x.what();}
            unused.push_back(row);
        }
        return {{"unapplied_observables",unused},{"energy_eV",E},{"scope","supported_conditional_comparison_subset"},{"complete_mechanism",false},
          {"one_keV_policy","outer3_selected_at_exact_1000_eV"},{"conditional_models",selected},{"missing_rates",holes},
          {"site_policy","retain_scattered_electron_at_primary_and_other_products_at_target; capture_requires_exact_3D_coincidence"}};
    }
    std::vector<pc::Channel> comparison_channels(Scope scope)const {
        pc::require(scope==Scope::CompareSupportedSubset,"complete water mechanism unavailable; see coverage ledger");
        std::vector<pc::Channel> out;for(const auto&s:specs_)out.push_back(channel(s));return out;
    }
    CandidateBatch candidates(const std::vector<PhysicalPacket>& electrons,const std::vector<PhysicalPacket>& water,
          const std::vector<pc::Cell>& cells,Scope scope,std::size_t pair_cap=1000000,std::size_t query_cap=1000000,
          std::size_t row_cap=100000,std::size_t coverage_cap=256)const {
        pc::require(electrons.size()<=100000&&water.size()<=100000,"generic population cap");
        pc::require(coverage_cap<=1024,"coverage evidence hard cap");
        std::size_t bytes=0;constexpr std::size_t byte_cap=64*1024*1024;
        auto validate=[&](const std::vector<PhysicalPacket>& population,bool electron){for(const auto&p:population){
            valid_packet(p,electron);
            for(const auto*text:{&p.packet.id,&p.packet.state}){pc::require(text->size()<=byte_cap-bytes,"generic input identity byte cap");bytes+=text->size();}
        }};validate(electrons,true);validate(water,false);
        auto channels=comparison_channels(scope);std::map<double,Json> ledger;
        // Every pair's energy is captured during the single indexed builder.
        // No second Cartesian traversal and no missing rate renormalization.
        for(auto&ch:channels){auto f=ch.sigma_m2;ch.sigma_m2=[&ledger,this,f,coverage_cap](double E){
            if(!ledger.contains(E)){pc::require(ledger.size()<coverage_cap,"coverage evidence cap");ledger.emplace(E,coverage(E));}return f(E);};}
        std::vector<pc::Packet> e,n;e.reserve(electrons.size());n.reserve(water.size());
        for(const auto&p:electrons){valid_packet(p,true);e.push_back(p.packet);}for(const auto&p:water){valid_packet(p,false);n.push_back(p.packet);}
        CandidateBatch out;pc::CandidateInputLimits caps;caps.projectiles=caps.targets=100000;
        out.rows=pc::candidates(e,n,cells,channels,pair_cap,query_cap,row_cap,caps,&out.work);
        out.coverage=Json::array();for(auto&[E,r]:ledger){(void)E;out.coverage.push_back(std::move(r));}return out;
    }
    Outcome sample(const pc::Candidate& event,const PhysicalPacket&e,const PhysicalPacket&n,const std::vector<double>& u)const {
        valid_packet(e,true);valid_packet(n,false);const auto&s=spec(event.channel_id);
        pc::require(u.size()==s.variates,"incorrect generic random operand count");for(double x:u)pc::require(std::isfinite(x)&&x>=0&&x<1,"random variate outside [0,1)");
        pc::require(!s.captures_primary||same(e.position_m,n.position_m),"capture requires exact 3D coincidence; no qualified nonlocal current/work path");
        const auto expected=pc::candidates({e.packet},{n.packet},{{event.cell,event.volume_m3}},{channel(s)},1,1,1);
        pc::require(expected.size()==1,"no positive supported selected event");const auto&a=expected[0];
        pc::require(event.projectile_id==a.projectile_id&&event.target_id==a.target_id&&event.channel_id==a.channel_id&&event.cell==a.cell&&event.data_sha256==a.data_sha256&&event.kernel_sha256==a.kernel_sha256&&
          event.event_weight==a.event_weight&&event.projectile_before==a.projectile_before&&event.target_before==a.target_before&&event.rate_s==a.rate_s&&event.energy_eV==a.energy_eV&&event.flux_m_s==a.flux_m_s&&event.sigma_m2==a.sigma_m2&&
          event.projectile_remainder.high==a.projectile_remainder.high&&event.projectile_remainder.low==a.projectile_remainder.low&&event.target_remainder.high==a.target_remainder.high&&event.target_remainder.low==a.target_remainder.low,"foreign/stale/inconsistent generic candidate");
        const auto beam=lorentz::incident_direction(e.packet,n.packet,a.energy_eV);const auto b=k::conditional_water_mass_basis(*data_,M_,H_);
        std::vector<std::string> names;std::vector<Vec> momenta;std::vector<double> masses;double sigma=0;int level=-1;
        if(s.model=="parent_ground_equal_pair_CM_v1"){
            auto p=k::conditional_h2o_parent_products(*data_,s.model,a.energy_eV,a.event_weight,M_,beam,isotropic(u[0],u[1]));names={"electron","electron","H2Op_ground"};momenta={p.primary_momentum_kg_m_s,p.secondary_momentum_kg_m_s,p.ion_momentum_kg_m_s};masses={me,me,p.ion_mass_kg};sigma=p.cross_section_m2;
        }else if(s.model=="Hminus_Ram65_factorized_window_v1"){
            auto p=k::conditional_hminus_Ram65_products(*data_,s.model,a.energy_eV,a.event_weight,M_,b.Hminus_kg,beam,u[0],u[1],u[2]);level=p.OH_vibrational_level;names={"Hm","OH_X_v"+std::to_string(level)};momenta={p.Hminus_momentum_kg_m_s,p.OH_momentum_kg_m_s};masses={b.Hminus_kg,p.OH_state_mass_kg};sigma=p.cross_section_m2;
        }else if(s.model=="OHminus_ground_isotropic_CM_v1"){
            auto p=k::conditional_ohminus_ground_products(*data_,s.model,a.energy_eV,a.event_weight,M_,H_,beam,isotropic(u[0],u[1]));names={"OHm_X_v0","H_1s"};momenta={p.anion_momentum_kg_m_s,p.neutral_momentum_kg_m_s};masses={p.anion_state_mass_kg,p.neutral_state_mass_kg};sigma=p.cross_section_m2;
        }else if(s.model=="A1B1_fixed_loss_dipole_ground_v1"){
            auto p=k::conditional_A1B1_dissociation_products(*data_,s.model,a.energy_eV,a.event_weight,M_,H_,beam,u[0],u[1],isotropic(u[2],u[3]));names={"electron","H_1s","OH_X_v0"};momenta={p.electron_momentum_kg_m_s,p.H_momentum_kg_m_s,p.OH_momentum_kg_m_s};masses={me,H_,p.OH_state_mass_kg};sigma=p.cross_section_m2;
        }else if(s.model=="IAM_NIST64_piecewise_mu_elastic_v1"){
            auto p=k::conditional_IAM_elastic_products(*data_,s.model,a.energy_eV,a.event_weight,M_,beam,u[0],u[1]);names={"electron","H2O"};momenta={p.electron_momentum_kg_m_s,p.H2O_momentum_kg_m_s};masses={me,p.H2O_mass_kg};sigma=p.cross_section_m2;
        }else{
            auto p=k::conditional_outer3_ionization_products(*data_,s.model,a.energy_eV,a.event_weight,M_,H_,beam,isotropic(u[0],u[1]),isotropic(u[2],u[3]));for(int i=0;i<p.product_count;++i){names.push_back(p.state_id[i]=="e_g"?"electron":p.state_id[i]);momenta.push_back(p.momentum_kg_m_s[i]);masses.push_back(p.mass_kg[i]);}sigma=p.cross_section_m2;
        }
        pc::require(std::find(s.allowed_product_states.begin(),s.allowed_product_states.end(),names)!=s.allowed_product_states.end(),"product state family differs");
        pc::require(std::abs(sigma-a.sigma_m2)<=2e-13*a.sigma_m2,"kernel/rate sigma mismatch");
        Outcome out;out.channel=s.channel;out.model=s.model;out.data_sha256=dataset_;out.kernel_sha256=kernel_;out.event_weight=a.event_weight;out.energy_eV=a.energy_eV;out.sigma_m2=a.sigma_m2;out.captures_primary=s.captures_primary;out.OH_level=level;out.random_variates=u.size();out.projectile_remainder=a.projectile_remainder;out.target_remainder=a.target_remainder;
        int nh=0,no=0,nq=0;long double before=lorentz::kinetic(me,e.packet.proper_velocity_m_s)+lorentz::kinetic(M_,n.packet.proper_velocity_m_s),after=0,pscale=0;std::array<long double,3> dp{};
        for(int d=0;d<3;++d){dp[d]=-static_cast<long double>(me)*e.packet.proper_velocity_m_s[d]-static_cast<long double>(M_)*n.packet.proper_velocity_m_s[d];pscale+=std::abs(static_cast<long double>(me)*e.packet.proper_velocity_m_s[d])+std::abs(static_cast<long double>(M_)*n.packet.proper_velocity_m_s[d]);}
        for(std::size_t i=0;i<names.size();++i){const auto& state=states_.at(names[i]);pc::require(std::abs(masses[i]-state.mass_kg)<=8*std::numeric_limits<double>::epsilon()*state.mass_kg,"kernel/common catalog mass mismatch");Vec proper{};for(int d=0;d<3;++d)proper[d]=momenta[i][d]/state.mass_kg;proper=lorentz::to_lab(proper,n.packet.proper_velocity_m_s);out.products.push_back({state,proper,!s.captures_primary&&i==0,a.event_weight});nh+=state.H;no+=state.O;nq+=state.q;after+=state.reference_J+state.internal_J+lorentz::kinetic(state.mass_kg,proper);for(int d=0;d<3;++d){const long double p=static_cast<long double>(state.mass_kg)*proper[d];dp[d]+=p;pscale+=std::abs(p);}}
        pc::require(nh==2&&no==1&&nq==-1,"unweighted nuclei/charge closure");
        const long double factor=256*std::numeric_limits<double>::epsilon(),ebound=factor*(std::abs(after)+std::abs(before)),pbound=factor*pscale;
        pc::require(std::abs(after-before)<=ebound,"unweighted generic lab energy closure");
        pc::require(ebound<=lorentz::maximum_lab_roundoff_fraction_of_reaction_energy*a.energy_eV*qe,"lab representation cannot resolve reaction energy");
        for(auto x:dp)pc::require(std::abs(x)<=pbound,"unweighted generic lab momentum closure");
        out.weighted_energy_defect_J=static_cast<double>((after-before)*a.event_weight);out.weighted_energy_bound_J=static_cast<double>(ebound*a.event_weight);out.weighted_momentum_bound=static_cast<double>(pbound*a.event_weight);
        pc::require(std::isfinite(out.weighted_energy_defect_J)&&std::isfinite(out.weighted_energy_bound_J)&&out.weighted_energy_bound_J>0&&std::isfinite(out.weighted_momentum_bound)&&out.weighted_momentum_bound>0,"weighted receipt overflow/positive underflow");
        for(int d=0;d<3;++d){out.weighted_momentum_defect[d]=static_cast<double>(dp[d]*a.event_weight);pc::require(std::isfinite(out.weighted_momentum_defect[d]),"weighted momentum overflow");}
        return out;
    }
    Outcome sample(const pc::Candidate&a,const PhysicalPacket&e,const PhysicalPacket&n,std::mt19937_64&rng)const {
        auto staged=rng;std::vector<double> u;for(std::size_t i=0;i<spec(a.channel_id).variates;++i)u.push_back(uniform(staged));
        auto out=sample(a,e,n,u);rng=staged;return out; // Successful private RNG step only.
    }
};
} // namespace raikou::generic_water
