#include "generic_water_bridge.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
namespace g=raikou::generic_water;namespace pc=raikou::pic_collisions;using J=g::Json;
std::size_t checks=0;
void check(bool yes,const char*s){++checks;if(!yes)throw std::runtime_error(s);}
template<class F>void refuses(F f,const char*s){bool threw=false;try{f();}catch(const std::exception&){threw=true;}check(threw,s);}
std::string state(const std::mt19937_64&r){std::ostringstream o;o<<r;return o.str();}
constexpr double amu=1.66053906660e-27,M=18.01056468403*amu,H=1.00782503223*amu;
std::pair<g::PhysicalPacket,g::PhysicalPacket> pair(double E,double boost=0,double ew=2,double nw=7){
    const double x=E*pc::qe/(g::me*pc::c*pc::c),u=pc::c*std::sqrt(x*(x+2));
    const g::Vec target_u{boost*pc::c,0,0};
    return {{{"e","electron",7,ew,g::me,raikou::pic_water::to_lab({.6*u,0,.8*u},target_u)},{.004,.004,.0125}},{{"n","H2O",7,nw,M,target_u},{.004,.004,.0125}}};
}
J serialize(const g::Outcome&o,const g::PhysicalPacket&e,const g::PhysicalPacket&n){J p=J::array();for(const auto&x:o.products)p.push_back({{"state",x.state.id},{"H",x.state.H},{"O",x.state.O},{"q",x.state.q},{"mass",x.state.mass_kg},{"internal",x.state.internal_J},{"reference",x.state.reference_J},{"u",x.proper_velocity_m_s},{"weight",x.weight},{"at_primary",x.at_primary}});return {{"channel",o.channel},{"model",o.model},{"weight",o.event_weight},{"capture",o.captures_primary},{"OH_level",o.OH_level},{"incoming_e_mass",e.packet.mass_kg},{"incoming_n_mass",n.packet.mass_kg},{"incoming_e_u",e.packet.proper_velocity_m_s},{"incoming_n_u",n.packet.proper_velocity_m_s},{"products",p},{"energy_bound",o.weighted_energy_bound_J},{"momentum_bound",o.weighted_momentum_bound}};}
int main(int argc,char**argv){try{
    if(argc!=4)return 2;
    const auto data=std::make_shared<const raikou::kinetics::StrictCollisionDataset>(raikou::kinetics::StrictCollisionDataset::read(argv[1]));J binding;{std::ifstream f(argv[2]);f>>binding;}
    g::Bridge b(data,binding.at("dataset_sha256"),binding.at("kernel_manifest_sha256"),M,H,g::SitePolicy::PreserveReactantSites);
    const std::string out=argv[3];J records=J::array(),coverage=J::array();std::set<int> levels;
    check(b.specs().size()==8&&b.specs()[1].allowed_product_states.size()==5,"one conditional family and eight models");
    check(b.OH_family_metadata().at("OH_vibrational_relative_intensities").size()==5&&b.OH_family_metadata().at("base_rate_multiplicity")==1,"conditional source family metadata");
    refuses([&]{b.comparison_channels(g::Scope::RequireCompleteMechanism);},"full mechanism accidentally admitted");
    std::size_t model_index=0;
    for(const auto&s:b.specs()){
        const double E=model_index==0?100:model_index==1||model_index==2?6.5:model_index==3?20:8084.470604111418;
        for(double boost:{0.,.0001})for(int seed=0;seed<32;++seed){
            auto[e,n]=pair(E,boost);auto batch=b.candidates({e},{n},{{7,1e-12}},g::Scope::CompareSupportedSubset);
            check(batch.work.channel_queries==8,"all eight channels queried together");
            const auto it=std::find_if(batch.rows.begin(),batch.rows.end(),[&](auto&a){return a.channel_id==s.channel;});check(it!=batch.rows.end(),"supported candidate absent");
            std::mt19937_64 rng(17+seed),repeat(17+seed);const auto result=b.sample(*it,e,n,rng),again=b.sample(*it,e,n,repeat);
            check(state(rng)==state(repeat)&&serialize(result,e,n)==serialize(again,e,n),"product/RNG determinism");
            check(result.captures_primary==s.captures_primary&&result.products.size()>=2&&result.products.size()<=4,"variable product contract");
            if(result.OH_level>=0)levels.insert(result.OH_level);
            records.push_back(serialize(result,e,n));
        }
        ++model_index;
    }
    check(levels.size()==5,"conditional OHv family coverage");
    for(double E:{0.,6.5,20.,100.,1000.,8084.470604111418,12000.})coverage.push_back(b.coverage(E));
    const auto high=b.coverage(8084.470604111418);check(high.at("missing_rates")[0].at("sigma_m2").get<double>()>0,"positive 2a1 hidden");
    check(!high.at("complete_mechanism").get<bool>(),"incomplete coverage mislabeled");
    const auto channels=b.comparison_channels(g::Scope::CompareSupportedSubset);
    check(channels[0].sigma_m2(1000)==0&&channels[4].sigma_m2(1000)>0,"1keV overlap double counted");
    check(channels[1].sigma_m2(6.5)==data->observable_value_SI(channels[1].id,6.5),"conditional OHv duplicated base rate");
    // One actual competing clock, including terminal no-event draws. Model
    // frequencies compare to their actual rate fractions, not equal weights.
    J ensembles=J::array();for(double E:{6.5,8084.470604111418}){
        auto[e,n]=pair(E);auto batch=b.candidates({e},{n},{{7,1e-12}},g::Scope::CompareSupportedSubset);long double rate=0;for(auto&a:batch.rows)rate+=a.rate_s;
        const double dt=static_cast<double>(1/rate);std::map<std::string,std::size_t> counts;std::size_t events=0;constexpr std::size_t N=8192;
        for(std::size_t seed=0;seed<N;++seed){pc::Scheduler scheduler(seed+109);pc::Clock clock{0,0,0,dt,std::string(64,'a'),std::string(64,'b'),std::string(64,'c'),0,0};const auto p=scheduler.prepare_exposure({clock,dt,dt},batch.rows);if(p.event()){++events;++counts[p.event()->channel_id];}check(scheduler.try_publish(p),"competing clock publish");}
        const long double probability=1-std::exp(-rate*dt);check(std::abs(static_cast<long double>(events)-N*probability)<6*std::sqrt(N*probability*(1-probability))+4,"actual total-rate clock frequency");
        J rows=J::array();for(const auto&a:batch.rows){const long double p=probability*a.rate_s/rate;check(std::abs(static_cast<long double>(counts[a.channel_id])-N*p)<6*std::sqrt(N*p*(1-p))+4,"actual competing channel frequency");const long double first=static_cast<long double>(a.event_weight)*a.rate_s,expected=static_cast<long double>(e.packet.weight)*n.packet.weight*a.sigma_m2*a.flux_m_s/a.volume_m3;check(std::abs(first-expected)<=8*std::numeric_limits<double>::epsilon()*expected,"unequal-weight first moment");rows.push_back({{"channel",a.channel_id},{"rate",a.rate_s},{"events",counts[a.channel_id]},{"expected",static_cast<double>(N*p)}});}
        ensembles.push_back({{"energy_eV",E},{"clocks",N},{"events",events},{"rows",rows}});
    }
    auto[e,n]=pair(6.5);auto batch=b.candidates({e},{n},{{7,1e-12}},g::Scope::CompareSupportedSubset);auto candidate=batch.rows.front();std::mt19937_64 rng(73);const auto before=state(rng);
    auto distant=n;distant.position_m[2]+=.001;refuses([&]{b.sample(candidate,e,distant,rng);},"nonlocal capture relocated");check(state(rng)==before,"refused capture RNG mutation");
    for(int kind=0;kind<5;++kind){auto a=candidate;if(kind==0)a.rate_s=std::nextafter(a.rate_s,INFINITY);if(kind==1)a.sigma_m2=std::nextafter(a.sigma_m2,INFINITY);if(kind==2)a.data_sha256[0]='f';if(kind==3)a.projectile_remainder.low=1;if(kind==4)a.event_weight*=.5;refuses([&]{b.sample(a,e,n,rng);},"changed candidate accepted");check(state(rng)==before,"refusal RNG mutation");}
    refuses([&]{b.candidates({e},{n},{{7,1e-12}},g::Scope::CompareSupportedSubset,1,7);},"channel query cap bypass");
    refuses([&]{b.candidates({e},{n},{{7,1e-12}},g::Scope::CompareSupportedSubset,1,8,100,0);},"coverage cap bypass");
    auto fragment=n;fragment.packet.state="OH_X_v0";refuses([&]{b.candidates({e},{fragment},{{7,1e-12}},g::Scope::CompareSupportedSubset);},"fragment continuation silently ignored");
    auto long_id=e;long_id.packet.id=std::string(4097,'x');refuses([&]{b.candidates({long_id},{n},{{7,1e-12}},g::Scope::CompareSupportedSubset);},"packet identity cap bypass");
    refuses([&]{b.comparison_channels(static_cast<g::Scope>(99));},"invalid scope admitted");
    {auto[high_e,high_n]=pair(12000);const auto unsupported=b.candidates({high_e},{high_n},{{7,1e-12}},g::Scope::CompareSupportedSubset);check(unsupported.rows.empty()&&unsupported.coverage.size()==1&&!unsupported.complete_mechanism,"outside-subset evidence hidden");}
    {auto staged=rng;const auto first=b.sample(candidate,e,n,staged);(void)first;std::istringstream in(state(staged));std::mt19937_64 restored;in>>restored;const auto left=b.sample(candidate,e,n,staged),right=b.sample(candidate,e,n,restored);check(serialize(left,e,n)==serialize(right,e,n)&&state(staged)==state(restored),"product RNG restart parity");}
    std::ofstream(out+"/outcome-records.json")<<records.dump(2)<<'\n';std::ofstream(out+"/coverage.json")<<coverage.dump(2)<<'\n';std::ofstream(out+"/clock-ensembles.json")<<ensembles.dump(2)<<'\n';
    std::ofstream(out+"/bridge-results.json")<<J{{"status","passed"},{"assertions",checks},{"actual_kernel_outcomes",records.size()},{"supported_models",8},{"OH_levels_observed",levels},{"competing_clocks",16384},{"complete_mechanism",false},{"physical_time_advanced",false}}.dump(2)<<'\n';std::cout<<checks<<" assertions; "<<records.size()<<" actual conditional outcomes\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
