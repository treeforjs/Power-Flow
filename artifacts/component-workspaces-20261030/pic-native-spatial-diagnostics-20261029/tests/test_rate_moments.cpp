#include "rate_moments.hpp"
#include "native_diagnostic_view.hpp"
#include <iostream>
using namespace native_diagnostics;
void check(bool x,const char*m){if(!x)throw std::runtime_error(m);}
template<class F>void rejects(F f){bool caught=false;try{f();}catch(const std::exception&){caught=true;}check(caught,"Expected closed rejection");}
AssignedRateEntry point(std::uint64_t id,const char*state,double w,double mass,double ux){ArenaRateEntry e;e.packet_id=std::string("packet/")+std::to_string(id);e.state=state;e.native.id=id;e.native.weight=w;e.native.mass_kg=mass;e.native.charge_C=std::string(state)=="electron"?-1.602176634e-19:0;e.native.proper_velocity_m_s={ux,0,0};return {e,1};}
int main(){try{
 std::string state(64,'a'),geometry(64,'b');RateQuery q{"H2O/test","electron","H2O",std::string(64,'c'),std::string(64,'d'),"manufactured_constant_sigma","exclusive_event",[](double){return SigmaCoverage{true,2e-20,"",true};}};
 const double me=9.1093837139e-31;std::vector<AssignedRateEntry> e{point(1,"electron",2,me,1e6),point(2,"H2O",3,3e-26,0)};std::vector<ReactionCell>cells{{1,2e-6}};
 auto j=rate_moments(e,cells,{q},state,geometry,10);const double v=1e6/std::sqrt(1+1e12/(pc::c*pc::c)),expected=6*2e-20*v/2e-6,actual=j["rows"][0]["supported_partial_physical_number_rate_s"]["value"].get<double>();check(std::abs(actual-expected)<=32*std::numeric_limits<double>::epsilon()*expected,"Actual WeWn sigma relative-flux /physical volume");
 e.push_back(point(3,"electron",1,me,2e6));q.evaluate=[](double E){return E>5?SigmaCoverage{true,2e-20,"",false}:SigmaCoverage{false,0,"explicit unsupported low-energy interval",false};};j=rate_moments(e,cells,{q},state,geometry,10);const auto&r=j["rows"][0];check(std::abs(r["unsupported_pair_fraction"]["value"].get<double>()-2./3)<1e-15,"Unsupported pair fraction");check(r["unsupported_projectile_weight"]["value"].get<double>()==2,"Unsupported positive electron weight");check(r["missing_exclusive_event_kernel_pair_weight"]["value"].get<double>()==3,"Observable support is not event-kernel support");
 const double v2=2e6/std::sqrt(1+4e12/(pc::c*pc::c));check(std::abs(r["supported_partial_sum_We_target_averaged_sigma_flux_m3_s"]["value"].get<double>()-2e-20*v2)<=64*std::numeric_limits<double>::epsilon()*2e-20*v2,"No unsupported-population renormalization");
 rejects([&]{rate_moments(e,cells,{q},state,geometry,1);});auto duplicate=e;duplicate.push_back(e[0]);rejects([&]{rate_moments(duplicate,cells,{q},state,geometry,10);});auto corrupt=e;corrupt[0].entry.native.proper_velocity_m_s[0]=NAN;rejects([&]{rate_moments(corrupt,cells,{q},state,geometry,10);});auto same=q;same.target_state=same.projectile_state;rejects([&]{rate_moments(e,cells,{same},state,geometry,10);});q.evaluate=[](double){return SigmaCoverage{false,0,"",false};};rejects([&]{rate_moments(e,cells,{q},state,geometry,10);});
 std::vector<AssignedRateEntry>only_e{e[0]};j=rate_moments(only_e,cells,{q},state,geometry,0);check(j["actual_pair_count"].get<int>()==0&&!j["rows"][0]["complete_sigma_coverage"].get<bool>(),"Missing target is explicit; no target velocity invented");
 std::cout<<"Actual-weight relative-flux moments, unsupported fractions, event distinction and caps passed; no arena/PIC construction\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
