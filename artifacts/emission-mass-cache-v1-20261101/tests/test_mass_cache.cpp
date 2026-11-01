#include "fast_emission.hpp"
#include "distribution_bytes.hpp"
#include "registry_mass_fixture.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <typeinfo>
namespace se=surface_emission;namespace fast=electrode_heat_v6_1::emission;
unsigned checks=0;
void need(bool x,const char*s){if(!x)throw std::runtime_error(s);++checks;}
struct Outcome{std::string bytes,type,error;};
template<class F>Outcome result(F f){try{return{distribution_bytes(f()),{}, {}};}catch(const std::exception&e){return{{},typeid(e).name(),e.what()};}}
void compare(double T,double F,const se::Config&c){auto a=result([&]{return fast::evaluate(T,F,c);});auto b=result([&]{return se::evaluate(T,F,c);});need(a.bytes==b.bytes && a.type==b.type && a.error==b.error,"actual scalar byte/refusal parity");}
int main(){try{
 se::Config legacy;legacy.work_function_eV=4.3;legacy.energy_panels=32;legacy.action_nodes=32;
 auto declared=legacy;declared.electron_mass_kg=std::nextafter(se::me,std::numeric_limits<double>::infinity());declared.electron_mass_identity="diagnostic/declared-nextafter-mass";
 auto doubled=declared;doubled.electron_mass_kg=2*se::me;doubled.electron_mass_identity="diagnostic/2me";
 auto provenance=legacy;provenance.electron_mass_identity="diagnostic/same-mass-distinct-policy";
 auto actual_registry=registry_mass_fixture(legacy);
 for(const auto&c:{legacy,declared,doubled,provenance,actual_registry})for(double T:{0.,1.,300.,1000.,3000.})for(double F:{-1e8,-0.,0.,1.,300.,2e6,3e7})compare(T,F,c);
 // A warm retarding distribution may be reused only within its complete key.
 fast::reset_statistics();auto a=fast::evaluate(877.,-3.,declared);auto b=fast::evaluate(877.,-4.,declared);
 auto st=fast::statistics();need(st.retarding_sealed==1 && st.retarding_reused==1,"declared mass reaches bounded retarding cache");
 auto identity=declared;identity.electron_mass_identity+="/other";compare(877.,-4.,identity);
 need(fast::statistics().retarding_sealed==2,"provenance change cannot reuse distribution metadata");
 auto mass=declared;mass.electron_mass_kg=std::nextafter(mass.electron_mass_kg,std::numeric_limits<double>::infinity());compare(877.,-4.,mass);
 need(fast::statistics().retarding_sealed==3,"distinct mass bits cannot reuse spectrum");
 auto swapped=a;swapped.config.electron_mass_kg=mass.electron_mass_kg;need(!fast::same_distribution(a,swapped),"mass is part of distribution equality");swapped=a;swapped.config.electron_mass_identity+="x";need(!fast::same_distribution(a,swapped),"provenance is part of distribution equality");
 fast::reset_statistics();compare(733.,2e6,declared);need(fast::statistics().extracting_fast==1,"declared mass exact extracting path is enabled after self-check");
 compare(733.,2e6,actual_registry);need(fast::statistics().extracting_fast==2,"actual qualified registry mass reaches exact extracting path");
 fast::set_fast_path(false);compare(877.,-4.,declared);compare(733.,2e6,declared);need(fast::statistics().sealed_calls>=2,"explicit fallback remains available");fast::set_fast_path(true);
 for(double m:{0.,-1.,std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity()}){auto c=declared;c.electron_mass_kg=m;compare(300.,0.,c);compare(-1.,1e10,c);}
 auto noid=declared;noid.electron_mass_identity.clear();compare(300.,0.,noid);auto longid=declared;longid.electron_mass_identity=std::string(4097,'x');compare(300.,0.,longid);
 for(double T:{-1.,.5,3001.,std::numeric_limits<double>::quiet_NaN()})compare(T,0.,declared);
 for(double F:{.5,5e9+1,std::numeric_limits<double>::quiet_NaN()})compare(300.,F,declared);
 // Default quadrature is independently exercised, beyond the bounded matrix.
 se::Config defaults;defaults.work_function_eV=4.3;defaults.electron_mass_kg=declared.electron_mass_kg;defaults.electron_mass_identity=declared.electron_mass_identity;compare(1000.,2e6,defaults);compare(1000.,-1.,defaults);
 se::Config default_legacy;default_legacy.work_function_eV=4.3;compare(1000.,2e6,default_legacy);
 auto default_registry=registry_mass_fixture(default_legacy);fast::reset_statistics();compare(1000.,2e6,default_registry);need(fast::statistics().extracting_fast==1,"actual registry default quadrature uses qualified exact fast path");
 std::cout<<"PASS "<<checks<<" independent byte/refusal/cache controls\n";
}catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<" after "<<checks<<" controls\n";return 1;}}
