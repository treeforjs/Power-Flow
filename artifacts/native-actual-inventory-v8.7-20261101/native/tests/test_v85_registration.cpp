#include "chem_v8_fixture.hpp"
#include <iostream>
int main(int argc,char **argv) {
  try {
    const bool baseline=argc==2 && std::string(argv[1])=="baseline";
    auto o=chemfix::options();auto c=chemfix::chemistry_configuration(o);
    chemfix::Population pop;pop.carriers=10;
    auto packets=chemfix::population(pop);auto A=chemfix::make_arena(o,c,packets,1000);
    unsigned refused=0,accepted=0;
    for(int kind=0;kind<2;++kind) {
      auto states=A.states;
      auto it=std::find_if(states.begin(),states.end(),[](const auto&s){return s.id=="pfc-v3/H2O";});
      if(it==states.end())throw std::runtime_error("missing water fixture");
      if(kind==0)it->internal_J=1e-19;else it->reference_J+=1e-19;
      try{auto wrong=raikou::mfem_backend::PicArena::create_full_time(o,states,A.channels,packets,A.g,A.cfg);++accepted;}
      catch(const raikou::mfem_backend::FullTimeChemistryRefusal&e) {
        if(e.kind()!=raikou::mfem_backend::FullTimeChemistryRefusal::Kind::configuration)throw;
        ++refused;
      }
    }
    if((baseline && accepted!=2)||(!baseline && refused!=2))throw std::runtime_error("energy registration outcome");
    std::cout<<(baseline?"BASELINE accepted ":"FIXED refused ")<<(baseline?accepted:refused)<<" mismatched species energy conventions; valid genesis control passed\n";
    return 0;
  }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}
}
