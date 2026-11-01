#include "chem_v8_fixture.hpp"
#include "chem_v3_native.hpp"
#include <iostream>
#include <stdexcept>
namespace c=raikou::mfem_backend::chem8;
void need(bool ok,const char*m){if(!ok)throw std::runtime_error(m);}
int main(){try{
  auto cfg=chemfix::chemistry_configuration(chemfix::options(45));
  auto setup=c::make_setup(cfg,chemfix::options(45));
  auto state=c::initial_state(*setup);
  const auto& b=*setup->reg->blob;
  std::uint32_t species=0;
  while(species<b.n_species && b.species_names[species]!="H2O") ++species;
  need(species<b.n_species,"water species");
  c::Key key{0,species}; auto& R=state->reservoirs[key];
  const double m=b.view.species[species].mass_eV;
  for(int i=0;i<2;++i){
    pfc::ReservoirEntry e{1.0,{1000.0*i,0,0},{.003+i*.00001,.004,0}};
    R.parcels.push_back(e);const auto u=pfc::v3(e.u[0],e.u[1],e.u[2]);const auto p=pfc::p_from_u(m,u);
    R.count.add(e.weight);R.momentum[0].add_product(e.weight,p.x);R.momentum[1].add_product(e.weight,p.y);
    R.momentum[2].add_product(e.weight,p.z);R.kinetic.add_product(e.weight,pfc::ke_from_u(m,u));
  }
  const auto bytes=c::encode_state(*state);
  auto restored=c::decode_state(*setup,bytes);
  need(c::state_digest(*restored)==c::state_digest(*state),"parcel codec identity");
  need(c::encode_state(*restored)==bytes,"parcel codec canonical bytes");
  std::map<pfc::InventoryKey,pfc::Reservoir> a{{{0,species},R}}, d{{{0,species},restored->reservoirs.at(key)}};
  const auto x=pfc::release_reservoirs(b,a,.75,1,0);const auto y=pfc::release_reservoirs(b,d,.75,999,999);
  need(x.size()==4 && y.size()==x.size(),"partial deterministic parcel split");
  for(std::size_t i=0;i<x.size();++i){
    need(x[i].weight==y[i].weight && x[i].u[0]==y[i].u[0] && x[i].x[0]==y[i].x[0],"restored release phase space");
  }
  need(a.begin()->second.count.is_zero() && a.begin()->second.kinetic.is_zero(),"empty reservoir energy");
  auto corrupt=c::initial_state(*setup);corrupt->reservoirs[key]=R;corrupt->reservoirs[key].kinetic.add(-1.0);
  bool refused=false;try{(void)c::decode_state(*setup,c::encode_state(*corrupt));}catch(const std::exception&){refused=true;}
  need(refused,"nonrealizable checkpoint stock refused");
  auto old=bytes;const std::string tag="native-v8.2/chemistry-state-codec-v3-parcels-support";
  const auto pos=old.find(tag);need(pos!=std::string::npos,"codec tag");old[pos]='X';
  refused=false;try{(void)c::decode_state(*setup,old);}catch(const std::exception&){refused=true;}
  need(refused,"incompatible codec refused");
  std::cout<<"audit parcels: persisted thermal spread, partial release, exact stock, corrupt codec refusal passed\n";
  return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
