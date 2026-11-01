#include "pf_gas_v1/gas_owner.hpp"
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
using namespace pf_gas_v1;
Configuration config(){ Configuration c;c.feedback=true;c.species={{"test-neutral",2.99e-26,0,true}};c.grid.low={0,0,0};c.grid.high={2,1,1};c.grid.n={2,1,1};c.grid.volume=1;c.processes={{0,2,pfg::kNone,1}};return c; }
void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
int main(int argc,char**argv){try{
 need(argc==3,"usage: baseline|fixed|refuse checkpoint");const std::string mode=argv[1];auto c=config();
 if(mode=="refuse"){
  std::ifstream in(argv[2],std::ios::binary);need(bool(in),"open legacy checkpoint");std::vector<std::uint8_t>b{std::istreambuf_iterator<char>(in),{}};
  bool rejected=false;try{auto restored=Owner::restore(b,c,"");}catch(const Refusal&e){rejected=e.kind()==Refusal::Kind::restore && std::string(e.what()).find("configuration identity differs")!=std::string::npos;}
  need(rejected,"legacy checkpoint must refuse on actual configuration identity");std::cout<<"old checkpoint refused by configuration identity\n";return 0;
 }
 need(mode=="baseline"||mode=="fixed","unknown mode");Owner o(c);o.begin_step(0);RecoilRecord r{};r.stage=1;auto&e=r.ev;e.weight=2;e.n_out=2;e.u_in[0]=3000;e.ut[0]=-4000;e.ut[1]=2500;e.ut[2]=1500;e.out_dest[0]=pfg::kDestSurvivor;e.out_dest[1]=pfg::kDestReservoir;for(unsigned a=0;a<3;++a){e.out_u[0][a]=e.u_in[a];e.out_u[1][a]=e.ut[a];}
 o.take(0,1,{&r,1},{});FlowRecord f{};f.kind=pfg::kFlowRelease;f.native_id=100;f.weight=2;for(unsigned a=0;a<3;++a)f.u[a]=e.ut[a];o.take(0,2,{},{&f,1});o.commit(o.candidate_sha256());
 Exact zero,old;old.add_product(2,pfg::kinetic(c.species[0].mass_kg,e.ut));need(o.global_digits(pfg::G_RES_E,false)==(mode=="baseline"?old.normalised():zero.normalised()),"actual full-release energy mirror");
 for(unsigned a=0;a<3;++a){Exact p;p.add_triple(2,c.species[0].mass_kg,e.ut[a]);need(o.global_digits(pfg::G_RES_P+a,false)==(mode=="baseline"?p.normalised():zero.normalised()),"actual full-release momentum mirror");}
 auto bytes=o.encode();std::ofstream out(argv[2],std::ios::binary);out.write(reinterpret_cast<const char*>(bytes.data()),std::streamsize(bytes.size()));need(bool(out),"checkpoint write");std::cout<<mode<<" mirror proved; configuration "<<configuration_sha256(c)<<"\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
