#include "raikou/exact_particle_count.h"
#include "dyadic.hpp"
#include <fstream>
#include <iostream>
using C=raikou::mfem_backend::ExactParticleCount;
using J=nlohmann::json;
using metal_transfer::D;
struct Two {double hi=0,lo=0;};
double value(const J&j){
 if(j.is_number())return j.get<double>();
 const std::string s=j.get<std::string>();
 if(s=="max")return std::numeric_limits<double>::max();
 if(s=="tiny")return std::numeric_limits<double>::denorm_min();
 if(s=="inf")return INFINITY;
 if(s=="-inf")return -INFINITY;
 if(s=="nan")return NAN;
 if(s=="-zero")return -0.;
 throw std::invalid_argument("unknown constant");
}
std::string units(const C&c){auto n=c.copy_units();char*p=BN_bn2hex(n.get());if(!p)throw std::bad_alloc();std::string s(p);OPENSSL_free(p);return s;}
J snapshot(const C&c,double atom,double R,double K){
 J j={{"units_hex",units(c)},{"carriers",c.positive_carriers()},{"zero",c.zero()}};
 try{j["rounded"]=c.rounded_for_rate();}catch(const std::exception&e){j["rounded_refusal"]=e.what();}
 try{auto out=c.require_two_word<Two>();j["two_word"]={out.hi,out.lo};j["two_word_reimport_equal"]=(C::from_two_word(out)==c);}catch(const std::exception&e){j["two_word_refusal"]=e.what();}
 D count;for(double w:c.positive_carriers())count=count+D(w);
 j["material_exact_atoms"]=count.json();j["material_mass_kg"]=(count*D(atom)).json();j["material_packet_energy_J"]=(count*(D(R)+D(K))).json();
 auto copy=c;auto external=c.copy_units();BN_zero(external.get());j["copy_units_cannot_mutate_value"]=(copy==c&&units(copy)==j.at("units_hex").get<std::string>());
 return j;
}
int main(int argc,char**argv){try{
 if(argc!=3)throw std::invalid_argument("input and output");
 J input;{std::ifstream f(argv[1]);f>>input;}J rows=J::array();double atom=input.at("atom_mass"),R=input.at("reference"),K=input.at("manufactured_K");
 for(const auto&test:input.at("cases")){
  J row={{"id",test.at("id")},{"input",test}};
  try{
   C c;
   if(test.at("kind")=="two")c=C::from_two_word(Two{value(test.at("hi")),value(test.at("lo"))});
   else{std::vector<double>w;for(auto&x:test.at("weights"))w.push_back(value(x));c=C::from_positive_carriers(w);}
   row["admitted"]=true;row["initial"]=snapshot(c,atom,R,K);row["operations"]=J::array();
   for(const auto&op:test.at("operations")){
    const auto before=c;const std::string old=units(c);J attempt={{"op",op},{"before",snapshot(c,atom,R,K)}};
    try{const double w=value(op.at("amount"));c=op.at("kind")=="debit"?c.debit(w):c.credited(w);attempt["admitted"]=true;attempt["after"]=snapshot(c,atom,R,K);}catch(const std::exception&e){attempt["admitted"]=false;attempt["reason"]=e.what();attempt["unchanged"]=(c==before&&units(c)==old);}
    row["operations"].push_back(attempt);
   }
  }catch(const std::exception&e){row["admitted"]=false;row["reason"]=e.what();}
  rows.push_back(row);
 }
 std::ofstream f(argv[2]);f<<J({{"status","completed"},{"cases",rows},{"maximum_bits",C::maximum_bits},{"maximum_carriers",C::maximum_carriers},{"atom_mass",atom},{"reference",R},{"manufactured_K",K}}).dump(2)<<'\n';if(!f)throw std::runtime_error("output write");return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
