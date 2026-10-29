#include "raikou/mfem_gap_pic.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <cmath>
#include <filesystem>
using namespace raikou::mfem_backend;
struct State {GapPicOptions o;GapPicDiagnostics d;std::map<std::uint64_t,GapParticle> p;};
State read(const char* directory) {
 std::filesystem::path dir(directory);std::ifstream meta(dir/"restart.meta");std::string name;meta>>name;
 std::ifstream f(dir/name,std::ios::binary);State s;std::uint64_t magic,id,n;
 const auto get=[&](auto& x){f.read(reinterpret_cast<char*>(&x),sizeof(x));};
 get(magic);get(s.o);get(s.d);get(id);get(n);
 if(!f || magic!=0x5241494b50494332ULL || n>2000000)throw std::runtime_error("Invalid checkpoint");
 for(std::uint64_t i=0;i<n;++i){GapParticle p;get(p);if(!s.p.emplace(p.id,p).second)throw std::runtime_error("Duplicate particle");}
 if(!f || f.peek()!=std::char_traits<char>::eof())throw std::runtime_error("Truncated checkpoint");
 return s;
}
int main(int argc,char** argv) {
 try {
  if(argc!=3)throw std::runtime_error("Expected two checkpoint output directories");
  auto a=read(argv[1]),b=read(argv[2]);double pos=0,vel=0,vscale=0,weight=0;bool ids=a.p.size()==b.p.size();
  for(const auto& [id,p]:a.p) {
   auto it=b.p.find(id);if(it==b.p.end()){ids=false;continue;}const auto& q=it->second;
   for(int k=0;k<3;++k){pos=std::max(pos,std::abs(p.position_m[k]-q.position_m[k]));vel=std::max(vel,std::abs(p.proper_velocity_m_s[k]-q.proper_velocity_m_s[k]));vscale=std::max(vscale,std::abs(p.proper_velocity_m_s[k]));}
   weight=std::max(weight,std::abs(p.weight-q.weight));
  }
  double relative=vel/std::max(vscale,1e-100);
  bool pass=ids && weight==0 && a.d.time_s==b.d.time_s && pos<1e-12 && relative<1e-10;
  std::cout<<std::setprecision(17)<<"{\"passed\":"<<(pass?"true":"false")<<",\"identical_particle_ids\":"<<(ids?"true":"false")<<",\"particles\":"<<a.p.size()<<",\"max_position_difference_m\":"<<pos<<",\"max_proper_velocity_difference_m_s\":"<<vel<<",\"global_scale_relative_velocity_difference\":"<<relative<<",\"max_weight_difference\":"<<weight<<",\"position_limit_m\":1e-12,\"relative_velocity_limit\":1e-10}\n";
  return pass?0:1;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 2;}
}
