#include "cudss_tuning_options.hpp"
#include <array>
#include <iostream>
int main(){
 try {
  using pic_cudss_tuning::Options;auto o=Options::parse(nullptr,nullptr,nullptr);if(o.amd || o.nd_levels!=-1 || o.hybrid_execute!=0)return 1;
  for(const char* level:{"128","256","512","1024"})if(Options::parse("default",level,"0").nd_levels<128)return 2;
  if(!Options::parse("amd","default","1").amd)return 3;
  int rejected=0;
  for(const auto& args:{std::array<const char*,3>{"amd","256","0"},{"other",nullptr,nullptr},{nullptr,"123",nullptr},{nullptr,nullptr,"2"},{"",nullptr,nullptr},{nullptr,"",nullptr}}) {
   try{(void)Options::parse(args[0],args[1],args[2]);}catch(const std::exception&){++rejected;}
  }
  if(rejected!=6)return 4;std::cout<<"Default, four ND values, AMD/hybrid, six invalid combinations: passed\n";return 0;
 }catch(...){return 5;}
}
