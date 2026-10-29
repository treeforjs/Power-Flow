#pragma once
#include <cstdlib>
#include <stdexcept>
#include <string>
namespace pic_cudss_tuning {
struct Options {
 bool amd=false;int nd_levels=-1,hybrid_execute=0;
 static Options parse(const char* reorder,const char* levels,const char* hybrid) {
  Options o;
  if(reorder){const std::string v(reorder);if(v=="amd")o.amd=true;else if(v!="default")throw std::runtime_error("CUDSS_REORDER must be default or amd");}
  if(levels && std::string(levels)!="default") {
   const std::string v(levels);
   if(v=="128")o.nd_levels=128;else if(v=="256")o.nd_levels=256;else if(v=="512")o.nd_levels=512;
   else if(v=="1024")o.nd_levels=1024;else throw std::runtime_error("CUDSS_ND_NLEVELS must be default, 128, 256, 512 or 1024");
  }
  if(o.amd && o.nd_levels!=-1)throw std::runtime_error("ND_NLEVELS applies only to default reordering");
  if(hybrid){const std::string v(hybrid);if(v=="1")o.hybrid_execute=1;else if(v!="0")throw std::runtime_error("CUDSS_HYBRID_EXECUTE must be 0 or 1");}
  return o;
 }
 static Options environment(){return parse(std::getenv("RAIKOU_PIC_CUDSS_REORDER"),std::getenv("RAIKOU_PIC_CUDSS_ND_NLEVELS"),std::getenv("RAIKOU_PIC_CUDSS_HYBRID_EXECUTE"));}
};
}
