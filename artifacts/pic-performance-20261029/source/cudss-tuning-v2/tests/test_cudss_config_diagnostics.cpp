#include "cudss_tuning_audit.hpp"
#include <iostream>
// CPU mocks only. No CUDA runtime or cuDSS library is linked into this test.
extern "C" cudssStatus_t cudssConfigSet(cudssConfig_t,cudssConfigParam_t,void*,size_t bytes) {
 return bytes==4?CUDSS_STATUS_SUCCESS:CUDSS_STATUS_INVALID_VALUE;
}
extern "C" cudssStatus_t cudssConfigGet(cudssConfig_t,cudssConfigParam_t param,void* value,size_t bytes,size_t* written) {
 if(bytes!=4)return CUDSS_STATUS_INVALID_VALUE;
 if(param==CUDSS_CONFIG_HOST_NTHREADS)return CUDSS_STATUS_NOT_SUPPORTED;
 int result=param==CUDSS_CONFIG_ND_NLEVELS?512:0;std::memcpy(value,&result,4);
 *written=param==CUDSS_CONFIG_ND_NLEVELS?0:4;return CUDSS_STATUS_SUCCESS;
}
extern "C" cudssStatus_t cudssDataGet(cudssHandle_t,cudssData_t,cudssDataParam_t,void*,size_t,size_t*) {
 return CUDSS_STATUS_NOT_SUPPORTED;
}
int main(int argc,char** argv) {
 try {
  if(argc!=2)return 1;::setenv("RAIKOU_PIC_CUDSS_AUDIT_PREFIX",argv[1],1);
  pic_cudss_tuning::ConfigDiagnostics d(1);
  if(d.set(nullptr,CUDSS_CONFIG_IR_N_STEPS,"CUDSS_CONFIG_IR_N_STEPS",2)!=CUDSS_STATUS_SUCCESS)return 2;
  if(d.set(nullptr,CUDSS_CONFIG_REORDERING_ALG,"CUDSS_CONFIG_REORDERING_ALG",CUDSS_ALG_3)!=CUDSS_STATUS_SUCCESS)return 3;
  const auto unknown=d.get<int>(nullptr,CUDSS_CONFIG_ND_NLEVELS,"CUDSS_CONFIG_ND_NLEVELS");
  if(unknown.verified() || unknown.raw!=512 || unknown.written!=0)return 4;
  const auto unsupported=d.get<int>(nullptr,CUDSS_CONFIG_HOST_NTHREADS,"CUDSS_CONFIG_HOST_NTHREADS");
  if(unsupported.verified() || unsupported.status!=CUDSS_STATUS_NOT_SUPPORTED)return 5;
  const auto known=d.get<cudssAlgType_t>(nullptr,CUDSS_CONFIG_SOLVE_ALG,"CUDSS_CONFIG_SOLVE_ALG");
  if(!known.verified() || known.raw!=0)return 6;
  pic_cudss_tuning::audit(nullptr,nullptr,nullptr,{0,1},{0},{1.},d);
  std::cout<<"CPU mocks: enum/int setters, unsupported getter, zero written count, known getter and nonfatal data queries passed\n";return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 7;}
}
