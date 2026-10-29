#pragma once
#include "cudss_tuning_options.hpp"
#include "cudss_matrix_audit.hpp"
#include <cudss.h>
#include <unistd.h>
namespace pic_cudss_tuning {
inline void audit(cudssHandle_t handle,cudssConfig_t config,cudssData_t data,
                  const std::vector<int>& rows,const std::vector<int>& cols,const std::vector<double>& values) {
 const char* root=std::getenv("RAIKOU_PIC_CUDSS_AUDIT_PREFIX");if(!root)return;
 if(!*root || *root!='/')throw std::runtime_error("CUDSS_AUDIT_PREFIX must be absolute and nonempty");
 static std::size_t next_id=0;const auto id=next_id++;const std::string prefix=std::string(root)+"."+std::to_string(id);
 if(std::filesystem::exists(prefix+".json"))throw std::runtime_error("Use a fresh CUDSS_AUDIT_PREFIX");
 const auto integer=[&](cudssConfigParam_t param) {
  int value=-1;size_t written=0;
  if(cudssConfigGet(config,param,&value,sizeof(value),&written)!=CUDSS_STATUS_SUCCESS || written!=sizeof(value))throw std::runtime_error("Cannot read effective cuDSS tuning configuration");
  return value;
 };
 const int reorder=integer(CUDSS_CONFIG_REORDERING_ALG),nd=integer(CUDSS_CONFIG_ND_NLEVELS);
 const int hybrid=integer(CUDSS_CONFIG_HYBRID_EXECUTE_MODE),ir=integer(CUDSS_CONFIG_IR_N_STEPS);
 if(ir!=2 || integer(CUDSS_CONFIG_HOST_NTHREADS)!=1 || integer(CUDSS_CONFIG_HYBRID_MODE)!=0)throw std::runtime_error("cuDSS frozen configuration guard failed");
 std::int64_t fill=-1,memory[16]{};size_t fill_bytes=0,memory_bytes=0;
 const auto fill_status=cudssDataGet(handle,data,CUDSS_DATA_LU_NNZ,&fill,sizeof(fill),&fill_bytes);
 const auto memory_status=cudssDataGet(handle,data,CUDSS_DATA_MEMORY_ESTIMATES,memory,sizeof(memory),&memory_bytes);
 Structure structure(rows,cols);std::ofstream out(prefix+".json");
 out<<"{\"diagnostics_only\":true,\"solver_id\":"<<id<<",\"pid\":"<<::getpid()<<",\"matrix_rows\":"<<rows.size()-1
    <<",\"stored_upper_nnz\":"<<cols.size()<<",\"symmetric_nonzero_nnz\":"<<2*cols.size()-(rows.size()-1)<<',';structure.json(out);
 out<<",\"reorder_enum\":"<<reorder<<",\"nd_nlevels\":"<<nd<<",\"hybrid_execute\":"<<hybrid<<",\"ir_steps\":"<<ir
    <<",\"host_threads\":1,\"hybrid_memory\":0,\"fp64_spd\":true,\"lu_nnz_status\":"<<int(fill_status)<<",\"lu_nnz_bytes\":"<<fill_bytes<<",\"lu_nnz\":"<<fill
    <<",\"memory_estimates_status\":"<<int(memory_status)<<",\"memory_estimates_bytes\":"<<memory_bytes<<",\"memory_estimates\":[";
 for(int i=0;i<16;++i){if(i)out<<',';out<<memory[i];}out<<']';
 if(const char* dump=std::getenv("RAIKOU_PIC_CUDSS_DUMP_MATRIX")) {
  if(std::string(dump)!="1")throw std::runtime_error("CUDSS_DUMP_MATRIX must be unset or 1");
  matrix_arrays(prefix,rows,cols,values);std::vector<int> permutation(rows.size()-1,-1);size_t written=0;
  const auto status=cudssDataGet(handle,data,CUDSS_DATA_PERM_REORDER_ROW,permutation.data(),permutation.size()*sizeof(int),&written);
  bool valid=status==CUDSS_STATUS_SUCCESS && written==permutation.size()*sizeof(int);
  if(valid) {
   std::vector<unsigned char> seen(permutation.size());
   for(int index:permutation)if(index<0 || std::size_t(index)>=seen.size() || seen[index]++){valid=false;break;}
   if(valid)binary(prefix+".permutation.i32",permutation);
  }
  out<<",\"permutation_status\":"<<int(status)<<",\"permutation_bytes\":"<<written<<",\"permutation_valid\":"<<(valid?"true":"false");
 }
 out<<"}\n";out.close();if(!out)throw std::runtime_error("Cannot write cuDSS audit record");
}
}
