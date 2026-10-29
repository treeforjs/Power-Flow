#pragma once
#include "cudss_tuning_options.hpp"
#include "cudss_matrix_audit.hpp"
#include "cudss_config_diagnostics.hpp"
#include <cudss.h>
#include <unistd.h>
namespace pic_cudss_tuning {
inline void audit(cudssHandle_t handle,cudssConfig_t config,cudssData_t data,
                  const std::vector<int>& rows,const std::vector<int>& cols,const std::vector<double>& values,
                  ConfigDiagnostics& diagnostics) {
 if(!diagnostics.enabled())return;
 const auto& prefix=diagnostics.prefix;const auto id=diagnostics.id;
 const auto reorder=diagnostics.get<cudssAlgType_t>(config,CUDSS_CONFIG_REORDERING_ALG,"CUDSS_CONFIG_REORDERING_ALG");
 const auto nd=diagnostics.get<int>(config,CUDSS_CONFIG_ND_NLEVELS,"CUDSS_CONFIG_ND_NLEVELS");
 const auto hybrid=diagnostics.get<int>(config,CUDSS_CONFIG_HYBRID_EXECUTE_MODE,"CUDSS_CONFIG_HYBRID_EXECUTE_MODE");
 const auto ir=diagnostics.get<int>(config,CUDSS_CONFIG_IR_N_STEPS,"CUDSS_CONFIG_IR_N_STEPS");
 const auto threads=diagnostics.get<int>(config,CUDSS_CONFIG_HOST_NTHREADS,"CUDSS_CONFIG_HOST_NTHREADS");
 const auto memory_mode=diagnostics.get<int>(config,CUDSS_CONFIG_HYBRID_MODE,"CUDSS_CONFIG_HYBRID_MODE");
 (void)diagnostics.get<cudssAlgType_t>(config,CUDSS_CONFIG_SOLVE_ALG,"CUDSS_CONFIG_SOLVE_ALG");
 (void)diagnostics.get<cudssAlgType_t>(config,CUDSS_CONFIG_FACTORIZATION_ALG,"CUDSS_CONFIG_FACTORIZATION_ALG");
 std::int64_t fill=-1,memory[16]{};size_t fill_bytes=0,memory_bytes=0;
 const auto fill_status=diagnostics.data_get(handle,data,CUDSS_DATA_LU_NNZ,"CUDSS_DATA_LU_NNZ",&fill,sizeof(fill),&fill_bytes);
 const auto memory_status=diagnostics.data_get(handle,data,CUDSS_DATA_MEMORY_ESTIMATES,"CUDSS_DATA_MEMORY_ESTIMATES",memory,sizeof(memory),&memory_bytes);
 Structure structure(rows,cols);std::ofstream out(prefix+".json");
 out<<"{\"diagnostics_only\":true,\"solver_id\":"<<id<<",\"pid\":"<<::getpid()<<",\"matrix_rows\":"<<rows.size()-1
    <<",\"stored_upper_nnz\":"<<cols.size()<<",\"symmetric_nonzero_nnz\":"<<2*cols.size()-(rows.size()-1)<<',';structure.json(out);
 out<<",\"reorder_enum\":";reorder.json(out);out<<",\"nd_nlevels\":";nd.json(out);
 out<<",\"hybrid_execute\":";hybrid.json(out);out<<",\"ir_steps\":";ir.json(out);
 out<<",\"host_threads\":";threads.json(out);out<<",\"hybrid_memory\":";memory_mode.json(out);
 out<<",\"fp64_spd\":true,\"lu_nnz_status\":"<<int(fill_status)<<",\"lu_nnz_bytes\":"<<fill_bytes<<",\"lu_nnz\":";
 if(fill_status==CUDSS_STATUS_SUCCESS && fill_bytes==sizeof(fill))out<<fill;else out<<"null";
 out<<",\"memory_estimates_status\":"<<int(memory_status)<<",\"memory_estimates_bytes\":"<<memory_bytes<<",\"memory_estimates\":";
 if(memory_status==CUDSS_STATUS_SUCCESS && memory_bytes==sizeof(memory)) {
  out<<'[';for(int i=0;i<16;++i){if(i)out<<',';out<<memory[i];}out<<']';
 } else out<<"null";
 if(const char* dump=std::getenv("RAIKOU_PIC_CUDSS_DUMP_MATRIX")) {
  if(std::string(dump)!="1")throw std::runtime_error("CUDSS_DUMP_MATRIX must be unset or 1");
  matrix_arrays(prefix,rows,cols,values);std::vector<int> permutation(rows.size()-1,-1);size_t written=0;
  const auto status=diagnostics.data_get(handle,data,CUDSS_DATA_PERM_REORDER_ROW,"CUDSS_DATA_PERM_REORDER_ROW",permutation.data(),permutation.size()*sizeof(int),&written);
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
