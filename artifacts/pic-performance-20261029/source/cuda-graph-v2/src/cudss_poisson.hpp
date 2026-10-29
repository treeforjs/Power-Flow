#pragma once
// NVIDIA license and disclaimer retained alongside unmodified SDK headers.
#include <cuda_runtime.h>
#include <cudss.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include "bench_profile.hpp"
#include "cudss_graph_record.hpp"

namespace raikou_pic_parallel {
class CudssPoisson {
 enum class GraphMode { Off, Required, Auto };
 int n; int *di=nullptr,*dj=nullptr;
 double *da=nullptr,*db=nullptr,*dx=nullptr;
 double *host_rhs=nullptr,*host_answer=nullptr;
 cudssHandle_t handle=nullptr; cudssConfig_t config=nullptr;
 cudssData_t data=nullptr; cudssMatrix_t A=nullptr,b=nullptr,x=nullptr;
 cudaStream_t stream=nullptr;
 cudaGraph_t graph=nullptr;
 cudaGraphExec_t graph_exec=nullptr;
 GraphMode graph_mode=GraphMode::Off;
 bool graph_qualified=false;
 bool fell_back=false;
 std::size_t record_id=std::size_t(-1);
 std::size_t graph_nodes=0,kernel_nodes=0,memcpy_nodes=0;

 static void cuda_ok(cudaError_t code) {
  if(code!=cudaSuccess)throw std::runtime_error(std::string("CUDA: ")+cudaGetErrorString(code));
 }
 static void ok(cudssStatus_t code) {
  if(code!=CUDSS_STATUS_SUCCESS)throw std::runtime_error("cuDSS status "+std::to_string(int(code)));
 }
 static GraphMode graph_setting() {
  const char* value=std::getenv("RAIKOU_PIC_CUDSS_GRAPH");
  if(!value || std::string(value)=="off")return GraphMode::Off;
  if(std::string(value)=="required")return GraphMode::Required;
  if(std::string(value)=="auto")return GraphMode::Auto;
  throw std::runtime_error("RAIKOU_PIC_CUDSS_GRAPH must be off, required, or auto");
 }
 static int device_alloc(void*,void** ptr,size_t bytes,cudaStream_t allocation_stream) {
  return int(cudaMallocAsync(ptr,bytes,allocation_stream));
 }
 static int device_free(void*,void* ptr,size_t,cudaStream_t allocation_stream) {
  return int(cudaFreeAsync(ptr,allocation_stream));
 }
 static bool documented_unsupported(cudaError_t code) {
  return code==cudaErrorStreamCaptureUnsupported || code==cudaErrorNotSupported;
 }
 static const char* mode_name(GraphMode mode) {
  return mode==GraphMode::Off?"off":mode==GraphMode::Required?"required":"auto";
 }
 void record(const char* status,const char* reason="") {
  std::cerr<<"cuDSS graph status="<<status<<" reason="<<reason
           <<" ir_steps=2 nodes="<<graph_nodes<<" kernels="<<kernel_nodes
           <<" memcpy="<<memcpy_nodes<<'\n';
  pic_graph::RunRecord::instance().update(record_id,status,reason,graph_nodes,kernel_nodes,memcpy_nodes,graph_qualified);
 }
 bool unsupported(const char* reason) {
  if(graph_mode==GraphMode::Required) {
   record("unsupported",reason);
   throw std::runtime_error(std::string("Required cuDSS graph capture unsupported: ")+reason);
  }
  graph_mode=GraphMode::Off;
  fell_back=true;
  record("fallback",reason);
  return false;
 }
 void release() noexcept {
  bool failed=false;
  const auto cuda_cleanup=[&](cudaError_t status){failed|=status!=cudaSuccess;};
  const auto cudss_cleanup=[&](cudssStatus_t status){failed|=status!=CUDSS_STATUS_SUCCESS;};
  if(stream)cuda_cleanup(cudaStreamSynchronize(stream));
  if(graph_exec)cuda_cleanup(cudaGraphExecDestroy(graph_exec));
  if(graph)cuda_cleanup(cudaGraphDestroy(graph));
  graph_exec=nullptr;graph=nullptr;
  if(A)cudss_cleanup(cudssMatrixDestroy(A));if(b)cudss_cleanup(cudssMatrixDestroy(b));if(x)cudss_cleanup(cudssMatrixDestroy(x));
  A=nullptr;b=nullptr;x=nullptr;
  if(data)cudss_cleanup(cudssDataDestroy(handle,data));if(config)cudss_cleanup(cudssConfigDestroy(config));
  if(handle)cudss_cleanup(cudssDestroy(handle));
  data=nullptr;config=nullptr;handle=nullptr;
  if(stream)cuda_cleanup(cudaStreamSynchronize(stream));
  if(di)cuda_cleanup(cudaFree(di));if(dj)cuda_cleanup(cudaFree(dj));if(da)cuda_cleanup(cudaFree(da));
  if(db)cuda_cleanup(cudaFree(db));if(dx)cuda_cleanup(cudaFree(dx));
  di=nullptr;dj=nullptr;da=nullptr;db=nullptr;dx=nullptr;
  if(host_rhs)cuda_cleanup(cudaFreeHost(host_rhs));if(host_answer)cuda_cleanup(cudaFreeHost(host_answer));
  host_rhs=nullptr;host_answer=nullptr;
  if(stream)cuda_cleanup(cudaStreamDestroy(stream));
  stream=nullptr;
  if(record_id!=std::size_t(-1)) {
   auto& run=pic_graph::RunRecord::instance();
   if(failed)run.fail(record_id,"CUDA/cuDSS cleanup failed");
   run.close(record_id);
  }
 }
 void direct_solve(const double* rhs,double* answer) {
  if(stream) {
   cuda_ok(cudaMemcpyAsync(db,rhs,n*sizeof(double),cudaMemcpyHostToDevice,stream));
   ok(cudssExecute(handle,CUDSS_PHASE_SOLVE,config,data,A,x,b));
   cuda_ok(cudaMemcpyAsync(answer,dx,n*sizeof(double),cudaMemcpyDeviceToHost,stream));
   cuda_ok(cudaStreamSynchronize(stream));
  } else {
   cuda_ok(cudaMemcpy(db,rhs,n*sizeof(double),cudaMemcpyHostToDevice));
   ok(cudssExecute(handle,CUDSS_PHASE_SOLVE,config,data,A,x,b));
   cuda_ok(cudaMemcpy(answer,dx,n*sizeof(double),cudaMemcpyDeviceToHost));
  }
  pic_graph::RunRecord::instance().direct_solve(record_id);
 }
 void inspect_graph(cudaGraph_t inspected) {
  size_t count=0;cuda_ok(cudaGraphGetNodes(inspected,nullptr,&count));
  std::vector<cudaGraphNode_t> nodes(count);
  cuda_ok(cudaGraphGetNodes(inspected,nodes.data(),&count));
  graph_nodes+=count;
  for(auto node:nodes) {
   cudaGraphNodeType type;cuda_ok(cudaGraphNodeGetType(node,&type));
   if(type==cudaGraphNodeTypeKernel)++kernel_nodes;
   if(type==cudaGraphNodeTypeMemcpy)++memcpy_nodes;
   if(type==cudaGraphNodeTypeGraph) {
    cudaGraph_t child=nullptr;cuda_ok(cudaGraphChildGraphNodeGetGraph(node,&child));
    inspect_graph(child);
   }
  }
 }
 bool capture() {
  pic_bench::Timer timer("cuda_graph_capture");
  cudaError_t begin=cudaStreamBeginCapture(stream,cudaStreamCaptureModeThreadLocal);
  if(documented_unsupported(begin))return unsupported("cudaStreamBeginCapture unsupported");
  cuda_ok(begin);
  cudaError_t first=cudaMemcpyAsync(db,host_rhs,n*sizeof(double),cudaMemcpyHostToDevice,stream);
  cudssStatus_t solve=CUDSS_STATUS_SUCCESS;
  cudaError_t last=cudaSuccess;
  if(first==cudaSuccess) {
   solve=cudssExecute(handle,CUDSS_PHASE_SOLVE,config,data,A,x,b);
   if(solve==CUDSS_STATUS_SUCCESS)
    last=cudaMemcpyAsync(host_answer,dx,n*sizeof(double),cudaMemcpyDeviceToHost,stream);
  }
  // End capture even after an error so an invalidated stream leaves capture mode.
  cudaError_t end=cudaStreamEndCapture(stream,&graph);
  const bool unsupported_call=documented_unsupported(first) || documented_unsupported(last) ||
                              documented_unsupported(end) || solve==CUDSS_STATUS_NOT_SUPPORTED;
  // Never let a secondary capture error hide an unrelated runtime failure.
  if(first!=cudaSuccess && !documented_unsupported(first))cuda_ok(first);
  if(solve!=CUDSS_STATUS_SUCCESS && solve!=CUDSS_STATUS_NOT_SUPPORTED)ok(solve);
  if(last!=cudaSuccess && !documented_unsupported(last))cuda_ok(last);
  if(end!=cudaSuccess && !documented_unsupported(end) &&
     !(end==cudaErrorStreamCaptureInvalidated && unsupported_call))cuda_ok(end);
  if(unsupported_call) {
   if(graph){cuda_ok(cudaGraphDestroy(graph));graph=nullptr;}
   cudaGetLastError();
   return unsupported("CUDA or cuDSS reported a documented unsupported feature during capture");
  }
  cuda_ok(first);ok(solve);cuda_ok(last);cuda_ok(end);
  inspect_graph(graph);
  if(kernel_nodes==0 || memcpy_nodes<2) {
   record("failed","Captured graph lacks solve kernels or transfer nodes");
   throw std::runtime_error("Captured cuDSS graph lacks solve kernels or transfer nodes");
  }
  cudaError_t instantiated=cudaGraphInstantiate(&graph_exec,graph,nullptr,nullptr,0);
  if(documented_unsupported(instantiated)) {
   cuda_ok(cudaGraphDestroy(graph));graph=nullptr;
   return unsupported("cudaGraphInstantiate unsupported");
  }
  cuda_ok(instantiated);
  return true;
 }
 void qualify(const mfem::Vector& rhs,mfem::Vector& answer) {
  pic_bench::Timer timer("cuda_graph_first_rhs_qualification");
  std::copy_n(rhs.GetData(),n,host_rhs);
  // Warm up the exact unchanged solve, including both refinement iterations.
  direct_solve(host_rhs,host_answer);
  std::vector<double> reference(host_answer,host_answer+n);
  if(!capture()) {
   std::copy(reference.begin(),reference.end(),answer.GetData());
   return;
  }
  cuda_ok(cudaGraphLaunch(graph_exec,stream));
  cuda_ok(cudaStreamSynchronize(stream));
  pic_graph::RunRecord::instance().graph_replay(record_id);
  if(std::memcmp(reference.data(),host_answer,n*sizeof(double))!=0) {
   record("failed","First RHS direct and graph answers differ bitwise");
   throw std::runtime_error("cuDSS graph qualification failed: first RHS differs bitwise");
  }
  graph_qualified=true;
  record("qualified","First RHS direct and graph answers are bitwise identical");
  std::copy_n(host_answer,n,answer.GetData());
 }
public:
 explicit CudssPoisson(const mfem::SparseMatrix& matrix):n(matrix.Height()) {
  pic_bench::Timer timer("cuda_factor_setup");
  try {
   graph_mode=graph_setting();
   record_id=pic_graph::RunRecord::instance().add(n,mode_name(graph_mode));
   cuda_ok(cudaSetDevice(0));ok(cudssCreate(&handle));
   if(graph_mode!=GraphMode::Off) {
    int pools=0;cuda_ok(cudaDeviceGetAttribute(&pools,cudaDevAttrMemoryPoolsSupported,0));
    if(!pools)unsupported("CUDA device does not support stream-ordered memory pools");
    if(graph_mode!=GraphMode::Off) {
     cuda_ok(cudaStreamCreateWithFlags(&stream,cudaStreamNonBlocking));
     ok(cudssSetStream(handle,stream));
     cudssDeviceMemHandler_t memory{};
     memory.device_alloc=&device_alloc;memory.device_free=&device_free;
     std::strcpy(memory.name,"stream_ordered_cuda_allocator");
     ok(cudssSetDeviceMemHandler(handle,&memory));
     cuda_ok(cudaMallocHost(reinterpret_cast<void**>(&host_rhs),n*sizeof(double)));
     cuda_ok(cudaMallocHost(reinterpret_cast<void**>(&host_answer),n*sizeof(double)));
    }
   }
   ok(cudssConfigCreate(&config));ok(cudssDataCreate(handle,&data));
   int one=1,zero=0,two=2;
   ok(cudssConfigSet(config,CUDSS_CONFIG_HOST_NTHREADS,&one,sizeof(one)));
   ok(cudssConfigSet(config,CUDSS_CONFIG_HYBRID_EXECUTE_MODE,&zero,sizeof(zero)));
   ok(cudssConfigSet(config,CUDSS_CONFIG_HYBRID_MODE,&zero,sizeof(zero)));
   ok(cudssConfigSet(config,CUDSS_CONFIG_IR_N_STEPS,&two,sizeof(two)));
   std::vector<int> rows(n+1),cols;std::vector<double> values;
   for(int i=0;i<n;++i) {
    std::vector<std::pair<int,double>> row;bool diagonal=false;
    for(int k=matrix.GetI()[i];k<matrix.GetI()[i+1];++k) {
     int j=matrix.GetJ()[k];double v=matrix.GetData()[k];
     if(std::abs(v-matrix(j,i))>1e-14*std::max({std::abs(v),std::abs(matrix(j,i)),1e-100}))
      throw std::runtime_error("Nonsymmetric SPD input");
     if(j==i){diagonal=true;if(!(v>0))throw std::runtime_error("Nonpositive Poisson diagonal");}
     if(j>=i && v!=0)row.emplace_back(j,v);
    }
    if(!diagonal)throw std::runtime_error("Missing Poisson diagonal");
    std::sort(row.begin(),row.end());
    for(auto [j,v]:row){cols.push_back(j);values.push_back(v);}
    rows[i+1]=int(cols.size());
   }
   cuda_ok(cudaMalloc(reinterpret_cast<void**>(&di),rows.size()*sizeof(int)));
   cuda_ok(cudaMalloc(reinterpret_cast<void**>(&dj),cols.size()*sizeof(int)));
   cuda_ok(cudaMalloc(reinterpret_cast<void**>(&da),values.size()*sizeof(double)));
   cuda_ok(cudaMalloc(reinterpret_cast<void**>(&db),n*sizeof(double)));
   cuda_ok(cudaMalloc(reinterpret_cast<void**>(&dx),n*sizeof(double)));
   cuda_ok(cudaMemcpy(di,rows.data(),rows.size()*sizeof(int),cudaMemcpyHostToDevice));
   cuda_ok(cudaMemcpy(dj,cols.data(),cols.size()*sizeof(int),cudaMemcpyHostToDevice));
   cuda_ok(cudaMemcpy(da,values.data(),values.size()*sizeof(double),cudaMemcpyHostToDevice));
   cuda_ok(cudaMemset(db,0,n*sizeof(double)));cuda_ok(cudaMemset(dx,0,n*sizeof(double)));
   // The cuDSS nonblocking stream does not inherit default-stream setup order.
   if(stream)cuda_ok(cudaDeviceSynchronize());
   ok(cudssMatrixCreateCsr(&A,n,n,cols.size(),di,nullptr,dj,da,CUDA_R_32I,CUDA_R_64F,
                           CUDSS_MTYPE_SPD,CUDSS_MVIEW_UPPER,CUDSS_BASE_ZERO));
   ok(cudssMatrixCreateDn(&b,n,1,n,db,CUDA_R_64F,CUDSS_LAYOUT_COL_MAJOR));
   ok(cudssMatrixCreateDn(&x,n,1,n,dx,CUDA_R_64F,CUDSS_LAYOUT_COL_MAJOR));
   ok(cudssExecute(handle,CUDSS_PHASE_ANALYSIS,config,data,A,x,b));
   ok(cudssExecute(handle,CUDSS_PHASE_FACTORIZATION,config,data,A,x,b));
   cuda_ok(cudaDeviceSynchronize());int info=-1;size_t written=0;
   ok(cudssDataGet(handle,data,CUDSS_DATA_INFO,&info,sizeof(info),&written));
   if(info!=0)throw std::runtime_error("cuDSS factorization info "+std::to_string(info));
   if(!fell_back)record(graph_mode==GraphMode::Off?"off":"pending");
  } catch(const std::exception& error) {
   pic_graph::RunRecord::instance().fail(record_id,error.what());release();throw;
  } catch(...) {
   pic_graph::RunRecord::instance().fail(record_id,"Unknown solver construction failure");release();throw;
  }
 }
 CudssPoisson(const CudssPoisson&)=delete;
 CudssPoisson& operator=(const CudssPoisson&)=delete;
 void Solve(const mfem::Vector& rhs,mfem::Vector& answer) {
  try {
  pic_bench::Timer timer("cuda_solve_including_transfers");answer.SetSize(n);
  if(rhs.Size()!=n)throw std::runtime_error("cuDSS Poisson right hand side size changed");
  if(graph_mode==GraphMode::Off) {direct_solve(rhs.GetData(),answer.GetData());return;}
  if(!graph_qualified) {qualify(rhs,answer);return;}
  std::copy_n(rhs.GetData(),n,host_rhs);
  cuda_ok(cudaGraphLaunch(graph_exec,stream));
  cuda_ok(cudaStreamSynchronize(stream));
  pic_graph::RunRecord::instance().graph_replay(record_id);
  std::copy_n(host_answer,n,answer.GetData());
  } catch(const std::exception& error) {
   pic_graph::RunRecord::instance().fail(record_id,error.what());throw;
  } catch(...) {
   pic_graph::RunRecord::instance().fail(record_id,"Unknown solver runtime failure");throw;
  }
 }
 ~CudssPoisson() {release();}
};
}
