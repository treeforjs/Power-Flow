#pragma once
// NVIDIA license and disclaimer retained alongside unmodified SDK headers.
#include <cuda_runtime.h>
#include <cudss.h>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>
#include "bench_profile.hpp"
namespace raikou_pic_parallel {
class CudssPoisson {
 int n; int *di=nullptr,*dj=nullptr;
 double *da=nullptr,*db=nullptr,*dx=nullptr;
 cudssHandle_t handle=nullptr; cudssConfig_t config=nullptr;
 cudssData_t data=nullptr; cudssMatrix_t A=nullptr,b=nullptr,x=nullptr;
 static void cuda_ok(cudaError_t code) {if(code!=cudaSuccess)throw std::runtime_error(std::string("CUDA: ")+cudaGetErrorString(code));}
 static void ok(cudssStatus_t code) {if(code!=CUDSS_STATUS_SUCCESS)throw std::runtime_error("cuDSS status "+std::to_string(int(code)));}
public:
 explicit CudssPoisson(const mfem::SparseMatrix& matrix):n(matrix.Height()) {
  pic_bench::Timer timer("cuda_factor_setup");
  pic_phase::Scope phase_factor(pic_phase::Phase::FactorSetup); // PHASE_TIMING
  cuda_ok(cudaSetDevice(0));ok(cudssCreate(&handle));
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
    if(std::abs(v-matrix(j,i))>1e-14*std::max({std::abs(v),std::abs(matrix(j,i)),1e-100}))throw std::runtime_error("Nonsymmetric SPD input");
    if(j==i){diagonal=true;if(!(v>0))throw std::runtime_error("Nonpositive Poisson diagonal");}
    if(j>=i && v!=0)row.emplace_back(j,v);
   }
   if(!diagonal)throw std::runtime_error("Missing Poisson diagonal");
   std::sort(row.begin(),row.end());
   for(auto [j,v]:row){cols.push_back(j);values.push_back(v);}rows[i+1]=int(cols.size());
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
  ok(cudssMatrixCreateCsr(&A,n,n,cols.size(),di,nullptr,dj,da,CUDA_R_32I,CUDA_R_64F,CUDSS_MTYPE_SPD,CUDSS_MVIEW_UPPER,CUDSS_BASE_ZERO));
  ok(cudssMatrixCreateDn(&b,n,1,n,db,CUDA_R_64F,CUDSS_LAYOUT_COL_MAJOR));
  ok(cudssMatrixCreateDn(&x,n,1,n,dx,CUDA_R_64F,CUDSS_LAYOUT_COL_MAJOR));
  ok(cudssExecute(handle,CUDSS_PHASE_ANALYSIS,config,data,A,x,b));
  ok(cudssExecute(handle,CUDSS_PHASE_FACTORIZATION,config,data,A,x,b));
  cuda_ok(cudaDeviceSynchronize());int info=-1;size_t written=0;
  ok(cudssDataGet(handle,data,CUDSS_DATA_INFO,&info,sizeof(info),&written));
  if(info!=0)throw std::runtime_error("cuDSS factorization info "+std::to_string(info));
 }
 void Solve(const mfem::Vector& rhs,mfem::Vector& answer) {
  pic_bench::Timer timer("cuda_solve_including_transfers");answer.SetSize(n);
  pic_phase::Scope phase_cuda(pic_phase::solve_phase(pic_phase::Phase::CudaFirst,pic_phase::Phase::CudaSecond,pic_phase::Phase::CudaOther)); // PHASE_TIMING
  cuda_ok(cudaMemcpy(db,rhs.GetData(),n*sizeof(double),cudaMemcpyHostToDevice));
  ok(cudssExecute(handle,CUDSS_PHASE_SOLVE,config,data,A,x,b));
  cuda_ok(cudaMemcpy(answer.GetData(),dx,n*sizeof(double),cudaMemcpyDeviceToHost));
 }
 ~CudssPoisson() {
  if(A)cudssMatrixDestroy(A);if(b)cudssMatrixDestroy(b);if(x)cudssMatrixDestroy(x);
  if(data)cudssDataDestroy(handle,data);if(config)cudssConfigDestroy(config);if(handle)cudssDestroy(handle);
  cudaFree(di);cudaFree(dj);cudaFree(da);cudaFree(db);cudaFree(dx);
 }
};
}
