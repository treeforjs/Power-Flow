#include <mfem.hpp>
#include "cudss_poisson.hpp"
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc,char** argv) {
 if(argc>1 && std::string(argv[1])=="--help") {
  std::cout<<"Qualify unchanged cuDSS solves against required CUDA Graph replay on 32 right hand sides.\n";
  return 0;
 }
 try {
  constexpr int n=33;
  mfem::SparseMatrix matrix(n);
  for(int i=0;i<n;++i) {
   matrix.Add(i,i,4.0+0.125*(i%5));
   if(i>0)matrix.Add(i,i-1,-1.0);
   if(i+1<n)matrix.Add(i,i+1,-1.0);
  }
  matrix.Finalize();
  ::setenv("RAIKOU_PIC_CUDSS_GRAPH","off",1);
  raikou_pic_parallel::CudssPoisson direct(matrix);
  ::setenv("RAIKOU_PIC_CUDSS_GRAPH","required",1);
  raikou_pic_parallel::CudssPoisson captured(matrix);
  for(int test=0;test<32;++test) {
   mfem::Vector rhs(n),reference,answer;
   for(int i=0;i<n;++i) {
    const double scale=test<8?1.0:std::ldexp(1.0,test-16);
    rhs[i]=test==1?0.0:scale*((i*7+test*11)%23-11);
   }
   direct.Solve(rhs,reference);captured.Solve(rhs,answer);
   if(std::memcmp(reference.GetData(),answer.GetData(),n*sizeof(double))!=0)
    throw std::runtime_error("Direct and graph answers differ bitwise on RHS "+std::to_string(test));
   mfem::Vector residual(n);matrix.Mult(answer,residual);residual-=rhs;
   if(!std::isfinite(residual.Norml2()) ||
      residual.Norml2()>std::max(1e-30,2e-13*rhs.Norml2()))
    throw std::runtime_error("Unchanged Poisson residual guard failed on RHS "+std::to_string(test));
  }
  std::cout<<"32 right hand sides: direct and graph answers bitwise identical; residual guards passed.\n";
  return 0;
 } catch(const std::exception& error) {
  std::cerr<<error.what()<<'\n';return 1;
 }
}
