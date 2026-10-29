// CPU-only RHS reconstruction using the frozen production MFEM assembly order.
// No linear solve or CUDA dependency is used.
#include "raikou/mfem_gap_pic.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <filesystem>
using namespace raikou::mfem_backend;
int main(int argc,char** argv){try{
 if(argc!=2)throw std::runtime_error("Expected checkpoint file");
 mfem::Device device("cpu");
 std::ifstream in(argv[1],std::ios::binary);std::uint64_t magic=0,next=0,n=0;
 GapPicOptions o;GapPicDiagnostics d;
 auto get=[&](auto& x){in.read(reinterpret_cast<char*>(&x),sizeof(x));};
 get(magic);get(o);get(d);get(next);get(n);
 if(!in || magic!=0x5241494b50494332ULL || n>2000000 || o.nx<1 || o.nx>4096 || o.ny<1 || o.ny>4096)throw std::runtime_error("Invalid checkpoint header");
 for(double v:{o.width_m,o.gap_m,o.depth_m})if(!std::isfinite(v)||v<=0)throw std::runtime_error("Invalid geometry");
 if(std::filesystem::file_size(argv[1])!=8+sizeof(o)+sizeof(d)+16+n*sizeof(GapParticle))throw std::runtime_error("Truncated or trailing checkpoint");
 auto mesh=mfem::Mesh::MakeCartesian2D(o.nx,o.ny,mfem::Element::QUADRILATERAL,true,o.width_m,o.gap_m,false);
 mfem::H1_FECollection collection(1,2);mfem::FiniteElementSpace space(&mesh,&collection);
 mfem::ConstantCoefficient epsilon(8.8541878188e-12);mfem::BilinearForm form(&space);
 form.AddDomainIntegrator(new mfem::DiffusionIntegrator(epsilon));form.Assemble();form.Finalize();
 mfem::LinearForm rhs(&space);rhs=0.;mfem::Array<int>dofs;mfem::Vector shape;
 for(std::uint64_t i=0;i<n;++i){
  GapParticle p;get(p);if(!in)throw std::runtime_error("Truncated particle");
  for(double v:p.position_m)if(!std::isfinite(v))throw std::runtime_error("Nonfinite particle position");
  if(p.position_m[0]<0||p.position_m[0]>o.width_m||p.position_m[1]<0||p.position_m[1]>o.gap_m||!std::isfinite(p.charge_C)||!std::isfinite(p.weight)||p.weight<=0)throw std::runtime_error("Invalid particle");
  double sx=p.position_m[0]/o.width_m*o.nx,sy=p.position_m[1]/o.gap_m*o.ny;
  int ix=std::min(static_cast<int>(sx),o.nx-1),iy=std::min(static_cast<int>(sy),o.ny-1),element=iy*o.nx+ix;
  mfem::IntegrationPoint point;point.Set2(sx-ix,sy-iy);
  const auto* fe=space.GetFE(element);shape.SetSize(fe->GetDof());fe->CalcShape(point,shape);
  space.GetElementDofs(element,dofs);rhs.AddElementVector(dofs,shape*=p.charge_C*p.weight/o.depth_m);
 }
 mfem::Array<int> marker(mesh.bdr_attributes.Max());marker=0;mfem::Array<int>vertices;
 for(int b=0;b<mesh.GetNBE();++b){mesh.GetBdrElementVertices(b,vertices);double y0=mesh.GetVertex(vertices[0])[1],y1=mesh.GetVertex(vertices[1])[1];if(y0==y1&&(y0==0||y0==o.gap_m))marker[mesh.GetBdrAttribute(b)-1]=1;}
 mfem::Array<int> essential;space.GetEssentialTrueDofs(marker,essential);
 mfem::FunctionCoefficient voltage([&](const mfem::Vector& x){return o.lower_potential_V+(o.upper_potential_V-o.lower_potential_V)*x[1]/o.gap_m;});
 mfem::GridFunction trial(&space);trial=0.;trial.ProjectBdrCoefficient(voltage,marker);
 if(o.finite_foils){
  mfem::GridFunction xfield(&space),yfield(&space);mfem::FunctionCoefficient xc([](const mfem::Vector& x){return x[0];}),yc([](const mfem::Vector& x){return x[1];});xfield.ProjectCoefficient(xc);yfield.ProjectCoefficient(yc);
  essential.SetSize(0);double tol=1e-12*std::max(o.width_m,o.gap_m);
  for(int i=0;i<trial.Size();++i){int who=0;for(const auto& f:{o.lower_foil,o.upper_foil}){if(xfield[i]>=f.left_m-tol&&xfield[i]<=f.right_m+tol&&yfield[i]>=f.bottom_m-tol&&yfield[i]<=f.top_m+tol){essential.Append(i);trial[i]=who==0?o.lower_potential_V:o.upper_potential_V;}++who;}}
 }
 mfem::OperatorPtr matrix;mfem::Vector x,b;form.FormLinearSystem(essential,trial,rhs,matrix,x,b);
 const double norm=b.Norml2(),limit=std::max(1e-30,2e-13*norm);
 if(!std::isfinite(norm)||!std::isfinite(d.poisson_residual_C_per_m)||d.poisson_residual_C_per_m<0)throw std::runtime_error("Nonfinite RHS/residual");
 std::cout<<std::setprecision(17)<<"{\"rhs_norm\":"<<norm<<",\"residual_norm\":"<<d.poisson_residual_C_per_m<<",\"guard_limit\":"<<limit<<",\"guard_passed\":"<<(d.poisson_residual_C_per_m<=limit?"true":"false")<<",\"ndofs\":"<<space.GetVSize()<<",\"essential_dofs\":"<<essential.Size()<<",\"particles\":"<<n<<",\"cpu_only\":true}\n";
 return d.poisson_residual_C_per_m<=limit?0:1;
}catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 2;}}
