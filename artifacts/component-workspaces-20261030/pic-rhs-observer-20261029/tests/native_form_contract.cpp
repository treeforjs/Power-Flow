#include <mfem.hpp>
#include "../source/src/rhs_observer.hpp"
#include <cassert>
#include <iostream>

int main(int argc,char** argv) {
 assert(argc==2);
 const std::filesystem::path output(argv[1]);
 auto mesh=mfem::Mesh::MakeCartesian2D(2,2,mfem::Element::QUADRILATERAL,true,1.,1.);
 mfem::H1_FECollection collection(1,2);
 mfem::FiniteElementSpace space(&mesh,&collection);
 mfem::BilinearForm form(&space);
 mfem::ConstantCoefficient one(1.),ten(10.);
 form.AddDomainIntegrator(new mfem::DiffusionIntegrator(one));form.Assemble();form.Finalize();
 mfem::Array<int> marker(mesh.bdr_attributes.Max()),essential;marker=1;
 space.GetEssentialTrueDofs(marker,essential);
 const int n=space.GetVSize();
 std::vector<bool> mask(n,false);for(int dof:essential)mask[dof]=true;
 int interior=-1;for(int i=0;i<n;++i)if(!mask[i])interior=i;
 assert(interior>=0);
 // Match the production ordering: boundary projection first, then a replaced
 // final essential list. Projected outer DOFs need not be final essential.
 essential.SetSize(1);essential[0]=interior;
 mfem::GridFunction trial(&space);trial=0.;
 rhs_observer::Recorder observer(output,8,1024*1024,512*1024,true);
 std::vector<double> first_projected;
 for(unsigned step=1;step<=8;++step) {
  observer.driver_step(step,step*1e-14,1e-14);
  observer.begin_advance(true,1,step*1e-14,1e-14,0,0);
  for(unsigned half=1;half<=2;++half) {
   observer.set_half(half);
   if(half==1)for(int i=0;i<n;++i)if(!mask[i])trial[i]=100.+step+i;
   trial.ProjectBdrCoefficient(ten,marker);trial[interior]=0.;
   mfem::LinearForm rhs(&space);rhs=0.;
   mfem::Vector full_rhs(rhs),x,b,error(n);error=0.;
   if(step==1 && half==1)first_projected.assign(trial.GetData(),trial.GetData()+n);
   auto token=observer.before_form(trial.GetData(),full_rhs.GetData(),n,
                                  essential.GetData(),essential.Size(),0x11);
   mfem::OperatorPtr matrix;
   // The exact frozen MFEM API, with its real default copy_interior=0.
   form.FormLinearSystem(essential,trial,rhs,matrix,x,b);
   auto* sparse=dynamic_cast<mfem::SparseMatrix*>(matrix.Ptr());assert(sparse);
   for(int i=0;i<n;++i)assert(x[i]==0. && b[i]==0.);
   observer.after_form(token,sparse->GetI(),sparse->GetJ(),sparse->GetData(),n,
                       sparse->GetI()[n],x.GetData(),b.GetData(),
                       reinterpret_cast<std::uintptr_t>(sparse));
   // Zero is the known exact answer and residual for this synthetic zero-RHS
   // contract. No linear solver, PIC engine, CUDA API, or GPU is called.
   observer.accepted(b.GetData(),x.GetData(),error.GetData(),n,0,0,9,0x22);
  }
  observer.end_advance();
 }
 observer.finish();
 std::vector<double> recorded(n);
 std::ifstream input(output/"call_0001_projected_previous_phi.f64",std::ios::binary);
 input.read(reinterpret_cast<char*>(recorded.data()),n*sizeof(double));assert(input);
 assert(recorded==first_projected);
 assert(std::any_of(recorded.begin(),recorded.end(),[](double value){return value!=0.;}));
 std::cout<<"Frozen native MFEM FormLinearSystem default overwrite and actual pre-form snapshot passed; no solver/PIC/GPU execution.\n";
}
