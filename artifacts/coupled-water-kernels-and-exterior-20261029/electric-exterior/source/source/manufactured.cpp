#include "electric_coordinates.hpp"
#include <iostream>
#include <iomanip>
#include <cstring>
using namespace mfem;
namespace ex=electric_exterior;
namespace io=snapshot_io;
double analytic(double x,double y){return 2+3*x-4*y+5*x*y;}
std::pair<double,double> smooth_poisson(int n,int margin){
 auto original=Mesh::MakeCartesian2D(n,n,Element::QUADRILATERAL,true,.009,.009,false);H1_FECollection collection(1,2);FiniteElementSpace oldspace(&original,&collection);GridFunction oldfield(&oldspace);oldfield=0.;
 std::vector<int>elements(original.GetNE());std::iota(elements.begin(),elements.end(),0);pic_batch::HostSnapshot snapshot;snapshot.parameters.self_nx=snapshot.parameters.self_ny=n;snapshot.parameters.width=snapshot.parameters.height=.009;snapshot.self_cells=io::export_grid(oldfield,elements);
 ex::Layout layout(snapshot,40,margin);auto mesh=layout.mesh();FiniteElementSpace space(&mesh,&collection);GridFunction field(&space);
 FunctionCoefficient exact([](const Vector&p){return std::exp(20*p[0])*std::cos(15*p[1]);});field.ProjectCoefficient(exact);
 ConstantCoefficient epsilon(8.8541878188e-12);BilinearForm form(&space);form.AddDomainIntegrator(new DiffusionIntegrator(epsilon));form.Assemble();form.Finalize();
 FunctionCoefficient rho([](const Vector&p){return -175*8.8541878188e-12*std::exp(20*p[0])*std::cos(15*p[1]);});LinearForm rhs(&space);rhs.AddDomainIntegrator(new DomainLFIntegrator(rho));rhs.Assemble();
 Array<int>boundary(mesh.bdr_attributes.Max());boundary=1;Array<int>essential;space.GetEssentialTrueDofs(boundary,essential);OperatorPtr op;Vector x,b;form.FormLinearSystem(essential,field,rhs,op,x,b);auto*A=dynamic_cast<SparseMatrix*>(op.Ptr());io::require(A,"No smooth Poisson matrix");GSSmoother preconditioner(*A);CGSolver solver;solver.SetOperator(*A);solver.SetPreconditioner(preconditioner);solver.SetRelTol(1e-14);solver.SetAbsTol(1e-30);solver.SetMaxIter(4000);solver.SetPrintLevel(-1);solver.Mult(b,x);io::require(solver.GetConverged(),"Smooth manufactured solve failed");Vector r(b.Size());A->Mult(x,r);r-=b;io::require(r.Norml2()<=std::max(1e-30,2e-13*b.Norml2()),"Smooth manufactured residual failed");form.RecoverFEMSolution(x,rhs,field);
 double value_error=0,gradient_error=0;
 for(int j=0;j<layout.y.cells();++j)for(int i=0;i<layout.x.cells();++i)for(double a:{.17,.53,.89}){
  double px=layout.x.nodes[i]+a*(layout.x.nodes[i+1]-layout.x.nodes[i]),py=layout.y.nodes[j]+.37*(layout.y.nodes[j+1]-layout.y.nodes[j]);auto sample=ex::gather(layout,field,{px,py});double f=std::exp(20*px)*std::cos(15*py);value_error=std::max(value_error,std::abs(sample.value-f));gradient_error=std::max({gradient_error,std::abs(sample.dx-20*f),std::abs(sample.dy+15*std::exp(20*px)*std::sin(15*py))});
 }
 return {value_error,gradient_error};
}
int main(){try{
 // A small manufactured fixture tests the coordinate/FE adapter, not PIC physics.
 auto original=Mesh::MakeCartesian2D(8,8,Element::QUADRILATERAL,true,.009,.009,false);
 H1_FECollection collection(1,2);FiniteElementSpace oldspace(&original,&collection);GridFunction oldfield(&oldspace);oldfield=0.;
 std::vector<int>elements(original.GetNE());std::iota(elements.begin(),elements.end(),0);
 pic_batch::HostSnapshot snapshot;snapshot.parameters.self_nx=snapshot.parameters.self_ny=8;snapshot.parameters.width=snapshot.parameters.height=.009;
 snapshot.self_cells=io::export_grid(oldfield,elements);
 ex::Layout uniform(snapshot,9,0);auto same=uniform.mesh();FiniteElementSpace same_space(&same,&collection);GridFunction same_field(&same_space);same_field=0.;
 const auto exported=io::export_grid(same_field,elements);io::require(std::memcmp(exported.data(),snapshot.self_cells.data(),exported.size()*sizeof(gather_push::Cell))==0,"Manufactured uniform topology mismatch");
 ex::Axis ratio_one(uniform.x.nodes,.009,8*(.009/8),8);io::require(ratio_one.growth==1,"Unit growth limit changed");
 ex::Axis nested_short(uniform.x.nodes,.009,.02,12,1.05),nested_long(uniform.x.nodes,.009,.04,24,1.05);
 for(int i=0;i<int(nested_short.nodes.size());++i)io::require(nested_short.nodes[i]==nested_long.nodes[i+12],"Nested interior coordinates changed");
 for(double x:{0.,.009,std::nextafter(.0045,0.),.0045,std::nextafter(.0045,1.)}){
  const auto u=uniform.x.locate(x);const double normalized=x/.009*8;const int index=std::min(int(normalized),7);
  io::require(u.cell==index&&u.xi==normalized-index,"Exact normalized locate changed");
 }
 bool invalid_rejected=false;try{ex::Axis invalid(uniform.x.nodes,.009,.00001,8);}catch(const std::exception&){invalid_rejected=true;}io::require(invalid_rejected,"Invalid exterior accepted");
 ex::Layout layout(snapshot,40,12);auto mesh=layout.mesh();layout.validate_mesh(mesh);
 FiniteElementSpace space(&mesh,&collection);GridFunction field(&space);
 FunctionCoefficient exact([](const Vector&p){return analytic(p[0],p[1]);});field.ProjectCoefficient(exact);
 Vector exact_nodes(field);ConstantCoefficient epsilon(8.8541878188e-12);BilinearForm form(&space);form.AddDomainIntegrator(new DiffusionIntegrator(epsilon));form.Assemble();form.Finalize();
 LinearForm rhs(&space);rhs=0.;Array<int>boundary(mesh.bdr_attributes.Max());boundary=1;Array<int>essential;space.GetEssentialTrueDofs(boundary,essential);
 OperatorPtr op;Vector x,b;form.FormLinearSystem(essential,field,rhs,op,x,b);auto*A=dynamic_cast<SparseMatrix*>(op.Ptr());io::require(A,"No sparse operator");
 GSSmoother preconditioner(*A);CGSolver solver;solver.SetOperator(*A);solver.SetPreconditioner(preconditioner);solver.SetRelTol(1e-14);solver.SetAbsTol(1e-30);solver.SetMaxIter(4000);solver.SetPrintLevel(-1);solver.Mult(b,x);
 io::require(solver.GetConverged(),"Manufactured solve did not converge");Vector residual(b.Size());A->Mult(x,residual);residual-=b;const double limit=std::max(1e-30,2e-13*b.Norml2());io::require(residual.Norml2()<=limit,"Manufactured original residual guard failed");
 form.RecoverFEMSolution(x,rhs,field);double max_node=0,max_value=0,max_gradient=0;
 for(int i=0;i<field.Size();++i)max_node=std::max(max_node,std::abs(field[i]-exact_nodes[i]));
 for(int j=0;j<layout.y.cells();++j)for(int i=0;i<layout.x.cells();++i)for(double a:{.17,.53,.89}){
  const double px=layout.x.nodes[i]+a*(layout.x.nodes[i+1]-layout.x.nodes[i]),py=layout.y.nodes[j]+.37*(layout.y.nodes[j+1]-layout.y.nodes[j]);
  const auto sample=ex::gather(layout,field,{px,py});max_value=std::max(max_value,std::abs(sample.value-analytic(px,py)));
  max_gradient=std::max({max_gradient,std::abs(sample.dx-(3+5*py)),std::abs(sample.dy-(-4+5*px))});
 }
 io::require(max_node<1e-10&&max_value<1e-10&&max_gradient<1e-7,"Manufactured bilinear field/gradient error");
 // Independent physical source work uses the analytic polynomial, not another
 // invocation of the candidate gather. Q1 should reproduce it on graded cells.
 LinearForm charge(&space);charge=0.;double direct=0,work=0,activity=0;
 for(int n=0;n<100;++n){double px=layout.x.nodes.front()+(n+.13)/100*(layout.x.nodes.back()-layout.x.nodes.front());double py=layout.y.nodes.front()+((37*n)%100+.27)/100*(layout.y.nodes.back()-layout.y.nodes.front());
  double q=(n%2?1:-1)*1.602176634e-19*(1+n%5);ex::deposit_point(layout,space,charge,{px,py,0},q);direct+=q;double value=q*analytic(px,py);work+=value;activity+=std::abs(value);
 }
 const double source_error=std::abs(charge.Sum()-direct),work_error=std::abs(exact_nodes*charge-work);
 io::require(source_error<1e-30&&work_error<=1024*std::numeric_limits<double>::epsilon()*activity,"Manufactured charge/work identity failed");
 const auto coarse=smooth_poisson(8,12),fine=smooth_poisson(16,24),finer=smooth_poisson(32,48);
 io::require(fine.first<.4*coarse.first&&finer.first<.4*fine.first&&fine.second<.7*coarse.second&&finer.second<.7*fine.second,"Smooth manufactured graded refinement failed");
 std::cout<<std::setprecision(17)<<"{\"smooth_Poisson_value_errors_V\":["<<coarse.first<<","<<fine.first<<","<<finer.first<<"],\"smooth_Poisson_gradient_errors_V_m\":["<<coarse.second<<","<<fine.second<<","<<finer.second<<"],\"smooth_Poisson_refinement_passed\":true,\"detail\":\"-epsilon Laplacian(exp(20x)cos(15y))=-175 epsilon exp(20x)cos(15y); exact outer Dirichlet; manufactured test only\"}\n";
 std::cout<<std::setprecision(17)<<"{\"status\":\"passed\",\"uniform_topology_bit_equal\":true,\"unit_growth_exact\":true,\"exact_core_locator_checked\":true,\"invalid_margin_rejected\":true,\"graded_nodes\":"<<space.GetVSize()<<",\"CG_iterations\":"<<solver.GetNumIterations()<<",\"residual\":"<<residual.Norml2()<<",\"original_residual_limit\":"<<limit<<",\"max_nodal_error_V\":"<<max_node<<",\"max_value_error_V\":"<<max_value<<",\"max_gradient_error_V_m\":"<<max_gradient<<",\"source_sum_error_C\":"<<source_error<<",\"source_work_error_J\":"<<work_error<<",\"GPU_or_PIC_evolution\":false}\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
