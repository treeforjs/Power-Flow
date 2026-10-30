#include "batch_snapshot.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <sys/resource.h>
#include "bench_profile.hpp"
#ifdef STATIC_DOMAIN_WITH_CUDSS
#include "cudss_poisson.hpp"
#endif
using namespace mfem;
using namespace raikou::mfem_backend;
namespace io=snapshot_io;
namespace gp=gather_push;
static_assert(sizeof(mfem::real_t)==8);
constexpr double epsilon0=8.8541878188e-12;
constexpr int max_nodes=1010000;
const double E_abs=1e-6, B_abs=1e-13, field_rel=5e-12;

struct FieldRecord {gp::Vec E,B;double phi,psi;};
static_assert(sizeof(FieldRecord)==64);
struct Residual {double absolute,relative,limit;bool passed;};
Residual residual(const SparseMatrix&A,const Vector&x,const Vector&b) {
 Vector r(b.Size());A.Mult(x,r);r-=b;
 double error=r.Norml2(),norm=b.Norml2(),limit=std::max(1e-30,2e-13*norm);
 return {error,norm>0?error/norm:error,limit,std::isfinite(error)&&error<=limit};
}
void vector_file(const std::filesystem::path&p,const Vector&x) {io::write(p,std::vector<double>(x.HostRead(),x.HostRead()+x.Size()));}
void csr_files(const std::filesystem::path&base,const SparseMatrix&A) {
 io::write(base.string()+"_I.bin",std::vector<int>(A.GetI(),A.GetI()+A.Height()+1));
 io::write(base.string()+"_J.bin",std::vector<int>(A.GetJ(),A.GetJ()+A.NumNonZeroElems()));
 io::write(base.string()+"_V.bin",std::vector<double>(A.GetData(),A.GetData()+A.NumNonZeroElems()));
}

struct StaticCase {
 GapPicOptions o;double shift;
 Mesh mesh;H1_FECollection collection;FiniteElementSpace space;
 GridFunction phi,psi,xfield,yfield;
 ConstantCoefficient epsilon;BilinearForm form;
 std::unique_ptr<SparseMatrix> full;
 LinearForm charge,unit_rhs;Vector raw_charge;Array<int> essential,outer_marker;
 OperatorPtr phi_operator,psi_operator;Vector phi_x,phi_b,psi_x,psi_b;
 double direct_charge=0,deposited_charge=0,independent_source_max_error=0;
 bool source_bitwise=false,topology_bitwise=false;
 std::vector<int> owner;

 static GapPicOptions options(const GapPicOptions&base,int mm) {
  io::require(mm==9||mm==13||mm==17||mm==25,"Only declared uniform domains are supported");
  auto o=base;const double width=mm==9?base.width_m:mm==13?.013:mm==17?.017:.025,offset=(width-base.width_m)/2;
  o.nx=o.ny=mm*40;o.width_m=o.gap_m=width;
  for(auto*f:{&o.lower_foil,&o.upper_foil}){f->left_m+=offset;f->right_m+=offset;f->bottom_m+=offset;f->top_m+=offset;}
  io::require((o.nx+1)*(o.ny+1)<=max_nodes,"Node ceiling exceeded");return o;
 }
 StaticCase(const GapPicOptions&base,int mm,const pic_batch::HostSnapshot&s):
  o(options(base,mm)),shift((o.width_m-base.width_m)/2),
  mesh(Mesh::MakeCartesian2D(o.nx,o.ny,Element::QUADRILATERAL,true,o.width_m,o.gap_m,false)),
  collection(1,2),space(&mesh,&collection),phi(&space),psi(&space),xfield(&space),yfield(&space),
  epsilon(epsilon0),form(&space),charge(&space),unit_rhs(&space) {
  pic_bench::Timer timer("geometry_assembly_and_RHS");
  phi=0;psi=0;charge=0;unit_rhs=0.;
  form.AddDomainIntegrator(new DiffusionIntegrator(epsilon));form.Assemble();form.Finalize();
  full=std::make_unique<SparseMatrix>(form.SpMat());
  FunctionCoefficient xc([](const Vector&x){return x[0];}),yc([](const Vector&x){return x[1];});
  xfield.ProjectCoefficient(xc);yfield.ProjectCoefficient(yc);
  Array<int>dofs;Vector shape;
  for(const auto&p:s.particles){IntegrationPoint ip;int e=locate(p,ip);const auto*fe=space.GetFE(e);shape.SetSize(fe->GetDof());fe->CalcShape(ip,shape);space.GetElementDofs(e,dofs);double q=p.charge_C*p.weight/o.depth_m;charge.AddElementVector(dofs,shape*=q);direct_charge+=q*o.depth_m;}
  deposited_charge=charge.Sum()*o.depth_m;
  raw_charge=charge;
  // Same initial marker, then finite-foil override, as the sealed solve_fields.
  outer_marker.SetSize(mesh.bdr_attributes.Max());outer_marker=0;Array<int>vertices;
  for(int b=0;b<mesh.GetNBE();++b){mesh.GetBdrElementVertices(b,vertices);double y0=mesh.GetVertex(vertices[0])[1],y1=mesh.GetVertex(vertices[1])[1];if(y0==y1&&(y0==0||y0==o.gap_m))outer_marker[mesh.GetBdrAttribute(b)-1]=1;}
  space.GetEssentialTrueDofs(outer_marker,essential);essential.SetSize(0);owner.assign(space.GetVSize(),-1);
  const double tol=1e-12*std::max(o.width_m,o.gap_m);
  for(int i=0;i<space.GetVSize();++i){int who=0;for(const auto&f:{o.lower_foil,o.upper_foil}){if(xfield[i]>=f.left_m-tol&&xfield[i]<=f.right_m+tol&&yfield[i]>=f.bottom_m-tol&&yfield[i]<=f.top_m+tol){owner[i]=who;essential.Append(i);}++who;}}
  io::require(essential.Size()==3618,"Foil Dirichlet geometry changed");
  trial(phi,o.upper_potential_V);trial(psi,1.);
  form.FormLinearSystem(essential,phi,charge,phi_operator,phi_x,phi_b);
  form.FormLinearSystem(essential,psi,unit_rhs,psi_operator,psi_x,psi_b);
  io::require(phi_operator.Ptr()==psi_operator.Ptr(),"Two RHSs do not share the identical constrained operator");
  if(mm==9) {
   std::vector<int>elements(o.nx*o.ny);std::iota(elements.begin(),elements.end(),0);
   const auto exported=io::export_grid(phi,elements);
   topology_bitwise=exported.size()==s.self_cells.size()&&std::memcmp(exported.data(),s.self_cells.data(),exported.size()*sizeof(gp::Cell))==0;
   io::require(topology_bitwise,"9mm topology/DOF bytes differ from capture");
   // Independent source route through the captured cell geometry and Q1 basis.
   Vector independent(charge.Size());independent=0.;
   std::vector<double>xs(o.nx+1),ys(o.ny+1);for(int i=0;i<o.nx;++i)xs[i]=s.self_cells[i].x[0];xs.back()=s.self_cells[o.nx-1].x[1];for(int j=0;j<o.ny;++j)ys[j]=s.self_cells[j*o.nx].y[0];ys.back()=s.self_cells[(o.ny-1)*o.nx].y[3];
   for(const auto&p:s.particles){int i=gp::bin(xs.data(),o.nx,p.position_m[0]),j=gp::bin(ys.data(),o.ny,p.position_m[1]);io::require(i>=0&&j>=0,"Independent particle outside grid");const auto&c=s.self_cells[j*o.nx+i];double x=(p.position_m[0]-xs[i])/(xs[i+1]-xs[i]),y=(p.position_m[1]-ys[j])/(ys[j+1]-ys[j]);double vx[2],vy[2];gp::basis1d(x,vx);gp::basis1d(y,vy);double weights[4]={vx[0]*vy[0],vx[1]*vy[0],vx[1]*vy[1],vx[0]*vy[1]};double q=p.charge_C*p.weight/o.depth_m;for(int k=0;k<4;++k)independent[c.dof[k]]+=weights[k]*q;}
   source_bitwise=std::memcmp(independent.HostRead(),raw_charge.HostRead(),charge.Size()*sizeof(double))==0;
   for(int i=0;i<charge.Size();++i)independent_source_max_error=std::max(independent_source_max_error,std::abs(independent[i]-raw_charge[i]));
  }
 }
 void trial(GridFunction&f,double voltage) {
  FunctionCoefficient bc([&](const Vector&x){return voltage*x[1]/o.gap_m;});f.ProjectBdrCoefficient(bc,outer_marker);
  for(int i=0;i<essential.Size();++i){int d=essential[i];f[d]=owner[d]==0?0:voltage;}
 }
 int locate(const GapParticle&p,IntegrationPoint&ip)const {
  const double x=p.position_m[0]+shift,y=p.position_m[1]+shift;
  io::require(x>=0&&x<=o.width_m&&y>=0&&y<=o.gap_m,"Translated particle outside static domain");
  double sx=x/o.width_m*o.nx,sy=y/o.gap_m*o.ny;int ix=std::min(int(sx),o.nx-1),iy=std::min(int(sy),o.ny-1);ip.Set2(sx-ix,sy-iy);return iy*o.nx+ix;
 }
 SparseMatrix&matrix(){auto*a=dynamic_cast<SparseMatrix*>(phi_operator.Ptr());io::require(a,"Expected sparse Poisson operator");return *a;}
 void load_capture(const pic_batch::HostSnapshot&s){io::require(phi.Size()==int(s.phi.size()),"Captured nodal vector mismatch");for(int i=0;i<phi.Size();++i){phi[i]=s.phi[i];psi[i]=s.psi[i];}}
 std::vector<FieldRecord> gather(const pic_batch::HostSnapshot&s) {
  std::vector<FieldRecord> result;result.reserve(s.particles.size());const auto&t=s.parameters;
  for(const auto&p:s.particles) {
   IntegrationPoint ip;const int e=locate(p,ip);auto*tr=mesh.GetElementTransformation(e);tr->SetIntPoint(&ip);Vector phigrad(2),psigrad(2);phi.GetGradient(*tr,phigrad);psi.GetGradient(*tr,psigrad);
   const double px=p.position_m[0]-.0045,py=p.position_m[1]-.0045;int ni=gp::bin(s.xs.data(),t.native_nx,px),nj=gp::bin(s.ys.data(),t.native_ny,py);io::require(ni>=0&&nj>=0,"Native coordinate mapping changed");int ne=nj*t.native_nx+ni;double nx=(px-s.xs[ni])/(s.xs[ni+1]-s.xs[ni]),ny=(py-s.ys[nj])/(s.ys[nj+1]-s.ys[nj]);gp::Grid grid{s.native_cells.data(),t.native_nx,t.native_ny};
   auto a0=gp::sample(grid,s.az0.data(),ne,nx,ny),a1=gp::sample(grid,s.az1.data(),ne,nx,ny),r0=gp::sample(grid,s.rate0.data(),ne,nx,ny),r1=gp::sample(grid,s.rate1.data(),ne,nx,ny);io::require(!a0.status&&!a1.status&&!r0.status&&!r1.status,"Invalid native gather");
   FieldRecord record{};record.phi=phi.GetValue(e,ip);record.psi=psi.GetValue(e,ip);record.E={{-phigrad[0],-phigrad[1],0}};record.E.v[0]+=0;record.E.v[1]+=0;record.E.v[2]+=(1-t.fraction)*(-r0.value)+t.fraction*(-r1.value);
   gp::Vec correction{{-psigrad[0],-psigrad[1],0}};double g=(p.position_m[2]-t.short_z)/t.depth;for(int k=0;k<2;++k)correction.v[k]*=t.voltage*(2*g-1);correction.v[2]=-2*t.voltage/t.depth*(record.psi-.5);for(int k=0;k<3;++k)record.E.v[k]+=correction.v[k];
   record.B={{(1-t.fraction)*a0.dy+t.fraction*a1.dy,(1-t.fraction)*(-a0.dx)+t.fraction*(-a1.dx),0}};
   io::require(gp::finite(record.E)&&gp::finite(record.B)&&std::isfinite(record.phi)&&std::isfinite(record.psi),"Nonfinite static gather");result.push_back(record);
  }
  return result;
 }
};

struct Compare {std::size_t failures=0;double maxE=0,maxB=0;};
Compare compare(const std::vector<FieldRecord>&a,const pic_batch::HostSnapshot&s) {
 Compare c;io::require(a.size()==s.reference.size(),"Comparison count changed");
 for(std::size_t i=0;i<a.size();++i)for(int k=0;k<3;++k){double e=std::abs(a[i].E.v[k]-s.reference[i].E.v[k]),b=std::abs(a[i].B.v[k]-s.reference[i].B.v[k]);c.maxE=std::max(c.maxE,e);c.maxB=std::max(c.maxB,b);if(e>E_abs+field_rel*std::abs(s.reference[i].E.v[k])||b>B_abs+field_rel*std::abs(s.reference[i].B.v[k]))++c.failures;}return c;
}

int main(int argc,char**argv) {
 try {
  io::require(argc==5,"Usage: static-domain cpu-preflight|gpu FIXTURE NEW_OUTPUT WIDTH_MM");
  std::string mode=argv[1];io::require(mode=="cpu-preflight"||mode=="gpu","Unknown mode");
  const std::filesystem::path fixture(argv[2]),out(argv[3]);int mm=std::stoi(argv[4]);
  io::require(mode!="cpu-preflight"||mm==9,"CPU preflight is the captured9mm control only");
  io::require(!std::filesystem::exists(out),"Output directory must be fresh");std::filesystem::create_directories(out);
  const auto snapshot=pic_batch::HostSnapshot::load(fixture,true);
  const auto original_particles=snapshot.particles;
  auto option_records=io::read<GapPicOptions>(fixture/"options.bin");io::require(option_records.size()==1,"Missing captured options");const auto&base=option_records[0];
  io::require(base.nx==360&&base.ny==360&&base.finite_foils&&base.depth_m==.025&&base.upper_potential_V==snapshot.parameters.voltage,"Unexpected fixture geometry");
  io::require(snapshot.particles.size()==46726,"Expected the bound actual46726-particle fixture");
  StaticCase c(base,mm,snapshot);Residual phi_r{},psi_r{};Compare control{};int factors=0,solves=0;
  if(mode=="cpu-preflight") {
   c.load_capture(snapshot);phi_r=residual(c.matrix(),c.phi,c.phi_b);psi_r=residual(c.matrix(),c.psi,c.psi_b);control=compare(c.gather(snapshot),snapshot);
   // CPU work validates the captured solution, never initializes CUDA or solves.
   csr_files(out/"full_operator",*c.full);csr_files(out/"constrained_operator",c.matrix());vector_file(out/"charge_rhs",c.raw_charge);vector_file(out/"phi_rhs",c.phi_b);vector_file(out/"psi_rhs",c.psi_b);
  } else {
#ifdef STATIC_DOMAIN_WITH_CUDSS
   const char*authorized=std::getenv("PIC_STATIC_DOMAIN_PARENT_GPU");io::require(authorized&&std::string(authorized)=="1","Parent GPU dispatch flag required");
   {raikou_pic_parallel::CudssPoisson solver(c.matrix());++factors;solver.Solve(c.phi_b,c.phi_x);++solves;phi_r=residual(c.matrix(),c.phi_x,c.phi_b);io::require(phi_r.passed,"Original FP64 residual guard failed for phi");solver.Solve(c.psi_b,c.psi_x);++solves;psi_r=residual(c.matrix(),c.psi_x,c.psi_b);io::require(psi_r.passed,"Original FP64 residual guard failed for psi");}
   c.form.RecoverFEMSolution(c.phi_x,c.charge,c.phi);c.form.RecoverFEMSolution(c.psi_x,c.unit_rhs,c.psi);
   auto fields=c.gather(snapshot);if(mm==9)control=compare(fields,snapshot);io::write(out/"particle_fields.bin",fields);vector_file(out/"phi_nodes.bin",c.phi);vector_file(out/"psi_nodes.bin",c.psi);
#else
   throw std::runtime_error("This CPU-only binary has no CUDA backend");
#endif
  }
  io::require(snapshot.particles.size()==original_particles.size()&&std::memcmp(snapshot.particles.data(),original_particles.data(),original_particles.size()*sizeof(GapParticle))==0,"Particle state was modified");
  const double charge_error=std::abs(c.deposited_charge-c.direct_charge);const bool charge_ok=charge_error<=1e-11*std::max(std::abs(c.direct_charge),pic_e);
  const double source_roundoff_limit=1024*std::numeric_limits<double>::epsilon()*std::abs(c.direct_charge)/base.depth_m;
  const bool source_ok=mm!=9||c.independent_source_max_error<=source_roundoff_limit;
  double max_phi_error=0,max_psi_error=0,phi_scale=0,psi_scale=0;
  if(mm==9)for(int i=0;i<c.phi.Size();++i){max_phi_error=std::max(max_phi_error,std::abs(c.phi[i]-snapshot.phi[i]));max_psi_error=std::max(max_psi_error,std::abs(c.psi[i]-snapshot.psi[i]));phi_scale=std::max(phi_scale,std::abs(snapshot.phi[i]));psi_scale=std::max(psi_scale,std::abs(snapshot.psi[i]));}
  const double phi_limit=(1e-10+256*std::numeric_limits<double>::epsilon())*phi_scale;
  const double psi_limit=(1e-10+256*std::numeric_limits<double>::epsilon())*psi_scale;
  const bool nodal_ok=mm!=9||(max_phi_error<=phi_limit&&max_psi_error<=psi_limit);
  bool pass=phi_r.passed&&psi_r.passed&&charge_ok&&source_ok&&nodal_ok&&(mm!=9||control.failures==0);
  std::ofstream result(out/"result.json");result<<std::setprecision(17)<<"{\"mode\":\""<<mode<<"\",\"case_full_width_mm\":"<<mm<<",\"nodes\":"<<c.space.GetVSize()<<",\"cells\":"<<c.mesh.GetNE()<<",\"nnz\":"<<c.matrix().NumNonZeroElems()<<",\"foil_dirichlet_nodes\":"<<c.essential.Size()<<",\"particles\":"<<snapshot.particles.size()<<",\"particle_state_unchanged\":true,\"source_topology_bitwise_equal_to_capture\":"<<(c.topology_bitwise?"true":"false")<<",\"independent_source_bitwise_equal\":"<<(c.source_bitwise?"true":"false")<<",\"independent_source_max_abs_error_C_per_m\":"<<c.independent_source_max_error<<",\"direct_charge_C\":"<<c.direct_charge<<",\"deposited_charge_C\":"<<c.deposited_charge<<",\"charge_error_C\":"<<charge_error<<",\"phi_residual_C_per_m\":"<<phi_r.absolute<<",\"phi_residual_limit_C_per_m\":"<<phi_r.limit<<",\"psi_residual_C_per_m\":"<<psi_r.absolute<<",\"psi_residual_limit_C_per_m\":"<<psi_r.limit<<",\"captured_field_failures\":"<<control.failures<<",\"max_E_difference_V_m\":"<<control.maxE<<",\"max_B_difference_T\":"<<control.maxB<<",\"factorizations\":"<<factors<<",\"field_solves\":"<<solves<<",\"same_operator_for_two_RHS\":true,\"graded_mesh\":false,\"no_particle_advance\":true,\"passed\":"<<(pass?"true":"false")<<"}\n";result.close();io::require(bool(result),"Result write failed");
  std::ofstream gates(out/"control-gates.json");
  gates<<std::setprecision(17)<<"{\"independent_source_checked\":"<<(mm==9?"true":"false")<<",\"source_roundoff_limit_C_per_m\":"<<source_roundoff_limit<<",\"source_gate_passed\":"<<(source_ok?"true":"false")<<",\"captured_nodal_comparison_performed\":"<<(mm==9?"true":"false")<<",\"max_phi_error_V\":"<<max_phi_error<<",\"phi_error_limit_V\":"<<phi_limit<<",\"max_psi_error\":"<<max_psi_error<<",\"psi_error_limit\":"<<psi_limit<<",\"nodal_gate_passed\":"<<(nodal_ok?"true":"false")<<",\"capture_field_comparison_performed\":"<<(mm==9?"true":"false")<<",\"original_E_abs_V_m\":1e-6,\"original_B_abs_T\":1e-13,\"original_field_relative\":5e-12}\n";
  gates.close();io::require(bool(gates),"Gate receipt write failed");
  pic_bench::resources((out/"resources.json").string());pic_bench::write((out/"timers.json").string());
  std::cout<<"Static field-only "<<mode<<" "<<mm<<"mm control "<<(pass?"passed":"failed")<<'\n';return pass?0:3;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}
}
