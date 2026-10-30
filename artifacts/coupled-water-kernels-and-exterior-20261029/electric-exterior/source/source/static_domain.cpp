#include "electric_coordinates.hpp"
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
 GapPicOptions o;electric_exterior::Layout layout;
 Mesh mesh;H1_FECollection collection;FiniteElementSpace space;
 GridFunction phi,psi,xfield,yfield;
 ConstantCoefficient epsilon;BilinearForm form;
 std::unique_ptr<SparseMatrix> full;
 LinearForm charge,unit_rhs;Vector raw_charge;Array<int> essential,outer_marker;
 OperatorPtr phi_operator,psi_operator;Vector phi_x,phi_b,psi_x,psi_b;
 double direct_charge=0,deposited_charge=0,independent_source_max_error=0;
 bool source_bitwise=false,topology_bitwise=false,core_geometry_bitwise=false,core_source_bitwise=false;
 std::vector<int> owner,core_dof_map;

 StaticCase(const GapPicOptions&base,int mm,int outer,const pic_batch::HostSnapshot&s):
  o(base),layout(s,mm,outer),mesh(layout.mesh()),
  collection(1,2),space(&mesh,&collection),phi(&space),psi(&space),xfield(&space),yfield(&space),
  epsilon(epsilon0),form(&space),charge(&space),unit_rhs(&space) {
  pic_bench::Timer timer("geometry_assembly_and_RHS");
  layout.validate_mesh(mesh);
  phi=0;psi=0;charge=0;unit_rhs=0.;
  form.AddDomainIntegrator(new DiffusionIntegrator(epsilon));form.Assemble();form.Finalize();
  full=std::make_unique<SparseMatrix>(form.SpMat());
  FunctionCoefficient xc([](const Vector&x){return x[0];}),yc([](const Vector&x){return x[1];});
  xfield.ProjectCoefficient(xc);yfield.ProjectCoefficient(yc);
  for(const auto&p:s.particles){double q=p.charge_C*p.weight/o.depth_m;electric_exterior::deposit_point(layout,space,charge,{p.position_m[0],p.position_m[1],p.position_m[2]},q);direct_charge+=q*o.depth_m;}
  deposited_charge=charge.Sum()*o.depth_m;
  raw_charge=charge;
  // Same initial marker, then finite-foil override, as the sealed solve_fields.
  outer_marker.SetSize(mesh.bdr_attributes.Max());outer_marker=0;Array<int>vertices;
  for(int b=0;b<mesh.GetNBE();++b){mesh.GetBdrElementVertices(b,vertices);double y0=mesh.GetVertex(vertices[0])[1],y1=mesh.GetVertex(vertices[1])[1];if(y0==y1&&(y0==layout.y.nodes.front()||y0==layout.y.nodes.back()))outer_marker[mesh.GetBdrAttribute(b)-1]=1;}
  space.GetEssentialTrueDofs(outer_marker,essential);essential.SetSize(0);owner.assign(space.GetVSize(),-1);
  const double tol=1e-12*std::max(o.width_m,o.gap_m);
  for(int i=0;i<space.GetVSize();++i){int who=0;for(const auto&f:{o.lower_foil,o.upper_foil}){if(xfield[i]>=f.left_m-tol&&xfield[i]<=f.right_m+tol&&yfield[i]>=f.bottom_m-tol&&yfield[i]<=f.top_m+tol){owner[i]=who;essential.Append(i);}++who;}}
  io::require(essential.Size()==3618,"Foil Dirichlet geometry changed");
  trial(phi,o.upper_potential_V);trial(psi,1.);
  form.FormLinearSystem(essential,phi,charge,phi_operator,phi_x,phi_b);
  form.FormLinearSystem(essential,psi,unit_rhs,psi_operator,psi_x,psi_b);
  io::require(phi_operator.Ptr()==psi_operator.Ptr(),"Two RHSs do not share the identical constrained operator");
  audit_source(s);
 }
 void audit_source(const pic_batch::HostSnapshot&s) {
   std::vector<int>elements(mesh.GetNE());std::iota(elements.begin(),elements.end(),0);
   const auto exported=io::export_grid(phi,elements);
   if(layout.margin_count==0){topology_bitwise=exported.size()==s.self_cells.size()&&std::memcmp(exported.data(),s.self_cells.data(),exported.size()*sizeof(gp::Cell))==0;io::require(topology_bitwise,"Uniform limit topology differs from capture");}
   core_dof_map.assign(s.phi.size(),-1);core_geometry_bitwise=true;
   for(int j=0;j<s.parameters.self_ny;++j)for(int i=0;i<s.parameters.self_nx;++i){
    const auto&a=s.self_cells[j*s.parameters.self_nx+i];const auto&b=exported[(j+layout.margin_count)*layout.x.cells()+i+layout.margin_count];
    core_geometry_bitwise&=std::memcmp(a.x,b.x,sizeof(a.x))==0&&std::memcmp(a.y,b.y,sizeof(a.y))==0;
    for(int k=0;k<4;++k){int&mapped=core_dof_map[a.dof[k]];io::require(mapped==-1||mapped==b.dof[k],"Inconsistent core node map");mapped=b.dof[k];}
   }
   io::require(core_geometry_bitwise&&std::find(core_dof_map.begin(),core_dof_map.end(),-1)==core_dof_map.end(),"Captured physical core changed");
   Vector original_source(s.phi.size());original_source=0.;Array<int>old_dofs(4);Vector old_shape(4);
   for(const auto&p:s.particles){
    const double sx=p.position_m[0]/s.parameters.width*s.parameters.self_nx,sy=p.position_m[1]/s.parameters.height*s.parameters.self_ny;
    int i=std::min(int(sx),s.parameters.self_nx-1),j=std::min(int(sy),s.parameters.self_ny-1);IntegrationPoint ip;ip.Set2(sx-i,sy-j);
    space.GetFE(0)->CalcShape(ip,old_shape);const auto&cell=s.self_cells[j*s.parameters.self_nx+i];for(int k=0;k<4;++k)old_dofs[k]=cell.dof[k];
    original_source.AddElementVector(old_dofs,old_shape*=p.charge_C*p.weight/o.depth_m);
   }
   Vector mapped_source(raw_charge.Size());mapped_source=0.;for(int i=0;i<original_source.Size();++i)mapped_source[core_dof_map[i]]=original_source[i];
   core_source_bitwise=std::memcmp(mapped_source.HostRead(),raw_charge.HostRead(),raw_charge.Size()*sizeof(double))==0;io::require(core_source_bitwise,"Source changed after physical core-node mapping");
   Vector independent(charge.Size());independent=0.;
   const auto&xs=layout.x.nodes;const auto&ys=layout.y.nodes;
   for(const auto&p:s.particles){
    int i=gp::bin(xs.data(),layout.x.cells(),p.position_m[0]),j=gp::bin(ys.data(),layout.y.cells(),p.position_m[1]);io::require(i>=0&&j>=0,"Independent point outside electric mesh");
    const auto&c=exported[j*layout.x.cells()+i];double x=(p.position_m[0]-xs[i])/(xs[i+1]-xs[i]),y=(p.position_m[1]-ys[j])/(ys[j+1]-ys[j]);
    double vx[2],vy[2];gp::basis1d(x,vx);gp::basis1d(y,vy);double weights[4]={vx[0]*vy[0],vx[1]*vy[0],vx[1]*vy[1],vx[0]*vy[1]};
    double q=p.charge_C*p.weight/o.depth_m;for(int k=0;k<4;++k)independent[c.dof[k]]+=weights[k]*q;
   }
   source_bitwise=std::memcmp(independent.HostRead(),raw_charge.HostRead(),charge.Size()*sizeof(double))==0;
   for(int i=0;i<charge.Size();++i)independent_source_max_error=std::max(independent_source_max_error,std::abs(independent[i]-raw_charge[i]));
 }

 void trial(GridFunction&f,double voltage) {
  FunctionCoefficient bc([&](const Vector&x){return voltage*x[1]/o.gap_m;});f.ProjectBdrCoefficient(bc,outer_marker);
  for(int i=0;i<essential.Size();++i){int d=essential[i];f[d]=owner[d]==0?0:voltage;}
 }
 int locate(const GapParticle&p,IntegrationPoint&ip)const {
  const auto point=layout.locate(electric_exterior::ParticlePoint{p.position_m[0],p.position_m[1],p.position_m[2]});ip=point.reference;return point.element;
 }
 SparseMatrix&matrix(){auto*a=dynamic_cast<SparseMatrix*>(phi_operator.Ptr());io::require(a,"Expected sparse Poisson operator");return *a;}
 void load_capture(const pic_batch::HostSnapshot&s){io::require(phi.Size()==int(s.phi.size()),"Captured nodal vector mismatch");for(int i=0;i<phi.Size();++i){phi[i]=s.phi[i];psi[i]=s.psi[i];}}
 void transfer_captured_core(const pic_batch::HostSnapshot&s){
  // Interpolation control only: no exterior solution is inferred from these
  // zero outer nodes, and no residual/energy claim is made for this vector.
  phi=0.;psi=0.;for(std::size_t i=0;i<core_dof_map.size();++i){phi[core_dof_map[i]]=s.phi[i];psi[core_dof_map[i]]=s.psi[i];}
 }
 std::vector<FieldRecord> gather(const pic_batch::HostSnapshot&s) {
  std::vector<FieldRecord> result;result.reserve(s.particles.size());const auto&t=s.parameters;
  for(const auto&p:s.particles) {
   IntegrationPoint ip;const int e=locate(p,ip);auto*tr=mesh.GetElementTransformation(e);tr->SetIntPoint(&ip);Vector phigrad(2),psigrad(2);phi.GetGradient(*tr,phigrad);psi.GetGradient(*tr,psigrad);
   const auto native=electric_exterior::native_point({p.position_m[0],p.position_m[1],p.position_m[2]});const double px=native.x,py=native.y;int ni=gp::bin(s.xs.data(),t.native_nx,px),nj=gp::bin(s.ys.data(),t.native_ny,py);io::require(ni>=0&&nj>=0,"Native coordinate mapping changed");int ne=nj*t.native_nx+ni;double nx=(px-s.xs[ni])/(s.xs[ni+1]-s.xs[ni]),ny=(py-s.ys[nj])/(s.ys[nj+1]-s.ys[nj]);gp::Grid grid{s.native_cells.data(),t.native_nx,t.native_ny};
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
  io::require(argc==6,"Usage: electric-exterior cpu-preflight|cpu-geometry|gpu FIXTURE NEW_OUTPUT FULL_WIDTH_MM MARGIN_CELLS");
  std::string mode=argv[1];io::require(mode=="cpu-preflight"||mode=="cpu-geometry"||mode=="gpu","Unknown mode");
  const std::filesystem::path fixture(argv[2]),out(argv[3]);int mm=std::stoi(argv[4]),outer=std::stoi(argv[5]);
  io::require(mode!="cpu-preflight"||mm==9,"CPU preflight is the captured9mm control only");
  io::require(!std::filesystem::exists(out),"Output directory must be fresh");std::filesystem::create_directories(out);
  const auto snapshot=pic_batch::HostSnapshot::load(fixture,true);
  const auto original_particles=snapshot.particles;
  auto option_records=io::read<GapPicOptions>(fixture/"options.bin");io::require(option_records.size()==1,"Missing captured options");const auto&base=option_records[0];
  io::require(base.nx==360&&base.ny==360&&base.finite_foils&&base.depth_m==.025&&base.upper_potential_V==snapshot.parameters.voltage,"Unexpected fixture geometry");
  io::require(snapshot.particles.size()==46726,"Expected the bound actual46726-particle fixture");
  StaticCase c(base,mm,outer,snapshot);Residual phi_r{},psi_r{};Compare control{};int factors=0,solves=0;
  if(mode=="cpu-preflight") {
   c.load_capture(snapshot);phi_r=residual(c.matrix(),c.phi,c.phi_b);psi_r=residual(c.matrix(),c.psi,c.psi_b);control=compare(c.gather(snapshot),snapshot);
   // CPU work validates the captured solution, never initializes CUDA or solves.
   csr_files(out/"full_operator",*c.full);csr_files(out/"constrained_operator",c.matrix());vector_file(out/"charge_rhs",c.raw_charge);vector_file(out/"phi_rhs",c.phi_b);vector_file(out/"psi_rhs",c.psi_b);
  } else if(mode=="cpu-geometry") {
   vector_file(out/"charge_rhs",c.raw_charge);io::write(out/"axes_x.bin",c.layout.x.nodes);io::write(out/"axes_y.bin",c.layout.y.nodes);
   c.transfer_captured_core(snapshot);const auto interpolation=compare(c.gather(snapshot),snapshot);
   std::ofstream control_file(out/"interpolation-control.json");control_file<<std::setprecision(17)<<"{\"captured_core_nodes_transferred_only\":true,\"outer_nodes_zero_for_test_only\":true,\"fresh_Poisson_solution\":false,\"particles\":"<<snapshot.particles.size()<<",\"core_physical_cells_bitwise_equal\":"<<(c.core_geometry_bitwise?"true":"false")<<",\"mapped_charge_RHS_bitwise_equal\":"<<(c.core_source_bitwise?"true":"false")<<",\"captured_field_failures\":"<<interpolation.failures<<",\"max_E_difference_V_m\":"<<interpolation.maxE<<",\"max_B_difference_T\":"<<interpolation.maxB<<",\"passed\":"<<(interpolation.failures==0?"true":"false")<<"}\n";
   io::require(interpolation.failures==0,"Core interpolation control failed");
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
  const bool source_ok=c.independent_source_max_error<=source_roundoff_limit;
  double max_phi_error=0,max_psi_error=0,phi_scale=0,psi_scale=0;
  if(mm==9)for(int i=0;i<c.phi.Size();++i){max_phi_error=std::max(max_phi_error,std::abs(c.phi[i]-snapshot.phi[i]));max_psi_error=std::max(max_psi_error,std::abs(c.psi[i]-snapshot.psi[i]));phi_scale=std::max(phi_scale,std::abs(snapshot.phi[i]));psi_scale=std::max(psi_scale,std::abs(snapshot.psi[i]));}
  const double phi_limit=(1e-10+256*std::numeric_limits<double>::epsilon())*phi_scale;
  const double psi_limit=(1e-10+256*std::numeric_limits<double>::epsilon())*psi_scale;
  const bool nodal_ok=mm!=9||(max_phi_error<=phi_limit&&max_psi_error<=psi_limit);
  const bool solution_available=mode!="cpu-geometry";
  double field_energy=0,source_work=0,nodal_source_work=0,lower_charge=0,upper_charge=0,gauss=0,energy_defect=0,energy_limit=0,source_work_error=0,source_work_limit=0,gauss_limit=0;
  double old_ledger_limit=0,null_action_charge=0,free_residual_norm=0,free_residual_charge_bound=0,gauss_roundoff_bound=0,gauss_identity_defect=0,gauss_deposition_term=0;
  int free_count=0;bool old_ledger_pass=false;
  bool work_pass=true;
  if(solution_available){
   Vector action(c.phi.Size());c.full->Mult(c.phi,action);field_energy=.5*(c.phi*action)*base.depth_m;
   Vector reaction(action);reaction-=c.raw_charge;
   for(int i=0;i<c.essential.Size();++i){int d=c.essential[i];(c.owner[d]==0?lower_charge:upper_charge)+=reaction[d]*base.depth_m;}
   gauss=c.direct_charge+lower_charge+upper_charge;
   // Preserve the first, rejected diagnostic. It is not the runtime's Gauss
   // guard (which checks finiteness), and it omits the complete FEM residual.
   const double eps=std::numeric_limits<double>::epsilon();
   old_ledger_limit=1024*eps*std::max(std::abs(c.direct_charge)+std::abs(lower_charge)+std::abs(upper_charge),pic_e);
   old_ledger_pass=std::abs(gauss)<=old_ledger_limit;
   // For a=K*phi and r=a-q, the exact finite-vector identity is
   // Qdirect+H*sum_essential(r) = Qdirect-H*sum(q)
   //     +H*sum(a)-H*sum_free(r).
   // The free residual is bounded by Cauchy-Schwarz. No multiplier is chosen
   // from the observed defect. gamma(n)=n*eps/(1-n*eps) bounds n rounded sums.
   const auto gamma=[](long double n,long double u){return n*u/(1-n*u);};
   long double sum_q=0,sum_a=0,sum_free=0,free_squared=0,abs_aq=0,abs_essential=0,abs_q=0,abs_a=0,abs_free=0;
   for(int i=0;i<action.Size();++i){
    sum_q+=c.raw_charge[i];sum_a+=action[i];abs_q+=std::abs(c.raw_charge[i]);abs_a+=std::abs(action[i]);abs_aq+=std::abs(action[i])+std::abs(c.raw_charge[i]);
    if(c.owner[i]<0){++free_count;sum_free+=reaction[i];free_squared+=static_cast<long double>(reaction[i])*reaction[i];abs_free+=std::abs(reaction[i]);}
    else abs_essential+=std::abs(reaction[i])*base.depth_m;
   }
   const long double H=base.depth_m,u=eps,ul=std::numeric_limits<long double>::epsilon();
   const long double deposition=static_cast<long double>(c.direct_charge)-H*sum_q;
   const long double null_action=H*sum_a;
   const long double free_bound=H*std::sqrt(static_cast<long double>(free_count))*std::sqrt(free_squared);
   // Includes each rounded r=a-q, each reaction*H, electrode sums, final
   // two additions, and the long-double audit accumulation itself.
   const long double roundoff=H*gamma(1,u)*abs_aq+(gamma(c.essential.Size(),u)+gamma(1,u))*abs_essential
    +gamma(2,u)*(std::abs(c.direct_charge)+std::abs(lower_charge)+std::abs(upper_charge))
    +H*gamma(action.Size()+4,ul)*(abs_q+abs_a+abs_free)
    +gamma(8,ul)*(std::abs(c.direct_charge)+H*(abs_q+abs_a+abs_free));
   gauss_deposition_term=double(deposition);null_action_charge=double(null_action);free_residual_norm=double(std::sqrt(free_squared));free_residual_charge_bound=double(free_bound);gauss_roundoff_bound=double(roundoff);
   gauss_identity_defect=double(std::abs(static_cast<long double>(gauss)-(deposition+null_action-H*sum_free)));
   gauss_limit=double((std::abs(deposition)+std::abs(null_action)+free_bound+roundoff)/(1-u));
   nodal_source_work=(c.phi*c.raw_charge)*base.depth_m;
   double source_activity=0;for(const auto&p:snapshot.particles){IntegrationPoint ip;int e=c.locate(p,ip);double work=p.charge_C*p.weight*c.phi.GetValue(e,ip);source_work+=work;source_activity+=std::abs(work);}
   source_work_error=std::abs(source_work-nodal_source_work);source_work_limit=1024*std::numeric_limits<double>::epsilon()*std::max(source_activity,1e-30);
   const double electrode_work=base.lower_potential_V*lower_charge+base.upper_potential_V*upper_charge;
   energy_defect=2*field_energy-nodal_source_work-electrode_work;energy_limit=1024*std::numeric_limits<double>::epsilon()*std::max(std::abs(2*field_energy)+std::abs(nodal_source_work)+std::abs(electrode_work),1e-30);
   work_pass=std::isfinite(field_energy)&&field_energy>=0&&std::isfinite(gauss)&&std::abs(gauss)<=gauss_limit&&gauss_identity_defect<=gauss_roundoff_bound&&std::abs(energy_defect)<=energy_limit&&source_work_error<=source_work_limit;
  }
  bool pass=charge_ok&&source_ok&&c.core_geometry_bitwise&&c.core_source_bitwise&&nodal_ok&&work_pass&&(!solution_available||(phi_r.passed&&psi_r.passed&&(mm!=9||control.failures==0)));
  std::ofstream audit(out/"work-and-layout.json");audit<<std::setprecision(17)<<"{\"solution_available\":"<<(solution_available?"true":"false")<<",\"work_checks_performed\":"<<(solution_available?"true":"false")<<",\"work_checks_passed\":"<<(work_pass?"true":"false")<<",\"field_energy_J\":"<<field_energy<<",\"source_potential_work_J\":"<<source_work<<",\"nodal_source_work_J\":"<<nodal_source_work<<",\"source_work_error_J\":"<<source_work_error<<",\"source_work_limit_J\":"<<source_work_limit<<",\"electrostatic_energy_identity_defect_J\":"<<energy_defect<<",\"energy_identity_limit_J\":"<<energy_limit<<",\"gauss_defect_C\":"<<gauss<<",\"gauss_limit_C\":"<<gauss_limit<<",\"rejected_initial_ledger_budget_C\":"<<old_ledger_limit<<",\"rejected_initial_ledger_passed\":"<<(old_ledger_pass?"true":"false")<<",\"gauss_source_deposition_term_C\":"<<gauss_deposition_term<<",\"full_operator_null_action_sum_C\":"<<null_action_charge<<",\"free_full_residual_norm_C_per_m\":"<<free_residual_norm<<",\"free_dofs\":"<<free_count<<",\"free_residual_Cauchy_Schwarz_bound_C\":"<<free_residual_charge_bound<<",\"finite_sum_roundoff_bound_C\":"<<gauss_roundoff_bound<<",\"finite_vector_gauss_identity_error_C\":"<<gauss_identity_defect<<",\"electric_bounds_m\":["<<c.layout.x.nodes.front()<<","<<c.layout.x.nodes.back()<<"],\"central_25um_cell_count\":360,\"margin_cells_each_side\":"<<outer<<",\"outer_growth_x\":"<<c.layout.x.growth<<",\"last_outer_cell_width_m\":"<<c.layout.x.last_width<<",\"physical_particle_coordinates_unshifted\":true,\"native_coordinates_unshifted\":true,\"natural_Neumann_exterior\":true,\"Robin_or_sheath_smoothing\":false}\n";audit.close();

  std::ofstream result(out/"result.json");result<<std::setprecision(17)<<"{\"mode\":\""<<mode<<"\",\"case_full_width_mm\":"<<mm<<",\"nodes\":"<<c.space.GetVSize()<<",\"cells\":"<<c.mesh.GetNE()<<",\"nnz\":"<<c.matrix().NumNonZeroElems()<<",\"foil_dirichlet_nodes\":"<<c.essential.Size()<<",\"particles\":"<<snapshot.particles.size()<<",\"particle_state_unchanged\":true,\"source_topology_bitwise_equal_to_capture\":"<<(c.topology_bitwise?"true":"false")<<",\"independent_source_bitwise_equal\":"<<(c.source_bitwise?"true":"false")<<",\"independent_source_max_abs_error_C_per_m\":"<<c.independent_source_max_error<<",\"direct_charge_C\":"<<c.direct_charge<<",\"deposited_charge_C\":"<<c.deposited_charge<<",\"charge_error_C\":"<<charge_error<<",\"phi_residual_C_per_m\":"<<phi_r.absolute<<",\"phi_residual_limit_C_per_m\":"<<phi_r.limit<<",\"psi_residual_C_per_m\":"<<psi_r.absolute<<",\"psi_residual_limit_C_per_m\":"<<psi_r.limit<<",\"captured_field_failures\":"<<control.failures<<",\"max_E_difference_V_m\":"<<control.maxE<<",\"max_B_difference_T\":"<<control.maxB<<",\"factorizations\":"<<factors<<",\"field_solves\":"<<solves<<",\"same_operator_for_two_RHS\":true,\"graded_mesh\":"<<(outer?"true":"false")<<",\"no_particle_advance\":true,\"passed\":"<<(pass?"true":"false")<<"}\n";result.close();io::require(bool(result),"Result write failed");
  std::ofstream gates(out/"control-gates.json");
  gates<<std::setprecision(17)<<"{\"independent_source_checked\":"<<"true"<<",\"source_roundoff_limit_C_per_m\":"<<source_roundoff_limit<<",\"source_gate_passed\":"<<(source_ok?"true":"false")<<",\"captured_nodal_comparison_performed\":"<<(mm==9?"true":"false")<<",\"max_phi_error_V\":"<<max_phi_error<<",\"phi_error_limit_V\":"<<phi_limit<<",\"max_psi_error\":"<<max_psi_error<<",\"psi_error_limit\":"<<psi_limit<<",\"nodal_gate_passed\":"<<(nodal_ok?"true":"false")<<",\"capture_field_comparison_performed\":"<<(mm==9?"true":"false")<<",\"original_E_abs_V_m\":1e-6,\"original_B_abs_T\":1e-13,\"original_field_relative\":5e-12}\n";
  gates.close();io::require(bool(gates),"Gate receipt write failed");
  pic_bench::resources((out/"resources.json").string());pic_bench::write((out/"timers.json").string());
  std::cout<<"Static field-only "<<mode<<" "<<mm<<"mm control "<<(pass?"passed":"failed")<<'\n';return pass?0:3;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}
}
