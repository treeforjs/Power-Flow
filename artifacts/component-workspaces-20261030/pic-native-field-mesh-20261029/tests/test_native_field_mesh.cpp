// CPU-only field geometry/operator/gather/checkpoint audit. This translation
// unit includes the owned native implementation; no CUDA backend is compiled.
#include "../source/src/mfem_gap_pic.cpp"
#include "raikou/json.hpp"
#include <iostream>
#include <numeric>
#include <sys/resource.h>
#include <sched.h>
using namespace mfem;
using namespace raikou::mfem_backend;
using Json=nlohmann::json;
namespace fs=std::filesystem;
namespace {
Json checks=Json::array();
void check(bool ok,const std::string& name){if(!ok)throw std::runtime_error(name);checks.push_back(name);}
template<class F>void rejects(F f,const char*name){bool rejected=false;try{f();}catch(const std::exception&){rejected=true;}check(rejected,name);}
template<class T>std::vector<T> read(const fs::path&p){std::ifstream f(p,std::ios::binary|std::ios::ate);if(!f)throw std::runtime_error("Missing input "+p.string());auto n=f.tellg();if(n<0||n%sizeof(T))throw std::runtime_error("Truncated input");std::vector<T>a(std::size_t(n)/sizeof(T));f.seekg(0);f.read(reinterpret_cast<char*>(a.data()),n);if(!f)throw std::runtime_error("Input read failed");return a;}
template<class T>void write(const fs::path&p,const std::vector<T>&a){std::ofstream f(p,std::ios::binary);f.write(reinterpret_cast<const char*>(a.data()),a.size()*sizeof(T));if(!f)throw std::runtime_error("Output write failed");}
void vector_file(const fs::path&p,const Vector&v){write(p,std::vector<double>(v.HostRead(),v.HostRead()+v.Size()));}
void csr(const fs::path&p,const SparseMatrix&a){write(p.string()+"_I.bin",std::vector<int>(a.GetI(),a.GetI()+a.Height()+1));write(p.string()+"_J.bin",std::vector<int>(a.GetJ(),a.GetJ()+a.NumNonZeroElems()));write(p.string()+"_V.bin",std::vector<double>(a.GetData(),a.GetData()+a.NumNonZeroElems()));}
bool same(double a,double b){return std::memcmp(&a,&b,sizeof(a))==0;}
GapPicOptions small(){GapPicOptions o;o.nx=o.ny=40;o.width_m=o.gap_m=.001;o.depth_m=.025;o.finite_foils=true;o.finite_axial_extent=true;o.lower_foil={.0002,.0008,.0002,.000225};o.upper_foil={.0002,.0008,.000775,.0008};o.upper_potential_V=-100;return o;}
}
namespace raikou::mfem_backend {
struct PicFieldMeshAudit {
 static void particles(GapPic&p,const std::vector<GapParticle>&a){p.particles_=a;}
 static Vector charge(GapPic&p,double&direct){LinearForm rhs(&p.space_);direct=p.assemble_field_charge(rhs);return Vector(rhs);}
 static void nodes(GapPic&p,const std::vector<double>&a){p.potential_=0.;for(std::size_t i=0;i<a.size();++i)p.potential_[p.core_field_dofs_.at(i)]=a[i];p.fields_current_=true;}
 static void drift(GapPic&p,double dt){p.drift(dt);}
 static bool transport_equal(const GapPic&a,const GapPic&b){return a.particles_.size()==b.particles_.size()&&std::memcmp(a.particles_.data(),b.particles_.data(),a.particles_.size()*sizeof(GapParticle))==0&&a.horizontal_transfers_==b.horizontal_transfers_&&a.vertical_transfers_==b.vertical_transfers_&&std::memcmp(&a.ledger_,&b.ledger_,sizeof(a.ledger_))==0;}
 static GapPicDiagnostics ledger(const GapPic&p){return p.ledger_;}
 static void move_first_vertex(GapPic&p,double delta){p.mesh_.GetVertex(0)[0]+=delta;}
 static void assemble_control(GapPic&p,const fs::path&out){
  fs::create_directories(out);p.poisson_=std::make_unique<PicPoissonState>(p.space_,p.potential_);auto&form=p.poisson_->form;
  LinearForm rhs(&p.space_);p.assemble_field_charge(rhs);Vector full_rhs(rhs);Array<int>essential;for(int i=0;i<int(p.field_owner_identity_.size());++i)if(p.field_owner_identity_[i]>=0)essential.Append(i);
  Array<int>marker(p.mesh_.bdr_attributes.Max());marker=0;Array<int>vertices;
  for(int i=0;i<p.mesh_.GetNBE();++i){p.mesh_.GetBdrElementVertices(i,vertices);double a=p.mesh_.GetVertex(vertices[0])[1],b=p.mesh_.GetVertex(vertices[1])[1];if(a==b&&(a==p.field_layout_.ys.front()||a==p.field_layout_.ys.back()))marker[p.mesh_.GetBdrAttribute(i)-1]=1;}
  GridFunction phi(&p.space_),psi(&p.space_);phi=0.;psi=0.;
  const auto trial=[&](GridFunction&g,double voltage){FunctionCoefficient bc([&](const Vector&x){return voltage*x[1]/p.options_.gap_m;});g.ProjectBdrCoefficient(bc,marker);for(int i=0;i<essential.Size();++i){int d=essential[i];g[d]=p.field_owner_identity_[d]==0?0:voltage;}};
  trial(phi,p.options_.upper_potential_V);trial(psi,1.);LinearForm unit(&p.space_);unit=0.;OperatorPtr A,B;Vector x,b,y,c;
  form.FormLinearSystem(essential,phi,rhs,A,x,b);form.FormLinearSystem(essential,psi,unit,B,y,c);
  check(A.Ptr()==B.Ptr(),"Native phi/psi use one constrained operator");
  csr(out/"full_operator",*p.poisson_->full_matrix);csr(out/"constrained_operator",*dynamic_cast<SparseMatrix*>(A.Ptr()));vector_file(out/"charge_rhs",full_rhs);vector_file(out/"phi_rhs",b);vector_file(out/"psi_rhs",c);write(out/"essential_owner.bin",p.field_owner_identity_);
 }
};
}
namespace {
void actual_core(const fs::path&fixture,const fs::path&out,Json&report){
 auto options=read<GapPicOptions>(fixture/"options.bin");auto particles=read<GapParticle>(fixture/"particles.bin");check(options.size()==1&&particles.size()==46726,"Bound actual captured source inputs");const auto original=particles;
 GapPic uniform(options[0]);PicFieldMeshOptions fm;fm.margin_cells=120;GapPic graded(options[0],fm);PicFieldMeshAudit::particles(uniform,particles);PicFieldMeshAudit::particles(graded,particles);
 PicFieldMeshAudit::assemble_control(uniform,out/"native9");PicFieldMeshAudit::assemble_control(graded,out/"native-graded120");
 write(out/"native-graded120/axes_x.bin",graded.field_layout().xs);write(out/"native-graded120/axes_y.bin",graded.field_layout().ys);
 double a=0,b=0;auto q=PicFieldMeshAudit::charge(uniform,a),r=PicFieldMeshAudit::charge(graded,b);bool source=true;Vector mapped(r.Size());mapped=0.;
 for(int i=0;i<q.Size();++i)mapped[graded.core_field_dofs()[i]]=q[i];for(int i=0;i<r.Size();++i)source&=same(mapped[i],r[i]);check(source&&same(a,b),"Graded native deposited RHS is bitwise old-core mapping");
 check(std::abs(q.Sum()*options[0].depth_m-a)<=1e-11*std::max(std::abs(a),pic_e),"Original source charge conservation budget");
 std::size_t samples=0;
 for(const auto*name:{"phi.bin","psi.bin"}){
  const auto nodes=read<double>(fixture/name);PicFieldMeshAudit::nodes(uniform,nodes);PicFieldMeshAudit::nodes(graded,nodes);
  if(std::string(name)=="phi.bin"){
   double point_work=0,activity=0;for(const auto&p:particles){const double term=p.charge_C*p.weight*graded.potential_value(p.position_m[0],p.position_m[1]);point_work+=term;activity+=std::abs(term);}
   const double nodal_work=(r*graded.potential())*options[0].depth_m;const double error=std::abs(point_work-nodal_work),limit=1024*std::numeric_limits<double>::epsilon()*std::max(activity,1e-30);
   check(error<=limit,"Native graded source work matches point work on actual source scale");report["source_work"]={{"point_J",point_work},{"nodal_J",nodal_work},{"error_J",error},{"limit_J",limit}};
  }
  for(const auto&p:particles){auto x=uniform.electric_field(p.position_m[0],p.position_m[1]),y=graded.electric_field(p.position_m[0],p.position_m[1]);for(int k=0;k<3;++k)if(!same(x[k],y[k]))throw std::runtime_error("Captured core gradient changed");if(!same(uniform.potential_value(p.position_m[0],p.position_m[1]),graded.potential_value(p.position_m[0],p.position_m[1])))throw std::runtime_error("Captured scalar gather changed");++samples;}
  for(int j=0;j<options[0].ny;++j)for(int i=0;i<options[0].nx;++i){PicKineticCellPoint point{i,j,.5,.5};IntegrationPoint ip;ip.Set2(.5,.5);const double old=uniform.potential().GetValue(j*options[0].nx+i,ip);if(!same(old,uniform.potential_value(point))||!same(old,graded.potential_value(point)))throw std::runtime_error("Snapshot cell-center gather changed");++samples;}
 }
 check(std::memcmp(particles.data(),original.data(),particles.size()*sizeof(GapParticle))==0,"Captured particle bytes unchanged");
 check(uniform.field_topology_identity()!=graded.field_topology_identity(),"Uniform and graded topology digests differ");
 report["actual_core"]={{"captured_particles",particles.size()},{"phi_psi_particle_and_snapshot_samples",samples},{"all_samples_bitwise_equal",true},{"source_bitwise_mapping",source},{"field_nodes",graded.potential().Size()},{"field_x_bounds_m",{graded.field_layout().xs.front(),graded.field_layout().xs.back()}},{"uniform_topology_sha256",uniform.field_topology_identity()},{"graded_topology_sha256",graded.field_topology_identity()},{"actual_fixture_solves",0}};
}
void transport_and_restart(const fs::path&out,Json&report){
 auto o=small();o.upper_potential_V=0;PicFieldMeshOptions fm;fm.margin_cells=8;GapPic a(o),b(o,fm);
 for(auto electrode:{GapElectrode::Lower,GapElectrode::Upper})for(auto face:{GapSurface::Gap,GapSurface::Back,GapSurface::Left,GapSurface::Right}){a.inject_thermal_electrons(electrode,2000,1,1e-12,7,123,face);b.inject_thermal_electrons(electrode,2000,1,1e-12,7,123,face);}
 check(PicFieldMeshAudit::transport_equal(a,b),"All surface source positions/IDs/weights/normals/ledgers unchanged");
 const auto scenario=[&](GapParticle p,double dt,const char*label){GapPic x(o),y(o,fm);x.add_particle(p);y.add_particle(p);PicFieldMeshAudit::drift(x,dt);PicFieldMeshAudit::drift(y,dt);check(PicFieldMeshAudit::transport_equal(x,y),label);return PicFieldMeshAudit::ledger(y);};
 GapParticle p;p.position_m={.00099,.0005,0};p.proper_velocity_m_s={2e7,0,0};p.weight=1e-8;auto escape=scenario(p,1e-12,"Original kinetic x escape and path/continuity remain unchanged");check(escape.escaped_charge_C<0,"Extended field does not extend kinetic escape");
 p.position_m={.0005,.0005,o.short_z_m+1e-6};p.proper_velocity_m_s={0,0,-2e7};auto short_end=scenario(p,1e-12,"Finite-z short footprint unchanged");check(short_end.short_absorbed_charge_C<0,"Short absorption still owns bridge footprint");
 p.position_m={.00005,.0005,o.short_z_m+1e-6};auto axial=scenario(p,1e-12,"Axial escape outside short footprint unchanged");check(axial.axial_escaped_charge_C<0&&axial.short_absorbed_charge_C==0,"Axial exterior sink distinct from short");
 p.position_m={-.00001,.0005,0};rejects([&]{b.add_particle(p);},"Field-only vacuum is not a kinetic insertion domain");
 auto charged=small();GapPic base(charged);p.position_m={.0005,.0005,0};p.proper_velocity_m_s={100,200,300};base.add_particle(p);base.solve_fields();const auto saved=out/"uniform.pic3";base.save_checkpoint(saved.string());GapPic same_layout(charged);same_layout.load_checkpoint(saved.string());check(same_layout.particles()[0].id==base.particles()[0].id,"PIC3 same-topology restart preserves identity");
 GapPic wrong(charged,fm);rejects([&]{wrong.load_checkpoint(saved.string());},"PIC3 rejects changed field topology before importing state");check(wrong.particles().empty(),"Failed field-topology restart leaves population untouched");
 auto bytes=read<char>(saved);const std::uint64_t legacy_magic=0x5241494b50494332ULL;std::memcpy(bytes.data(),&legacy_magic,sizeof(legacy_magic));bytes.erase(bytes.begin()+8+sizeof(GapPicOptions),bytes.begin()+8+sizeof(GapPicOptions)+64);auto legacy=out/"manufactured-legacy.pic2";write(legacy,bytes);GapPic migration(charged,fm);
 rejects([&]{migration.load_checkpoint(legacy.string());},"Legacy checkpoint cannot silently change field model");rejects([&]{migration.import_legacy_kinetic_checkpoint(legacy.string(),std::string(64,'0'));},"Legacy kinetic import requires exact file SHA256");migration.import_legacy_kinetic_checkpoint(legacy.string(),checkpoint_sha256(legacy.string()));check(migration.particles().size()==base.particles().size(),"Explicit hash-bound legacy kinetic-state import succeeds");
 const auto graded_saved=out/"graded.pic3";migration.save_checkpoint(graded_saved.string());PicFieldMeshOptions other=fm;other.growth_ratio=1.06;GapPic same_size_other_axes(charged,other);rejects([&]{same_size_other_axes.load_checkpoint(graded_saved.string());},"Equal node count cannot hide changed field axes");
 auto bias=charged;bias.upper_potential_V=1;GapPic unit(bias,fm);check(unit.field_topology_identity()==migration.field_topology_identity(),"Main phi and unit-bias psi have identical geometry identity");
 auto corrupt=read<char>(graded_saved);corrupt[8+sizeof(GapPicOptions)]^=1;write(out/"corrupt-topology.pic3",corrupt);GapPic corrupted_target(charged,fm);rejects([&]{corrupted_target.load_checkpoint((out/"corrupt-topology.pic3").string());},"Corrupt topology digest is rejected");
 PicFieldMeshAudit::move_first_vertex(migration,1e-10);rejects([&]{migration.save_checkpoint((out/"mutated-topology.pic3").string());},"Unversioned mesh-coordinate mutation cannot masquerade as old topology");PicFieldMeshAudit::move_first_vertex(migration,-1e-10);
 report["checkpoint"]={{"new_magic","PIC3"},{"topology_sha256",migration.field_topology_identity()},{"legacy_import_sha256",checkpoint_sha256(legacy.string())},{"kinetic_options_bytes",sizeof(GapPicOptions)},{"particle_bytes",sizeof(GapParticle)},{"diagnostics_bytes",sizeof(GapPicDiagnostics)},{"graded_CPU_residual_C_per_m",migration.diagnostics().poisson_residual_C_per_m}};
}
void manufactured(Json&report){
 PicFieldMeshOptions fm;fm.margin_cells=12;fm.full_width_m=fm.full_height_m=.04;PicFieldLayout layout(8,8,.009,.009,fm);auto mesh=layout.make_mesh();H1_FECollection collection(1,2);FiniteElementSpace space(&mesh,&collection);GridFunction field(&space);FunctionCoefficient exact([](const Vector&x){return 2+3*x[0]-4*x[1]+5*x[0]*x[1];});field.ProjectCoefficient(exact);Vector exact_nodes(field);ConstantCoefficient epsilon(8.8541878188e-12);BilinearForm form(&space);form.AddDomainIntegrator(new DiffusionIntegrator(epsilon));form.Assemble();form.Finalize();LinearForm rhs(&space);rhs=0.;Array<int>marker(mesh.bdr_attributes.Max());marker=1;Array<int>essential;space.GetEssentialTrueDofs(marker,essential);OperatorPtr A;Vector x,b;form.FormLinearSystem(essential,field,rhs,A,x,b);auto*sparse=dynamic_cast<SparseMatrix*>(A.Ptr());GSSmoother prec(*sparse);CGSolver solver;solver.SetOperator(*A);solver.SetPreconditioner(prec);solver.SetRelTol(1e-14);solver.SetAbsTol(1e-30);solver.SetMaxIter(4000);solver.SetPrintLevel(-1);solver.Mult(b,x);check(solver.GetConverged(),"Manufactured graded harmonic Poisson solve converges");form.RecoverFEMSolution(x,rhs,field);Q1ScalarGather gather(field);double maximum=0,gradient=0;
 for(int i=0;i<field.Size();++i)maximum=std::max(maximum,std::abs(field[i]-exact_nodes[i]));
 for(int j=0;j<layout.ny();++j)for(int i=0;i<layout.nx();++i){double px=layout.xs[i]+.23*(layout.xs[i+1]-layout.xs[i]),py=layout.ys[j]+.61*(layout.ys[j+1]-layout.ys[j]);IntegrationPoint ip;int e=layout.locate({px,py},ip);auto g=gather.gradient(e,ip);gradient=std::max({gradient,std::abs(g[0]-(3+5*py)),std::abs(g[1]-(-4+5*px))});}
 check(maximum<1e-10&&gradient<1e-7,"Manufactured graded scalar/Jacobian gradients match analytic field");
 report["manufactured"]={{"nodes",field.Size()},{"max_nodal_error_V",maximum},{"max_gradient_error_V_m",gradient},{"CPU_only",true}};
}
}
int main(int argc,char**argv){try{
 if(argc!=3)throw std::runtime_error("Expected READONLY_FIXTURE NEW_OUTPUT");cpu_set_t cpu;CPU_ZERO(&cpu);sched_getaffinity(0,sizeof(cpu),&cpu);check(CPU_COUNT(&cpu)==1&&CPU_ISSET(7,&cpu)&&getpriority(PRIO_PROCESS,0)>=15,"CPU7/nice15 contract");
 fs::path out=argv[2];check(!fs::exists(out),"Fresh test output");fs::create_directories(out);Json report;
 actual_core(argv[1],out,report);transport_and_restart(out,report);manufactured(report);
 report["status"]="passed";report["checks"]=checks;report["GPU_execution"]=false;report["production_actions"]=false;report["actual_fixture_evolution"]=false;report["field_extension_is_not_kinetic_or_sheath_or_3D_convergence"]=true;struct rusage usage;getrusage(RUSAGE_SELF,&usage);report["max_RSS_KiB"]=usage.ru_maxrss;
 std::ofstream(out/"native-field-tests.json")<<report.dump(2)<<'\n';std::cout<<"Passed "<<checks.size()<<" native field-mesh CPU checks\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
