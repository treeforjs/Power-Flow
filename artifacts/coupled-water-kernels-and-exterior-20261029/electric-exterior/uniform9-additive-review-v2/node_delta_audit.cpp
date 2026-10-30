// CPU-only old-fresh to new-fresh nodal propagation audit; no GPU or solve.
// CPU-only counterfactual audit. Include the sealed CPU assembly/gather without
// altering it or executing its entry point. STATIC_DOMAIN_WITH_CUDSS stays off.
#define main sealed_static_domain_entry_not_called
#include "/home/treefor/projects/runtimes/pic-static-domain-20261029/source/static_domain.cpp"
#undef main
#include <cfloat>
#include <sched.h>

struct AuditMetric {
 std::size_t failures=0,worst_particle=0;
 double maximum=0,maximum_limit=0,maximum_fraction=0;
 void add(double difference,double limit,std::size_t i) {
  io::require(std::isfinite(difference)&&std::isfinite(limit)&&limit>=0,"Nonfinite audit metric");
  if(difference>maximum){maximum=difference;worst_particle=i;}
  maximum_limit=std::max(maximum_limit,limit);
  maximum_fraction=std::max(maximum_fraction,limit>0?difference/limit:(difference?std::numeric_limits<double>::max():0.));
  if(difference>limit)++failures;
 }
};
void metric(std::ostream& f,const AuditMetric& m) {
 f<<"{\"failures\":"<<m.failures<<",\"maximum\":"<<m.maximum<<",\"maximum_limit\":"<<m.maximum_limit
  <<",\"maximum_fraction_of_limit\":"<<m.maximum_fraction<<",\"worst_particle_index\":"<<m.worst_particle<<'}';
}
template<class T>void load_grid(GridFunction& field,const std::vector<T>& values) {
 io::require(values.size()==std::size_t(field.Size()),"Nodal vector shape mismatch");
 for(int i=0;i<field.Size();++i){io::require(std::isfinite(values[i]),"Nonfinite nodal vector");field[i]=values[i];}
}
int main(int argc,char** argv) {
 try {
  io::require(argc==5,"Expected ACTUAL_FIXTURE OLD_FRESH NEW_FRESH NEW_AUDIT_OUTPUT");
  cpu_set_t cpus;CPU_ZERO(&cpus);
  io::require(sched_getaffinity(0,sizeof(cpus),&cpus)==0&&CPU_COUNT(&cpus)==1&&CPU_ISSET(7,&cpus)&&getpriority(PRIO_PROCESS,0)>=15,"Audit requires CPU7/nice15");
  const std::filesystem::path fixture=argv[1],baseline=argv[2],gpu=argv[3],out=argv[4];
  io::require(!std::filesystem::exists(out),"Audit output must be fresh");std::filesystem::create_directories(out);
  auto s=pic_batch::HostSnapshot::load(fixture,true);const auto particles_before=s.particles;
  auto options=io::read<GapPicOptions>(fixture/"options.bin");io::require(options.size()==1,"Expected one geometry");
  StaticCase c(options[0],9,s);const auto baseline_phi=io::read<double>(baseline/"phi_nodes.bin"),baseline_psi=io::read<double>(baseline/"psi_nodes.bin");load_grid(c.phi,baseline_phi);load_grid(c.psi,baseline_psi);
  const auto capture_phi_residual=residual(c.matrix(),c.phi,c.phi_b),capture_psi_residual=residual(c.matrix(),c.psi,c.psi_b);
  const auto captured=c.gather(s);const auto original_capture_gate=compare(captured,s);const auto baseline_stored=io::read<FieldRecord>(baseline/"particle_fields.bin");const bool identical_baseline_fields=baseline_stored.size()==captured.size()&&std::memcmp(baseline_stored.data(),captured.data(),captured.size()*sizeof(FieldRecord))==0;
  const auto fresh_phi=io::read<double>(gpu/"phi_nodes.bin"),fresh_psi=io::read<double>(gpu/"psi_nodes.bin");
  const auto stored=io::read<FieldRecord>(gpu/"particle_fields.bin");
  load_grid(c.phi,fresh_phi);load_grid(c.psi,fresh_psi);
  const auto fresh_phi_residual=residual(c.matrix(),c.phi,c.phi_b),fresh_psi_residual=residual(c.matrix(),c.psi,c.psi_b);
  const auto recomputed=c.gather(s);const auto original_fresh_gate=compare(recomputed,s);
  io::require(stored.size()==recomputed.size(),"Saved field count mismatch");
  const bool identical_fresh_fields=std::memcmp(stored.data(),recomputed.data(),stored.size()*sizeof(FieldRecord))==0;
  const double eps=std::numeric_limits<double>::epsilon(),gamma64=64*eps/(1-64*eps),h=c.o.width_m/c.o.nx;
  double phi_scale=0,psi_scale=0,max_dphi=0,max_dpsi=0;std::array<double,3>Escale{};
  for(std::size_t i=0;i<baseline_phi.size();++i){phi_scale=std::max({phi_scale,std::abs(baseline_phi[i]),std::abs(fresh_phi[i])});psi_scale=std::max({psi_scale,std::abs(baseline_psi[i]),std::abs(fresh_psi[i])});max_dphi=std::max(max_dphi,std::abs(fresh_phi[i]-baseline_phi[i]));max_dpsi=std::max(max_dpsi,std::abs(fresh_psi[i]-baseline_psi[i]));}
  for(std::size_t i=0;i<captured.size();++i)for(int k=0;k<3;++k)Escale[k]=std::max({Escale[k],std::abs(captured[i].E.v[k]),std::abs(recomputed[i].E.v[k])});
  std::array<AuditMetric,3>strict,regression,propagation;std::array<AuditMetric,2>potential;
  std::array<double,3>max_delta_prediction{},max_self_delta{},max_axial_delta{},max_remainder{};
  std::vector<gp::Vec> predicted,observed,round_bounds;predicted.reserve(s.particles.size());observed.reserve(s.particles.size());round_bounds.reserve(s.particles.size());
  std::ofstream failures(out/"original-strict-failures.csv");failures<<std::setprecision(17)<<"id,component,baseline_E,new_E,observed_delta,self_delta,axial_delta,linear_prediction,floating_remainder,a_priori_round_bound,original_strict_limit,existing_regression_limit\n";
  for(std::size_t i=0;i<s.particles.size();++i) {
   const auto& p=s.particles[i];IntegrationPoint ip;const int e=c.locate(p,ip);auto* tr=c.mesh.GetElementTransformation(e);tr->SetIntPoint(&ip);
   auto* fe=c.space.GetFE(e);Vector shape(4);DenseMatrix ds(4,2);fe->CalcShape(ip,shape);fe->CalcDShape(ip,ds);const auto& inv=tr->InverseJacobian();
   const auto& cell=s.self_cells[e];long double grad_phi[2]{},grad_psi[2]{},abs_grad_phi[2]{},abs_grad_psi[2]{},delta_psi_value=0,abs_psi_value=0;
   for(int j=0;j<4;++j) {
    const int d=cell.dof[j];const long double dp=static_cast<long double>(fresh_phi[d])-baseline_phi[d],du=static_cast<long double>(fresh_psi[d])-baseline_psi[d];
    delta_psi_value+=static_cast<long double>(shape[j])*du;abs_psi_value+=std::abs(static_cast<long double>(shape[j]))*(std::abs(fresh_psi[d])+std::abs(baseline_psi[d]));
    for(int k=0;k<2;++k) {
     const long double a=static_cast<long double>(ds(j,0))*inv(0,k),b=static_cast<long double>(ds(j,1))*inv(1,k),coefficient=a+b;
     grad_phi[k]+=coefficient*dp;grad_psi[k]+=coefficient*du;
     abs_grad_phi[k]+=(std::abs(a)+std::abs(b))*(std::abs(fresh_phi[d])+std::abs(baseline_phi[d]));
     abs_grad_psi[k]+=(std::abs(a)+std::abs(b))*(std::abs(fresh_psi[d])+std::abs(baseline_psi[d]));
    }
   }
   const auto& t=s.parameters;const double g=(p.position_m[2]-t.short_z)/t.depth,alpha=t.voltage*(2*g-1),beta=-2*t.voltage/t.depth;
   // Native terms are identical; reconstruct the unchanged native Ez solely
   // for its absolute operation scale in the floating-point forward bound.
   const double x=p.position_m[0]-.0045,y=p.position_m[1]-.0045;const int ni=gp::bin(s.xs.data(),t.native_nx,x),nj=gp::bin(s.ys.data(),t.native_ny,y),ne=nj*t.native_nx+ni;
   const double nx=(x-s.xs[ni])/(s.xs[ni+1]-s.xs[ni]),ny=(y-s.ys[nj])/(s.ys[nj+1]-s.ys[nj]);gp::Grid grid{s.native_cells.data(),t.native_nx,t.native_ny};
   auto r0=gp::sample(grid,s.rate0.data(),ne,nx,ny),r1=gp::sample(grid,s.rate1.data(),ne,nx,ny);io::require(!r0.status&&!r1.status,"Native sample invalid");
   const double native_Ez=(1-t.fraction)*(-r0.value)+t.fraction*(-r1.value);
   gp::Vec pred{},obs{},bound{};
   for(int k=0;k<3;++k) {
    const long double self=k<2?-grad_phi[k]:0,axial=k<2?-static_cast<long double>(alpha)*grad_psi[k]:static_cast<long double>(beta)*delta_psi_value;
    const long double prediction=self+axial,difference=static_cast<long double>(recomputed[i].E.v[k])-captured[i].E.v[k],remainder=difference-prediction;
    const long double absolute_scale=k<2?abs_grad_phi[k]+std::abs(alpha)*abs_grad_psi[k]:2*std::abs(native_Ez)+std::abs(beta)*(abs_psi_value+1);
    const double round_bound=static_cast<double>(gamma64*absolute_scale+64*LDBL_EPSILON*(std::abs(self)+std::abs(axial)));
    const double strict_limit=E_abs+field_rel*std::abs(s.reference[i].E.v[k]);
    // Preserve original Ez check. Existing checkpoint field policy defines
    // transverse gradient regression, not a new independent Ez tolerance.
    const double regression_limit=k<2?1e-10*Escale[k]+256*eps*phi_scale/h:strict_limit;
    strict[k].add(std::abs(static_cast<double>(difference)),strict_limit,i);regression[k].add(std::abs(static_cast<double>(difference)),regression_limit,i);
    propagation[k].add(std::abs(static_cast<double>(remainder)),round_bound,i);
    pred.v[k]=static_cast<double>(prediction);obs.v[k]=static_cast<double>(difference);bound.v[k]=round_bound;
    max_delta_prediction[k]=std::max(max_delta_prediction[k],std::abs(pred.v[k]));max_self_delta[k]=std::max(max_self_delta[k],std::abs(static_cast<double>(self)));max_axial_delta[k]=std::max(max_axial_delta[k],std::abs(static_cast<double>(axial)));max_remainder[k]=std::max(max_remainder[k],std::abs(static_cast<double>(remainder)));
    if(std::abs(static_cast<double>(difference))>strict_limit)failures<<p.id<<','<<k<<','<<captured[i].E.v[k]<<','<<recomputed[i].E.v[k]<<','<<obs.v[k]<<','<<self<<','<<axial<<','<<pred.v[k]<<','<<remainder<<','<<round_bound<<','<<strict_limit<<','<<regression_limit<<'\n';
   }
   predicted.push_back(pred);observed.push_back(obs);round_bounds.push_back(bound);
   potential[0].add(std::abs(recomputed[i].phi-captured[i].phi),(1e-10+256*eps)*phi_scale,i);
   potential[1].add(std::abs(recomputed[i].psi-captured[i].psi),(1e-10+256*eps)*psi_scale,i);
  }
  failures.close();io::require(bool(failures),"Audit failure CSV write failed");
  csr_files(out/"full_operator",*c.full);csr_files(out/"constrained_operator",c.matrix());vector_file(out/"charge_rhs",c.raw_charge);vector_file(out/"phi_rhs",c.phi_b);vector_file(out/"psi_rhs",c.psi_b);
  io::write(out/"captured_nodes_recomputed_fields.bin",captured);io::write(out/"fresh_nodes_recomputed_fields.bin",recomputed);io::write(out/"linear_predicted_delta_E.bin",predicted);io::write(out/"observed_delta_E.bin",observed);io::write(out/"a_priori_round_bounds.bin",round_bounds);
  bool pass=identical_baseline_fields&&identical_fresh_fields&&capture_phi_residual.passed&&capture_psi_residual.passed&&fresh_phi_residual.passed&&fresh_psi_residual.passed;
  for(auto&m:regression)pass&=!m.failures;for(auto&m:propagation)pass&=!m.failures;for(auto&m:potential)pass&=!m.failures;
  io::require(std::memcmp(s.particles.data(),particles_before.data(),s.particles.size()*sizeof(GapParticle))==0,"Particle input changed");
  std::ofstream f(out/"cpu-control-audit.json");f<<std::setprecision(17)<<std::boolalpha
   <<"{\"passed\":"<<pass<<",\"CPU_only\":true,\"GPU_executed\":false,\"factorizations\":0,\"solves\":0,\"particles\":"<<s.particles.size()
   <<",\"baseline_vs_original_strict_failures\":"<<original_capture_gate.failures<<",\"baseline_vs_original_max_E\":"<<original_capture_gate.maxE
   <<",\"baseline_nodal_CPU_gather_byte_equal_old_GPU_fields\":"<<identical_baseline_fields<<",\"fresh_nodal_CPU_gather_byte_equal_saved_GPU_fields\":"<<identical_fresh_fields<<",\"fresh_original_strict_failures\":"<<original_fresh_gate.failures
   <<",\"fresh_original_max_E_difference\":"<<original_fresh_gate.maxE<<",\"fresh_original_max_B_difference\":"<<original_fresh_gate.maxB
   <<",\"max_nodal_phi_delta\":"<<max_dphi<<",\"max_nodal_psi_delta\":"<<max_dpsi<<",\"phi_scale\":"<<phi_scale<<",\"psi_scale\":"<<psi_scale<<",\"cell_h_m\":"<<h
   <<",\"capture_phi_residual\":"<<capture_phi_residual.absolute<<",\"capture_psi_residual\":"<<capture_psi_residual.absolute
   <<",\"fresh_phi_residual\":"<<fresh_phi_residual.absolute<<",\"fresh_psi_residual\":"<<fresh_psi_residual.absolute
   <<",\"original_phi_residual_limit\":"<<fresh_phi_residual.limit<<",\"original_psi_residual_limit\":"<<fresh_psi_residual.limit
   <<",\"topology_bytes_equal\":"<<c.topology_bitwise<<",\"independent_source_max_error\":"<<c.independent_source_max_error
   <<",\"direct_charge_C\":"<<c.direct_charge<<",\"deposited_charge_C\":"<<c.deposited_charge<<",\"components\":[";
  for(int k=0;k<3;++k){if(k)f<<',';f<<"{\"component\":"<<k<<",\"original_strict\":";metric(f,strict[k]);f<<",\"existing_regression\":";metric(f,regression[k]);f<<",\"linear_propagation_remainder\":";metric(f,propagation[k]);f<<",\"max_predicted_delta\":"<<max_delta_prediction[k]<<",\"max_self_delta\":"<<max_self_delta[k]<<",\"max_axial_delta\":"<<max_axial_delta[k]<<",\"max_floating_remainder\":"<<max_remainder[k]<<'}';}
  f<<"],\"sampled_phi_regression\":";metric(f,potential[0]);f<<",\"sampled_psi_regression\":";metric(f,potential[1]);f<<"}\n";f.close();io::require(bool(f),"Audit receipt write failed");
  std::cout<<"CPU control audit "<<(pass?"passed":"failed")<<"; baseline-vs-original strict failures="<<original_capture_gate.failures<<"; fresh strict failures="<<original_fresh_gate.failures<<'\n';return pass?0:1;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 2;}
}
