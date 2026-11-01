#!/usr/bin/env python3
from pathlib import Path
import ctypes as C
import fcntl,hashlib,json,math,os,random,shutil,subprocess,sys,time
from decimal import Decimal as D,localcontext
S=Path(__file__).resolve().parent
O=Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs/ss304-thermal-properties-v2-20261029')
V1=O.parent/'ss304-thermal-properties-20261029'
LOCK=O.parent/'pic-profile/testing-core7.lock'
EPS=sys.float_info.epsilon
class Values(C.Structure):
 _fields_=[(n,C.c_double) for n in ['temperature_k','h_j_kg','u_model_j_kg','cp_sensible_j_kg_k','apparent_heat_capacity_j_kg_k','k_w_m_k','resistivity_ohm_m','sigma_s_m','liquid_fraction','energy_j_m3_reference','capacity_j_m3_reference_k']]+[('phase',C.c_int),('synthesized_high_temperature',C.c_int)]
class Gas(C.Structure):
 _fields_=[(n,C.c_double) for n in ['saturation_pressure_pa','kim_pressure_pa','vapor_h_j_kg','vapor_u_j_kg','latent_h_j_kg','gas_cp_j_kg_k','gas_cv_j_kg_k','atomic_reference_u0_j_kg']]
class State(C.Structure):
 _fields_=[('mass_kg',C.c_double),('energy_j',C.c_double),('reference_volume_m3',C.c_double),('revision',C.c_uint64),('last_event',C.c_uint64)]
class Trial(C.Structure):
 _fields_=[('next',State)]+[(n,C.c_double) for n in ['temperature_k','emitted_reference_plus_kinetic_j','signed_energy_residual_j','caloric_projection_residual_j','remaining_mass_fraction','removed_reference_volume_m3']]
class Support(C.Structure):
 _fields_=[(n,C.c_double) for n in ['minimum_j','maximum_j','cold_margin_j','hot_margin_j']]+[('range_status',C.c_int)]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def hold():
 p=Path('/proc/278902');raw=(p/'stat').read_text();v=raw[raw.rfind(')')+2:].split();e=dict(line.split() for line in Path('/sys/fs/cgroup/powerflow-pic-gpu-reserved/cgroup.events').read_text().splitlines())
 assert v[0]=='T' and v[19]=='5288339' and e=={'populated':'1','frozen':'1'}
 assert str((p/'exe').resolve())=='/home/treefor/projects/runtimes/mfem-gpu-pic-20261029/qualified/heating-history-pic-gpu'
 return {'pid':278902,'state':v[0],'start_ticks':v[19],'events':e}
def verify():
 start=time.time();before=hold();source_before={str(p.relative_to(S)):sha(p) for p in S.rglob('*') if p.is_file()}
 seal=json.loads((V1/'source-manifest.json').read_text());root=Path(seal['source_root'])
 assert all(sha(root/k)==v for k,v in seal['files'].items())
 assert sha(S/'include/sealed1800.hpp')==seal['files']['include/ss304_thermal.hpp']
 assert sha(seal['compiled_library']['path'])==seal['compiled_library']['sha256']
 B=S.parent/'build';B.mkdir(exist_ok=True);compiler=Path(shutil.which('c++')).resolve()
 cmd=[str(compiler),'-std=c++17','-O2','-fno-fast-math','-Wall','-Wextra','-Werror','-pedantic','-fPIC','-shared',str(S/'c_api.cpp'),'-o',str(B/'libss304v2.so')]
 built=subprocess.run(cmd,capture_output=True,text=True);(O/'build.stdout.txt').write_text(built.stdout);(O/'build.stderr.txt').write_text(built.stderr);assert built.returncode==0,built.stderr
 lib=C.CDLL(str(B/'libss304v2.so'));old=C.CDLL(seal['compiled_library']['path'])
 lib.ss304v2_error.restype=C.c_char_p
 lib.ss304v2_evaluate.argtypes=[C.c_int,C.c_double,C.c_double,C.POINTER(Values)]
 old.ss304_evaluate.argtypes=lib.ss304v2_evaluate.argtypes
 lib.ss304v2_inverse.argtypes=[C.c_int,C.c_double,C.c_double,C.c_int,C.POINTER(C.c_double)]
 lib.ss304v2_vapor.argtypes=[C.c_double,C.POINTER(Gas)]
 lib.ss304v2_flux.argtypes=[C.c_double,C.c_double,C.c_double,C.POINTER(C.c_double)]
 lib.ss304v2_gas_k.argtypes=[C.c_double,C.c_double,C.c_double,C.POINTER(C.c_double)]
 lib.ss304v2_transfer.argtypes=[C.c_double,C.POINTER(State),C.c_uint64,C.c_uint64,C.c_double,C.c_double,C.c_double,C.POINTER(Trial)]
 lib.ss304v2_energy_support.argtypes=[C.c_int,C.c_double,C.c_double,C.c_double,C.POINTER(Support)]
 lib.ss304v2_evaluate_inventory.argtypes=[C.c_int,C.c_double,C.c_double,C.c_double,C.c_double,C.POINTER(Values)]
 counts={};maxima={}
 def check(ok,label):
  if not ok:raise AssertionError(label)
  counts[label]=counts.get(label,0)+1
 def near(a,b,label,rel=64*EPS,absolute=0):
  error=abs(a-b);limit=max(absolute,rel*max(abs(a),abs(b)));check(error<=limit,label)
  maxima[label]=max(maxima.get(label,0),error)
 def rc(code):
  if code:raise ValueError(lib.ss304v2_error().decode())
 def ev(t,q=2):
  v=Values();rc(lib.ss304v2_evaluate(q,8000,t,C.byref(v)));return v
 def inv(e,q=2,vol=0):
  v=C.c_double();rc(lib.ss304v2_inverse(q,8000,e,vol,C.byref(v)));return v.value
 def gas(t):
  v=Gas();rc(lib.ss304v2_vapor(t,C.byref(v)));return v
 def flux(t,p=0,a=1):
  v=C.c_double();rc(lib.ss304v2_flux(t,p,a,C.byref(v)));return v.value
 def gk(t,d,ion=0):
  v=C.c_double();rc(lib.ss304v2_gas_k(t,d,ion,C.byref(v)));return v.value
 def rejected(f,label):
  try:f()
  except ValueError:check(True,label)
  else:check(False,label)
 def transfer(b,dm,k,q=0,rev=None,event=None):
  v=Trial();rc(lib.ss304v2_transfer(8000,C.byref(b),b.revision if rev is None else rev,b.last_event+1 if event is None else event,dm,k,q,C.byref(v)));return v
 rng=random.Random(20261030)
 with localcontext() as context:
  context.prec=60
  for t in [298.2,math.nextafter(298.2,math.inf),299.,math.nextafter(300,-math.inf),300.]:
   v=ev(t);T=D.from_float(t)
   h=(T-300)*(D(472)+D('.068')*(T+300)-D(2820000)/(300*T))
   near(v.h_j_kg,float(h),'cold_enthalpy_decimal',absolute=1e-28)
   near(v.cp_sensible_j_kg_k,float(D(472)+D('.136')*T-D(2820000)/(T*T)),'cold_cp_decimal')
   if t<300:check(v.h_j_kg<0,'negative_cold_energy_reference')
   check(v.h_j_kg==v.u_model_j_kg,'same_reference_h_u')
   near(inv(v.h_j_kg),t,'cold_inverse',rel=4*EPS)
  for t in [300.,773.15,1500.,1672.,1700.,1727.,1800.]+[300+1500*rng.random() for _ in range(1500)]:
   a=ev(t);b=Values();check(old.ss304_evaluate(1,8000,t,C.byref(b))==0,'sealed_v1_evaluation')
   check(bytes(a)==bytes(b),'sealed_v1_bitwise_values_300_to_1800')
  a=ev(1800);previous=None
  for t in sorted([298.2,300.,1500.,1672.,1727.,1800.,3000.]+[298.2+(3000-298.2)*rng.random() for _ in range(2500)]):
   v=ev(t)
   near(inv(v.h_j_kg),t,'full_caloric_inverse',rel=4*EPS)
   near(inv(v.energy_j_m3_reference,vol=1),t,'full_reference_volume_inverse',rel=4*EPS)
   check(v.k_w_m_k>0 and v.sigma_s_m>0 and v.apparent_heat_capacity_j_kg_k>0,'positive_condensed_properties')
   if previous:check(v.h_j_kg>previous,'strictly_increasing_h')
   previous=v.h_j_kg
   if t>1800:
    near(v.h_j_kg,a.h_j_kg+float(D('.190')*4184)*(t-1800),'source_liquid_energy_increment')
    near(v.cp_sensible_j_kg_k,794.96,'source_liquid_cp')
    near(v.k_w_m_k,28.99*(.1229+.00003248*t)/(.1229+.00003248*1800),'source_conditioned_k_shape')
    near(v.k_w_m_k*v.resistivity_ohm_m/t,28.99*1.517e-6/1800,'constant_anchored_lorenz')
  near(ev(math.nextafter(1800,math.inf)).h_j_kg,a.h_j_kg,'liquid_join_energy_continuity',rel=4*EPS)
  for t in [298.3,299.5,1800.1,2200.,2700.,2999.]:
   dt=1e-4;derivative=(ev(t+dt).h_j_kg-ev(t-dt).h_j_kg)/(2*dt)
   near(derivative,ev(t).apparent_heat_capacity_j_kg_k,'dH_dT',rel=1e-7)
  for q,max_t in [(0,1500.),(1,1800.),(2,3000.)]:
   for t in [math.nextafter(298.2,-math.inf),math.nextafter(max_t,math.inf),math.nan,math.inf]:rejected(lambda t=t:ev(t,q),'temperature_range_rejection')
   for vol in [0,1]:
    lo=ev(298.2,q).energy_j_m3_reference if vol else ev(298.2,q).h_j_kg
    hi=ev(max_t,q).energy_j_m3_reference if vol else ev(max_t,q).h_j_kg
    for x in [math.nextafter(lo,-math.inf),math.nextafter(hi,math.inf),math.nan]:rejected(lambda x=x:inv(x,q,vol),'energy_range_rejection')
 R=8.31446261815324/(.713*.055845+.202*.0519961+.085*.0586934)
 for t in [1800.,2000.,2400.,2800.,3000.]:
  v=gas(t);near(v.kim_pressure_pa,101325*10**(6.121-18836/t),'kim_pressure_source')
  near(v.vapor_h_j_kg-ev(t).h_j_kg,v.latent_h_j_kg,'coexistence_latent_energy')
  near(v.vapor_h_j_kg-v.vapor_u_j_kg,R*t,'ideal_gas_pressure_work')
  near(v.vapor_u_j_kg, v.atomic_reference_u0_j_kg+1.5*R*t,'constant_atomic_reference')
  near(flux(t,v.saturation_pressure_pa),0,'equilibrium_net_flux',absolute=1e-30)
  near(flux(t,2*v.saturation_pressure_pa),-flux(t),'signed_condensation_flux')
  check(v.latent_h_j_kg>0 and v.saturation_pressure_pa>0,'positive_vapor_branch')
  if 1800<t<3000:
   dt=.01;derivative=(math.log(gas(t+dt).saturation_pressure_pa)-math.log(gas(t-dt).saturation_pressure_pa))/(2*dt)
   near(R*t*t*derivative,v.latent_h_j_kg,'clausius_clapeyron_consistency',rel=1e-8)
   derivative=(gas(t+dt).latent_h_j_kg-gas(t-dt).latent_h_j_kg)/(2*dt)
   near(derivative,v.gas_cp_j_kg_k-794.96,'kirchhoff_consistency',rel=1e-8)
 near(gas(3000).latent_h_j_kg,7405680,'kim_latent_anchor')
 near(gas(3000).saturation_pressure_pa,gas(3000).kim_pressure_pa,'kim_pressure_anchor')
 for t in [1799.,math.nextafter(3000,math.inf),math.nan]:rejected(lambda:gas(t),'vapor_domain_rejection')
 for t,p,alpha in [(2000,-1,1),(2000,0,-1),(2000,0,1.1),(2000,math.nan,1)]:rejected(lambda:flux(t,p,alpha),'flux_input_rejection')
 # No measured collision diameter is supplied. Test algebra using a diameter
 # derived from the declared nominal reference atomic volume, labeled a model.
 M=.713*.055845+.202*.0519961+.085*.0586934;N=6.02214076e23;kb=1.380649e-23
 diameter=(6*M/(math.pi*8000*N))**(1/3)
 near(gk(3000,diameter),75*kb/(64*diameter**2)*math.sqrt(kb*3000/(math.pi*(M/N))),'neutral_hard_sphere_coefficient')
 near(gk(3000,2*diameter),gk(3000,diameter)/4,'gas_diameter_sensitivity')
 rejected(lambda:gk(3000,diameter,1e-20),'no_ionized_gas_transport_substitution')
 rejected(lambda:gk(3000,0),'no_invented_collision_diameter')
 # Signed metal transfers use actual packet kinetic energy, never mean2RT plus
 # actual kinetic energy. Reference energy remains invariant across T changes.
 T=2500.;mass=2e-8;u=ev(T).u_model_j_kg;b=State(mass,mass*u,mass/8000,7,100);original=bytes(b)
 dm=mass*1e-5;kinetic=dm*2*R*T
 tr=transfer(b,dm,kinetic)
 check(bytes(b)==original,'immutable_accepted_metal_state')
 near(tr.emitted_reference_plus_kinetic_j,dm*gas(T).atomic_reference_u0_j_kg+kinetic,'actual_packet_energy_once')
 near(tr.next.mass_kg+dm,b.mass_kg,'mass_inventory_closure')
 near(tr.next.energy_j+tr.emitted_reference_plus_kinetic_j,b.energy_j,'energy_inventory_closure')
 near(tr.removed_reference_volume_m3,dm/8000,'reference_removed_volume')
 near(tr.remaining_mass_fraction,tr.next.mass_kg/(8000*b.reference_volume_m3),'remaining_reference_mass_fraction')
 depleted=Values();rc(lib.ss304v2_evaluate_inventory(2,8000,2500,mass/2,b.reference_volume_m3,C.byref(depleted)))
 near(depleted.energy_j_m3_reference,ev(2500).energy_j_m3_reference/2,'depleted_mass_reference_energy')
 near(depleted.capacity_j_m3_reference_k,ev(2500).capacity_j_m3_reference_k/2,'depleted_mass_reference_capacity')
 check(depleted.k_w_m_k==ev(2500).k_w_m_k and depleted.sigma_s_m==ev(2500).sigma_s_m,'no_invented_effective_medium_transport')
 near(tr.next.mass_kg*ev(tr.temperature_k).u_model_j_kg,tr.next.energy_j,'caloric_projection',rel=8*EPS)
 check(tr.next.revision==8 and tr.next.last_event==101,'revision_watermark_advance')
 # Recondensing the SAME emitted mass/reference+actual kinetic packet reverses
 # the transfer. The gas owner must remove that packet exactly once.
 reverse=transfer(tr.next,-dm,-kinetic)
 near(reverse.next.mass_kg,b.mass_kg,'condensation_mass_roundtrip')
 near(reverse.next.energy_j,b.energy_j,'condensation_energy_roundtrip')
 for f,label in [
  (lambda:transfer(b,dm,kinetic,rev=6),'stale_revision_rejected'),
  (lambda:transfer(tr.next,dm,kinetic,event=101),'duplicate_event_rejected'),
  (lambda:transfer(b,mass,kinetic),'empty_cell_requires_geometry'),
  (lambda:transfer(b,dm,-kinetic),'wrong_signed_kinetic_rejected'),
  (lambda:transfer(b,dm,kinetic,q=-1e10),'cooling_outside_support_rejected'),
  (lambda:transfer(b,0,0,q=1e10),'heating_outside_support_rejected')]:rejected(f,label)
 sentinel=Trial();C.memset(C.byref(sentinel),0x2a,C.sizeof(sentinel));raw=bytes(sentinel)
 check(lib.ss304v2_transfer(8000,C.byref(b),6,101,dm,kinetic,0,C.byref(sentinel))!=0 and bytes(sentinel)==raw,'failed_transfer_output_atomicity')
 # Cold patch can lose real heat through299K without changing initial300K.
 cold=ev(300);target=ev(299);cb=State(mass,mass*cold.u_model_j_kg,mass/8000,0,0)
 cooled=transfer(cb,0,0,q=mass*(target.u_model_j_kg-cold.u_model_j_kg))
 near(cooled.temperature_k,299,'physical_cooling_below300_without_clamp',rel=4*EPS)
 for energy,status in [(mass*ev(298.2).h_j_kg-1e-6,-1),(mass*ev(3000).h_j_kg+1e-6,1),(0,0)]:
  budget=Support();rc(lib.ss304v2_energy_support(2,8000,mass,energy,C.byref(budget)))
  check(budget.range_status==status,'explicit_unsupported_energy_budget')
  near(budget.cold_margin_j,energy-mass*ev(298.2).h_j_kg,'cold_heat_margin')
  near(budget.hot_margin_j,mass*ev(3000).h_j_kg-energy,'hot_heat_margin')
 after=hold();assert all(sha(root/k)==v for k,v in seal['files'].items())
 assert source_before=={str(p.relative_to(S)):sha(p) for p in S.rglob('*') if p.is_file()}
 linked=subprocess.check_output(['ldd',str(B/'libss304v2.so')],text=True);assert not any(x in linked.lower() for x in ['cuda','cudss','mfem'])
 points=[]
 for t in [298.2,299.,300.,1500.,1800.,2000.,2400.,2800.,3000.]:
  v=ev(t);point={n:getattr(v,n) for n,_ in Values._fields_}
  if t>=1800:
   g=gas(t);point['vapor']={n:getattr(g,n) for n,_ in Gas._fields_};point['vacuum_flux_alpha1_kg_m2_s']=flux(t)
   point['mass_per_m2_over60ns_at_constantT']=flux(t)*60e-9
   point['mean_effusive_energy_cooling_per_m2_over60ns']=flux(t)*60e-9*(g.latent_h_j_kg-.5*R*t)
   point['raw_kim_liquid_k_w_m_k']=100*(.1229+.00003248*t)
   point['density_model_gas_k_w_m_k']=gk(t,diameter)
  points.append(point)
 result={'status':'passed_scoped_source_conditioned_material_checks_not_coupled60ns_validation','checks':sum(counts.values()),'groups':counts,'maximum_absolute_numeric_errors':maxima,'source_files_sha256':source_before,'library_sha256':sha(B/'libss304v2.so'),'compiler':str(compiler),'compiler_sha256':sha(compiler),'build_command':cmd,'linked_dependencies':linked,'v1_source_and_library_unchanged':True,'v1_overlap_bitwise_identical':True,'cold_bound_K':298.2,'energy_zero_K':300,'cold_energy_headroom_J_kg':-ev(298.2).u_model_j_kg,'source_pressure_implied_latent_J_kg':R*18836*math.log(10),'source_latent_J_kg':7405680,'atomic_mean_molar_mass_kg_mol':M,'gas_collision_diameter_model_m':diameter,'collision_diameter_is_not_measured':True,'kim_k_normalization_factor':28.99/(100*(.1229+.00003248*1800)),'points':points,'production_before':before,'production_after':after,'cpu_affinity':sorted(os.sched_getaffinity(0)),'nice':os.getpriority(os.PRIO_PROCESS,0),'lease':str(LOCK),'wall_s':time.time()-start,'gpu_or_production_signals':False,'existing_production_and_mechanical_guards_changed':False}
 (O/'verification.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({k:result[k] for k in ['status','checks','wall_s','cold_energy_headroom_J_kg','kim_k_normalization_factor']}))
if __name__=='__main__':
 os.sched_setaffinity(0,{7});priority=os.getpriority(os.PRIO_PROCESS,0)
 if priority<15:os.nice(15-priority)
 assert os.getpriority(os.PRIO_PROCESS,0)==15
 for n in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[n]='1'
 os.environ['CUDA_VISIBLE_DEVICES']=''
 with LOCK.open('a') as f:
  fcntl.flock(f,fcntl.LOCK_EX);verify()
