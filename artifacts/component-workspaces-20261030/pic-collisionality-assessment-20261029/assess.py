#!/usr/bin/env python3
"""Frozen electron-field/collision scale assessment; no solve or evolution."""
import argparse,csv,hashlib,json,math,os,pathlib,resource,subprocess,time
import numpy as np
QE=1.602176634e-19;ME=9.1093837139e-31;C=299792458.;KB=1.380649e-23;EPS0=8.8541878128e-12;H=6.62607015e-34
BASE=pathlib.Path('/home/treefor/sim-workspace/powerflow/outputs')
FIX=BASE/'pic-profile/cuda-gather-investigation/actual-midpoint-fixture'
FD=BASE/'cathode-emission-audit-20261029/fields'
DATA=pathlib.Path('/home/treefor/projects/runtimes/pic-native-spatial-diagnostics-20261029/arena-v5/inputs/h2o-gas-20261029.json')
OWNER=pathlib.Path('/home/treefor/projects/runtimes/pic-electrode-heat-feedback-20261029/native/mfem_electrode_heat_feedback_v2')
QUERY=BASE/'h2o-collision-data-20261029/attempt-003/h2o-collision-query'
DIAG=BASE/'pic-native-spatial-diagnostics-20261029/qualified/preserved-midpoint/diagnostics.json'
CHANNELS=['H2O:ion_yield:H2O+','H2O:ion_yield:total','H2O:elastic_momentum','H2O:munoz2007:recommended_total','H2O:nist_BEB_valence','H2O:excite_v010']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def norm(x):return np.sqrt(np.sum(x*x,axis=-1))
def evaluate(channel,E):
 E=np.asarray(E,dtype=float);assert np.all(np.isfinite(E)) and np.all(E>=0)
 lo,hi=channel['support_eV'];zero=channel['zero_below_eV'];known=(E>=lo)&(E<=hi);sigma=np.full(E.shape,np.nan)
 z=(E<=zero)&(E<=hi) if zero is not None else np.zeros(E.shape,dtype=bool);known|=z;sigma[z]=0;active=known&~z
 if channel['model']=='linear_table':
  x=np.array(channel['energy_eV']);y=np.array(channel['values']);k=np.minimum(np.searchsorted(x,E[active],side='right')-1,len(x)-2);q=(E[active]-x[k])/(x[k+1]-x[k]);sigma[active]=((1-q)*y[k]+q*y[k+1])*channel['ordinate_scale_to_SI']
 elif channel['model']=='nist_BEB':
  result=np.zeros(E.shape)
  for orbital in channel['orbitals']:
   B=orbital['B_eV'];mask=active&(E>B);t=E[mask]/B;u=orbital['U_eV']/B;Q=orbital['Q'];N=orbital['N'];a=channel['bohr_radius_m'];R=channel['rydberg_eV'];S=4*math.pi*a*a*N*(R/B)*(R/B);l=np.log(t);result[mask]+=S/(t+u+1)*(Q*l/2*(1-1/(t*t))+(2-Q)*(1-1/t-l/(t+1)))
  sigma[active]=result[active]
 else:raise ValueError('Unsupported observable model')
 assert np.all(np.isfinite(sigma[known])) and np.all(sigma[known]>=0);return known,sigma
def csvwrite(p,rows):
 if not rows:return
 with p.open('w') as stream:
  writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--out',type=pathlib.Path,required=True);a=parser.parse_args();assert not a.out.exists();assert os.sched_getaffinity(0)=={7} and os.getpriority(os.PRIO_PROCESS,0)>=15
 resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));a.out.mkdir();start=time.monotonic()
 pins={FIX/name:expected for name,expected in json.loads((FIX/'sha256.json').read_text()).items()};pins.update({p:sha(p) for p in [DATA,OWNER/'water-scenario-v2.json',OWNER/'water_scenario.hpp',QUERY,FD/'probes.csv',FD/'field-geometry.json',DIAG,pathlib.Path(__file__)]});assert all(sha(p)==h for p,h in pins.items())
 manifest=json.loads((FIX/'manifest.json').read_text());source=json.loads((OWNER/'water-scenario-v2.json').read_text());dataset=json.loads(DATA.read_text());channels={c['id']:c for c in dataset['channels']};diag=json.loads(DIAG.read_text())
 pt=np.dtype([('id','<u8'),('r','<f8',(3,)),('u','<f8',(3,)),('q','<f8'),('m','<f8'),('w','<f8')]);ft=np.dtype([('E','<f8',(3,)),('B','<f8',(3,)),('post_u','<f8',(3,)),('status','<u8')]);p=np.fromfile(FIX/'particles.bin',dtype=pt);fields=np.fromfile(FIX/'reference.bin',dtype=ft);assert len(p)==len(fields)==46726 and np.all(fields['status']==0)
 r,u,w=p['r'],p['u'],p['w'];assert np.all(w>0) and np.all(p['m']==ME) and np.all(p['q']==-QE) and all(np.all(np.isfinite(x)) for x in [r,u,w,fields['E'],fields['B']]);W=math.fsum(w)
 n2=(u[:,0]*u[:,0]+u[:,1]*u[:,1])+u[:,2]*u[:,2];gamma=np.sqrt(1+n2/(C*C));v=u/gamma[:,None];speed=norm(v);energy=p['m']*n2/(gamma+1)/QE
 E,B=fields['E'],fields['B'];bm=norm(B);assert np.all(bm>0);b=B/bm[:,None];ep=np.sum(E*b,axis=1);et=norm(E-ep[:,None]*b);up=np.sum(u*b,axis=1);rho=p['m']*norm(u-up[:,None]*b)/(QE*bm);omega=QE*bm/(gamma*ME);dt=manifest['dt'];beta=et/(C*bm);assert np.all(beta<1)
 coldK=2*ME*C*C/QE*beta**2/(1-beta**2);coldx=2*ME*et/(QE*bm*bm)/(1-beta**2);parallel_threshold_length=np.divide(12.621,np.abs(ep),out=np.full_like(ep,np.inf),where=ep!=0)
 distance=np.minimum(*[np.hypot(np.maximum(np.maximum(x0-r[:,0],0),r[:,0]-x1),np.maximum(np.maximum(y0-r[:,1],0),r[:,1]-y1)) for x0,x1,y0,y1 in [(.002,.007,.0018,.002),(.002,.007,.007,.0072)]])
 regions={'all':np.ones(len(p),bool),'above_upper_back_over_width':(r[:,1]>.0072)&(r[:,0]>=.002)&(r[:,0]<=.007),'gap_over_width':(r[:,1]>=.002)&(r[:,1]<=.007)&(r[:,0]>=.002)&(r[:,0]<=.007),'near_any_foil_25um':distance<=25e-6,'near_any_foil_100um':distance<=100e-6}
 def fraction(mask,region=None):
  if region is None:region=np.ones(len(w),bool)
  denom=math.fsum(w[region]);return math.fsum(w[region&mask])/denom if denom else None
 def stats(values,mask=None):
  if mask is None:mask=np.ones(len(p),bool)
  ww=w[mask];xx=values[mask];assert np.all(np.isfinite(xx));total=math.fsum(ww)
  if not total:return {'positive_weight':0,'count':0}
  order=np.argsort(xx,kind='stable');cum=np.cumsum(ww[order]);return dict(count=len(xx),positive_weight=total,min=float(xx.min()),weighted_p05_median_p95=[float(xx[order[min(int(np.searchsorted(cum,q*total)),len(xx)-1)]]) for q in [.05,.5,.95]],max=float(xx.max()),weighted_mean=math.fsum(ww*xx)/total)
 field_summary={name:dict(population_fraction=fraction(mask),energy_eV=stats(energy,mask),B_T=stats(bm,mask),signed_Eparallel_V_m=stats(ep,mask),abs_Eparallel_V_m=stats(abs(ep),mask),Eperp_V_m=stats(et,mask),abs_Eparallel_over_Emag=stats(abs(ep)/norm(E),mask),cold_uniform_crossed_Kmax_eV=stats(coldK,mask),cold_uniform_crossed_excursion_m=stats(coldx,mask),fraction_cold_crossed_Kmax_above_parent_threshold=fraction(coldK>=12.621,mask),parallel_threshold_distance_m=stats(parallel_threshold_length,mask),rho_m=stats(rho,mask),dx_over_rho=stats(25e-6/rho,mask),omega_c_dt=stats(omega*dt,mask),gyroperiod_s=stats(2*math.pi/omega,mask)) for name,mask in regions.items()}
 # Direct parity against the sealed observational C++ query, including gaps.
 parity=[]
 for cid in CHANNELS:
  channel=channels[cid];lo,hi=channel['support_eV'];samples=sorted(set([0.,float(lo),float(hi),float(hi+1),float((lo+hi)/2)]+([float(channel['zero_below_eV']),float((channel['zero_below_eV']+lo)/2)] if channel['zero_below_eV'] is not None else [])))
  for sample in samples:
   covered,values=evaluate(channel,[sample]);run=subprocess.run([str(QUERY),str(DATA),cid,format(sample,'.17g')],text=True,capture_output=True,timeout=10,check=False);assert run.returncode==(0 if covered[0] else 78),(cid,sample,run.stderr)
   relative_error=None
   if covered[0]:
    expected=json.loads(run.stdout)['value'];error=abs(values[0]-expected);assert error<=128*np.finfo(float).eps*max(expected,abs(values[0]),1e-300),(cid,sample,error);relative_error=float(error/max(expected,1e-300))
   parity.append(dict(channel=cid,energy_eV=sample,covered=bool(covered[0]),query_returncode=run.returncode,relative_error=relative_error))
 # Stationary ground-state H2O hypothesis. T below sets ideal-gas density only.
 collision={};coefficients={};region_rates=[];sigma_arrays={}
 for cid in CHANNELS:
  channel=channels[cid];known,sigma=evaluate(channel,energy);sigma_arrays[cid]=(known,sigma);coefficients[cid]=math.fsum(w[known]*sigma[known]*speed[known])/W
  collision[cid]=dict(metadata={k:channel.get(k) for k in ['observable','model','support_eV','zero_below_eV','source','source_locator','limitations','use']},supported_population_weight=math.fsum(w[known]),unsupported_population_weight=math.fsum(w[~known]),unsupported_population_fraction=fraction(~known),supported_positive_sigma_population_fraction=fraction(known&(np.nan_to_num(sigma)>0)),known_support_coefficient_m3_s=coefficients[cid],known_support_mean_sigma_m2=math.fsum(w[known]*sigma[known])/W,normalization='ALL electron weight; unsupported contribution remains unknown, not zero physical rate',complete_stationary_target_observable=bool(np.all(known)),attenuation_interpretation='Known-support single-scattering probability under declared straight/frozen stationary-target assumptions' if channel['observable']=='total_scattering' else 'Formal 1-exp(-tau_observable) scale only; yield/momentum/model observables are not admitted exclusive collision hazards',exclusive_event_or_full_mechanism=False)
  for name,mask in regions.items():
   denominator=math.fsum(w[mask]);region_rates.append(dict(region=name,channel=cid,population_weight=denominator,unsupported_fraction=fraction(~known,mask),known_support_coefficient_m3_s=math.fsum(w[mask&known]*sigma[mask&known]*speed[mask&known])/denominator if denominator else None))
 background=[]
 for pressure in [3e-5,5e-5,1e-4]:
  for T in [250.,300.,350.]:
   density=pressure*133.32236842105263/(KB*T)
   for cid,(known,sigma) in sigma_arrays.items():
    row=dict(pressure_Torr=pressure,T_K=T,H2O_fraction_assumed=1.,density_m3=density,density_cm3=density/1e6,channel=cid,partial_known_support_mean_observable_frequency_s=density*coefficients[cid],unsupported_population_fraction=fraction(~known),attenuation_interpretation=collision[cid]['attenuation_interpretation'])
    for L in [25e-6,100e-6,.001,.005]:row['partial_formal_attenuation_straight_fixed_energy_L_'+format(L,'.7g')+'_m']=math.fsum(w[known]*(-np.expm1(-density*sigma[known]*L)))/W
    for t in [1e-11,1e-9,manifest['absolute_native_time_s']]:row['partial_formal_attenuation_frozen_energy_t_'+format(t,'.7g')+'_s']=math.fsum(w[known]*(-np.expm1(-density*sigma[known]*speed[known]*t)))/W
    background.append(row)
 # Explicit finite-column scenarios, not a reconstruction of missing neutrals.
 columns=[];column_resolution=[];ML=source['provenance']['literature_facts']['monolayer_molecules_m2'];theta0=source['provenance']['selected_model_assumptions']['initial_coverage']
 for theta in [.06,theta0,.67,1.]:
  for released in [.01,.1,1.]:
   column=theta*ML*released
   for cid,(known,sigma) in sigma_arrays.items():columns.append(dict(theta0_assumed=theta,released_fraction_assumed=released,column_m2=column,channel=cid,unsupported_population_fraction=fraction(~known),partial_known_support_mean_optical_depth=column*collision[cid]['known_support_mean_sigma_m2'],partial_known_support_formal_single_normal_crossing_attenuation=math.fsum(w[known]*(-np.expm1(-column*sigma[known])))/W,attenuation_interpretation=collision[cid]['attenuation_interpretation']))
   for delta in [.1e-6,1e-6,12.5e-6,25e-6,100e-6,.001]:column_resolution.append(dict(theta0_assumed=theta,released_fraction_assumed=released,column_m2=column,layer_thickness_m_assumed=delta,neutral_density_cm3=column/delta/1e6,layer_cells_at_25um=delta/25e-6,peak_density_dilution_if_smeared_to_25um=25e-6/delta if delta<25e-6 else 1.,geometric_actual_electron_weight_fraction_within_layer=fraction(distance<=delta),known_support_total_scattering_mean_nu_dt=(column/delta)*coefficients['H2O:munoz2007:recommended_total']*dt,known_support_BEB_valence_mean_observable_nu_dt=(column/delta)*coefficients['H2O:nist_BEB_valence']*dt,inverse_known_support_total_mean_opacity_m=delta/(column*collision['H2O:munoz2007:recommended_total']['known_support_mean_sigma_m2'])))
 # Pointwise source-law evaluation; this does not step a reservoir.
 water_mass=18.0106*1.66053906660e-27;law=[]
 for T in [300.,500.,700.,1000.,1300.]:
  barrier=.87*(1-.5*theta0);hazard=KB*T/H*math.exp(-barrier*QE/(KB*T));vn=math.sqrt(math.pi*KB*T/(2*water_mass));law.append(dict(scenario='source-v2 ideal-oxide Reviewed300To1300 theta0=0.33',T_K=T,theta_start_assumed=theta0,barrier_eV=barrier,rate_s=hazard,rounded_frozen_start_released_fraction_over_1ps=-math.expm1(-hazard*1e-12),positive_remaining_fraction_over_1ps=math.exp(-hazard*1e-12),rounded_frozen_start_released_fraction_over_1ns=-math.expm1(-hazard*1e-9),positive_remaining_fraction_over_1ns=math.exp(-hazard*1e-9),neutral_mass_assumed_kg=water_mass,effusive_mean_outward_normal_speed_m_s=vn,normal_flight_scale_over_1ns_m=vn*1e-9,normal_flight_scale_over_capture_age_m=vn*manifest['absolute_native_time_s']))
 # Local second-central-velocity-moment proxies: no Maxwellian closure inferred.
 geometry=diag['diagnostic_geometry'];xs=np.array(geometry['xs_m']);ys=np.array(geometry['ys_m']);nx=len(xs)-1;ix=np.minimum(np.searchsorted(xs,r[:,0],side='right')-1,nx-1);iy=np.minimum(np.searchsorted(ys,r[:,1],side='right')-1,len(ys)-2);assert np.all(ix>=0) and np.all(iy>=0);cell=iy*nx+ix;volume_by_cell={z['cell_id']:z['physical_reaction_volume_m3'] for z in diag['species_cell_densities']};density_by_cell={z['cell_id']:z['density_cm3']*1e6 for z in diag['species_cell_densities']};cellrows=[]
 for cellid in np.unique(cell):
  selected=cell==cellid;ww=w[selected];vv=v[selected];total=math.fsum(ww);mean=np.array([math.fsum(ww*vv[:,j])/total for j in range(3)]);central=vv-mean;Tcomponents=[ME*math.fsum(ww*central[:,j]*central[:,j])/total/QE for j in range(3)];Tproxy=math.fsum(Tcomponents)/3;ne=density_by_cell[int(cellid)];op=math.sqrt(ne*QE*QE/(EPS0*ME));lam=math.sqrt(EPS0*Tproxy*QE/(ne*QE*QE)) if Tproxy>0 else None
  cellrows.append(dict(cell_id=int(cellid),carrier_records=int(selected.sum()),physical_electron_count=total,ne_cm3=ne/1e6,central_velocity_Tx_proxy_eV=Tcomponents[0],central_velocity_Ty_proxy_eV=Tcomponents[1],central_velocity_Tz_proxy_eV=Tcomponents[2],scalar_T_second_moment_proxy_eV=Tproxy,lambda_D_Maxwellian_proxy_m=lam,lambda_D_proxy_over_dx=lam/25e-6 if lam else None,omega_pe_classical_s=op,omega_pe_dt=op*dt,plasma_period_classical_s=2*math.pi/op,Debye_sphere_electron_count_proxy=(4*math.pi/3)*ne*lam**3 if lam else None))
 plasma_scenarios=[]
 for delta in [1e-6,25e-6,100e-6]:
  nn=theta0*ML/delta
  for fi in [1e-6,1e-3,1.]:
   ne=nn*fi
   for T in [1.,10.,100.]:
    lam=math.sqrt(EPS0*T/(ne*QE));op=math.sqrt(ne*QE*QE/(EPS0*ME));plasma_scenarios.append(dict(scenario='source-v2 theta0=.33 fully released uniform layer; hypothetical ionized fraction and Maxwellian T',layer_thickness_m=delta,ionized_fraction_assumed=fi,T_eV_assumed=T,ne_cm3=ne/1e6,lambda_D_m=lam,lambda_D_over_25um=lam/25e-6,omega_pe_s=op,omega_pe_dt=op*dt,plasma_period_s=2*math.pi/op,Debye_sphere_count=(4*math.pi/3)*ne*lam**3))
 # The old static/cold planar control is a comparison only; actual Eparallel !=0.
 voltage=abs(manifest['voltage']);gap=.005;BH=math.sqrt(2*ME*voltage/QE+voltage*voltage/(C*C))/gap
 probes=[]
 with (FD/'probes.csv').open() as stream:
  for row in csv.DictReader(stream):
   ee=np.array([float(row[k]) for k in ['Ex_V_m','Ey_V_m','Ez_V_m']]);bb=np.array([float(row[k]) for k in ['Bx_T','By_T','Bz_T']]);Bmag=float(norm(bb));eparallel=float(np.dot(ee,bb)/Bmag);eperp=float(norm(ee-eparallel*bb/Bmag));betaD=eperp/(C*Bmag)
   probes.append(dict(name=row['name'],x_m=float(row['x_m']),y_m=float(row['y_m']),z_m=float(row['z_m']),vacuum_offset_m=float(row['offset_m']),tangent_edge_distance_m=float(row['tangent_edge_distance_m']),B_T=Bmag,Eparallel_V_m=eparallel,Eperp_V_m=eperp,cold_crossed_Kmax_eV=2*ME*C*C/QE*betaD**2/(1-betaD**2) if betaD<1 else None,parent_threshold_parallel_distance_m=12.621/abs(eparallel) if eparallel else None,native_dx_m=float(row['native_dx_m']),native_dy_m=float(row['native_dy_m'])))
 # Basic analytic controls and count conservation; all scientific source files pinned.
 assert math.isclose(W,diag['species_counts']['electron']['rounded'],rel_tol=3e-15);assert np.all(energy>=0) and np.all(np.isfinite(coldK));assert math.isclose(BH*BH*gap*gap,2*ME*voltage/QE+voltage*voltage/C**2,rel_tol=2e-15);assert all(sha(p)==h for p,h in pins.items())
 for name,rows in [('background-scenarios.csv',background),('finite-column-optical-depth.csv',columns),('finite-layer-resolution.csv',column_resolution),('source-v2-pointwise-law.csv',law),('local-cell-resolution-proxies.csv',cellrows),('hypothetical-plasma-resolution.csv',plasma_scenarios),('region-rate-support.csv',region_rates),('controlled-probe-scales.csv',probes)]:csvwrite(a.out/name,rows)
 result=dict(schema='frozen_electron_collisionality_resolution/v1',status='conditional_readonly_assessment',time_s=manifest['absolute_native_time_s'],fixture_scope=manifest['capture_stage'],physical_electron_weight=W,energy_distribution=dict(min_eV=float(energy.min()),max_eV=float(energy.max()),fraction_below_parent_threshold=fraction(energy<=12.621),fraction_above_1keV=fraction(energy>1000)),collision_observables=collision,particle_field_regions=field_summary,background_assumption=dict(central_P_Torr=5e-5,central_T_K=300,H2O_fraction=1.,water_partial_pressure_fraction_unknown=True,RGA_infers_no_surface_coverage=True,stationary_ground_state_target_approximation=True,thermal_target_velocity_not_convolved=True),source_v2_scenario=dict(file=str(OWNER/'water-scenario-v2.json'),sha256=sha(OWNER/'water-scenario-v2.json'),initial_coverage_assumed=theta0,monolayer_m2=ML,range='Reviewed300To1300',legacy_fixture_contains_no_water=True,finite_columns_are_hypotheses=True),conditional_planar_Hull=dict(gap_m=gap,voltage_magnitude_V=voltage,B_H_T=BH,B_over_Hull=stats(bm/BH),scope='Uniform static planar cold E perpendicular B only. Actual finite-edge Eparallel and pulse invalidate use as a crossing switch.'),resolution_summary=dict(kinetic_field_cell_m=25e-6,actual_dt_s=dt,local_cell_count=len(cellrows),max_ne_cm3=max(z['ne_cm3'] for z in cellrows),positive_variance_cell_count=sum(z['scalar_T_second_moment_proxy_eV']>0 for z in cellrows),min_positive_lambda_D_Maxwellian_proxy_m=min(z['lambda_D_Maxwellian_proxy_m'] for z in cellrows if z['lambda_D_Maxwellian_proxy_m']),max_omega_pe_dt=max(z['omega_pe_dt'] for z in cellrows),nonneutral_nonMaxwellian=True,proxy_not_actual_shielding=True),data_sources=dataset['sources'],limits=['No neutral/ion distribution or trajectory history exists in the legacy electron fixture','Known-support rates/optical depths retain all electron weight in denominator; unsupported rate contribution is unknown','Observables overlap and must not be summed into an exclusive event frequency; total ion yield is not event multiplicity','Stationary targets are an explicitly conditional diagnostic, not the thermal moving-target collision kernel','Column normal-crossing optical depth is thickness-independent only for the declared straight constant-energy full-column crossing','Cell averaging can dilute near-surface peak neutral density while preserving column inventory; reaction overlap correlations and electron energy must still be resolved','Actual 25um field and coarse native magnetic cells do not gain resolution from 1um face probes','No forced shielding, complete mechanism, unique halo diagnosis or evolving/domain-converged claim'],execution=dict(affinity=sorted(os.sched_getaffinity(0)),nice=os.getpriority(os.PRIO_PROCESS,0),wall_s=time.monotonic()-start,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,address_space_cap_bytes=4*1024**3,GPU=False,Poisson=False,particle_evolution=False),inputs={str(p):h for p,h in pins.items()})
 write(a.out/'assessment.json',result);write(a.out/'strict-query-parity.json',dict(status='PASS',samples=parity,policy='Same supported zeros and rejected gaps; <=128 epsilon arithmetic parity',all_inputs_unchanged=True));print(json.dumps(dict(status=result['status'],energy=result['energy_distribution'],resolution=result['resolution_summary'],channels={c:dict(unsupported_fraction=collision[c]['unsupported_population_fraction'],known_support_k_m3_s=collision[c]['known_support_coefficient_m3_s']) for c in CHANNELS},execution=result['execution']),indent=2))
if __name__=='__main__':main()
