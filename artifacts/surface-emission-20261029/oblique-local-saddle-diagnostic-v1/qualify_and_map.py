import csv,hashlib,json,math,os,pathlib,time
import saddle_modes as sm

C=pathlib.Path(__file__).resolve().parent
O=pathlib.Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-oblique-emission-20261029')
O.mkdir(parents=True,exist_ok=True)
assert os.sched_getaffinity(0)=={7} and os.getpriority(os.PRIO_PROCESS,0)>=15
checks=0
def need(b,text):
 global checks
 assert b,text
 checks+=1

max_symplectic=max_diagonal=max_lorentz=0.
for bn in [-3.,-1.,-.1,-.001,.001,.1,1.,3.]:
 for bt in [0.,.2,.8,1.2,5.]:
  m=sm.modes(bn,bt)
  max_symplectic=max(max_symplectic,m['symplectic_residual'])
  max_diagonal=max(max_diagonal,m['diagonal_energy_residual'])
  max_lorentz=max(max_lorentz,m['independent_Lorentz_residual'])
  need(m['symplectic_residual']<1e-10,'Symplectic branch failed')
  need(m['diagonal_energy_residual']<1e-10,'Energy normal form failed')
  need(m['independent_Lorentz_residual']<1e-13,'Independent physical Lorentz generator failed')
  need(abs(m['unstable']**2*m['stable']**2-bn*bn)<1e-10*max(1,bn*bn),'Mode product invariant failed')
for b in [.01,.3,.9,2,5]:
 m=sm.modes(b,0);need(abs(m['unstable']-1)<1e-13 and abs(m['stable']-b)<1e-13,'Normal-only mode limit')
for b in [.1,.9]:
 m=sm.modes(0,b);need(abs(m['unstable']-math.sqrt(1-b*b))<1e-13 and m['stable']==0,'Weak tangent limit')
for b in [1.1,5.]:
 m=sm.modes(0,b);need(m['unstable']==0 and abs(m['stable']-math.sqrt(b*b-1))<1e-13,'Strong tangent stability limit')
for b in [.5,2.]:
 m=sm.modes(1e-6,b)
 target=math.sqrt(1-b*b) if b<1 else 1e-6/math.sqrt(b*b-1)
 need(abs(m['unstable']/target-1)<1e-8,'Small normal-B limiting branch')

# Check the analytic activated FD bound, then independently integrate a
# manufactured lower barrier where the FD/MB difference is resolvable.
for T in [300.,773.15,1490.]:
 r=sm.spectrum(T,97724679.40458213,27,14)
 need(r['FD_relative_current_bound']<1e-8,'Activated FD bound too large for high-field diagnostic')
 need(r['physical_current'] is None and not r['particle_phase_mapping_qualified'],'Rate diagnostic was relabeled as physical source')

def direct_fd(T,F,Bn,Bt,phi,panels):
 r=sm.spectrum(T,F,Bn,Bt,phi)
 t=sm.kb*T/sm.qe;Us=r['saddle_energy_eV'];omega=r['omega0_s-1'];lam=r['unstable_frequency_over_omega0']
 gamma=sm.hbar*omega*lam/(2*math.pi*sm.qe);delta=sm.hbar*omega*r['stable_frequency_over_omega0']/sm.qe
 count=math.ceil(48*t/delta);energies=[delta*(n+.5) for n in range(count)]
 lo=-Us-energies[-1]-48*gamma;hi=48*t
 current=heat=0.
 for i in range(panels+1):
  Eu=lo+(hi-lo)*i/panels;x=Eu/gamma
  D=1/(1+math.exp(-x)) if abs(x)<700 else (1. if x>0 else 0.)
  w=1 if i in [0,panels] else (4 if i%2 else 2)
  for stable in energies:
   total=Us+Eu+stable;y=total/t
   f=1/(1+math.exp(y)) if abs(y)<700 else (0. if y>0 else 1.)
   current+=w*D*f;heat+=w*D*f*total
 factor=sm.C*lam*delta*(hi-lo)/(3*panels)
 return current*factor,heat*factor

T=1490.;F=97724679.40458213;Bn=700.;Bt=14.;t=sm.kb*T/sm.qe
phi=2*math.sqrt(sm.c*F)+3*t
manufactured=sm.spectrum(T,F,Bn,Bt,phi)
coarse=direct_fd(T,F,Bn,Bt,phi,8192);fine=direct_fd(T,F,Bn,Bt,phi,16384)
need(abs(coarse[0]/fine[0]-1)<1e-6,'Direct FD quadrature refinement failed')
need(0<=manufactured['activated_MB_model_current_A_m2']-fine[0]<=manufactured['FD_current_difference_upper_A_m2'],'Direct FD current falls outside derived MB difference bound')
need(abs(-manufactured['activated_MB_complete_heat_W_m2']-fine[1])<=manufactured['FD_absolute_heat_difference_upper_W_m2'],'Direct FD complete heat falls outside derived bound')
fd_check={'manufactured_not_SS304':True,'T_K':T,'F_V_m':F,'B_normal_T':Bn,'B_tangent_T':Bt,'phi_eV':phi,
          'direct_FD_current_A_m2':fine[0],'direct_FD_energy_flux_W_m2':fine[1],
          'MB_model_current_A_m2':manufactured['activated_MB_model_current_A_m2'],
          'current_difference_upper_A_m2':manufactured['FD_current_difference_upper_A_m2'],
          'heat_difference_upper_W_m2':manufactured['FD_absolute_heat_difference_upper_W_m2'],
          'direct_quadrature_relative_current_difference':abs(coarse[0]/fine[0]-1)}

joint=pathlib.Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-magnetized-emission-20261029/source-backup-tangential-onepass-v1/evidence')
comparison=[]
for T in [298.2,300.,773.15,1490.]:
 full0=json.loads((joint/f'joint-{int(T)}-0-2.json').read_text())
 fullB=json.loads((joint/f'joint-{int(T)}-30-2.json').read_text())
 p0=sm.spectrum(T,97724679.40458213,0,0)
 pB=sm.spectrum(T,97724679.40458213,0,30)
 full_ratio=fullB['known_component_current_A_m2']/full0['known_component_current_A_m2']
 parabolic_ratio=pB['activated_MB_model_current_A_m2']/p0['activated_MB_model_current_A_m2']
 comparison.append({'T_K':T,'full_SN_tangent_ratio':full_ratio,'quadratic_tangent_ratio':parabolic_ratio,
  'magnetic_ratio_disagreement':abs(full_ratio-parabolic_ratio),
  'quadratic_B0_relative_current_error_vs_full_SN':p0['activated_MB_model_current_A_m2']/full0['known_component_current_A_m2']-1,
  'quadratic_mean_energy_B0_error_eV':p0['activated_MB_model_mean_Fermi_energy_eV']-full0['known_component_mean_metal_energy_eV'],
  'quadratic_is_not_full_SN_error_bound':True})

probe=pathlib.Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs/cathode-emission-audit-20261029/fields/probes.csv')
rows=list(csv.DictReader(probe.open()));maps=[]
for index,row in enumerate(rows):
 n=[float(row[x]) for x in ['nx','ny','nz']]
 E=[float(row[x]) for x in ['Ex_V_m','Ey_V_m','Ez_V_m']]
 B=[float(row[x]) for x in ['Bx_T','By_T','Bz_T']]
 F=-sum(x*y for x,y in zip(n,E));Bn=sum(x*y for x,y in zip(n,B));Bt=math.sqrt(max(0,sum(x*x for x in B)-Bn*Bn))
 for T in [298.2,300.,773.15,1490.]:
  if F<=0:record={'supported_local_unstable_channel':False,'reason':'Retarding surface has no SN saddle; thermal supply still required','physical_current':None,'physical_foil_model_admitted':False}
  else:record=sm.spectrum(T,F,Bn,Bt)
  record.update({'legacy_row':index,'face':row['name'],'T_K':T,'signed_F_V_m':F,'B_normal_T':Bn,'B_tangent_T':Bt,
                 'actual_E_tangent_over_abs_F':math.sqrt(sum((e+F*nn)**2 for e,nn in zip(E,n)))/max(abs(F),1),
                 'legacy_field_is_finite_offset_not_qualified_native_boundary':True})
  maps.append(record)
need(len(maps)==4224,'Legacy regime coverage failed')
(O/'local-saddle-regime-map.jsonl').write_text(''.join(json.dumps(x)+'\n' for x in maps))
report={'schema':'oblique-local-quadratic-mode-diagnostic-v1','status':'DIAGNOSTIC_ONLY_NOT_EMISSION_SOURCE',
 'created_unix_s':time.time(),'checks':checks,'maximum_symplectic_residual':max_symplectic,
 'maximum_diagonal_energy_residual':max_diagonal,'maximum_independent_Lorentz_residual':max_lorentz,
 'branch_tests_include_both_signs_of_1_minus_cyclotron_squared':True,
 'independent_manufactured_FD_integral_check':fd_check,
 'tangent_full_SN_comparison':comparison,'legacy_regime_states':len(maps),'probe_sha256':hashlib.sha256(probe.read_bytes()).hexdigest(),
 'pure_tangent_strong_field_no_unstable_mode_is_not_zero_thermal_surface_supply':True,
 'current_scope':'Activated local parabolic mode spectrum, orbital two-degenerate-spin-label model with no Zeeman splitting; FD difference bounded within that model',
 'spin_or_magnetic_band_error_bound':None,'full_SN_anharmonic_error_bound':None,'joint_physical_phase_projection_qualified':False,
 'inverse_matching_qualified':False,'native_midpoint_qualified':False,'physical_foil_model_admitted':False,'GPU_executed':False,'production_resumed':False}
(O/'qualification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
