from pathlib import Path
import json,math,sys
HERE=Path(__file__).resolve().parent
run=Path(sys.argv[1]);rows=json.loads((run/'results.json').read_text());limits=json.loads((HERE/'controls.json').read_text())['limits'];checks=[]
def check(name,passed,**details):checks.append(dict(name=name,passed=bool(passed),**details))
def bound(name,value,key):check(name,math.isfinite(value) and value<=limits[key],measured=value,limit=limits[key])
def order(name,a,b,factor,key):
    value=math.log(a/b)/math.log(factor)
    check(name,value>=limits[key],measured=value,minimum=limits[key]);return value

t=next(r for r in rows if r['family']=='transactions')
for k in ['invalid_second_product_rollback','field_cache_invalidated','duplicate_event_rejected','wrong_species_restart_rejected','step_cap_rejected','corrupt_checkpoint_rejected']:check(k,t[k])
check('all products enter one field solve',t['collocated_pair_field_J']==0 and t['charge_C']==0 and t['single_product_spurious_field_J']>0)
check('failed batch preserves ID allocation',t['first_ID']==1 and t['second_ID']==2)
w=next(r for r in rows if r['family']=='signed_walls')
check('signed species wall collection',w['electron_lower_charge_C']<0<w['ion_upper_charge_C'] and w['resident_count']==0 and w['charge_defect_C']==0)
bound('wall kinetic throughput',w['relative_kinetic_throughput_defect'],'kick_work_relative')
for sign in [-1,1]:
    own=sorted([r for r in rows if r['family']=='uniform_drift' and r['B_T']*sign>0],key=lambda r:r['steps'])
    bound(f'drift coarse sign {sign}',abs(own[0]['relative_total_current_error']),'drift_coarse_relative')
    bound(f'drift fine sign {sign}',abs(own[-1]['relative_total_current_error']),'drift_fine_relative')
    order(f'drift time refinement sign {sign}',abs(own[0]['relative_total_current_error']),abs(own[-1]['relative_total_current_error']),4,'drift_min_order')
    for r in own:
        check(f'signed currents {sign} {r["steps"]}',sign*r['Jy_A_m2'][0]>0>sign*r['Jy_A_m2'][1])
        for s in range(2):bound(f'drift work {sign} {r["steps"]} species {s}',abs(r['kick_work_defect_J_per_particle'][s])/r['work_activity_J_per_particle'][s],'kick_work_relative')
common=next(r for r in rows if r['family']=='balanced_common_drift')
check('common balanced drift carries no relative current',common['velocity_error_m_s']<1e-10)
native=[r for r in rows if r['family']=='native_wave']
for r in native:
    name=r['name']
    for k,lim in [('max_energy_relative_defect','native_energy_relative'),('max_charge_relative_defect','native_charge_relative'),('max_gauss_relative_defect','native_gauss_relative'),('max_species_continuity_relative_defect','native_species_continuity_relative'),('max_native_path_parity_relative_defect','native_path_parity_relative'),('max_poisson_relative_residual','native_Poisson_residual_relative'),('kick_work_relative_defect','kick_work_relative')]:bound(name+' '+k,r[k],lim)
    check(name+' both species move',all(x!=0 for x in r['final_xi_m']))
    check(name+' gyro/plasma/displacement/10fs guards active',r['dt_s']<=1e-14)
    if r['restart_tested']:
        bound('restart position',r['restart_position_error_m'],'restart_position_m')
        bound('restart proper velocity',r['restart_velocity_error_m_s'],'restart_velocity_m_s')
        bound('restart work ledger',r['restart_work_error_J']/r['U0_J'],'restart_work_relative')
mesh=sorted([r for r in native if r['steps']==400 and r['ppc_y']==4 and r['mass_ratio']==25],key=lambda r:r['ny'])
bound('fine Poisson potential',mesh[-1]['initial_phi_relative_L2_error'],'native_phi_fine_relative_L2')
order('Poisson spatial order',mesh[0]['initial_phi_relative_L2_error'],mesh[-1]['initial_phi_relative_L2_error'],4,'native_min_mesh_order')
for s in range(2):
    bound(f'fine species {s} displacement all samples',mesh[-1]['max_displacement_error_scaled'][s],'native_species_displacement_fine_scaled')
    bound(f'fine species {s} velocity all samples',mesh[-1]['max_velocity_error_scaled'][s],'native_species_displacement_fine_scaled')
    order(f'species {s} displacement spatial order',mesh[0]['max_displacement_error_scaled'][s],mesh[-1]['max_displacement_error_scaled'][s],4,'native_min_mesh_order')
time=sorted([r for r in native if r['ny']==64 and r['ppc_y']==4],key=lambda r:r['steps'])
order('native energy time order',time[0]['max_energy_relative_defect'],time[-1]['max_energy_relative_defect'],4,'native_min_energy_time_order')
ppc=sorted([r for r in native if r['ny']==32 and r['mass_ratio']==25],key=lambda r:r['ppc_y'])
for s in range(2):order(f'sampling convergence toward fixed-grid limit species {s}',abs(ppc[0]['final_xi_m'][s]-ppc[1]['final_xi_m'][s]),abs(ppc[1]['final_xi_m'][s]-ppc[2]['final_xi_m'][s]),2,'native_min_mesh_order')
em=[r for r in rows if r['family']=='em_wave']
for r in em:
    name=r['name'];check(name+' bounded perturbation',r['delta_B_over_B0']<1e-3)
    if r['current_feedback'] and r['faraday_feedback']:
        bound(name+' energy',r['max_energy_relative_defect'],'EM_energy_relative')
        bound(name+' work',r['max_midpoint_work_relative_defect'],'EM_work_relative')
        bound(name+' analytic midpoint frequency',r['midpoint_relative_error'],'EM_midpoint_frequency_relative')
    else:check(name+' feedback ablation changes mode',r['endpoint_E_relative_mismatch']>=limits['EM_ablation_E_mismatch_min'],measured=r['endpoint_E_relative_mismatch'])
meshem=sorted([r for r in em if r['helicity']==1 and r['dt_omega_pe']==.01 and r['current_feedback'] and r['faraday_feedback']],key=lambda r:r['N'])
bound('finest Maxwell dispersion',meshem[-1]['dispersion_relative_error'],'EM_finest_dispersion_relative')
order('Maxwell spatial order',meshem[0]['dispersion_relative_error'],meshem[-1]['dispersion_relative_error'],8,'EM_min_mesh_order')
timeem=sorted([r for r in em if r['N']==64 and r['helicity']==1 and r['current_feedback'] and r['faraday_feedback']],key=lambda r:-r['dt_omega_pe'])
order('Maxwell time order',timeem[0]['temporal_relative_error'],timeem[-1]['temporal_relative_error'],4,'EM_min_time_order')
reverse=next(r for r in em if r['helicity']==-1)
check('B reversal exchanges circular polarization',abs(reverse['omega_measured_over_pe']-timeem[0]['omega_measured_over_pe'])<1e-13)
hall=next(r for r in rows if r['family']=='hall_limit')['rows']
bound('small-k Hall-MHD limit',hall[0]['relative_difference'],'Hall_limit_small_k_relative')
check('Hall-MHD discrepancy decreases with electron-scale separation',all(a['relative_difference']<b['relative_difference'] for a,b in zip(hall,hall[1:])))
execution=next(r for r in rows if r['family']=='execution');bound('bounded benchmark wall time',execution['elapsed_s'],'max_benchmark_wall_s')
commands=json.loads((run/'commands.json').read_text())
check('frozen CPU6 nice15 one-thread commands',all(c['returncode']==0 and c['cpu_affinity']==[6] and c['nice']==15 and c['threads']==1 and c['shared_lease_exclusion_verified'] for c in commands))
check('native foundation regression',any(c['label']=='foundation' and c['returncode']==0 for c in commands))
failures=[c for c in checks if not c['passed']]
result={'status':'passed' if not failures else 'failed','checks':checks,'failures':failures,'count':len(checks),'limits_changed_after_runs':False}
(run/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
print(result['status'],len(checks),'checks',len(failures),'failures')
for c in failures:print(c)
if failures:raise SystemExit(1)
