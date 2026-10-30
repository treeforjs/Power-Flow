import csv,hashlib,json,os,pathlib,shutil,time
C=pathlib.Path('/home/treefor/projects/runtimes/pic-magnetized-emission-20261029')
O=pathlib.Path('<workspace>/outputs/pic-magnetized-emission-20261029')
E=O/'tangential-onepass'
assert os.sched_getaffinity(0)=={7} and os.getpriority(os.PRIO_PROCESS,0)>=15
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
api=pathlib.Path('/home/treefor/projects/runtimes/pic-thermal-field-emission-20261029/snapshots/usable-api-v1')
assert sha(api/'SOURCE-SHA256.json')=='9bd0fa5771b8cba405cc0b63b99c0e32a24da8d735eef240606c36b1fb49f0ca'
for path,value in json.loads((api/'SOURCE-SHA256.json').read_text()).items():assert sha(api/path)==value
reverse=pathlib.Path('<workspace>/outputs/pic-thermal-field-emission-20261029/source-backup-epoch-reentry-v2/inward_layer_adapter.hpp')
assert sha(reverse)=='f47c5978e1ada5a8f74342df41b7b591f3839f06157c0796f912f07adfa5e294'
cases=[json.loads(p.read_text()) for p in E.glob('joint-*.json')]
assert len(cases)==16 and all(x['passed'] and not x['physical_foil_model_admitted'] and not x['exact_magnetized_solid_supply_solved'] for x in cases)
groups={}
for row in cases:groups.setdefault((row['T_K'],row['B_tangent_T']),{})[row['refinement']]=row
refinement=[]
for key,value in sorted(groups.items()):
 a,b=value[1],value[2]
 refinement.append({'T_K':key[0],'B_tangent_T':key[1],
  'relative_current_refinement':abs(a['known_component_current_A_m2']/b['known_component_current_A_m2']-1),
  'relative_heat_refinement':abs(a['known_component_complete_heat_W_m2']/b['known_component_complete_heat_W_m2']-1),
  'fine_known_current_A_m2':b['known_component_current_A_m2'],
  'fine_nominal_unmagnetized_current_A_m2':b['nominal_unmagnetized_current_A_m2'],
  'unresolved_current_upper_A_m2':b['unresolved_current_upper_A_m2'],
  'unresolved_current_relative_to_known':b['unresolved_current_upper_A_m2']/b['known_component_current_A_m2'],
  'magnetic_reservoir_energy_scale_over_kT':b['magnetic_reservoir_energy_scale_over_kT']})
assert max(x['relative_current_refinement'] for x in refinement)<2e-5
assert max(x['relative_heat_refinement'] for x in refinement)<2e-5
cycles=[]
for p in E.glob('inverse-cycle-*.csv'):
 rows=list(csv.DictReader(p.open()));assert len(rows)==12
 for row in rows:row['case']=p.stem
 cycles+=rows
assert len(cycles)==96
complete=[r for r in cycles if r['completed_material_return']=='1']
retained=[r for r in cycles if r['completed_material_return']=='0']
assert len(complete)>=72
assert all(int(r['quantum_absorptions'])==1 for r in complete)
assert all(int(r['quantum_absorptions'])==0 for r in retained)
assert sum(int(r['quantum_reflections']) for r in complete)>0
admission=json.loads((E/'admission-and-matching.json').read_text());assert admission['passed'] and admission['unresolved_case_sampler_refused']
assert (E/'final-joint-cases.log').read_text().count('Passed joint prototype')==16
report={'schema':'tangential-onepass-joint-component-qualification-v1','status':'SEALED_CONDITIONAL_COMPONENT',
 'created_unix_s':time.time(),'source_API_v1_manifest_sha256':sha(api/'SOURCE-SHA256.json'),
 'source_API_v1_all_files_verified_unchanged':True,'reverse_adapter_sha256':sha(reverse),
 'nominal_API_v1_guard':.2,'original_guard_widened':False,
 'conditional_model':['continuum free-electron FD reservoir','stationary locally planar canonical tangential-B one-pass Kemble/JWKB barrier','incoherent reciprocal quantum return'],
 'numerical_test_temperatures_K':[298.2,300,773.15,1490],'numerical_test_F_V_m':97724679.40458213,'quadrature_cases':16,
 'refinement':refinement,'fixed_clock_trajectory_cases':len(cycles),'completed_material_return':len(complete),'retained_in_one_inventory_at_4ps':len(retained),
 'quantum_reflections_total':sum(int(r['quantum_reflections']) for r in cycles),
 'maximum_phase_reconstruction_error_eV':max(float(r['phase_reconstruction_error_eV']) for r in cycles),
 'maximum_ledger_residual_over_transported_energy':max(abs(float(r['complete_residual_J']))/float(r['transported_energy_scale_J']) for r in cycles),
 'maximum_abs_material_heat_over_transported_energy_for_closed_cycles':max(abs(float(r['net_wall_heat_J']))/float(r['transported_energy_scale_J']) for r in complete),
 'matching_and_unsupported_contribution_checks':admission,
 'quadrature_error_is_separate_from_unresolved_tail_bounds':True,
 'actual_magnetized_solid_band_error_bound':None,'oblique_B_error_bound':None,'coherent_quantum_matching_error_bound':None,
 'strict_tangential_replacement_admits_legacy_probe_states':0,
 'actual_foil_current':None,'physical_foil_model_admitted':False,'native_midpoint_Maxwell_heat_qualification':False,
 'native_return_counter_rollback_restart_qualified':False,'GPU_executed':False,'production_resumed':False,
 'cpu':7,'nice':15,'one_job':True,'frozen_runtime':'arch-base-20261028/run-frozen',
 'root_admission_requirements':['actual one-sided E/B and evolved supported T at converged common physical midpoint','native total-field/current/work owner and one inventory','qualified originating source-model tag plus quantum-return counter in restart/journal','quantified actual conductor-image self-response overlap','separate material/oblique/coherence applicability; no projection or unknown-to-zero fallback']}
(E/'QUALIFICATION.json').write_text(json.dumps(report,indent=2)+'\n')
port={'schema':'tangential-joint-known-component-source-port-v1','status':'sealed_conditional_component_not_native_admission',
 'header':str(C/'tangential-onepass/tangential_joint.hpp'),'header_sha256':sha(C/'tangential-onepass/tangential_joint.hpp'),
 'functions':['tangential_onepass::evaluate(Stage,Config,LayerPolicy)','sample_known_component(Distribution,LayerPolicy,weight,joint_uniform,event_id)','resolve_quantum_return(Stage,Config,LayerPolicy,provisional_collected_trial,return_uniform)'],
 'input_source_model_flags':['continuum_free_electron_reservoir_selected','incoherent_onepass_matching_selected'],
 'input_domain':{'T_K':[298.2,1500],'phi_eV':[2,6],'F_positive_V_m':[5e7,1e8],'B_tangent_T':[0,40],'B_normal':'strictly zero within floating equality tolerance','E_tangent':'strictly zero within floating equality tolerance'},
 'rate_owner':'known_component_current only; unresolved upper bound retained separately; no old nominal current renormalized on accepted phases',
 'heat_owner':'same joint weights; sampled COMPLETE wall heat once; mean diagnostic never also debited',
 'transport_owner':'same LayerState and one native arena inventory, Active or ordinary PIC',
 'return_owner':'v1 recollection is PROVISIONAL; dispatch reciprocal resolve_quantum_return BEFORE publish; reflection removes that material bridge and keeps Active state',
 'caller_restart_fields':['originating source model/seal identity','same LayerState','distinct quantum return event counter and caller-owned draw identity'],
 'current_path_owner':'real micro/PIC paths and explicit quantum bridges integrated once over accepted dt; reflected quantum bridge is removed atomically',
 'native_arena_observed':'runtimes/pic-reaction-transactions-20261029/native-source-port SurfaceSourceBatch/source_layers',
 'new_native_dispatch_qualified':False,'physical_foil_model_admitted':False,'qualification':str(E/'QUALIFICATION.json')}
(O/'interfaces/tangential-joint-source-port-v1.json').write_text(json.dumps(port,indent=2)+'\n')
P=O/'source-backup-tangential-onepass-v1';assert not P.exists()
shutil.copytree(C,P,ignore=shutil.ignore_patterns('build','snapshots','__pycache__'))
shutil.copytree(E,P/'evidence')
shutil.copy(O/'interfaces/tangential-joint-source-port-v1.json',P/'tangential-joint-source-port-v1.json')
shutil.copy(reverse,P/'inward_layer_adapter.hpp')
(P/'DEPENDENCIES.json').write_text(json.dumps({'immutable_API_v1_manifest_sha256':sha(api/'SOURCE-SHA256.json'),'complete_API_v1_outside_backup':'<workspace>/outputs/pic-thermal-field-emission-20261029/source-backup-usable-api-v1','reverse_adapter_sha256':sha(reverse)},indent=2)+'\n')
manifest={str(p.relative_to(P)):sha(p) for p in sorted(P.rglob('*')) if p.is_file()}
(P/'SOURCE-SHA256.json').write_text(json.dumps(manifest,indent=2)+'\n')
seal={'status':'SEALED_TANGENTIAL_ONEPASS_COMPONENT_V1','manifest_sha256':sha(P/'SOURCE-SHA256.json'),'source_snapshot':str(P),'physical_foil_model_admitted':False,'oblique_B_supported':False,'new_native_dispatch_qualified':False}
(P/'SEALED.json').write_text(json.dumps(seal,indent=2)+'\n');shutil.make_archive(str(P),'gztar',P.parent,P.name)
items=[{'path':'artifacts/surface-emission-20261029/tangential-onepass-v1/'+str(p.relative_to(P)),'mode':'100644','type':'blob','content':p.read_text()} for p in sorted(P.rglob('*')) if p.is_file()]
batches=[];batch=[];size=0
for item in items:
 n=len(json.dumps(item).encode())
 if batch and size+n>180000:batches.append(batch);batch=[];size=0
 batch.append(item);size+=n
if batch:batches.append(batch)
for i,b in enumerate(batches):(O/('publish-joint-batch-'+str(i)+'.json')).write_text(json.dumps(b))
(O/'TANGENTIAL-ONEPASS-V1.json').write_text(json.dumps({**seal,'publication_batches':len(batches),'publication_files':len(items)},indent=2)+'\n')
print(json.dumps({'seal':seal,'publication_batches':len(batches),'publication_files':len(items),'completed':len(complete),'retained':len(retained),'max_current_refinement':max(x['relative_current_refinement'] for x in refinement),'max_phase_error_eV':report['maximum_phase_reconstruction_error_eV']}))
