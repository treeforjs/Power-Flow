import collections, csv, hashlib, json, os, pathlib, shutil, time

C = pathlib.Path('/home/treefor/projects/runtimes/pic-thermal-field-emission-20261029')
O = pathlib.Path('<workspace>/outputs/pic-thermal-field-emission-20261029')
S = C/'snapshots/usable-api-v1'
Q = C/'qualification-v2'
E = O/'qualification-v2'
assert os.sched_getaffinity(0)=={7} and os.getpriority(os.PRIO_PROCESS,0)>=15
digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
api_manifest=json.loads((S/'SOURCE-SHA256.json').read_text())
assert digest(S/'SOURCE-SHA256.json')=='9bd0fa5771b8cba405cc0b63b99c0e32a24da8d735eef240606c36b1fb49f0ca'
for name,expected in api_manifest.items(): assert digest(S/name)==expected,name
material_origin=C.parent/'ss304-thermal-properties-20261029/source'
for local,origin in [('ss304_thermal.hpp','include/ss304_thermal.hpp'),('scenario.json','data/scenario.json'),('README.md','README.md')]:
    assert (Q/'dependencies/ss304'/local).read_bytes()==(material_origin/origin).read_bytes()
read_csv=lambda name:list(csv.DictReader((E/name).open()))
base=json.loads((E/'integration-qualification.json').read_text())
assert base['passed'] and base['exact_field_epoch_restart'] and base['prepare_and_failure_rollback']
assert '"passed":true' in (E/'explicit-transition-validation.log').read_text()
cases=read_csv('matching-substep-sensitivity.csv')
returns=read_csv('magnetic-return-sensitivity.csv')
assert len(cases)==9 and len(returns)==9
def relative_range(rows,name):
    values=[float(row[name]) for row in rows]
    return (max(values)-min(values))/max(abs(v) for v in values)
for rows in [cases,returns]:
    assert all(float(x['common_dt_s'])==1e-14 for x in rows)
    assert relative_range(rows,'mean_net_outward_current_A_m2')<1e-12
assert all(float(x['returned_lower'])>0 for x in returns)
assert all(int(x['image_reentries'])>0 for x in returns if float(x['eta_energy'])>=.03)
assert all(int(x['image_reentries'])==0 and int(x['outgoing_handoffs'])==0 for x in returns if float(x['eta_energy'])==.01)
failures=[x for x in read_csv('nominal-barrier-admission.csv') if x['admitted']=='0']
assert len(failures)==46
assert all('magnetic energy criterion' in x['reason'] for x in failures)
assert (E/'build-and-run-v3.log').exists() and 'did not resolve an extracting-side return' in (E/'build-and-run-v3.log').read_text()
report={
 'schema':'emission-epoch-reentry-qualification-v2', 'status':'SEALED_COMPONENT_QUALIFICATION',
 'created_unix_s':time.time(), 'immutable_API_v1_manifest_sha256':digest(S/'SOURCE-SHA256.json'),
 'immutable_API_v1_all_manifest_files_verified':True, 'API_v1_source_checks_preserved':80, 'API_v1_micro_cases_preserved':4,
 'common_clock_pilot':base,
 'resident_phase_space_vs_ledger_checked_each_epoch':True,
 'complete_wall_heat_Joule_caloric_carry_checked_each_epoch':True,
 'reverse_adapter':'inward_layer_adapter.hpp', 'reverse_adapter_source_sha256':digest(Q/'inward_layer_adapter.hpp'),
 'reverse_duplicate_restore_rejected':True,
 'matching_substep_variants':len(cases), 'magnetic_return_variants':len(returns),
 'magnetic_return_stages':200, 'magnetic_return_initial_F_V_m':1e7, 'magnetic_return_initial_B_tangent_T':25,
 'magnetic_return_initial_B_normal_T':0, 'nominal_quantum_magnetic_fraction':0.2,
 'magnetic_return_emitted_lower':float(returns[0]['emitted_lower']),
 'magnetic_return_returned_lower':float(returns[0]['returned_lower']),
 'magnetic_return_image_reentries_each_variant':[int(x['image_reentries']) for x in returns],
 'matching_field_work_relative_span':relative_range(cases,'field_work_J'),
 'magnetic_return_field_work_relative_span':relative_range(returns,'field_work_J'),
 'magnetic_return_wall_heat_relative_span':relative_range(returns,'wall_heat_J'),
 'maximum_absolute_ledger_residual_J':max(abs(float(x['energy_residual_J'])) for x in cases+returns),
 'nominal_scan_failures_by_face_T':dict(collections.Counter(x['face']+'@'+x['T_K'] for x in failures)),
 'nominal_scan_failed_F_range_V_m':[min(float(x['signed_F_V_m']) for x in failures),max(float(x['signed_F_V_m']) for x in failures)],
 'canonical_barrier_diagnostic':json.loads((E/'canonical-barrier-diagnostic.json').read_text()),
 'native_midpoint_Maxwell_heat_conduction_qualification':False,
 'native_particle_self_image_overlap_qualified':False,
 'physical_foil_model_admitted':False,'GPU_executed':False,'production_resumed':False,
 'runtime':'arch-base-20261028/run-frozen','cpu':7,'nice':15,'exclusive_lock':str(O.parent/'pic-profile/testing-core7.lock'),
 'new_followup_package':'/home/treefor/projects/runtimes/pic-magnetized-emission-20261029',
 'remaining_admission':'Retain nominal guard failures; qualify native reverse matching and a real common physical midpoint. Develop rate-weighted rejected-phase map and justified finite-T magnetized joint supply/transmission separately.'
}
(E/'QUALIFICATION.json').write_text(json.dumps(report,indent=2)+'\n')
port={
 'schema':'boundary-layer-inward-matching-port-v2','status':'sealed_component_helper',
 'immutable_API_dependency':str(S),'API_v1_manifest_sha256':digest(S/'SOURCE-SHA256.json'),
 'source':str(Q/'inward_layer_adapter.hpp'),'source_sha256':digest(Q/'inward_layer_adapter.hpp'),
 'function':'surface_emission_integration::reenter_image_layer(LayerState,Stage,LayerPolicy)->LayerTrial',
 'caller_event':'Resolve first inward crossing of originating LayerState.handoff_distance_m in ordinary PIC; retain LayerState provenance after outgoing handoff.',
 'transaction':'K+=|Uimage|; image_energy+=Uimage<0; preserve tangential proper momentum; normal proper momentum inward; status Active; physical position/time unchanged.',
 'zero_exchange':['number','particle_charge','electrode_charge','wall_heat','electron_reference','field_work','displacement_current'],
 'next':'Advance Active state over remaining same accepted interval; do not also advance it in ordinary PIC.',
 'guards':['PICHandedOff only; duplicate re-entry rejected','same field/thermal accepted-state, field time, mesh/material/patch and normal','actual first inward handoff-plane crossing in admitted local geometry','no separately already applied analytic microscopic image force'],
 'ledger_ownership':'LayerTrial.exchange incremental only; Packet.ledger remains birth diagnostic; commit once with inverse status transition.',
 'native_integration_admitted':False,'physical_model_admitted':False,'nominal_magnetic_guard':0.2,
 'evidence':str(E/'QUALIFICATION.json')
}
(O/'interfaces/boundary-layer-inward-matching-port-v2.json').write_text(json.dumps(port,indent=2)+'\n')
package=O/'source-backup-epoch-reentry-v2'
assert not package.exists(),'Refuse overwriting sealed package'
shutil.copytree(Q,package,ignore=shutil.ignore_patterns('build','__pycache__'))
shutil.copytree(E,package/'evidence')
shutil.copy(O/'interfaces/boundary-layer-inward-matching-port-v2.json',package/'boundary-layer-inward-matching-port-v2.json')
manifest={str(p.relative_to(package)):digest(p) for p in sorted(package.rglob('*')) if p.is_file()}
(package/'SOURCE-SHA256.json').write_text(json.dumps(manifest,indent=2)+'\n')
seal={'status':'SEALED_EPOCH_REENTRY_V2','source_manifest_sha256':digest(package/'SOURCE-SHA256.json'),
      'source_snapshot':str(package),'source_API_unchanged':True,'physical_model_admitted':False,'production_resume_authorized':False}
(package/'SEALED.json').write_text(json.dumps(seal,indent=2)+'\n')
shutil.make_archive(str(package),'gztar',package.parent,package.name)
(O/'EPOCH-REENTRY-V2.json').write_text(json.dumps(seal,indent=2)+'\n')
publish_files=[]
for p in sorted(package.rglob('*')):
    if p.is_file() and p.suffix!='.bin':
        publish_files.append({'path':'artifacts/surface-emission-20261029/qualification-v2/'+str(p.relative_to(package)),
                              'mode':'100644','type':'blob','content':p.read_text()})
payload={'repository_full_name':'treeforjs/Power-Flow','base_tree_sha':'b864e19df6c6f68e80b9e71a7dae81e49454db7f','tree_elements':publish_files}
(O/'publish-epoch-reentry-v2.json').write_text(json.dumps(payload))
print(json.dumps({'sealed':seal,'publish_files':len(publish_files),'maximum_ledger_residual_J':report['maximum_absolute_ledger_residual_J'],
                  'magnetic_reentries_each_variant':report['magnetic_return_image_reentries_each_variant']}))
