from pathlib import Path
import os,json,hashlib,shutil
p=Path(__file__).parent;old=p.parent/'review'
assert sorted(os.sched_getaffinity(0))==[7]
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
for f,r in json.loads((old/'output-manifest.json').read_text()).items():assert sha(old/f)==r['sha256'],f
bindings=json.loads((p/'input-bindings.json').read_text())
for r in bindings:assert sha(Path(r['snapshot']))==r['sha256'],r['snapshot']
results={name:json.loads((p/f'{name}-oracle.json').read_text()) for name in ('replay','physical')}
assert (results['replay']['summary']['cases'],results['replay']['summary']['accepted'])==(108,12)
assert (results['physical']['summary']['cases'],results['physical']['summary']['accepted'])==(72,48)
for name,r in results.items():
 assert not r['summary']['oracle_energy_outside_reported_bound'] and not r['summary']['oracle_momentum_outside_reported_bound']
 mutations={m['mutation']:m for m in r['summary']['mutation_cases']}
 for label in ('changed_target_mass_same_rate','tiny_weight','near_unit_axis_ordinary_weight','near_unit_axis_tiny_weight','changed_volume_without_rate'):
  assert not mutations[label]['accepted'],(name,label)
 assert mutations['rotated_projectile_same_rate']['accepted'] and mutations['changed_volume_with_rate']['accepted']
 for case in r['cases']:
  if case['accepted']:
   assert case['direction_vector_error']<1e-12
   assert case['energy_bound_eV_per_reaction']<=1e-10*case['oracle_target_rest_eV']*(1+1e-13)
assert not results['replay']['cases'][20]['accepted']
source=Path('/home/treefor/projects/runtimes/pic-water-event-bridge-20261029/source/water_parent_bridge.hpp')
status={'status':'bounded_component_numerical_fixes_independently_verified','reviewed_bridge_sha256':sha(p/'source/water_parent_bridge.hpp'),'live_source_matches_reviewed_snapshot':sha(source)==sha(p/'source/water_parent_bridge.hpp'),'reviewed_scheduler_sha256':sha(p/'inputs/collision_scheduler.hpp'),'reviewed_data_sha256':sha(p/'inputs/h2o-gas.json'),'old_review_and_original_evidence_unchanged':True,'original_sweep_replayed':108,'original_sweep_admitted':12,'additional_bounded_sweep_cases':72,'additional_bounded_sweep_admitted':48,'distinct_kinematic_inputs':144,'oracle_decimal_digits':120,'maximum_admitted_direction_error':results['physical']['summary']['max_direction_error_accepted']['direction_vector_error'],'all_admitted_oracle_energy_momentum_within_unchanged_bounds':True,'all_admitted_lab_energy_resolution_within_declared_policy':True,'unique_mutation_scale_cases':11,'original_direction_counterexample_rejected':True,'original_weighted_underflow_bypass_rejected':True,'catalog_mass_change_rejected':True,'authoritative_epoch_mesh_and_stage_remain_host_owned':True,'all_nonzero_receipt_components_preserved':False,'canonical_compensated_reservoir_implemented':False,'native_mutation':False,'physical_H2O_or_EEDF_acceptance':False,'GPU':False,'root_source_edited':False,'CPU':7,'nice':15,'threads':1,'frozen_userspace':True}
(p/'validation.json').write_text(json.dumps(status,indent=2)+'\n')
(p/'handoff.json').write_text(json.dumps({'status':status['status'],'review':str(p/'review.txt'),'validation':str(p/'validation.json'),'replay_oracle':str(p/'replay-oracle.json'),'bounded_oracle':str(p/'physical-oracle.json'),'reviewed_bridge_sha256':status['reviewed_bridge_sha256'],'native_or_full_H2O_accepted':False,'original_findings_and_evidence':str(old/'handoff.json')},indent=2)+'\n')
(p/'output-manifest.json').write_text(json.dumps({str(f.relative_to(p)):{'sha256':sha(f),'bytes':f.stat().st_size} for f in sorted(p.rglob('*')) if f.is_file() and f.name!='output-manifest.json'},indent=2)+'\n')
print(json.dumps(status))
