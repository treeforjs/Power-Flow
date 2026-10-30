from pathlib import Path
import hashlib,json,os,shutil
p=Path(__file__).parent
assert sorted(os.sched_getaffinity(0))==[7]
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
records=json.loads((p/'input-bindings.json').read_text())
for r in records:assert sha(Path(r['snapshot']))==r['sha256'],r['snapshot']
d=json.loads((p/'independent-oracle.json').read_text())['summary']
assert (d['cases'],d['accepted'],d['rejected'])==(108,46,62)
assert not d['oracle_energy_outside_reported_bound'] and not d['oracle_momentum_outside_reported_bound']
mut={r['mutation']:r for r in d['mutation_cases']}
assert not mut['near_unit_axis_ordinary_weight']['accepted'] and mut['near_unit_axis_tiny_weight']['accepted']
assert mut['changed_target_mass_same_rate']['accepted']
status={'status':'two_actionable_numerical_findings_returned_not_accepted','reviewed_bridge_sha256':sha(p/'source/water_parent_bridge.hpp'),'reviewed_scheduler_sha256':sha(p/'inputs/collision_scheduler.hpp'),'reviewed_data_sha256':sha(p/'inputs/h2o-gas.json'),'independent_kinematic_cases':108,'admitted':46,'rejected':62,'mutation_and_scale_probes':len(mut),'oracle_decimal_digits':120,'original_reported_lab_bounds_passed_for_all_admitted':True,'reaction_scale_accuracy_not_qualified_by_original_bounds':True,'findings':[{'id':'P2-direction-conditioning','file':'source/water_parent_bridge.hpp','line':43,'counterexample_case':20},{'id':'P2-weighted-conservation-underflow','file':'source/water_parent_bridge.hpp','line':105,'ordinary_rejection_tiny_weight_acceptance':True}],'native_or_physical_H2O_acceptance':False,'GPU':False,'native_mutation':False,'root_source_edited':False,'CPU':7,'nice':15,'threads':1,'frozen_userspace':True,'input_snapshots_verified':len(records),'scheduler_sparse_work_budget_recommendation_included':True}
(p/'validation.json').write_text(json.dumps(status,indent=2)+'\n')
(p/'handoff.json').write_text(json.dumps({'status':status['status'],'review':str(p/'review.txt'),'validation':str(p/'validation.json'),'oracle':str(p/'independent-oracle.json'),'probes':str(p/'probes.jsonl'),'reviewed_bridge_sha256':status['reviewed_bridge_sha256'],'native_compensated_reservoir_accepted':False},indent=2)+'\n')
(p/'output-manifest.json').write_text(json.dumps({str(f.relative_to(p)):{'sha256':sha(f),'bytes':f.stat().st_size} for f in sorted(p.rglob('*')) if f.is_file() and f.name!='output-manifest.json'},indent=2)+'\n')
print(json.dumps(status))
# New review target is a separate immutable snapshot, preserving this v1 review.
r=Path('/home/treefor/projects/runtimes/pic-water-event-bridge-20261029');v=p.parent/'review-v2'
v.mkdir(exist_ok=False)
for name in ('source','inputs'):shutil.copytree(r/name,v/name)
for name in ('input-sha256.json','run_checks.sh'):shutil.copy2(r/name,v/name)
shutil.copy2(p.parent/'checks.json',v/'owner-checks.json')
bindings=[]
for f in sorted(v.rglob('*')):
 if f.is_file():
  original=r/f.relative_to(v)
  if original.exists():assert sha(f)==sha(original)
  bindings.append({'snapshot':str(f),'original':str(original),'sha256':sha(f),'bytes':f.stat().st_size})
(v/'input-bindings.json').write_text(json.dumps(bindings,indent=2)+'\n')
print('New v2 snapshot bridge',sha(v/'source/water_parent_bridge.hpp'))
