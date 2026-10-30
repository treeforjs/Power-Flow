from pathlib import Path
import json,hashlib,os
p=Path(__file__).parent;assert sorted(os.sched_getaffinity(0))==[7]
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
inputs=json.loads((p/'input-bindings.json').read_text())
for x in inputs:assert sha(Path(x['snapshot']))==x['sha256']
for name in ['review','review-v2']:
 old=p.parent/name
 for file,rec in json.loads((old/'output-manifest.json').read_text()).items():assert sha(old/file)==rec['sha256'],str(old/file)
kin=json.loads((p/'high-precision-oracle.json').read_text());oracle=kin['summary'];rem=json.loads((p/'remainder-oracle.json').read_text())
assert oracle['all_admitted_cases']==54 and oracle['all_admitted_max_energy_ULP_difference']==oracle['all_admitted_max_flux_ULP_difference']==0
assert len(oracle['expected_conversion_underflow_rejections'])==3 and len(oracle['ordinary_overflow_rejections'])==3 and not oracle['unexpected_rejections']
assert not oracle['positive_energy_underflow_accepted'] and not oracle['positive_flux_underflow_accepted']
assert rem['status']=='passed' and rem['cases']==110
assert 'checks=149' in (p/'adversarial.txt').read_text() and '14 prepare-size and stage checks passed' in (p/'guard-checks.txt').read_text()
status={'status':'bounded_v3_component_review_passed_with_prepare_guard_patch','reviewed_live_header_sha256':inputs[0]['sha256'],'patched_review_header_sha256':sha(p/'patched-collision_scheduler.hpp'),'proposed_patch_sha256':sha(p/'prepare-size-guard.patch'),'root_source_edited':False,'owner_live_patch_application_verified':False,'independent_adversarial_assertions':149,'prepare_guard_and_stage_assertions':14,'high_precision_kinematic_cases':60,'admitted_cases_correctly_rounded_energy_and_flux':54,'expected_energy_overflow_rejections':3,'expected_conversion_underflow_rejections':3,'unexpected_rejections':0,'silent_underflow_cases':0,'exact_two_limb_difference_cases':110,'nonzero_low_conversion_rejections':69,'v1_v2_reviews_preserved':True,'remaining_integration_gaps':['Canonical compensated native/neutral storage and iterative compensated-weight input are not implemented','One-double conversion of a nonzero low limb must fail before common physical/RNG publication','Actual stage/state/revision/native deposition/product/energy integration remains caller-owned','Manufactured clock tests do not qualify physical H2O cross sections or product kernels'],'CPU':7,'nice':15,'threads':1,'frozen_userspace':True,'GPU':False,'native_mutation':False,'physical_H2O_acceptance':False}
(p/'validation.json').write_text(json.dumps(status,indent=2)+'\n')
(p/'handoff.json').write_text(json.dumps({'status':status['status'],'review':str(p/'review.txt'),'validation':str(p/'validation.json'),'patch':str(p/'prepare-size-guard.patch'),'patched_header':str(p/'patched-collision_scheduler.hpp'),'live_base_header_sha256':inputs[0]['sha256'],'patched_header_sha256':status['patched_review_header_sha256'],'owner_source_edited':False,'native_compensated_reservoir_accepted':False},indent=2)+'\n')
(p/'output-manifest.json').write_text(json.dumps({f.name:{'sha256':sha(f),'bytes':f.stat().st_size} for f in p.iterdir() if f.is_file() and f.name!='output-manifest.json'},indent=2)+'\n')
print(json.dumps(status))
