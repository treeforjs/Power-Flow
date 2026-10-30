from pathlib import Path
import hashlib,json,os
p=Path(__file__).parent;assert sorted(os.sched_getaffinity(0))==[7]
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
original=json.loads((p/'input-bindings.json').read_text())
for x in original:assert sha(Path(x['snapshot']))==x['sha256']
v1=p.parent/'review';old=json.loads((v1/'output-manifest.json').read_text())
assert all(sha(v1/name)==rec['sha256'] for name,rec in old.items())
status={'status':'five_structural_fixes_verified_with_remaining_requests','reviewed_header_sha256':original[0]['sha256'],'independent_assertions':149,'high_precision_kinematic_cases':60,'ordinary_cases_correctly_rounded_energy_and_flux':48,'expected_energy_overflow_rejections':3,'unexpected_kinematic_rejections':0,'exact_remainder_cases':110,'nonzero_low_conversion_rejections':69,'original_v1_review_preserved':True,'remaining_requests':[{'priority':2,'issue':'Nonzero energy/flux conversion underflow must reject or remain scaled through propensity','evidence':'underflow.txt; high-precision-oracle.json boost8keV_7; representable rate4.462882568429323e-213/s disappears'},{'priority':2,'issue':'Declare/qualify nearly comoving radial/angular precision domain','evidence':'one-ULP parallel/skew oracle cases; relative errors up to8.592e-4 at tiny energy'},{'priority':2,'issue':'Guard direct prepare candidate copy/sort size or enforce bounded provenance','evidence':'Reviewed prepare has no input-size/work parameter or precheck'}],'integration_limits':['Packet still takes one-double counts; no iterative compensated-neutral storage','Native owner must preserve both limbs or fail before physical/RNG commit','Physical frozen-stage state/revision and elapsed progress are caller-owned','No product kinematics, native mutation or physical H2O validation'],'CPU':7,'nice':15,'threads':1,'frozen_userspace':True,'source_edits':False,'GPU':False}
(p/'validation.json').write_text(json.dumps(status,indent=2)+'\n')
(p/'handoff.json').write_text(json.dumps({'status':status['status'],'review':str(p/'review.txt'),'validation':str(p/'validation.json'),'oracle':str(p/'high-precision-oracle.json'),'remainders':str(p/'remainder-oracle.json'),'reviewed_header_sha256':original[0]['sha256'],'original_v1_preserved':True,'owner_source_modified':False},indent=2)+'\n')
(p/'output-manifest.json').write_text(json.dumps({f.name:{'sha256':sha(f),'bytes':f.stat().st_size} for f in p.iterdir() if f.is_file() and f.name!='output-manifest.json'},indent=2)+'\n')
print(json.dumps(status))
