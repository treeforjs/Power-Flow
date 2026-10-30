from pathlib import Path
from fractions import Fraction as F
import os,json,hashlib,math,struct,shutil
p=Path(__file__).parent;s=p/'snapshot';o=Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-native-subflow-v6_1-review-20261029');assert sorted(os.sched_getaffinity(0))==[7]
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
assert sha(p/'driver-source.tar.gz')=='f7e4ceddc69ad4eed43721b24415691c4d668b340b3b06ec3b0380a8b8da25fc'
assert sha(s/'source/native_parent_subflow.hpp')=='ad499c5549a78bae371c5b2f72ee66085c3b141126d1931d1ec654d8a9e25675'
for rec in json.loads((o/'input-bindings.json').read_text()):assert sha(Path(rec['snapshot']))==rec['sha256']
for name in ('pic-native-arena-review-20261029','pic-native-arena-review-v5-20261029','pic-water-native-subflow-review-20261029'):
 old=o.parent/name
 for n,rec in json.loads((old/'output-manifest.json').read_text()).items():assert sha(old/n)==rec['sha256'],str(old/n)
d=json.loads((o/'probes.json').read_text());assert d['assertions']==39787 and len(d['records'])==8
records={r['name']:r for r in d['records']};r=records['compensated_carrier_first_moment'];e=sum(map(F,r['electron_count']),F());n=sum(map(F,r['target_count']),F())
assert sum(map(F,r['electron_carrier_weights']),F())==e and sum(map(F,r['target_carrier_weights']),F())==n
rows=r['rows'];assert len(rows)==6
coefficient=F(rows[0]['sigma_m2'])*F(rows[0]['flux_m_s'])/F(rows[0]['volume_m3']);expected=e*n*coefficient;actual=sum((F(x['event_weight'])*F(x['rate_s']) for x in rows),F());defect=float((actual-expected)/expected);assert abs(defect)<=8*math.ulp(1.)
old_records=json.loads((o.parent/'pic-water-native-subflow-review-20261029/probes.json').read_text());old=next(r for r in old_records if r['name']=='compensated_carrier_first_moment')
for a,b in zip(old['rows'],rows):
 for key in a:assert struct.pack('!d',a[key])==struct.pack('!d',b[key]),key
oracle={'method':'Fraction from exact binary64 leaves/rates plus numeric bit comparison to archived428 control','electron_carriers':3,'target_carriers':2,'rows':6,'carrier_sums_exact':True,'candidate_numeric_fields_bitwise_equal_archived_control':True,'weighted_event_rate_expected_per_s':float(expected),'weighted_event_rate_sum_per_s':float(actual),'weighted_event_rate_relative_rounding_defect':defect}
(o/'weighted-rate-oracle.json').write_text(json.dumps(oracle,indent=2)+'\n')
ports=json.loads((o/'port-bounds-probes.json').read_text());assert len(ports)==3 and all(x['prepared'] and x['published'] and x['save_restore_admitted'] for x in ports)
libs=json.loads((o.parent/'pic-native-arena-review-v5-20261029/native-library-hashes.json').read_text());libs={k:v for k,v in libs.items() if k.endswith('.a')}
for k,v in libs.items():assert sha(Path(k))==v
(o/'frozen-library-bindings.json').write_text(json.dumps(libs,indent=2)+'\n')
status={'status':'focused_driver_fix_family_passed_separate_standalone_port_findings','driver_sha256':sha(s/'source/native_parent_subflow.hpp'),'archive_sha256':sha(p/'driver-source.tar.gz'),'native_source_identity':'fd318a6a9ae081bdaf57656c7e47da3e528c9ea2b532224de15bab37b56f3dc6','original_normal_checks':442,'independent_assertions':39787,'focused_records':8,'zeroevent_mixed_age_rejected':True,'retagged_old_rng_and_foreign_product_rng_rejected':True,'fresh_restored_driver_and_stale_handles_passed':True,'one_real_event_then_six_events_restart_exact':True,'zeroevent_and_actual_reaction_completion_rollback_restart_passed':True,'vacuum_cells':9875,'solid_cells':250,'vacuum_volume_m3':records['aligned_mask_physical_volumes']['physical_vacuum_volume_m3'],'cutcell_surface_solid_refused':True,'six_leaf_numeric_rates_bitwise_equal_archived_control':True,'weighted_rate_relative_rounding_defect':defect,'standalone_api_SSA_cross_record_clock_variants_admitted':3,'standalone_variants_use_driver_Scheduler':False,'earlier_review_manifests_unchanged':True,'live_source_changed':False,'GPU':False,'physical_time_or_full_water_qualified':False,'CPU':7,'nice':15,'threads':1,'frozen_userspace':True}
(o/'validation.json').write_text(json.dumps(status,indent=2)+'\n');(o/'handoff-review.json').write_text(json.dumps({'status':status['status'],'review':str(o/'review.txt'),'validation':str(o/'validation.json'),'driver_controls':str(o/'probes.json'),'rate_oracle':str(o/'weighted-rate-oracle.json'),'standalone_port_probes':str(o/'port-bounds-probes.json'),'standalone_repro':str(p/'port-elapsed.cpp'),'scope':'ad499/fd318 only; driver control gates passed, standalone port semantic gap remains'},indent=2)+'\n')
(o/'independent-source-manifest.json').write_text(json.dumps({str(f.relative_to(p)):{'sha256':sha(f),'bytes':f.stat().st_size} for f in sorted(p.rglob('*')) if f.is_file() and 'build' not in f.relative_to(p).parts},indent=2)+'\n');(o/'output-manifest.json').write_text(json.dumps({str(f.relative_to(o)):{'sha256':sha(f),'bytes':f.stat().st_size} for f in sorted(o.rglob('*')) if f.is_file() and f.name!='output-manifest.json'},indent=2)+'\n');print(json.dumps(status))
v6=o.parent/'pic-native-completion-review-v6-20261029'
(v6/'independent-status.json').write_text(json.dumps({'status':'immutable_v6_source_inspection_only_superseded_before_execution','standalone_execution_performed':False,'interval_and_common_clock_port_limits_found_by_source_inspection':True,'shared_compound_core_and_clock_fixes_executed_in':str(o/'handoff-review.json')},indent=2)+'\n')
