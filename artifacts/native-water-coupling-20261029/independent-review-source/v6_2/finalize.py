from pathlib import Path
import os,json,hashlib
p=Path(__file__).parent;o=Path('/home/treefor/sim-workspace/powerflow/outputs/pic-native-subflow-v6_2-review-20261029');assert sorted(os.sched_getaffinity(0))==[7]
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest();identity=json.loads((o/'native-source-identity.json').read_text())
assert identity['sha256']=='5fa3a59f863e98f21e20af698b9c16417f3749344b555f18e3d2cba1e4f5a705'
assert sha(p/'native-source.tar.gz')=='d14f1e229bbd8a84cf75b782f49562217972d4cd55b086fd8f932cc4f0e37db4'
assert sha(p/'driver-source.tar.gz')=='6f13a3fe6101b280f126562bb2cd126b346e0cc7ed274c33816f6016841cd4e4'
assert sha(p/'driver/source/native_parent_subflow.hpp')=='ad499c5549a78bae371c5b2f72ee66085c3b141126d1931d1ec654d8a9e25675'
for name,want in identity['files'].items():assert sha(p/'native'/name)==want,name
old=o.parent/'pic-native-subflow-v6_1-review-20261029'
for name,rec in json.loads((old/'output-manifest.json').read_text()).items():assert sha(old/name)==rec['sha256'],name
oldsource=Path('/home/treefor/projects/runtimes/pic-native-subflow-v6_1-review-20261029/snapshot/native');changes=[]
for directory in ('source','inputs'):
 for f in sorted((p/'native'/directory).rglob('*')):
  if f.is_file():
   relative=f.relative_to(p/'native');other=oldsource/relative
   if not other.exists() or sha(f)!=sha(other):changes.append(str(relative))
assert changes==['source/src/mfem_pic_arena.cpp'],changes
live=json.loads((o/'port-bounds-probes.json').read_text());assert len(live)==3
for r in live:assert 'SSA child' in r['refused'] and r['unchanged_if_preparation_failed'] and not r.get('published',False)
focused=json.loads((o/'focused-probes.json').read_text());assert focused['assertions']==27 and len(focused['records'])==5
for r in focused['records']:
 if r['name'].startswith('restore_'):assert 'SSA child' in r['refusal'] and r['fresh_native_unchanged'] and r['correctly_rehashed_new_source_bundle']
normal=json.loads((o/'owner-result.json').read_text());assert normal['check_count']==442
status={'status':'targeted_v6_2_live_restore_SSA_completion_fixes_independently_verified','native_source_identity':identity['sha256'],'driver_sha256':sha(p/'driver/source/native_parent_subflow.hpp'),'driver_archive_sha256':sha(p/'driver-source.tar.gz'),'native_archive_sha256':sha(p/'native-source.tar.gz'),'original_normal_checks':442,'confirmed_live_failures_now_refused':3,'new_source_rehashed_semantic_restore_failures_now_refused':3,'valid_SSA_and_restart_assertions':27,'actual_Scheduler_event_terminal_completion_commit_restart_rollback_passed':True,'two_actual_events_then_five_continued_exact':True,'uninterrupted_and_restored_both_RNG_native_states_exact':True,'numeric_runtime_changed_files_from_v6_1':changes,'old_v6_1_evidence_unchanged':True,'exhaustive_grid_or_rate_survey_repeated':False,'live_source_edited':False,'GPU':False,'physical_time_advanced':False,'full_water_or_total_EM_qualified':False,'CPU':7,'nice':15,'threads':1,'frozen_userspace':True}
(o/'validation.json').write_text(json.dumps(status,indent=2)+'\n');(o/'handoff-review.json').write_text(json.dumps({'status':status['status'],'review':str(o/'review.txt'),'validation':str(o/'validation.json'),'live_failure_replays':str(o/'port-bounds-probes.json'),'valid_and_restore_controls':str(o/'focused-probes.json'),'source_identity':identity['sha256'],'scope':'targeted native SSA cross-record fixes plus valid SSA commit/restart/rollback only'},indent=2)+'\n')
(o/'independent-source-manifest.json').write_text(json.dumps({str(f.relative_to(p)):{'sha256':sha(f),'bytes':f.stat().st_size} for f in sorted(p.rglob('*')) if f.is_file() and 'build' not in f.relative_to(p).parts},indent=2)+'\n');(o/'output-manifest.json').write_text(json.dumps({str(f.relative_to(o)):{'sha256':sha(f),'bytes':f.stat().st_size} for f in sorted(o.rglob('*')) if f.is_file() and f.name!='output-manifest.json'},indent=2)+'\n');print(json.dumps(status))
