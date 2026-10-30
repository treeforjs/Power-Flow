from pathlib import Path
import os,json,hashlib
p=Path(__file__).parent;s=p/'snapshot';o=Path('<workspace>/outputs/pic-native-arena-review-v5-20261029')
assert sorted(os.sched_getaffinity(0))==[7]
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
identity=json.loads((o/'source-identity.json').read_text());assert identity['sha256']=='e892ca6803a0646e75a33a169525a416fc2fb18526c9f5f1a46822cdf1cc8898'
assert sha(p/'source-v5.tar.gz')=='eeb47dfcc4631fc9d0d01cca271b2520a7b9331fa10e7619b12c7439e33f73f1'
for name,want in identity['files'].items():assert sha(s/name)==want,name
old=Path('<workspace>/outputs/pic-native-arena-review-20261029')
for name,rec in json.loads((old/'output-manifest.json').read_text()).items():assert sha(old/name)==rec['sha256'],name
probes={r['name']:r for r in json.loads((o/'probes.json').read_text())};oracle=json.loads((o/'independent-oracle.json').read_text())
for name in ('sink_count_conflicts_with_journal','sink_species_conflicts_with_journal','sink_reference_conflicts_with_catalog','arbitrary_sink_object','journal_event_weight_conflicts_with_counts','journal_negative_remainder','journal_unsupported_selection_method','journal_negative_elapsed_clock','native_mass_tamper','changed_ion_catalog_mass','subnormal_charged_carrier','catalog_INT_MIN_charge','catalog_stoichiometry_int_overflow'):assert not probes[name]['accepted'],name
for name in ('exclusive_owner_type_traits','copy_move_deletion'):assert not probes[name]['copy_constructible'] and not probes[name]['move_constructible']
assert probes['actual_parent_publish_rollback']['accepted']
controls=json.loads((o/'restore-controls.json').read_text());assert len(controls)==2 and all(x['accepted'] for x in controls)
counts=oracle['count_summary'];assert counts['count_cases']==377 and counts['carrier_admissions']==353 and counts['exact_debit_admissions']==215
status={'status':'focused_v5_fix_family_independently_verified','reviewed_source_identity':identity['sha256'],'archive_sha256':sha(p/'source-v5.tar.gz'),'original_v3_evidence_unchanged':True,'original_eight_restart_failure_variants_rejected':True,'untouched_restart_positive_controls':2,'copy_and_move_deleted':True,'stoichiometry_overflow_rejected':True,'charge_INT_MIN_rejected_without_abs_UB':True,'charged_subnormal_conversion_rejected':True,'mass_tamper_rejected':True,'count_input_executions':377,'exact_carrier_admissions':353,'exact_debit_admissions':215,'admitted_count_errors':0,'maximum_carriers':40,'actual_parent_exact_charge_H_O_and_120_digit_energy_passed':True,'native_and_arena_independently_recompiled':True,'CPU':7,'nice':15,'threads':1,'frozen_userspace':True,'GPU':False,'live_source_edited':False,'production_or_paused_history_changed':False,'v6_exposure_API_or_finite_rate_driver_qualified':False,'full_water_or_coupled_physics_qualified':False}
(o/'validation.json').write_text(json.dumps(status,indent=2)+'\n');(o/'handoff-review.json').write_text(json.dumps({'status':status['status'],'review':str(o/'review.txt'),'validation':str(o/'validation.json'),'oracle':str(o/'independent-oracle.json'),'probes':str(o/'probes.json'),'positive_controls':str(o/'restore-controls.json'),'source_identity':identity['sha256'],'scope':'exact v3 failure-family fixes only; subsequent completion/driver APIs remain separate'},indent=2)+'\n')
(o/'independent-source-manifest.json').write_text(json.dumps({str(f.relative_to(p)):{'sha256':sha(f),'bytes':f.stat().st_size} for f in sorted(p.rglob('*')) if f.is_file() and 'build' not in f.relative_to(p).parts},indent=2)+'\n')
(o/'output-manifest.json').write_text(json.dumps({str(f.relative_to(o)):{'sha256':sha(f),'bytes':f.stat().st_size} for f in sorted(o.rglob('*')) if f.is_file() and f.name!='output-manifest.json'},indent=2)+'\n');print(json.dumps(status))
