from pathlib import Path
import os,json,hashlib
p=Path(__file__).parent;s=p/'snapshot';o=Path('<workspace>/outputs/pic-native-arena-review-20261029')
assert sorted(os.sched_getaffinity(0))==[7]
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
identity=json.loads((o/'source-identity.json').read_text())
assert sha(p/'source-v3.tar.gz')=='52b758294d157b510757dd0447048c43d00d802d62c738beef45bc6cf60084c8'
for name,want in identity['files'].items():assert sha(s/name)==want,name
oracle=json.loads((o/'independent-oracle.json').read_text());c=oracle['count_summary'];r={x['name']:x for x in json.loads((o/'probes.json').read_text())}
assert c['count_cases']==377 and c['carrier_admissions']==353 and c['exact_debit_admissions']==215 and c['maximum_carriers']==40
for key in ('overdraw_admissions','carrier_sum_errors','nonpositive_or_nonfinite_carriers','debit_errors'):assert not c[key]
for key in ('copied_owner_forks_native_history','catalog_stoichiometry_int_overflow','sink_count_conflicts_with_journal','sink_species_conflicts_with_journal','sink_reference_conflicts_with_catalog','arbitrary_sink_object','journal_event_weight_conflicts_with_counts','journal_negative_remainder','journal_unsupported_selection_method','journal_negative_elapsed_clock'):assert r[key]['accepted']
status={'status':'material_native_arena_v3_findings_returned_not_accepted','reviewed_source_identity':identity['sha256'],'source_archive_sha256':sha(p/'source-v3.tar.gz'),'original_normal_checks':190,'independent_probe_records':len(r),'count_input_executions':377,'exact_carrier_admissions':353,'exact_debit_admissions':215,'accepted_count_errors':0,'maximum_positive_carriers':40,'actual_parent_charge_H_O_exact':True,'actual_parent_120_digit_energy_within_original_bound':True,'copy_owner_history_fork_reproduced':True,'sink_semantic_restart_corruptions_admitted':4,'reaction_semantic_restart_corruptions_admitted':4,'stoichiometry_integer_overflow_admitted':True,'subnormal_charged_carrier_deposition_underflow_admitted':True,'native_mass_and_changed_product_catalog_mass_rejected':True,'live_source_edited':False,'GPU':False,'production_or_paused_history_changed':False,'CPU':7,'nice':15,'threads':1,'frozen_userspace':True,'address_space_cap_KiB':1572864,'CPU_time_cap_s':180,'native_and_arena_recompiled_from_sealed_sources':True,'physical_H2O_or_native_driver_acceptance':False}
(o/'validation.json').write_text(json.dumps(status,indent=2)+'\n')
(o/'handoff-review.json').write_text(json.dumps({'status':status['status'],'review':str(o/'review.txt'),'validation':str(o/'validation.json'),'oracle':str(o/'independent-oracle.json'),'probes':str(o/'probes.json'),'source':str(p/'probe.cpp'),'runner':str(p/'run-review.py'),'reviewed_source_identity':identity['sha256'],'scope':'immutable v3 only; later fixes/compound/finite-rate driver not reviewed by this receipt'},indent=2)+'\n')
(o/'independent-source-manifest.json').write_text(json.dumps({str(f.relative_to(p)):{'sha256':sha(f),'bytes':f.stat().st_size} for f in sorted(p.rglob('*')) if f.is_file() and 'build' not in f.relative_to(p).parts},indent=2)+'\n')
(o/'output-manifest.json').write_text(json.dumps({str(f.relative_to(o)):{'sha256':sha(f),'bytes':f.stat().st_size} for f in sorted(o.rglob('*')) if f.is_file() and f.name!='output-manifest.json'},indent=2)+'\n')
print(json.dumps(status))
