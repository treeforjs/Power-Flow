from pathlib import Path
from fractions import Fraction as F
import os,json,hashlib,math
p=Path(__file__).parent;s=p/'snapshot';o=Path('/home/treefor/sim-workspace/powerflow/outputs/pic-water-native-subflow-review-20261029')
assert sorted(os.sched_getaffinity(0))==[7]
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
for rec in json.loads((o/'input-bindings.json').read_text()):assert sha(Path(rec['snapshot']))==rec['sha256']
assert sha(s/'source/native_parent_subflow.hpp')=='96e9b3a5f5ddf5b571a84314cc7c6fe9944bbc83ea6604e58b6cf4588663415b'
r={x['name']:x for x in json.loads((o/'probes.json').read_text())};rate=r['compensated_carrier_first_moment'];assert rate['accepted']
ne=sum(map(F,rate['electron_count']),F());nn=sum(map(F,rate['target_count']),F())
assert sum(map(F,rate['electron_carrier_weights']),F())==ne and sum(map(F,rate['target_carrier_weights']),F())==nn
rows=rate['rows'];assert len(rows)==len(rate['electron_carrier_weights'])*len(rate['target_carrier_weights'])==6
coefficient=F(rows[0]['sigma_m2'])*F(rows[0]['flux_m_s'])/F(rows[0]['volume_m3'])
observed=sum((F(x['event_weight'])*F(x['rate_s']) for x in rows),F());expected=ne*nn*coefficient
for x in rows:assert F(x['event_weight'])==min(F(x['projectile_before']),F(x['target_before']))
relative=float((observed-expected)/expected);assert abs(relative)<=8*math.ulp(1.)
oracle={'method':'exact Fraction from stored binary64 native leaf counts, weights, cross section, flux, volume and rounded rates','logical_electron_count_exact':str(ne),'logical_target_count_exact':str(nn),'electron_carriers':len(rate['electron_carrier_weights']),'target_carriers':len(rate['target_carrier_weights']),'rows':len(rows),'both_carrier_sum_counts_exact':True,'weighted_event_rate_relative_rounding_defect':relative,'weighted_event_rate_expected_per_s':float(expected),'weighted_event_rate_sum_per_s':float(observed),'relation':'sum min(w_e,w_n)*max(w_e,w_n)*sigma*g/V = exact logical_Ne*logical_Nn*sigma*g/V up to individual propensity binary64 rounding'}
(o/'weighted-rate-oracle.json').write_text(json.dumps(oracle,indent=2)+'\n')
assert r['old_checkpoint_after_stochastic_noevent']['old_checkpoint_restored_exactly'] and r['old_checkpoint_after_stochastic_noevent']['arena_prefix_unchanged'] and r['old_checkpoint_after_stochastic_noevent']['clock_rng_advanced']
assert r['constructor_wrong_binding_file_noevent']['accepted'] and not r['constructor_wrong_binding_file_event']['accepted']
volume=r['cut_foil_collision_cell'];ratio=volume['grid_volume_m3']/volume['vacuum_intersection_m3']
status={'status':'archived_428_driver_review_findings_returned_not_accepted','reviewed_driver_sha256':sha(s/'source/native_parent_subflow.hpp'),'reviewed_native_identity':json.loads((s/'native-handoff.json').read_text())['source_identity'],'archived_commit':'f2d84f3acf3f60cc23deff809d527279cfe2f970','original_checks':428,'independent_probe_records':len(r),'noevent_old_unmodified_checkpoint_accepted':True,'foreign_binding_file_accepted_through_noevent_publish':True,'foreign_binding_file_refused_once_event_selected':True,'exhausted_interval_accepted_at_restore':True,'compensated_positive_leaf_first_moment_correct':True,'weighted_rate_relative_rounding_defect':relative,'cut_foil_grid_allowed':True,'cut_cell_full_over_vacuum_volume':ratio,'aligned_45x45x5_volume_correct':True,'root_live_source_or_paused_run_changed':False,'CPU':7,'nice':15,'threads':1,'frozen_userspace':True,'GPU':False,'v5_or_v6_driver_rebase_qualified':False,'full_water_or_physical_time_qualified':False}
(o/'validation.json').write_text(json.dumps(status,indent=2)+'\n')
(o/'handoff-review.json').write_text(json.dumps({'status':status['status'],'review':str(o/'review.txt'),'validation':str(o/'validation.json'),'probes':str(o/'probes.json'),'weighted_rate_oracle':str(o/'weighted-rate-oracle.json'),'runtime':str(p),'reviewed_header_sha256':status['reviewed_driver_sha256'],'scope':'archived428/96e9 and arena v4 only; current live edits/v5/v6 not covered'},indent=2)+'\n')
(o/'independent-source-manifest.json').write_text(json.dumps({str(f.relative_to(p)):{'sha256':sha(f),'bytes':f.stat().st_size} for f in sorted(p.rglob('*')) if f.is_file() and 'build' not in f.relative_to(p).parts},indent=2)+'\n')
(o/'output-manifest.json').write_text(json.dumps({str(f.relative_to(o)):{'sha256':sha(f),'bytes':f.stat().st_size} for f in sorted(o.rglob('*')) if f.is_file() and f.name!='output-manifest.json'},indent=2)+'\n');print(json.dumps(status));print(json.dumps(oracle))
