import collections,csv,hashlib,json,math,os,pathlib,statistics
base=pathlib.Path('<workspace>/outputs/pic-magnetized-emission-20261029')
assert os.sched_getaffinity(0)=={7} and os.getpriority(os.PRIO_PROCESS,0)>=15
rows=[]
for p in (base/'rate-map').glob('guard-rate-map-*.csv'):rows.extend(csv.DictReader(p.open()))
rows.sort(key=lambda x:(int(x['legacy_row']),float(x['T_K'])))
assert len(rows)==3168
assert len({(int(x['legacy_row']),float(x['T_K'])) for x in rows})==3168
assert all(x['actual_magnetized_current_known']=='0' for x in rows)
summary=[]
for key in sorted({(x['face'],float(x['T_K'])) for x in rows}):
 r=[x for x in rows if (x['face'],float(x['T_K']))==key]
 active=r[0]['SN_guard_applies']=='1'
 item={'face':key[0],'T_K':key[1],'legacy_spatial_probes':len(r),'SN_guard_applies':active,
       'signed_F_range_V_m':[min(float(x['signed_F_V_m']) for x in r),max(float(x['signed_F_V_m']) for x in r)],
       'nominal_supply_current_range_A_m2':[min(float(x['nominal_MG_current_A_m2']) for x in r),max(float(x['nominal_MG_current_A_m2']) for x in r)],
       'actual_magnetized_current':None,'physical_replacement_admitted':False,
       'strict_tangential_B_probe_count':sum(x['geometry_class']=='purely_tangential' for x in r)}
 if active:
  fractions=[float(x['nominal_rejected_fraction']) for x in r]
  numerator=sum(float(x['nominal_rejected_current_A_m2']) for x in r)
  denominator=sum(float(x['nominal_MG_current_A_m2']) for x in r)
  item.update({'nominal_guard_rejected_fraction_range':[min(fractions),max(fractions)],
       'nominal_guard_rejected_fraction_median_over_spatial_probes':statistics.median(fractions),
       'nominal_rejected_current_range_A_m2':[min(float(x['nominal_rejected_current_A_m2']) for x in r),max(float(x['nominal_rejected_current_A_m2']) for x in r)],
       'equal_probe_nominal_current_weighted_rejected_fraction':numerator/denominator,
       'equal_probe_weighting_is_not_surface_area_integration':True})
 summary.append(item)
with (base/'rate-weighted-admission-map.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
probes=pathlib.Path('<workspace>/outputs/cathode-emission-audit-20261029/fields/probes.csv')
report={'schema':'nominal-MG-rate-weighted-magnetic-applicability-map-v1',
 'coverage':{'legacy_field_probes':1056,'temperatures_K':[300,773.15,1490],'states':3168,'all_rows_retained':True},
 'source_phi_eV':4.3,'source_phi_status':'clean SS304 benchmark; no temperature/surface calibration',
 'nominal_guard':0.2,'image_energy_fraction':0.03,'analytic_transverse_integration':True,
 'native_predicate_and_refinement_qualification':json.loads((base/'guard-assessment-qualification.json').read_text()),
 'probe_sha256':hashlib.sha256(probes.read_bytes()).hexdigest(),
 'probe_scope':'legacy finite-offset magnitudes; not qualified one-sided native patch fields',
 'rate_scope':'rejected contribution to nominal MG supply only; actual magnetized supply remains unknown',
 'retarding_scope':'thermal supply retained; SN predicate bypass is not proof of no magnetic surface correction',
 'oblique_scope':'strict tangential replacement does not cover any nonzero B_normal without a separate justified bound/model',
 'zero_unknown_contributions':False,'physical_foil_model_admitted':False,'production_resumed':False,'GPU_executed':False,
 'summary':summary}
(base/'rate-weighted-admission-map.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'states':len(rows),'summary':summary},indent=2))
