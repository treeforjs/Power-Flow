from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal as D,localcontext
import os,json,hashlib,csv
p=Path(__file__).parent;o=Path('/home/treefor/sim-workspace/powerflow/outputs/pic-generic-products-native-audit-20261029');assert sorted(os.sched_getaffinity(0))==[7]
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest();assert sha(p/'products-source.tar.gz')=='7c7feb6e315f0433cfb6b1635ded4640c34c3191e9fe0bd385c0b752347f8596'
d=json.loads((o/'probes.json').read_text());spec={x['id']:x for x in d['catalog']};rows=d['records'];assert len(rows)==17
accepted=[r for r in rows if r['accepted']];assert len(accepted)==13
def totals(pop):
 out={'H':F(),'O':F(),'q_e':F(),'electron':F()}
 for packet in pop:
  n=sum(map(F,packet['count']),F());s=spec[packet['state']]
  for key in ('H','O','q_e'):out[key]+=n*s[key]
  if packet['state']=='electron':out['electron']+=n
 return out
proofs=[]
with localcontext() as ctx:
 ctx.prec=120;c=D.from_float(299792458.);eps=D.from_float(2**-52)
 def kinetic(m,u):return m*sum(x*x for x in u)/((1+sum(x*x for x in u)/(c*c)).sqrt()+1)
 for row in accepted:
  before=totals(row['before_population']);after=totals(row['after_population']);assert all(before[k]==after[k] for k in ('H','O','q_e'))
  expected=F(row['event_weight'])*(-1 if row['capture'] else 0 if row['model']=='A1B1_fixed_loss_dipole_ground_v1' else 1)
  assert after['electron']-before['electron']==expected
  audit=row['audit'];u0=list(map(D.from_float,audit['primary_u']));ut=list(map(D.from_float,audit['target_u']));me=D.from_float(spec['electron']['mass_kg']);M=D.from_float(spec['H2O']['mass_kg']);before_energy=kinetic(me,u0)+kinetic(M,ut);after_energy=D(0);dp=[-me*u0[i]-M*ut[i] for i in range(3)];activity=sum(abs(me*x) for x in u0)+sum(abs(M*x) for x in ut)
  for product in audit['products']:
   s=spec[product['state']];m=D.from_float(s['mass_kg']);u=list(map(D.from_float,product['u']));after_energy+=kinetic(m,u)+D.from_float(s['internal_J'])+D.from_float(s['reference_J'])
   for i in range(3):dp[i]+=m*u[i];activity+=abs(m*u[i])
  W=D.from_float(row['event_weight']);energy=(after_energy-before_energy)*W;bound=D.from_float(audit['energy_roundoff_bound_J']);assert abs(energy)<=bound
  assert max(map(abs,dp))<=256*eps*activity
  assert row['accepted_native_node_delta_Linf_C']==0 and row['horizontal_transfers_zero'] and row['vertical_transfers_zero'] and row['rollback_exact']
  proofs.append({'model':row['model'],'OH_level':row['OH_level'],'weight':row['event_weight'],'H_exact':str(after['H']),'O_exact':str(after['O']),'charge_e_exact':str(after['q_e']),'electron_change_exact':str(expected),'weighted_energy_defect_J':float(energy),'original_native_energy_bound_J':float(bound),'energy_bound_fraction':float(abs(energy)/bound),'momentum_bound_fraction':float(max(map(abs,dp))/(256*eps*activity))})
refused=[r for r in rows if not r['accepted']];assert len(refused)==4 and all(r['refusal_unchanged'] for r in refused)
assert sum('nonlocal charge' in r['refusal'] for r in refused)==3
assert any(r['wrong_product_contract'] and r['refusal']=='channel products mismatch' for r in refused)
assert any(r.get('capture_restart_population_exact') for r in accepted)
site=json.loads((o/'site-capability.json').read_text());assert site['admitted'] and site['native_node_delta_Linf_C']==0 and not site['physical_collision_site_model_qualified']
oracle={'accepted_count':13,'exact_species_atom_charge_and_electron_multiplicity_checks':True,'independent_energy_momentum_decimal_digits':120,'maximum_energy_bound_fraction':max(x['energy_bound_fraction'] for x in proofs),'maximum_momentum_bound_fraction':max(x['momentum_bound_fraction'] for x in proofs),'maximum_kernel_catalog_mass_ULP':max(r['kernel_catalog_mass_max_ULP'] for r in rows),'all_tested_accepted_native_Q1_nodal_deltas_zero':True,'proofs':proofs}
(o/'independent-oracle.json').write_text(json.dumps(oracle,indent=2)+'\n')
with (o/'admission-matrix.csv').open('w') as f:
 writer=csv.writer(f);writer.writerow(['model','OH_level','E_eV','capture','placement_0_colocated_1_separate_xy_2_separate_z','electron_weight','water_weight','accepted','reason','kernel_catalog_mass_ULP','attempted_node_delta_Linf_C','accepted_node_delta_Linf_C'])
 for r in rows:writer.writerow([r['model'],r['OH_level'],r['energy_eV'],r['capture'],r['placement'],r['electron_weight'],r['water_weight'],r['accepted'],r.get('refusal',''),r['kernel_catalog_mass_max_ULP'],r['attempted_cohort_node_delta_Linf_C'],r.get('accepted_native_node_delta_Linf_C','')])
requirements={'state_catalog':'one explicitly bound common mass/reference/internal-state basis; e_g→electron mapping and H/O/q are explicit','outcome':'variable-length products and scattered/captured primary; convert every product momentum to proper velocity with catalog mass; target-rest→lab boost required for moving targets','OH_v_family':'base Hminus event rate sampled once, OHv0..4 conditional state draw once; fixed ArenaChannel product multiset needs allowed product sets or explicit provenance-linked branch contract; never repeat the full base rate for each state','completion':'current completion binds one channel; generic competing SSA needs an immutable mechanism/channel-set binding and per-child registered membership/data/model validation, preserving strict clock checks','spatial':'explicit channel collision-site contract; target-site DEA with separated positions remains unsupported without qualified reaction current/work; co-location is supported; all-primary placement is API-capable but physically unqualified here','missing_rates':'do not filter unsupported nonlocal attachment pairs or set missing rates to zero; full mechanism must reject positive/unknown holes, with only explicit subset comparison eligible to omit application while retaining incompleteness','one_keV_boundary':'choose one low empirical or high outer3 model at exact1keV, never add overlapping prescriptions','no_new_current_path_or_field_solver':'none implemented or assumed','target_frame_test_scope':'stationary targets only; existing qualified Lorentz helper can supply future moving-target wrapper, but this audit does not qualify that integration'}
(o/'adapter-requirements.json').write_text(json.dumps(requirements,indent=2)+'\n')
status={'status':'narrow_actual_v5_native_product_admission_audit_complete','products_archive_sha256':sha(p/'products-source.tar.gz'),'products_commit':'c9f6fd294ad32bfa7317070665a7803d9b49cba1','dataset_sha256':sha(p/'products/data/collisions/h2o-gas-20261029-outer3.json'),'native_source_identity':'5fa3a59f863e98f21e20af698b9c16417f3749344b555f18e3d2cba1e4f5a705','models_tested':7,'main_records':17,'native_admitted':13,'nonlocal_target_site_capture_refused_unchanged':3,'fixed_OHv_state_mismatch_refused_unchanged':1,'captured_primary_and_partial_remainders_supported':True,'capture_restart_and_all_accepted_rollbacks_passed':True,'native_Q1_node_deltas_zero_in_all_accepted_fixtures':True,'edge_transfer_arrays_zero_in_all_accepted_fixtures':True,'nodewise_shape_source':'qualified native locate/FE/element-dof/charge assembly via read-only test inspection','separate_all_primary_site_capability_admitted':True,'all_primary_site_physical_model_qualified':False,'conditional_OH_mass_rounding_max_ULP':oracle['maximum_kernel_catalog_mass_ULP'],'single_channel_completion_not_generic_mechanism':True,'unmapped_2a1_m2':d['unmapped_2a1_m2'],'mapped_four_valence_fraction':d['mapped_four_valence_fraction'],'full_mechanism_or_physical_history_qualified':False,'moving_target_integration_qualified':False,'no_rate_omission_or_current_path_implemented':True,'live_source_changed':False,'GPU':False,'physical_time_advanced':False,'CPU':7,'nice':15,'threads':1,'frozen_userspace':True}
(o/'validation.json').write_text(json.dumps(status,indent=2)+'\n');(o/'handoff-review.json').write_text(json.dumps({'status':status['status'],'review':str(o/'review.txt'),'validation':str(o/'validation.json'),'admission_matrix':str(o/'admission-matrix.csv'),'requirements':str(o/'adapter-requirements.json'),'oracle':str(o/'independent-oracle.json'),'main_probes':str(o/'probes.json'),'site_capability_only':str(o/'site-capability.json')},indent=2)+'\n')
(o/'source-manifest.json').write_text(json.dumps({str(f.relative_to(p)):{'sha256':sha(f),'bytes':f.stat().st_size} for f in sorted(p.rglob('*')) if f.is_file() and 'build' not in f.relative_to(p).parts},indent=2)+'\n');(o/'output-manifest.json').write_text(json.dumps({str(f.relative_to(o)):{'sha256':sha(f),'bytes':f.stat().st_size} for f in sorted(o.rglob('*')) if f.is_file() and f.name!='output-manifest.json'},indent=2)+'\n');print(json.dumps(status))
