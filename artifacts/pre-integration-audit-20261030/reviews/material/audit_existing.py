from pathlib import Path
from fractions import Fraction as F
import json,hashlib,os,fcntl,time,math
O=Path('/home/treefor/sim-workspace/powerflow/outputs')
R=Path(__file__).resolve().parent
V=O/'ss304-thermal-expansion-20261029/foil-support-v4'
M=O/'ss304-thermal-properties-v2-20261029'
S=O/'ss304-spatial-metal-transfer-20261029'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text())
F64=lambda x:F.from_float(float(x))
def flat(v):
 if isinstance(v,list):return [y for x in v for y in flat(x)]
 return [v]
def sumf(v):return sum(map(F64,flat(v)),F())
def hold():
 p=Path('/proc/278902');x=(p/'stat').read_text().split(') ',1)[1].split()
 assert x[0]=='T'and x[19]=='5288339'and'frozen 1'in Path('/sys/fs/cgroup/powerflow-pic-gpu-reserved/cgroup.events').read_text()
 return {'pid':278902,'start_ticks':x[19],'state':x[0],'frozen':True,'exe':os.readlink(p/'exe'),'exe_sha256':sha(p/'exe')}
os.sched_setaffinity(0,{7});os.nice(max(0,15-os.getpriority(os.PRIO_PROCESS,0)))
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
with (O/'pic-profile/testing-core7.lock').open('a')as lease:
 fcntl.flock(lease,fcntl.LOCK_EX);start=time.time();before=hold();checks=[];manifest_checks=[]
 def check(ok,name):
  checks.append({'name':name,'passed':bool(ok)})
  if not ok:raise AssertionError(name)
 vm=read(V/'source-manifest.json');mm=read(M/'source-manifest.json');sm=read(S/'source-manifest.json')
 for label,manifest,key in [('expansion-v4',vm,'source_files_sha256'),('material-v2',mm,'files')]:
  root=Path(manifest['source_root']);mismatches=[rel for rel,digest in manifest[key].items()if sha(root/rel)!=digest]
  check(not mismatches,label+'/sealed source hashes');manifest_checks.append({'component':label,'files_verified':len(manifest[key]),'mismatches':mismatches})
 check(sha(M/'source-manifest.json')=='1a1845e04f1f6fda33338d9bd284f5e1ead9cbc4b7747d0c6355e32cb5661c2d','material v2 owner seal identity')
 spatialroot=Path('/home/treefor/projects/runtimes/ss304-spatial-metal-transfer-20261029/source')
 check(all(sha(spatialroot/rel)==digest for rel,digest in sm['files'].items()),'spatial source hashes')
 manifest_checks.append({'component':'spatial-mass-energy','files_verified':len(sm['files']),'mismatches':[]})
 for label,record in [('material-library',mm['compiled_library']),*vm['native_executables'].items(),('frozen-MFEM',vm['frozen_mfem_archive'])]:check(sha(record['path'])==record['sha256'],label+'/binary hash')
 rows=[];source_bound_factor=F64(1e-10+256*math.ulp(1.));jointmax=0.;nodalmax=0.;pathmax=0.;zero_bridges=0
 for case in ['expansion','compression']:
  history=read(V/'qualified-final'/f'{case}-history.json')
  for row in history:
   a,b=[row[k]for k in ['native_source_state0','native_source_state1']];budget=row['native_budget'];mech=row['mechanical_exchange'];new=b['gas_packets'][len(a['gas_packets']):]
   reservoirs=[]
   for i,(p0,p1)in enumerate(zip(a['reservoirs'],b['reservoirs'])):
    if not a['patches'][i]['water_enabled']:continue
    old=F64(p0['adsorbed_N_hi'])+F64(p0['adsorbed_N_lo']);now=F64(p1['adsorbed_N_hi'])+F64(p1['adsorbed_N_lo']);packet=sum((F64(p['weight'])for p in new if p['patch_index']==i),F());gas=sum((F64(p['weight'])for p in b['gas_packets']if p['patch_index']==i),F());stored=F64(p1['gas_N_hi'])+F64(p1['gas_N_lo'])
    check(old-now==packet,f'{case}/{b["revision"]}/exact adsorbed delta equals actual packet weights')
    check(stored==gas,f'{case}/{b["revision"]}/exact component gas Count2 equals actual packets')
    reservoirs.append({'patch':i,'source_debit_minus_native_packet_weights':str(old-now-packet),'gas_Count2_minus_actual_packet_weights':str(stored-gas)})
   check(all(p['reference_J']==0 and p['charge_C']==0 for p in new),f'{case}/{b["revision"]}/neutral water reference and charge')
   check(F64(budget['particle_charge_change_C'])+F64(budget['electrode_charge_change_C'])==0,f'{case}/{b["revision"]}/source electrode charge closure')
   paths=b['paths'][len(a['paths']):];E=[F64(row['source_stage'][f'E_{d}_V_m'])for d in ['x','y','z']];pathdetails=[]
   for path in paths:
    calculated=F64(path['charge_C'])*sum((E[i]*(F64(path[f'to_{d}_m'])-F64(path[f'from_{d}_m']))for i,d in enumerate(['x','y','z'])),F());stored=F64(path['field_work_J']);bound=source_bound_factor*(abs(calculated)+abs(stored));error=abs(calculated-stored);check(error<=bound,f'{case}/{b["revision"]}/weighted source path qE dx')
    pathmax=max(pathmax,float(error));zero_bridges+=path['start_time_s']==path['end_time_s']and any(path[f'to_{d}_m']!=path[f'from_{d}_m']for d in ['x','y','z'])
   pathsum=sum((F64(p['field_work_J'])for p in paths),F());pathdelta=pathsum-F64(budget['source_field_work_J']);check(abs(pathdelta)<=F64(row['source_partition_bound_J']),f'{case}/{b["revision"]}/path work matches incremental source work')
   # Reconstruct the published subsystem equation from separate owner operands.
   joint=F64(budget['thermal_delta_J'])+sumf(mech['cell_elastic_internal_energy_change_J'])+F64(mech['kinetic_energy_change_J'])+sum((F64(budget[k])for k in ['source_particle_K_change_J','source_electron_reference_change_J','source_image_energy_change_J','neutral_K_J','neutral_internal_J','adsorbed_reference_change_J']),F())-F64(budget['source_field_work_J'])
   check(abs(joint)<=F64(row['energy_tolerance_J']),f'{case}/{b["revision"]}/joint equation within existing owner bound');jointmax=max(jointmax,float(abs(joint)))
   node=budget['nodal_ledger'];delta=flat(node['delta_H_J']);channels={k:flat(v)for k,v in node.items()};eq=[]
   for i in range(len(delta)):
    residual=F64(delta[i])-F64(channels['mechanical_load_J'][i])-F64(channels['electron_load_J'][i])-F64(channels['water_load_J'][i])+F64(channels['conduction_out_J'][i]);eq.append(float(abs(residual)))
   nodalmax=max(nodalmax,max(eq));check(max(eq)<=budget['roundoff_bound_J'],f'{case}/{b["revision"]}/native nodal equation within existing native bound')
   rows.append({'case':case,'revision':b['revision'],'reservoir_checks':reservoirs,'new_paths':len(paths),'new_water_packets':len(new),'joint_residual_exact_operand_sum_J':float(joint),'joint_existing_bound_J':row['energy_tolerance_J'],'path_sum_minus_incremental_work_J':float(pathdelta),'maximum_nodal_equation_residual_J':max(eq)})
 pilot=read(O/'ss304-material-feedback-pilot-20261029/results.json');geom=read(V/'actual-pilot-geometry-summary.json');component=read(V/'qualified-final/summary.json');v4checkpoint=read(V/'qualified-final/expansion-checkpoint/joint.json')
 evidence=[M/'HANDOFF.md',M/'source-manifest.json',M/'verification.json',M/'thermal-material-port-v2.json',S/'HANDOFF.md',S/'source-manifest.json',S/'disposition.json',S/'spatial-material-port-v1.txt',O/'pic-profile/exact-count-independent-review-20261029/REVIEW.md',O/'pic-profile/exact-count-independent-review-20261029/disposition.json',V/'HANDOFF.md',V/'handoff.json',V/'source-manifest.json',V/'compatibility.json',V/'independent-thermodynamic-review.md',V/'native-independent-review.md',V/'qualified-final/summary.json',V/'qualified-final/expansion-history.json',V/'qualified-final/compression-history.json',V/'qualified-final/expansion-checkpoint/joint.json',V/'qualified-final/compression-checkpoint/joint.json',V/'actual-pilot-geometry-summary.json',V/'pilot-audit/mesh32-audit-001-supported-prefix.json',O/'ss304-material-feedback-pilot-20261029/HANDOFF.txt',O/'ss304-material-feedback-pilot-20261029/results.json',O/'ss304-material-feedback-pilot-20261029/component-contract.json',O/'ss304-material-feedback-pilot-20261029/restart-proof.json']
 result={'status':'existing_evidence_and_source_audit_complete_not_full_driver_validation','fresh_simulations_or_component_test_runs':0,'existing_history_intervals_reconstructed':len(rows),'assertions_on_existing_evidence':len(checks),'checks':checks,'manifest_checks':manifest_checks,'rows':rows,'maximum_reconstructed_joint_residual_J':jointmax,'maximum_reconstructed_nodal_residual_J':nodalmax,'maximum_path_operand_difference_J':pathmax,'retained_zero_duration_nonzero_displacement_bridges':zero_bridges,'density_scenarios':{'joint_v4_component_kg_m3':v4checkpoint['config']['caloric']['reference_density_kg_m3'],'spatial_and_material_pilot_kg_m3':8000,'actual_pilot_geometry_scenario':geom['density'],'no_implicit_equivalence_claim':True},'mechanics_actual_prefix':{'accepted_intervals':geom['accepted_intervals'],'accepted_time_interval_s':geom['accepted_time_interval_s'],'final_native_temperature_range_K':geom['final_native_temperature_range_K'],'next_stage_stop':geom['prefix_termination_reason'],'formal_source_build_run_seal':geom['formal_pilot_source_to_binary_build_seal_verified']},'solid_only_refinement':pilot['comparison'],'evidence_sha256':{str(p):sha(p)for p in evidence},'CPU':[7],'nice':15,'threads':1,'lease_s':time.time()-start,'production_before':before,'production_after':hold(),'GPU':False,'source_changes':False,'guard_changes':False,'new_model_assumptions':False}
 (R/'evidence-audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k]for k in ['status','fresh_simulations_or_component_test_runs','assertions_on_existing_evidence','manifest_checks','maximum_reconstructed_joint_residual_J','maximum_reconstructed_nodal_residual_J','retained_zero_duration_nonzero_displacement_bridges','lease_s']},indent=2))
