from pathlib import Path
import fcntl,hashlib,json,os,subprocess,math
from fractions import Fraction
P=Path('/home/treefor/projects/runtimes/pic-reaction-transactions-20261029/native-exact-count');R=Path('/home/treefor/sim-workspace/powerflow/outputs/pic-profile/reaction-transactions-20261029/native-exact-count');frozen='/home/treefor/projects/runtimes/arch-base-20261028/run-frozen';kernel=P/'control-inputs/source-kernel';B=P/'build'
assert os.sched_getaffinity(0)=={7}and os.getpriority(os.PRIO_PROCESS,0)==15
with(R.parent.parent/'testing-core7.lock').open('a')as h:
 try:fcntl.flock(h,fcntl.LOCK_EX|fcntl.LOCK_NB)
 except BlockingIOError:pass
 else:raise RuntimeError('CPU7 lease required')
cmd=[frozen,'/usr/bin/g++','-std=c++20','-O2','-g','-Wall','-Wextra','-Wpedantic','-Werror','-I'+str(P/'source/include'),'-isystem',str(kernel),'-I'+str(P.parent/'source/vendor'),str(P/'tests/test_exact_count_controls.cpp'),str(kernel/'libsurface_emission.a'),'-lcrypto','-o',str(B/'test_exact_count_controls')]
with(R/'controls-build.log').open('w')as log:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=180)
inputs=P/'control-inputs/count2-cases.json';run=[frozen,str(B/'test_exact_count_controls'),str(inputs),str(R/'exact-count-controls.json')]
with(R/'controls-run.log').open('w')as log:subprocess.run(run,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=90)
j=json.loads((R/'exact-count-controls.json').read_text());F=lambda x:Fraction.from_float(float(x));sumf=lambda a:sum((F(x)for x in a),Fraction());tot={'count2_inputs':len(j['count2_cases']),'count2_admissions':0,'exact_debit_admissions':0,'count2_export_refusals_after_debit':0,'admitted_source_controls':0,'refused_source_controls':0,'actual_source_weight_debits':0,'source_export_refusals':0,'maximum_positive_carriers':0}
for row in j['count2_cases']:
 if not row['accepted']:continue
 tot['count2_admissions']+=1;value=F(row['hi'])+F(row['lo']);assert value==sumf(row['carriers'])==sumf(row['export_two_limbs']);assert float(value)==row['rounded'];tot['maximum_positive_carriers']=max(tot['maximum_positive_carriers'],len(row['carriers']))
 if 'debit'in row:
  valid=math.isfinite(row['debit'])and row['debit']>0 and F(row['debit'])<=value
  assert row['debit_accepted']==valid
  if not valid:continue
  tot['exact_debit_admissions']+=1;after=value-F(row['debit']);assert after==sumf(row['after_carriers']);assert float(after)==row['round_after'];tot['maximum_positive_carriers']=max(tot['maximum_positive_carriers'],len(row['after_carriers']))
  if 'after_two_limbs'in row:assert after==sumf(row['after_two_limbs'])
  else:
   tot['count2_export_refusals_after_debit']+=1;hi=float(after);residual=after-F(hi);assert F(float(residual))!=residual
for case in j['source_controls']:
 if not case['admitted_source']:
  assert 'magnetic energy criterion'in case['reason'];tot['refused_source_controls']+=1;continue
 tot['admitted_source_controls']+=1
 for sequence in case['sequences']:
  value=F(sequence['initial'])
  assert len(sequence['records'])==16
  for row in sequence['records']:
   assert value==sumf(row['before']);value-=F(row['weight']);assert value>=0 and value==sumf(row['after']);tot['actual_source_weight_debits']+=1;tot['maximum_positive_carriers']=max(tot['maximum_positive_carriers'],len(row['after']))
   if 'two_limb_export'in row:assert value==sumf(row['two_limb_export'])
   else:
    tot['source_export_refusals']+=1;hi=float(value);residual=value-F(hi);assert F(float(residual))!=residual
  assert value==sumf(sequence['final'])
assert tot['count2_admissions']==353 and tot['admitted_source_controls']==20 and tot['refused_source_controls']==4 and tot['actual_source_weight_debits']==5120
# New count encoding must leave the previously supported native math unchanged.
comparisons=[]
old=R.parent/'native-completion-scope/tests';new=R/'tests'
for name in ['case-1.json','case-2.json','case-3.json','case-4.json','neutral.json','compound.json','completion-zero.json','completion-reaction.json','scope-ssa.json','scope-forced-comparison.json']:
 a=json.loads((old/name).read_text())['payload'];b=json.loads((new/name).read_text())['payload'];assert a['native_bytes_hex']==b['native_bytes_hex']and a['saved_potential_V']==b['saved_potential_V'];comparisons.append({'bundle':name,'native_bytes_identical':True,'saved_potential_identical':True})
(R/'controls-oracle.json').write_text(json.dumps({'status':'passed','totals':tot,'independent_oracle':'exact Python Fraction and IEEE float rounding/export checks','supported_Count2_native_comparison':comparisons,'limits':'conditional source controls and assumed inventory scales; no sampling policy/rate/force/transport change or actual C19 pilot claim'},indent=2)+'\n')
(R/'controls-commands.json').write_text(json.dumps({'build':cmd,'run':run,'count2_cases_sha256':hashlib.sha256(inputs.read_bytes()).hexdigest(),'source_kernel_sha256':hashlib.sha256((kernel/'libsurface_emission.a').read_bytes()).hexdigest()},indent=2)+'\n');print(json.dumps(tot))
