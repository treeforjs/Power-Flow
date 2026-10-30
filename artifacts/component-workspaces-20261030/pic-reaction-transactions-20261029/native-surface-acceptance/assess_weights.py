from pathlib import Path
import fcntl,hashlib,json,os,subprocess
from fractions import Fraction
P=Path('/home/treefor/projects/runtimes/pic-reaction-transactions-20261029/native-surface-acceptance');R=Path('/home/treefor/sim-workspace/powerflow/outputs/pic-profile/reaction-transactions-20261029/native-surface-acceptance');base=P.parent/'native-completion-clock';evidence=R.parent/'native-completion-clock'
assert os.sched_getaffinity(0)=={7}and os.getpriority(os.PRIO_PROCESS,0)==15
with(R.parent.parent/'testing-core7.lock').open('a')as h:
 try:fcntl.flock(h,fcntl.LOCK_EX|fcntl.LOCK_NB)
 except BlockingIOError:pass
 else:raise RuntimeError('CPU7 lease required')
C=json.loads((evidence/'commands.json').read_text());B=P/'build';B.mkdir(exist_ok=True)
compile=list(C['build'][2]);compile[compile.index(str(base/'tests/test_arena.cpp'))]=str(P/'tests/assess_source_weights.cpp');compile[-1]=str(B/'assess_source_weights.o');compile[compile.index('-c'):compile.index('-c')]=['-isystem',str(P/'inputs/source-kernel-v1')]
with(R/'source-weight-build.log').open('w')as log:subprocess.run(compile,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=180)
link=[str(B/'assess_source_weights.o')if x==str(base/'build/tests.o')else x for x in C['build'][-1]];link[-1]=str(B/'assess_source_weights');link[-2:-2]=[str(P/'inputs/source-kernel-v1/libsurface_emission.a')]
with(R/'source-weight-link.log').open('w')as log:subprocess.run(link,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=60)
run=[C['run'][0],link[-1],str(R/'admitted-source-weight-controls.json')]
with(R/'source-weight-run.log').open('w')as log:subprocess.run(run,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=60)
j=json.loads((R/'admitted-source-weight-controls.json').read_text());F=lambda x:Fraction.from_float(float(x));count=lambda a:F(a[0])+F(a[1]);totals={'source_controls':len(j['cases']),'admitted_source_controls':0,'refused_source_controls':0,'stresses_with_early_count_refusal':0,'accepted_exact_debits':0,'max_positive_native_carriers':0}
for case in j['cases']:
 if not case['admitted_source']:totals['refused_source_controls']+=1;continue
 totals['admitted_source_controls']+=1
 for test in case['count_stresses']:
  initial=F(test['initial_neutral_macro_count']);used=Fraction()
  for receipt in test['receipts']:
   assert count(receipt['before'])-F(receipt['event_weight'])==count(receipt['after']);assert sum((F(x)for x in receipt['positive_native_carriers']),Fraction())==count(receipt['after']);used+=F(receipt['event_weight']);totals['accepted_exact_debits']+=1;totals['max_positive_native_carriers']=max(totals['max_positive_native_carriers'],len(receipt['positive_native_carriers']))
  assert initial-used==count(test['final_count']);totals['stresses_with_early_count_refusal']+=int(test['stopped']and'Count'in test['reason'])
(R/'source-weight-oracle.json').write_text(json.dumps({'status':'passed','totals':totals,'exact_fraction_native_counts':True,'limits':'Conditional source controls and forced domain stress, not actual C19 fields/neutral state, sampled collision probabilities, or a coupled physical pilot. Cannot infer full-run admitted weight range.'},indent=2)+'\n')
(R/'source-weight-commands.json').write_text(json.dumps({'compile':compile,'link':link,'run':run,'base_archive_sha256':'441e3c0047d67f74046733d3cb0bf7a424d530d71536e4bc05a31537f0b797e0','source_kernel_library_sha256':hashlib.sha256((P/'inputs/source-kernel-v1/libsurface_emission.a').read_bytes()).hexdigest()},indent=2)+'\n')
print(json.dumps(totals))
