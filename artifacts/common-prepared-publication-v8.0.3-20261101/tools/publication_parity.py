#!/usr/bin/python -I
"""Targeted accepted-trajectory parity for publication changes; no speed claim."""
from pathlib import Path
import importlib.util,json,os,subprocess,sys
R=Path(__file__).resolve().parents[1]
def module(name,p):
 s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
m=module('runner',R/'run_v8.py');cmp=module('compare',R/'tools/compare_v8.py')
BASE=R.parent/'pic-water-common-commit-v8.0.2-cpu-20261031'
BO=Path('/home/treefor/projects/.workspace/outputs')/BASE.name
label=sys.argv[1];tag=sys.argv[2]
if not m.inside_frozen() or not os.getenv('SIMQ_JOB_ID'):raise RuntimeError('frozen scheduled allocation required')
rec={'admission':m.admission('publication-parity'),'passed':False,'speedup_qualified':False,'comparisons':[],'commands':[]}
D=m.OUT/tag;D.mkdir(exist_ok=False)
source=m.own_sources();pins=m.sealed_pins()
baseline=json.loads((BO/'build-c2/build-receipt.json').read_text())
if baseline['binary_sha256']!='e9bccca0428cb3a3065141cd233dce2921a27423833837cb97d153124543a7ea':raise RuntimeError('baseline binary pin changed')
if baseline['status']!='passed' or m.sha(baseline['binary'])!=baseline['binary_sha256']:raise RuntimeError('actual qualified v8.0.2 baseline missing or changed')
rec['baseline_binary_sha256']=baseline['binary_sha256']
try:
 for name,args in [('populated20',['--fixture','qual','--steps','20','--frame-every','10']),('adaptive',['--fixture','dense','--blob-weight','8e10','--dt-mode','adaptive','--sync-every','5e-14','--window','5e-13','--progress-every','50'])]:
  runs=[]
  for which,root,b in [('baseline',BASE,'c2'),('candidate',R,label)]:
   out=D/(name+'-'+which);command=[sys.executable,'-I',str(root/'run_v8.py'),'run',b,'run',str(out),*args]
   p=subprocess.run(command,cwd=root,capture_output=True,text=True,timeout=180)
   (D/(name+'-'+which+'.log')).write_text(p.stdout+p.stderr)
   rr=json.loads((out/'run-receipt.json').read_text());summary=json.loads((out/'run-summary.json').read_text())
   if p.returncode or rr['status']!='passed' or not all(rr[k] for k in ('sealed_unchanged_after_run','own_sources_unchanged_after_run','binary_unchanged_after_run')):raise RuntimeError('refused trajectory or provenance: '+name+which)
   rec['commands'].append({'argv':command,'returncode':p.returncode,'receipt':str(out/'run-receipt.json'),'receipt_sha256':m.sha(out/'run-receipt.json')})
   runs.append((out,summary))
  a,b=runs
  if any(a[1][k]!=b[1][k] for k in ('steps_completed','stop','final_time_s')) or a[1]['steps_completed']==0:raise RuntimeError('accepted trajectory differs')
  comparison=cmp.identity(a[0],b[0],a[1]['steps_completed']);comparison['fixture']=name;rec['comparisons'].append(comparison)
  if not comparison['bit_identical']:raise RuntimeError('accepted physical/ledger identity differs')
 rec['inputs_unchanged']=m.own_sources()==source and m.sealed_pins()==pins
 rec['passed']=rec['inputs_unchanged']
 if not rec['passed']:raise RuntimeError('inputs changed')
finally:(D/'comparison.json').write_text(json.dumps(rec,indent=2)+'\n')
print(json.dumps({'passed':True,'comparisons':len(rec['comparisons']),'speedup_qualified':False}))
