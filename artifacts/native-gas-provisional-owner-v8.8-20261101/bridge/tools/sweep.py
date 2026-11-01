"""Persistent allocation refusal at every measured final-seal allocation site.
Runs the real populated coupled binary; source/header/binary pins are unchanged.
"""
from pathlib import Path
import os,sys,json,hashlib,subprocess,time
if not __debug__:raise RuntimeError('optimized Python guards refused')
R=Path(__file__).resolve().parents[1];F=R.parent/'arch-base-20261028'
O=Path('/home/treefor/projects/.workspace/outputs/native-gas-common-stage-v1.2-provisional-owner-20261101')
assert len(sys.argv)==3
base=O/sys.argv[1];out=O/sys.argv[2];out.mkdir(exist_ok=False)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert os.getenv('SIMQ_JOB_ID','').isdigit() and os.getenv('SIMQ_LAUNCH_TOKEN')
aff=sorted(os.sched_getaffinity(0));assert len(aff)==1 and os.getenv('SIMQ_CPU_LIST')==str(aff[0])
assert os.getenv('SIMQ_PARTITION')=='powerflow' and os.environ['SIMQ_UNIT'] in Path('/proc/self/cgroup').read_text()
assert int(os.environ['SIMQ_MEM'])>=3*1024**3
assert all(os.getenv(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'))
assert os.getenv('CUDA_VISIBLE_DEVICES')=='' and not os.getenv('PF_BUILD_ISA_FLAGS','')
a=Path('/').stat();b=(F/'root').stat();assert (a.st_dev,a.st_ino)==(b.st_dev,b.st_ino),'frozen root required'
d=json.loads((base/'receipt.json').read_text());assert d['passed'] and d['mode']=='normal'
binary=base/'test_common';assert sha(binary)==d['binary_sha256']
pins={**d['source_inputs'],**d['compiler_headers'],str(binary):d['binary_sha256'],str(Path(__file__)):sha(__file__)}
assert all(sha(p)==v for p,v in pins.items())
count=d['results']['allocation-count']['seal_allocations'];assert isinstance(count,int) and 0<count<=10000
receipt={'passed':False,'job':os.environ['SIMQ_JOB_ID'],'affinity':aff,'frozen_root_identity':[a.st_dev,a.st_ino],
         'baseline_receipt':str(base/'receipt.json'),'baseline_receipt_sha256':sha(base/'receipt.json'),
         'binary_sha256':sha(binary),'source_pins':pins,'persistent_failure':True,'measured_allocations':count,
         'checks':[],'production_qualified':False}
start=time.monotonic()
try:
 for k in range(count):
  label='allocation-seal-'+str(k);dest=out/label
  with (out/(label+'.log')).open('w')as f:
   p=subprocess.run([str(binary),str(dest),label],stdout=f,stderr=subprocess.STDOUT,timeout=30)
  assert p.returncode==0,(k,p.returncode,(out/(label+'.log')).read_text()[-4000:])
  r=json.loads((dest/'result.json').read_text());assert r['passed'] and r['refused'] and r['retry_accepted'] and not r['common_accepted_publication'] and r['native_source']==d['native_source']
  receipt['checks'].append({'offset':k,'result_sha256':sha(dest/'result.json'),'log_sha256':sha(out/(label+'.log')),'checks':r['checks']})
  if k%128==0:print('passed offsets through',k,flush=True)
 assert all(sha(p)==v for p,v in pins.items())
 receipt['passed']=True
finally:
 receipt['wall_s']=time.monotonic()-start
 (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
