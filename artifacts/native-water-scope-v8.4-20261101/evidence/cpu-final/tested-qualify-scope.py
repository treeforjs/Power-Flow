#!/usr/bin/python3
import sys,os,json,subprocess,hashlib,time
from pathlib import Path
from tooling import *
if not __debug__:raise RuntimeError('optimized Python guards refused')
resource_guard(build=True)
need(bool(os.environ.get('SIMQ_JOB_ID')) and len(os.sched_getaffinity(0))==1,'one-core SimQueue allocation required')
need(int(os.environ.get('SIMQ_MEM','0'))>=3*1024**3,'memory claim required')
need(not os.environ.get('PF_BUILD_ISA_FLAGS',''),'unqualified ISA flags refused')
mode=sys.argv[1];need(mode in ('normal','sanitized'),'mode')
label=sys.argv[2];out=OUT/label;out.mkdir(parents=True,exist_ok=False)
def pins():return {str(f.relative_to(P)):sha256(f) for top in ['source','tests','forms','inputs'] for f in sorted((P/top).rglob('*')) if f.is_file()} | {n:sha256(P/n) for n in ['tooling.py','build-components.py','qualify-scope.py']}
before=pins();deps=sealed_inputs();receipt={'passed':False,'mode':mode,'job':os.environ['SIMQ_JOB_ID'],'affinity':sorted(os.sched_getaffinity(0)),'source_inputs':before,'dependencies':deps,'runs':[],'production_qualified':False,'network_complete':False,'device_census_qualified':False}
def run(cmd,name):
 env=dict(os.environ,PF_V8_CHEM='v3.1',PF_BUILD_ISA_FLAGS='',CUDA_VISIBLE_DEVICES='')
 if mode=='sanitized':env.update(ASAN_ENV)
 start=time.monotonic()
 with (out/(name+'.log')).open('w') as f:r=subprocess.run(list(map(str,cmd)),stdout=f,stderr=subprocess.STDOUT,env=env,timeout=900)
 receipt['runs'].append({'cmd':list(map(str,cmd)),'rc':r.returncode,'wall_s':time.monotonic()-start,'log':name+'.log','log_sha256':sha256(out/(name+'.log'))})
 print(name,r.returncode,flush=True)
 if r.returncode:raise RuntimeError((out/(name+'.log')).read_text()[-7000:])
try:
 run(['/usr/bin/python3',P/'build-components.py',mode,'core','hall','v84-scope','audit-target-guard'],'build')
 tag='sanitized-' if mode=='sanitized' else ''
 run([P/'build'/(tag+'test_v84_scope'),out/'controls'],'scope-controls')
 run([P/'build'/(tag+'test_audit_target_guard'),out/'target-energy-controls'],'target-energy-controls')
 need(pins()==before and sealed_inputs()==deps,'inputs changed')
 receipt['binary_sha256']={n:sha256(P/'build'/(tag+n)) for n in ['libpic_native_v8.a','test_v84_scope','test_audit_target_guard']}
 receipt['controls']=json.loads((out/'controls/result.json').read_text());need(receipt['controls']['passed'],'control result')
 receipt['result_sha256']=sha256(out/'controls/result.json');receipt['passed']=True
finally:(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
