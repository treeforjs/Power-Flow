#!/usr/bin/python3
"""Build the adapter against the actual pinned native archive and gas source.
Numerical work only in a one-core SimQueue allocation and frozen userspace.
Sanitizers cover adapter/test/gas source, not the prebuilt native archive.
"""
import hashlib, importlib.util, json, os, pathlib, shutil, subprocess, sys, time
if not __debug__: raise RuntimeError('optimized Python guards refused')
R=pathlib.Path(__file__).resolve().parents[1]
RT=R.parent
N=RT/'pic-reaction-transactions-20261029/native-full-time-v8.3-admitted-export-20261101'
G=RT/'gas-target-owner-v1.4-binding-20261101'
F=RT/'arch-base-20261028'
OUT=pathlib.Path('/home/treefor/projects/.workspace/outputs/native-gas-adapter-v1.1-lifetime-20261101')
MODE=sys.argv[1];assert MODE in ('normal','sanitized')
D=OUT/sys.argv[2];D.mkdir(parents=True,exist_ok=False)
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
assert os.getenv('SIMQ_JOB_ID','').isdigit() and os.getenv('SIMQ_LAUNCH_TOKEN')
aff=sorted(os.sched_getaffinity(0));assert len(aff)==1 and os.getenv('SIMQ_CPU_LIST')==str(aff[0])
assert os.getenv('SIMQ_PARTITION')=='powerflow' and os.environ['SIMQ_UNIT'] in pathlib.Path('/proc/self/cgroup').read_text()
assert aff[0] in (6,7) or (aff[0] in (4,5) and os.getenv('SIMQ_BORROW')=='1')
assert all(os.getenv(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'))
assert os.getenv('CUDA_VISIBLE_DEVICES')==''
assert int(os.getenv('SIMQ_MEM','0'))>=3*1024**3
memory={l.split(':')[0]:int(l.split()[1])*1024 for l in pathlib.Path('/proc/meminfo').read_text().splitlines()}
assert memory['MemAvailable']>=int(os.environ['SIMQ_MEM'])+2*1024**3
assert pathlib.Path('/proc/sys/kernel/osrelease').read_text()==(F/'kernel.release').read_text()
for p in ('usr/bin/c++','usr/bin/python3'):assert sha('/'+p)==sha(F/'root'/p)
spec=importlib.util.spec_from_file_location('pinned_native_tooling',N/'tooling.py')
native=importlib.util.module_from_spec(spec);spec.loader.exec_module(native)
native.resource_guard(build=True)
source=json.loads((N/'SOURCE-SHA256.json').read_text())
assert source['source_identity']=='64b327a4023ccef60488d5f1ae54b5d723d1b1441a40dec8ab50d3b0ad7ae1db'
assert sha(N/'build/libpic_native_v8.a')=='dacf80ce52738432d376eaecdef073b8da6af78244f72d0407036971ecd7fd75'
for name,digest in source['files'].items():assert sha(N/name)==digest,name
header_pins=json.loads((R/'NATIVE-HEADERS.json').read_text())['files']
gas_pins=json.loads((R/'GAS-BINDING.json').read_text())['files']
assert len(gas_pins)==2
for name,digest in gas_pins.items():assert sha(name)==digest
assert len(header_pins)==2
for name,digest in header_pins.items():assert sha(name)==digest
def pins():
    result={str(N/name):sha(N/name) for name in source['files']}
    result.update({name:sha(name) for name in header_pins})
    result.update({name:sha(name) for name in gas_pins})
    for root,dirs in ((R,('include','tests','tools')),(G,('include','src'))):
        for d in dirs:
            for p in (root/d).rglob('*'):
                if p.is_file() and '__pycache__' not in p.parts:result[str(p)]=sha(p)
    for p in (R/'build-template.json',R/'NATIVE-HEADERS.json',R/'GAS-BINDING.json',N/'tooling.py',N/'tests/chem_v8_fixture.hpp',N/'source/src/chem_v3_native.hpp',N/'build/libpic_native_v8.a',N/'build/kinetic-boundary.o'):
        result[str(p)]=sha(p)
    return result
before=pins();deps=native.sealed_inputs()
for path,digest in before.items():
    p=pathlib.Path(path)
    if p.is_relative_to(R):
        dst=D/'source'/p.relative_to(R);dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dst)
        assert sha(dst)==digest

receipt={'passed':False,'mode':MODE,'job':os.environ['SIMQ_JOB_ID'],'affinity':aff,'inputs':before,'dependencies':deps,
 'source_identity':source['source_identity'],'runs':[],
 'allocation':{k:v for k,v in os.environ.items() if k.startswith('SIMQ_') and k!='SIMQ_LAUNCH_TOKEN'},
 'sanitizer_scope':'adapter, fixture and gas source only; native archive is prebuilt' if MODE=='sanitized' else None}
def run(cmd,name,env=None):
    cmd=list(map(str,cmd));start=time.monotonic()
    with (D/(name+'.log')).open('w') as f:p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=env,timeout=240)
    receipt['runs'].append({'command':cmd,'rc':p.returncode,'wall_s':time.monotonic()-start,'log':name+'.log','log_sha256':sha(D/(name+'.log'))})
    print(name,p.returncode,flush=True)
    if p.returncode:raise RuntimeError((D/(name+'.log')).read_text()[-5000:])
try:
    original=next(x['command'] for x in json.loads((R/'build-template.json').read_text()) if x['log']=='v83-export-compile.log')
    cc=original[original.index('/usr/bin/mpic++'):original.index('-c')]
    cc+=['-I'+str(R/'include'),'-I'+str(G/'include'),'-I'+str(N/'forms/include/pf_v7')]
    san=['-fsanitize=address,undefined','-fno-omit-frame-pointer','-fno-sanitize-recover=all'] if MODE=='sanitized' else []
    run(cc+san+['-MMD','-MF',D/'adapter.d','-c',R/'tests/test_adapter.cpp','-o',D/'adapter.o'],'compile-adapter')
    run(['/usr/bin/c++','-std=c++20','-O2','-Wall','-Wextra','-Werror','-ffp-contract=off','-fno-fast-math',
      '-I'+str(G/'include'),'-I'+str(N/'forms/include/pf_v7'),*san,'-MMD','-MF',D/'gas.d','-c',G/'src/gas_owner.cpp','-o',D/'gas.o'],'compile-gas')
    actual_headers={}
    for dep in (D/'adapter.d',D/'gas.d'):
        for name in dep.read_text().replace('\\\n',' ').split(':',1)[1].split():
            actual_headers[name]=sha(name)
    receipt['compiler_headers']=actual_headers
    run(['/usr/bin/mpic++',*san,'-Wl,--gc-sections',D/'adapter.o',D/'gas.o',N/'build/kinetic-boundary.o',N/'build/libpic_native_v8.a',*native.LIBS,'-o',D/'test_adapter'],'link')
    env=dict(os.environ,PF_V8_CHEM='v3.1',ASAN_OPTIONS='detect_leaks=1:abort_on_error=1:hard_rss_limit_mb=3072',UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1')
    binary=sha(D/'test_adapter');run([D/'test_adapter',D/'result.json'],'controls',env)
    assert sha(D/'test_adapter')==binary and pins()==before and native.sealed_inputs()==deps
    for name,digest in actual_headers.items():assert sha(name)==digest,name
    receipt['binary_sha256']=binary;receipt['result_sha256']=sha(D/'result.json');receipt['passed']=True
finally:
    (D/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
