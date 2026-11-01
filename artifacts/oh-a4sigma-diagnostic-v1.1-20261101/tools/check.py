import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time
ROOT=Path(__file__).resolve().parents[1]
CHEM=Path('/home/treefor/projects/runtimes/pic-water-chemistry-v3.1.1-20261031')
ARR=Path('/home/treefor/projects/runtimes/pic-arrival-owner-v1-20261030/include/pf_arrival')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify(pins):
    sources={str(p.relative_to(ROOT)) for d in ('src','tests','tools') for p in (ROOT/d).rglob('*') if p.is_file() and p.suffix in ('.py','.hpp','.cpp')}
    assert sources==set(pins['source_files'])
    for section in ('source_files','dependencies'):
        for n,h in pins[section].items():assert sha(ROOT/n)==h,'changed input: '+n
    assert sha(Path(sys.executable).resolve())==pins['frozen_python_sha256']
def main():
    job=os.environ.get('SIMQ_JOB_ID');assert job and job.isdigit()
    assert len(os.sched_getaffinity(0))==1
    assert not os.environ.get('PF_BUILD_ISA_FLAGS') and not os.environ.get('CUDA_VISIBLE_DEVICES')
    pins=json.loads((ROOT/'INPUTS.json').read_text());verify(pins)
    out=ROOT/'evidence'/('job-'+job);out.mkdir(exist_ok=False)
    start=time.monotonic();commands=[]
    def run(args,name,timeout=60):
        result=subprocess.run(args,cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=timeout)
        (out/(name+'.stdout')).write_bytes(result.stdout);(out/(name+'.stderr')).write_bytes(result.stderr)
        commands.append({'argv':args,'returncode':result.returncode})
        assert result.returncode==0,(name,result.stderr.decode()[-3000:])
        return result.stdout
    run([sys.executable,'-I',str(ROOT/'tools/generate.py')],'generate')
    a=(ROOT/'build/generated_channel.hpp').read_bytes()
    run([sys.executable,'-I',str(ROOT/'tools/generate.py')],'generate-repeat')
    assert a==(ROOT/'build/generated_channel.hpp').read_bytes()
    compile_args=['/usr/bin/c++','-std=c++20','-O2','-Wall','-Wextra','-Wpedantic','-Werror','-ffp-contract=off','-fno-fast-math',
                  '-I'+str(CHEM/'include/pf_chem_v3'),'-I'+str(ARR),'-I'+str(ROOT/'src'),'-I'+str(ROOT/'build'),
                  str(ROOT/'tests/fixture.cpp'),'-o',str(ROOT/'build/fixture')]
    run(compile_args,'compile',120)
    first=run([str(ROOT/'build/fixture')],'events')
    second=run([str(ROOT/'build/fixture')],'events-repeat')
    assert first==second,'fresh fixture differs'
    report=run([sys.executable,'-I',str(ROOT/'tools/verify.py'),str(out/'events.stdout')],'independent',120)
    parsed=json.loads(report);assert parsed['passed'] and parsed['events']==2048
    # Same immutable fixture with UB instrumentation; no scientific tolerance change.
    san=compile_args.copy();san[san.index('-O2')]='-O1'
    san[2:2]=['-fsanitize=undefined','-fno-sanitize-recover=all']
    san[-1]=str(ROOT/'build/fixture-ubsan');run(san,'compile-ubsan',120)
    assert run([str(ROOT/'build/fixture-ubsan')],'events-ubsan')==first
    verify(pins)
    (out/'model.json').write_bytes((ROOT/'build/model.json').read_bytes())
    receipt={'schema':'oh-a4sigma-diagnostic-receipt-v1','job':job,'source_identity':pins['source_identity'],
             'inputs_sha256':sha(ROOT/'INPUTS.json'),'allocated_cpu':sorted(os.sched_getaffinity(0)),
             'elapsed_seconds':time.monotonic()-start,'maximum_rss_kib':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
             'generated_header_sha256':sha(ROOT/'build/generated_channel.hpp'),
             'binary_sha256':sha(ROOT/'build/fixture'),'ubsan_binary_sha256':sha(ROOT/'build/fixture-ubsan'),
             'frozen_python_sha256':pins['frozen_python_sha256'],'commands':commands,
             'files':{p.name:sha(p) for p in sorted(out.iterdir())},'passed':True,
             'production_qualified':False,'physical_model_validated':False,'network_complete':False,'native_adopted':False}
    (out/'receipt.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n')
    print('PASS 2048 events; independent source rates, energy/momentum/fragment partition, threshold, refusal, replay and UBSan controls',flush=True)
if __name__=='__main__':main()
