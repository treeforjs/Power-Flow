"""Run only as a bounded, one-core frozen SimQueue job."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify(pins):
    found={str(p.relative_to(ROOT)) for d in ('src','tests','tools') for p in (ROOT/d).rglob('*')
           if p.is_file() and '__pycache__' not in p.parts}
    assert found==set(pins['source_files']),'source set differs'
    for section in ('source_files','dependencies'):
        for n,h in pins[section].items():assert sha(ROOT/n)==h,'input changed: '+n
    assert sha(Path(sys.executable).resolve())==pins['frozen_python_sha256'],'interpreter is not pinned frozen Python'

def main():
    job=os.environ.get('SIMQ_JOB_ID')
    assert job and job.isdigit(),'submit through SimQueue'
    assert len(os.sched_getaffinity(0))==1,'one allocated CPU required'
    assert not os.environ.get('PF_BUILD_ISA_FLAGS'),'unqualified build flags'
    assert not os.environ.get('CUDA_VISIBLE_DEVICES'),'CPU-only source component'
    for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
    pins=json.loads((ROOT/'INPUTS.json').read_text());verify(pins)
    out=ROOT/'evidence'/('job-'+job);out.mkdir(exist_ok=False)
    start=time.monotonic()
    p=subprocess.run([sys.executable,'-I',str(ROOT/'tests/test_oh_primary.py')],cwd=ROOT,
                     text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=45)
    (out/'tests.log').write_text(p.stdout);print(p.stdout,flush=True)
    assert p.returncode==0,'controls failed; preserved tests.log'
    reports=[]
    for i in range(2):
        b=subprocess.check_output([sys.executable,'-I',str(ROOT/'tools/report.py')],cwd=ROOT,timeout=15)
        (out/f'report-{i+1}.json').write_bytes(b);reports.append(b)
    assert reports[0]==reports[1],'fresh-process numerical report differs'
    verify(pins)
    result={'schema':'oh-primary-supplement-qualification-v1','simq_job_id':job,
            'allocated_cpu':sorted(os.sched_getaffinity(0)),
            'memory_claim':os.environ.get('SIMQ_MEM'),'elapsed_seconds':time.monotonic()-start,
            'python':sys.version,'frozen_python_sha256':pins['frozen_python_sha256'],
            'inputs_sha256':sha(ROOT/'INPUTS.json'),'source_identity':pins['source_identity'],
            'tests_returncode':p.returncode,'fresh_process_reports_identical':True,
            'files':{p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()},
            'passed':True,'production_qualified':False,'physical_model_validated':False,
            'network_complete':False,'native_adopted':False}
    (out/'receipt.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print('PASS pinned OH workbook controls and two identical fresh-process reports',flush=True)

if __name__=='__main__':main()
