import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
ROOT=Path(__file__).resolve().parents[1]
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def verify(pins):
    assert {str(p.relative_to(ROOT)) for p in (ROOT/'tools').rglob('*.py')}==set(pins['source_files'])
    for section in ('source_files','dependencies'):
        for n,h in pins[section].items():assert sha(ROOT/n)==h,'input changed: '+n
    assert sha(Path(sys.executable).resolve())==pins['frozen_python_sha256']
def main():
    job=os.environ.get('SIMQ_JOB_ID');assert job and job.isdigit()
    assert len(os.sched_getaffinity(0))==1
    assert not os.environ.get('PF_BUILD_ISA_FLAGS') and not os.environ.get('CUDA_VISIBLE_DEVICES')
    pins=json.loads((ROOT/'INPUTS.json').read_text());verify(pins)
    out=ROOT/'evidence'/('job-'+job);out.mkdir(exist_ok=False)
    start=time.monotonic();reports=[]
    for i in range(2):
        p=subprocess.run([sys.executable,'-I',str(ROOT/'tools/assess.py')],cwd=ROOT,
                         stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=30)
        (out/f'run-{i+1}.stderr').write_bytes(p.stderr)
        assert p.returncode==0,p.stderr.decode()
        (out/f'report-{i+1}.json').write_bytes(p.stdout);reports.append(p.stdout)
    assert reports[0]==reports[1],'fresh process output differs'
    report=json.loads(reports[0]);assert len(report['controls'])==8
    verify(pins)
    result={'schema':'water-fragment-energy-audit-receipt-v1','simq_job_id':job,
            'allocated_cpu':sorted(os.sched_getaffinity(0)),'memory_claim':os.environ.get('SIMQ_MEM'),
            'elapsed_seconds':time.monotonic()-start,'frozen_python_sha256':pins['frozen_python_sha256'],
            'source_identity':pins['source_identity'],'inputs_sha256':sha(ROOT/'INPUTS.json'),
            'source_controls_passed':8,'fresh_reports_identical':True,
            'files':{p.name:sha(p) for p in sorted(out.iterdir())},
            'passed':True,'production_qualified':False,'physical_model_validated':False,
            'network_complete':False,'native_adopted':False}
    (out/'receipt.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print('PASS eight source/invariant checks and identical reports',flush=True)
if __name__=='__main__':main()
