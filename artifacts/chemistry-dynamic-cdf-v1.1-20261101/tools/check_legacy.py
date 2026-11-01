import hashlib,json,os,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE=Path('/home/treefor/projects/runtimes/pic-water-chemistry-v3.1.1-20261031')
RNG=Path('/home/treefor/projects/runtimes/chemistry-rng-layout-v1-20261101')
ARR=Path('/home/treefor/projects/runtimes/pic-arrival-owner-v1-20261030/include/pf_arrival')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    j=json.loads((ROOT/'LEGACY-INPUTS.json').read_text())
    def verify():
        for p,h in j['files'].items():assert sha(ROOT/p)==h,p
        assert sha(Path(sys.executable).resolve())==j['frozen_python_sha256']
    job=os.environ.get('SIMQ_JOB_ID');assert job and job.isdigit() and len(os.sched_getaffinity(0))==1 and not os.environ.get('PF_BUILD_ISA_FLAGS') and not os.environ.get('CUDA_VISIBLE_DEVICES');verify()
    out=ROOT/'evidence'/('legacy-job-'+job);out.mkdir(exist_ok=False)
    def run(args,label):
        p=subprocess.run(args,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=90);(out/(label+'.stdout')).write_bytes(p.stdout);(out/(label+'.stderr')).write_bytes(p.stderr);assert p.returncode==0,(label,p.stderr.decode());return p.stdout
    outputs=[]
    for name,inc in [('RNG2',RNG/'include/pf_chem_v3'),('CDF-extension',ROOT/'include/pf_chem_v3')]:
        binary=ROOT/'build'/('legacy-all-'+name);run(['/usr/bin/c++','-std=c++20','-O2','-Wall','-Wextra','-Werror','-ffp-contract=off','-fno-fast-math','-I'+str(inc),'-I'+str(ARR),str(ROOT/'tests/legacy_kinds.cpp'),'-lcrypto','-o',str(binary)],'compile-'+name)
        outputs.append(run([str(binary),str(BASE/'data/tables/chemistry-v3.1.blob'),sha(BASE/'data/tables/chemistry-v3.1.blob')],name))
    assert outputs[0]==outputs[1];verify()
    receipt={'schema':'dynamic-cdf-legacy-kinds-v1','job':job,'allocated_cpu':sorted(os.sched_getaffinity(0)),'numerical_source_identity':j['numerical_source_identity'],'inputs_sha256':sha(ROOT/'LEGACY-INPUTS.json'),'summary':json.loads(outputs[0]),'files_sha256':{p.name:sha(p) for p in out.iterdir()},'passed':True,'production_qualified':False,'physical_model_validated':False,'network_complete':False,'native_adopted':False};(out/'receipt.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n');print('PASS all six legacy process kinds: byte-identical records and ledgers with same RNG2')
if __name__=='__main__':main()
