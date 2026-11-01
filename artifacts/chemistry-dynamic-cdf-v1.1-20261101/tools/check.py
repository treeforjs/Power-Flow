import hashlib,json,os,struct,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE=Path('/home/treefor/projects/runtimes/pic-water-chemistry-v3.1.1-20261031')
RNG=Path('/home/treefor/projects/runtimes/chemistry-rng-layout-v1-20261101')
MODEL=Path('/home/treefor/projects/runtimes/h2-b3sigma-cdf-diagnostic-v1-20261101')
ARR=Path('/home/treefor/projects/runtimes/pic-arrival-owner-v1-20261030/include/pf_arrival')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def verify(pins):
    for sec in ('source_files','dependencies'):
        for p,h in pins[sec].items():assert sha(ROOT/p)==h,p
    assert sha(Path(sys.executable).resolve())==pins['frozen_python_sha256']
def main():
    job=os.environ.get('SIMQ_JOB_ID');assert job and job.isdigit() and len(os.sched_getaffinity(0))==1
    assert not os.environ.get('PF_BUILD_ISA_FLAGS') and not os.environ.get('CUDA_VISIBLE_DEVICES')
    pins=json.loads((ROOT/'INPUTS.json').read_text());verify(pins);out=ROOT/'evidence'/('job-'+job);out.mkdir(exist_ok=False);start=time.monotonic()
    def run(args,label,good=True):
        p=subprocess.run(args,cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=100)
        (out/(label+'.stdout')).write_bytes(p.stdout);(out/(label+'.stderr')).write_bytes(p.stderr)
        assert (p.returncode==0)==good,(label,p.stderr.decode()[-2000:]);return p.stdout
    run([sys.executable,'-I',str(ROOT/'tools/generate.py')],'generate');blob=ROOT/'data/chemistry-cdf.blob';a=blob.read_bytes()
    run([sys.executable,'-I',str(ROOT/'tools/generate.py')],'generate-repeat');assert a==blob.read_bytes()
    def compile(source,name,include,extra=()):
        run(['/usr/bin/c++','-std=c++20','-O2','-Wall','-Wextra','-Wpedantic','-Werror','-ffp-contract=off','-fno-fast-math','-I'+str(include),'-I'+str(ARR),'-I'+str(MODEL/'src'),'-I'+str(MODEL/'build'),str(source),'-lcrypto',*extra,'-o',str(ROOT/'build'/name)],'compile-'+name)
    inc=ROOT/'include/pf_chem_v3';fixture=ROOT/'build/fixture';compile(ROOT/'tests/fixture.cpp','fixture',inc)
    first=run([str(fixture),str(blob),sha(blob)],'fixture');assert first==run([str(fixture),str(blob),sha(blob)],'fixture-repeat')
    compile(ROOT/'tests/fixture.cpp','fixture-ubsan',inc,['-O1','-fsanitize=undefined','-fno-sanitize-recover=all']);assert first==run([str(ROOT/'build/fixture-ubsan'),str(blob),sha(blob)],'fixture-ubsan')
    # Same seeded ordinary events before/after the ABI extension, same RNG2.
    for name,headers in [('legacy',RNG/'include/pf_chem_v3'),('extended',inc)]:
        compile(RNG/'tests/late_rejection.cpp',name,headers)
    assert run([str(ROOT/'build/legacy')],'legacy')==run([str(ROOT/'build/extended')],'extended')
    compile(BASE/'tests/test_audit_regressions.cpp','audit-regressions',inc);run([str(ROOT/'build/audit-regressions'),str(BASE/'data/tables/chemistry-v3.1.blob')],'audit-regressions')
    run([str(fixture),str(BASE/'data/tables/chemistry-v3.1.blob'),sha(BASE/'data/tables/chemistry-v3.1.blob')],'v1-loaded')
    run([str(fixture),str(blob),''],'missing-sha',False)
    run([str(fixture),str(blob),'0'*64],'wrong-sha',False)
    sections={}
    for i in range(struct.unpack_from('<I',a,12)[0]):
        n,off,size=struct.unpack_from('<24sQQ',a,16+i*40);sections[n.rstrip(b'\0').decode()]=(off,size)
    descriptor=sections['cdf_desc'][0];q=sections['cdf_quantile'][0];frag=sections['cdf_fragment_eV'][0]
    c=json.loads((ROOT/'data/contract.json').read_text());pi=sections['procs'][0]+512*c['process_index']
    changes=[('descriptor-schema',descriptor,'I',2),('model',descriptor+4,'I',2),('target-v',descriptor+20,'I',1),('RNG-layout',descriptor+24,'I',1),('RNG-overlap',descriptor+28,'I',8),('point-count',descriptor+36,'I',99),('point-overflow',descriptor+32,'I',0xffffffff),('quantile-support',descriptor+64,'d',0.),('omitted-probability',descriptor+80,'d',0.),('reference-cost',descriptor+88,'d',0.),('event-threshold',descriptor+96,'d',4.5),('model-sha',descriptor+176,'I',0),('source-sha',descriptor+112,'I',0),('rate-sha',descriptor+144,'I',0),('NaN-fragment',frag+40,'d',float('nan')),('quantile-duplicate',q+8,'d',.001),('CDF-index',pi+88,'I',1),('CDF-fixed-loss',pi+96,'d',1.),('wrong-kind',pi,'I',7),('v1-new-sections',8,'I',1)]
    for label,offset,fmt,value in changes:
        bad=bytearray(a);struct.pack_into('<'+fmt,bad,offset,value);path=ROOT/'build'/('bad-'+label+'.blob');path.write_bytes(bad)
        run([str(fixture),str(path),sha(path)],label,False)
    verify(pins);(out/'contract.json').write_bytes((ROOT/'data/contract.json').read_bytes())
    receipt={'schema':'dynamic-CDF-receipt-v1','job':job,'source_identity':pins['source_identity'],'inputs_sha256':sha(ROOT/'INPUTS.json'),'allocated_cpu':sorted(os.sched_getaffinity(0)),'elapsed_seconds':time.monotonic()-start,'structural_fault_controls':len(changes),'identity_fault_controls':2,'base_fixed_seed_RNG2_parity':True,'base_audit_regressions':11,'summary':json.loads(first),'private_generated_sha256':{p.name:sha(p) for p in (ROOT/'data').iterdir() if p.is_file()},'files_sha256':{p.name:sha(p) for p in out.iterdir() if p.is_file()},'passed':True,'GPU_qualified':False,'native_adopted':False,'production_qualified':False,'physical_model_validated':False,'network_complete':False}
    (out/'receipt.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n');print('PASS explicit CDF descriptor, actual MCC/refusals, legacy parity, tamper, replay, UBSan')
if __name__=='__main__':main()
