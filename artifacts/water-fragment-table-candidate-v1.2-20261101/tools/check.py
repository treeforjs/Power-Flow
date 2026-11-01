import hashlib,json,os,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE=Path('/home/treefor/projects/runtimes/pic-water-chemistry-v3.1.1-20261031')
ARR=Path('/home/treefor/projects/runtimes/pic-arrival-owner-v1-20261030/include/pf_arrival')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def verify(pins):
 for section in ('source_files','dependencies'):
  for p,h in pins[section].items():assert sha(ROOT/p)==h,p
 assert sha(Path(sys.executable).resolve())==pins['frozen_python_sha256']
def main():
 job=os.environ.get('SIMQ_JOB_ID');assert job and job.isdigit() and len(os.sched_getaffinity(0))==1
 assert not os.environ.get('PF_BUILD_ISA_FLAGS') and not os.environ.get('CUDA_VISIBLE_DEVICES')
 pins=json.loads((ROOT/'INPUTS.json').read_text());verify(pins);out=ROOT/'evidence'/('job-'+job);out.mkdir(exist_ok=False);start=time.monotonic()
 def run(args,label):
  p=subprocess.run(args,cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=150)
  (out/(label+'.stdout')).write_bytes(p.stdout);(out/(label+'.stderr')).write_bytes(p.stderr);assert p.returncode==0,p.stderr.decode();return p.stdout
 run([sys.executable,'-I',str(ROOT/'tools/generate.py')],'generate')
 first={p.name:sha(p) for p in (ROOT/'data').iterdir() if p.is_file()}
 run([sys.executable,'-I',str(ROOT/'tools/generate.py')],'generate-repeat')
 assert first=={p.name:sha(p) for p in (ROOT/'data').iterdir() if p.is_file()},'nonreproducible candidate'
 contract=json.loads((ROOT/'data/contract.json').read_text());pol=json.loads((ROOT/'data/candidate-coverage-policy.json').read_text());man=json.loads((ROOT/'data/candidate-manifest.json').read_text());base=json.loads((BASE/'data/tables/chemistry-v3.1-manifest.json').read_text())
 assert pol['model_mode']=='strict-production' and pol['admitted_unbounded_regions']==[]
 assert len(man['groups'])==len(base['groups'])+1 and len(man['species'])==22
 assert 'e+H2' not in {g['name'] for g in man['groups']}
 for name in ('e+OH',):
  regs=[r for r in man['region_info'] if r['group']==name];assert regs and all(not r['bounded'] for r in regs)
  assert {r['index'] for r in regs}<={r['index'] for r in pol['refusing_unbounded_regions']}
 assert contract['H2']['source_window_eV']==[4.5,150.0] and contract['H2']['rate_bank_only'] and not contract['H2']['event_implementation_qualified_component_only'] and contract['H2']['native_group_absent']
 assert contract['OH']['source_window_eV']==[7.3605,10.0]
 run(['/usr/bin/c++','-std=c++20','-O2','-Wall','-Wextra','-Wpedantic','-Werror','-ffp-contract=off','-fno-fast-math','-I'+str(BASE/'include/pf_chem_v3'),'-I'+str(ARR),str(ROOT/'tests/fixture.cpp'),'-lcrypto','-o',str(ROOT/'build/fixture')],'compile')
 command=[str(ROOT/'build/fixture'),str(ROOT/'data/candidate.blob'),contract['candidate_blob_sha256']]
 result=run(command,'fixture');assert result==run(command,'fixture-repeat')
 p=subprocess.run([str(ROOT/'build/fixture'),str(ROOT/'data/candidate.blob'),'0'*64],cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=10);assert p.returncode!=0
 (out/'wrong-sha.stderr').write_bytes(p.stderr)
 # No restricted cross-section arrays or full generated blob/check tables in public evidence.
 (out/'contract.json').write_text(json.dumps(contract,sort_keys=True,indent=2)+'\n');verify(pins)
 receipt={'schema':'water-fragment-table-candidate-receipt-v1.1','job':job,'allocated_cpu':sorted(os.sched_getaffinity(0)),'source_identity':pins['source_identity'],'inputs_sha256':sha(ROOT/'INPUTS.json'),'binary_sha256':sha(ROOT/'build/fixture'),'elapsed_seconds':time.monotonic()-start,'generated_private_files_sha256':first,'files_sha256':{p.name:sha(p) for p in sorted(out.iterdir())},'passed':True,'production_qualified':False,'physical_model_validated':False,'network_complete':False,'native_adopted':False}
 (out/'receipt.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n');print('PASS candidate table generation, incomplete-group refusal support and OH blob-event controls; H2 differential admission remains closed')
if __name__=='__main__':main()
