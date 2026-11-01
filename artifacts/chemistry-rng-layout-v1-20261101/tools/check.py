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
  p=subprocess.run(args,cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=100)
  (out/(label+'.stdout')).write_bytes(p.stdout);(out/(label+'.stderr')).write_bytes(p.stderr);assert p.returncode==0,p.stderr.decode();return p.stdout
 def compile(source,binary,include,extra=()):
  run(['/usr/bin/c++','-std=c++20','-O2','-Wall','-Wextra','-Wpedantic','-Werror','-ffp-contract=off','-fno-fast-math','-I'+str(include),'-I'+str(ARR),str(source),'-lcrypto',*extra,'-o',str(ROOT/'build'/binary)],'compile-'+binary)
 source=ROOT/'tests/late_rejection.cpp';compile(source,'baseline',BASE/'include/pf_chem_v3');compile(source,'corrected',ROOT/'include/pf_chem_v3')
 raw={}
 for name in ('baseline','corrected'):
  raw[name]=run([str(ROOT/'build'/name)],name);assert raw[name]==run([str(ROOT/'build'/name)],name+'-repeat')
 old,new=map(json.loads,[raw['baseline'],raw['corrected']]);assert old['first_event_target_identity']==new['first_event_target_identity'] and old['ordinary_first_word']==8 and new['ordinary_first_word']==32
 compile(BASE/'tests/test_audit_regressions.cpp','audit-regressions',ROOT/'include/pf_chem_v3')
 run([str(ROOT/'build/audit-regressions'),str(BASE/'data/tables/chemistry-v3.1.blob')],'audit-regressions')
 compile(source,'corrected-ubsan',ROOT/'include/pf_chem_v3',['-O1','-fsanitize=undefined','-fno-sanitize-recover=all'])
 assert run([str(ROOT/'build/corrected-ubsan')],'corrected-ubsan')==raw['corrected']
 verify(pins)
 report={'schema':'chemistry-rng-layout-v1','source_identity':pins['source_identity'],'base_blob_sha256':pins['dependencies'][str(BASE/'data/tables/chemistry-v3.1.blob')],'old_overlap_reproduced':True,'actual_first_event_target_identity_equal':True,'old_fixed_seed_trajectory_parity_claimed':False,'layout_version':2,'half_open_ranges':{'candidate_time':[0,1],'target_species':[1,2],'channel':[2,3],'target_Maxwellian_rejection':[4,20],'ordinary_kinematics':[32,96],'reserved_H2_partition':[96,352]},'H2_component_implemented':False,'baseline':old,'corrected':new,'production_qualified':False,'physical_model_validated':False,'network_complete':False,'native_adopted':False,'GPU_qualified':False,'passed':True}
 (out/'report.json').write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
 receipt={'schema':'chemistry-rng-layout-receipt-v1','job':job,'allocated_cpu':sorted(os.sched_getaffinity(0)),'elapsed_seconds':time.monotonic()-start,'source_identity':pins['source_identity'],'inputs_sha256':sha(ROOT/'INPUTS.json'),'files_sha256':{p.name:sha(p) for p in sorted(out.iterdir())},'binaries_sha256':{p.name:sha(p) for p in sorted((ROOT/'build').iterdir()) if p.is_file()},'production_qualified':False,'native_adopted':False,'GPU_qualified':False,'passed':True}
 (out/'receipt.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n');print('PASS disjoint RNG layout; original late-target/angular correlation reproduced and corrected')
if __name__=='__main__':main()
