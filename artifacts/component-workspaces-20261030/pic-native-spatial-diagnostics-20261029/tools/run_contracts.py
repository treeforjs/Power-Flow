#!/usr/bin/env python3
"""Run bounded CPU contracts only; no arena construction or simulation."""
import argparse,hashlib,json,os,pathlib,resource,subprocess,sys,time,runpy
P=pathlib.Path(__file__).resolve().parents[1]
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--out-dir',type=pathlib.Path,required=True);a=p.parse_args()
 assert not a.out_dir.exists() and os.sched_getaffinity(0)=={7} and os.getpriority(os.PRIO_PROCESS,0)>=15
 resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));a.out_dir.mkdir();records=[]
 for name,command in [('spatial-energy',[sys.executable,str(P/'tests/test_diagnostics.py')]),('rate-moments',[str(P/'build/test_rate_moments')])]:
  t=time.monotonic();r=subprocess.run(command,capture_output=True,text=True,timeout=90,check=False);log=a.out_dir/(name+'.log');log.write_text(r.stdout+r.stderr)
  records.append(dict(name=name,command=command,returncode=r.returncode,wall_s=time.monotonic()-t,log=str(log),log_sha256=sha(log)));assert r.returncode==0
 fixture=runpy.run_path(str(P/'tests/test_diagnostics.py'))['fixture']();sys.path.insert(0,str(P/'source'));from diagnostics import diagnose
 for name,value in [('analytic-input.json',fixture),('analytic-diagnostics.json',diagnose(fixture))]:(a.out_dir/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
 receipt=dict(schema='native_diagnostics_cpu_contracts/v1',status='PASS',tests=records,affinity=sorted(os.sched_getaffinity(0)),nice=os.getpriority(os.PRIO_PROCESS,0),address_space_cap_bytes=4*1024**3,per_test_timeout_s=90,analytic_species=['electron','ion'],native_view_port='compiled with immutable v5 headers; no live Arena constructed',new_evolution=False,GPU=False,full_mechanism_qualification=False,files={p.name:sha(p) for p in sorted(a.out_dir.iterdir()) if p.is_file()})
 (a.out_dir/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print('CPU diagnostic contracts PASS; no arena/PIC construction or evolution')
if __name__=='__main__':main()
