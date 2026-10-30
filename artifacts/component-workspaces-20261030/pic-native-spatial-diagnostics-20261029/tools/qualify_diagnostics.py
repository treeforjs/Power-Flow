#!/usr/bin/env python3
"""Reproducible bounded read-only CPU qualification; launch under shared lease."""
import argparse,hashlib,json,os,pathlib,resource,subprocess,time
P=pathlib.Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--out-dir',type=pathlib.Path,required=True);a=p.parse_args();assert not a.out_dir.exists()
 assert os.sched_getaffinity(0)=={7} and os.getpriority(os.PRIO_PROCESS,0)>=15
 resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));a.out_dir.mkdir();records=[]
 fixture='/home/treefor/sim-workspace/powerflow/outputs/pic-profile/cuda-gather-investigation/actual-midpoint-fixture'
 commands=[('contracts',['python',str(P/'tools/run_contracts.py'),'--out-dir',str(a.out_dir/'contracts')]),('preserved-midpoint',['python',str(P/'tools/diagnose_preserved_fixture.py'),'--fixture',fixture,'--out-dir',str(a.out_dir/'preserved-midpoint')]),('plots',['/home/treefor/projects/Power-Flow/.venv/bin/python',str(P/'tools/plot_diagnostics.py'),'--diagnostics',str(a.out_dir/'preserved-midpoint/diagnostics.json'),'--out-dir',str(a.out_dir/'white-plots')])]
 for name,command in commands:
  t=time.monotonic();r=subprocess.run(command,text=True,capture_output=True,timeout=120,check=False);log=a.out_dir/(name+'.log');log.write_text(r.stdout+r.stderr);records.append(dict(name=name,command=command,returncode=r.returncode,wall_s=time.monotonic()-t,log_sha256=sha(log),cumulative_children_peak_RSS_KiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss));assert r.returncode==0;print(name+' PASS',flush=True)
 receipt=dict(schema='native_diagnostics_readonly_qualification/v1',status='PASS',records=records,affinity=sorted(os.sched_getaffinity(0)),nice=os.getpriority(os.PRIO_PROCESS,0),address_space_cap_bytes=4*1024**3,per_command_timeout_s=120,source_files={str(p.relative_to(P)):sha(p) for d in ['source','tools','tests'] for p in sorted((P/d).glob('*')) if p.is_file()},new_evolution=False,GPU=False,full_mechanism_qualification=False)
 (a.out_dir/'qualification-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
if __name__=='__main__':main()
