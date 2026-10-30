#!/usr/bin/env python3
import argparse,hashlib,json,os,pathlib,resource,subprocess,time
P=pathlib.Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--input',type=pathlib.Path,required=True);p.add_argument('--out',type=pathlib.Path,required=True);a=p.parse_args();assert not a.out.exists() and os.sched_getaffinity(0)=={7} and os.getpriority(os.PRIO_PROCESS,0)>=15;resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));a.out.mkdir()
 input_receipt=json.loads((a.input/'input-receipt.json').read_text());assert all(sha(a.input/name)==h for name,h in input_receipt['files'].items());origin=pathlib.Path('/home/treefor/projects/runtimes/pic-water-native-subflow-20261029/inputs/collision_scheduler.hpp');assert sha(origin)==sha(P/'inputs/collision_scheduler.hpp')
 commands=[('build',['c++','-std=c++17','-O2','-ffp-contract=off','-Wall','-Wextra','-I'+str(P/'inputs'),str(P/'test_sparse_candidates.cpp'),'-o',str(P/'build/test_sparse_candidates')]),('parity',[str(P/'build/test_sparse_candidates'),str(a.input)])];records=[]
 for name,command in commands:
  start=time.monotonic();run=subprocess.run(command,capture_output=True,text=True,timeout=90,check=False);log=a.out/(name+'.log');log.write_text(run.stdout+run.stderr);records.append(dict(name=name,command=command,returncode=run.returncode,wall_s=time.monotonic()-start,log_sha256=sha(log)));assert run.returncode==0,run.stderr
 rows=[json.loads(line) for line in (a.out/'parity.log').read_text().splitlines()];assert all(r.get('status','PASS')=='PASS' for r in rows);assert all(sha(a.input/name)==h for name,h in input_receipt['files'].items()) and sha(origin)==sha(P/'inputs/collision_scheduler.hpp')
 result=dict(schema='occupied_cell_candidate_parity/v1',status='PASS',cases=rows,input=str(a.input),input_receipt_sha256=sha(a.input/'input-receipt.json'),header_sha256=sha(origin),executable_sha256=sha(P/'build/test_sparse_candidates'),source_sha256=sha(P/'test_sparse_candidates.cpp'),source_files={p.name:sha(p) for p in [P/'microcase.py',P/'prepare_sparse_case.py',P/'run_sparse_test.py']},records=records,physical_rate_claim=False,event_selected_or_committed=False,evolution=False,GPU=False,execution=dict(affinity=sorted(os.sched_getaffinity(0)),nice=os.getpriority(os.PRIO_PROCESS,0),address_space_cap_bytes=4*1024**3,per_command_timeout_s=90,one_compile_job=True))
 (a.out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
