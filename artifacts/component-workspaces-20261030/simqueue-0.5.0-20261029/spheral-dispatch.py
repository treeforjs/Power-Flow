from pathlib import Path
import json,os,subprocess,time,fcntl
os.umask(0o077)
STATE=Path('/home/treefor/.local/state/simqueue')
HANDOFF=Path('/home/treefor/projects/runs/rz-spheral-setup-20261029/launch-handoff.json')
SOURCE=Path('/home/treefor/projects/runs/rz-spheral-setup-20261029/source')
PREFIX=Path('/home/treefor/projects/runtimes/rz-spheral-20261029')
ROOT=Path('/home/treefor/projects/runs/rz-spheral-setup-20261029')
STATUS=STATE/'spheral-dispatch.json'
def write(value):
 temp=STATUS.with_suffix('.json.tmp');temp.write_text(json.dumps(value,indent=2)+'\n');temp.replace(STATUS)
def run():
 previous=json.loads(STATUS.read_text()) if STATUS.exists() else {}
 if previous.get('job_id'):return
 h=json.loads(HANDOFF.read_text())
 if h.get('readiness') is not True or h.get('native_setup_checks_passed') is not True:
  write({'status':'waiting_for_native_readiness','checked_unix_s':time.time(),'handoff':str(HANDOFF),'why_not_ready':h.get('why_not_ready',[]),'queue_core_count':4})
  return
 if h.get('schema')!='electrothermal.rz_spheral_launch_handoff.v1' or h.get('coordinator_thread')!='01a18622-368c-7142-a554-80bd598d96e6':raise RuntimeError('Unexpected handoff owner or schema')
 if h.get('physical_cores')!=4 or h.get('source')!=str(SOURCE) or h.get('runtime_prefix')!=str(PREFIX):raise RuntimeError('Handoff changed reserved runtime or four-core allocation')
 command=h.get('future_command_argv')
 if not isinstance(command,list) or not command or any(not isinstance(v,str) or '\0' in v for v in command):raise RuntimeError('Qualified structured launch command required')
 expected=[str(ROOT/'run-rz-frozen'),str(ROOT/'python/cpython-3.12.16-linux-x86_64-gnu/bin/python3.12'),str(SOURCE/'tools/run_zapdos_rz_workstation.py')]
 if command[:3]!=expected:raise RuntimeError('Command is not the qualified frozen workstation launcher')
 def option(key):
  if command.count(key)!=1:raise RuntimeError('Expected exactly one '+key)
  i=command.index(key)
  if i+1>=len(command):raise RuntimeError('Missing '+key+' value')
  return command[i+1]
 for key,value in [('--physical-cores','4'),('--mode','run'),('--prefix',str(PREFIX)),('--host-source',str(SOURCE)),('--runs-root',str(ROOT/'prepared-runs'))]:
  if option(key)!=value:raise RuntimeError('Unqualified launch option '+key)
 output=Path(option('--output'))
 if output.parent!=ROOT/'prepared-runs' or output.exists():raise RuntimeError('Fresh prepared-runs stage required')
 # This environment is stable across repeated path events and service restarts,
 # so replaying a submission after a crash retains its idempotency identity.
 environment={'HOME':'/home/treefor','USER':'treefor','LOGNAME':'treefor','LANG':'C.UTF-8','PATH':'/home/treefor/.local/bin:/usr/local/bin:/usr/bin','XDG_RUNTIME_DIR':'/run/user/1000','DBUS_SESSION_BUS_ADDRESS':'unix:path=/run/user/1000/bus','SIMQUEUE_PCORES':'0-3','SIMQUEUE_PROJECT_ROOT':'/home/treefor/projects','SIMQUEUE_STATE_DIR':str(STATE),'PYTHONDONTWRITEBYTECODE':'1'}
 key='zapdos-raikou-spheral-arch-20261029-four-core-first-run'
 queued_argv=['/home/treefor/projects/runtimes/simqueue-0.5.0-20261029/python-frozen','/home/treefor/projects/runtimes/simqueue-resource-20261029/bounded-job.py','--',*command]
 args=['/home/treefor/.local/bin/simq','submit','--json','--idempotency-key',key,'-n','zapdos-raikou-spheral','--cores','4','-C',str(ROOT),'--',*queued_argv]
 result=subprocess.run(args,env=environment,text=True,capture_output=True,check=True)
 receipt=json.loads(result.stdout)
 if receipt.get('requested_cores')!=4 or not isinstance(receipt.get('job_id'),int):raise RuntimeError('Invalid four-core queue receipt')
 write({'status':'submitted','checked_unix_s':time.time(),'job_id':receipt['job_id'],'queue_receipt':receipt,'command':command,'queued_argv':queued_argv,'memory_scope':{'RAM_bytes':19327352832,'swap_bytes':1073741824,'controller_outside_job_scope':True},'source_commit':h.get('source_commit'),'handoff':str(HANDOFF),'output':str(output)})
with (STATE/'spheral-dispatch.lock').open('w') as lock:
 fcntl.flock(lock,fcntl.LOCK_EX)
 try:run()
 except Exception as exc:
  write({'status':'dispatch_failed_closed','checked_unix_s':time.time(),'error':str(exc),'handoff':str(HANDOFF)})
  raise
