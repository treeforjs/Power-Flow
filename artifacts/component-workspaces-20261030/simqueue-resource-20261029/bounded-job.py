from pathlib import Path
import argparse,json,os,subprocess,time,signal,re,sys
os.umask(0o077)
ROOT=Path('/home/treefor/projects/runtimes/simqueue-resource-20261029')
PYTHON='/home/treefor/projects/runtimes/simqueue-0.5.0-20261029/python-frozen'
p=argparse.ArgumentParser();p.add_argument('--validation',action='store_true');p.add_argument('--memory-bytes',type=int);p.add_argument('--swap-bytes',type=int);p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args()
command=a.command[1:] if a.command[:1]==['--'] else a.command
memory=18*1024**3 if a.memory_bytes is None else a.memory_bytes
swap=1024**3 if a.swap_bytes is None else a.swap_bytes
if not command or not 32*1024**2<=memory<=18*1024**3 or not 0<=swap<=1024**3:raise RuntimeError('Invalid bounded job request')
if not a.validation and (memory!=18*1024**3 or swap!=1024**3):raise RuntimeError('Production job must retain18GiBRAM/1GiBswap')
job=os.environ.get('SIMQ_JOB_ID','')
token=os.environ.get('SIMQ_LAUNCH_TOKEN','')
if not job.isdigit() or not re.fullmatch('[0-9a-fA-F-]{16,128}',token):raise RuntimeError('A recorded queue allocation is required')
unit='simq-job-'+job+'-'+re.sub('[^0-9a-fA-F]','',token)[:12]+'.scope'
records=Path('/home/treefor/.local/state/simqueue/resource-scopes');records.mkdir(mode=0o700,parents=True,exist_ok=True)
entry=records/(job+'-entry.json');result=records/(job+'.json')
record={'job_id':int(job),'unit':unit,'started_unix_s':time.time(),'controller_cgroup':Path('/proc/self/cgroup').read_text().strip(),'memory_max_bytes':memory,'swap_max_bytes':swap,'validation':a.validation,'command':command}
def save():
 t=result.with_suffix('.json.tmp');t.write_text(json.dumps(record,indent=2)+'\n');t.replace(result)
def properties():
 proc=subprocess.run(['systemctl','--user','show',unit,'-p','LoadState','-p','ControlGroup','-p','ActiveState','-p','SubState','-p','Result','-p','MemoryPeak','-p','MemorySwapPeak','-p','MemoryMax','-p','MemorySwapMax','-p','OOMPolicy'],capture_output=True,text=True)
 return dict(line.split('=',1) for line in proc.stdout.splitlines() if '=' in line)
def sample(cg):
 return {key:(cg/key).read_text().strip() for key in ['memory.current','memory.peak','memory.swap.current','memory.swap.peak','memory.events','cgroup.events','cgroup.procs'] if (cg/key).exists()}
args=['systemd-run','--user','--scope','--unit='+unit,'--property=MemoryMax='+str(memory),'--property=MemorySwapMax='+str(swap),'--property=OOMPolicy=kill',PYTHON,str(ROOT/'scope-entry.py'),'--memory-bytes',str(memory),'--swap-bytes',str(swap),'--unit',unit,'--receipt',str(entry),'--',*command]
child=None;cancelled=False;cg=None
def stop(sig,frame):
 global cancelled
 if cancelled:return
 cancelled=True
 subprocess.run(['systemctl','--user','stop',unit],check=False,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)
save()
try:
 child=subprocess.Popen(args)
 while child.poll() is None:
  if cg is None:
   props=properties()
   if props.get('ControlGroup'):cg=Path('/sys/fs/cgroup')/props['ControlGroup'].lstrip('/')
  if cg is not None:record['last_cgroup_sample']=sample(cg)
  time.sleep(.1)
 code=child.returncode
 if entry.exists():
  verified=json.loads(entry.read_text());record['entry_verification']=verified;cg=Path('/sys/fs/cgroup')/verified['cgroup'].lstrip('/')
  if verified['cgroup'] in record['controller_cgroup']:raise RuntimeError('Queue resource controller entered capped job scope')
 record['scope_properties']=properties()
 if cg is not None and cg.exists():record['last_cgroup_sample']=sample(cg)
 record['finished_unix_s']=time.time();record['process_returncode']=code;record['cancelled']=cancelled
 record['scope_empty']=cg is not None and (not cg.exists() or not (cg/'cgroup.procs').read_text().strip())
 record['memory_policy_verified']=entry.exists()
 record['oom_killed']=record['scope_properties'].get('Result')=='oom-kill' or 'oom_kill 0' not in record.get('last_cgroup_sample',{}).get('memory.events','oom_kill 0')
 if not entry.exists() and code==0:raise RuntimeError('Scoped command did not verify resource policy')
 if not record['scope_empty']:
  subprocess.run(['systemctl','--user','stop',unit],check=False)
  record['scope_empty']=cg is not None and (not cg.exists() or not (cg/'cgroup.procs').read_text().strip())
 if not record['scope_empty']:raise RuntimeError('Job scope cleanup is incomplete')
 save()
 if code<0:code=128-code
 sys.exit(code)
except Exception as exc:
 record['error']=str(exc);record['finished_unix_s']=time.time()
 subprocess.run(['systemctl','--user','stop',unit],check=False,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
 save();raise
