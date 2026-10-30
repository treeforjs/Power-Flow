#!/usr/bin/python
# Root helper: reserves only physical core7 (logical7,15), starts the GPU job
# as UID1000, and releases the reservation on exit. No package/boot changes.
from pathlib import Path
import os,sys,json,time,signal,pwd
CG=Path('/sys/fs/cgroup');GROUP=CG/'powerflow-pic-gpu-reserved'
SUPERVISOR='/home/treefor/projects/runtimes/mfem-gpu-pic-20261029/qualified/full60-gpu-supervisor.py'
def write(path,value):path.write_text(value+'\n')
def parse(value):
 out=set()
 for item in value.strip().split(','):
  if item:
   a,*b=item.split('-');out.update(range(int(a),int(b[0])+1) if b else [int(a)])
 return out
def verify():
 if (GROUP/'cpuset.cpus.partition').read_text().strip()!='isolated':raise RuntimeError('Kernel rejected isolated partition')
 if parse((GROUP/'cpuset.cpus.exclusive.effective').read_text())!={7,15}:raise RuntimeError('Exclusive set differs from7,15')
 if not {7,15}<=parse((CG/'cpuset.cpus.isolated').read_text()):raise RuntimeError('Kernel isolated CPU set is missing7,15')
 for name in ['user.slice','system.slice','init.scope']:
  p=CG/name/'cpuset.cpus.effective'
  if p.exists() and parse(p.read_text())&{7,15}:raise RuntimeError('Reserved CPUs remain in '+name)
def main():
 if os.geteuid()!=0:raise RuntimeError('Administrative authentication required for exclusive CPU reservation')
 assert Path('/sys/devices/system/cpu/cpu7/topology/thread_siblings_list').read_text().strip()=='7,15'
 if GROUP.exists():raise RuntimeError('Reservation group already exists; inspect it before retrying')
 # Enable the cpuset controller at the hierarchy root only. Existing processes
 # inherit their normal CPU sets except for CPUs granted to the isolated partition.
 if 'cpuset' not in (CG/'cgroup.subtree_control').read_text().split():write(CG/'cgroup.subtree_control','+cpuset')
 GROUP.mkdir()
 child=None
 try:
  write(GROUP/'cpuset.mems',(CG/'cpuset.mems.effective').read_text().strip())
  write(GROUP/'cpuset.cpus','7,15')
  if (GROUP/'cpuset.cpus.exclusive').exists():write(GROUP/'cpuset.cpus.exclusive','7,15')
  write(GROUP/'cpuset.cpus.partition','isolated')
  verify()
  account=pwd.getpwnam('treefor');assert account.pw_uid==1000
  child=os.fork()
  if child==0:
   try:
    os.setsid();write(GROUP/'cgroup.procs',str(os.getpid()))
    os.sched_setaffinity(0,{7})
    os.initgroups(account.pw_name,account.pw_gid);os.setgid(account.pw_gid);os.setuid(account.pw_uid)
    os.chdir('/home/treefor/sim-workspace/powerflow')
    frozen='/home/treefor/projects/runtimes/arch-base-20261028/run-frozen'
    os.execve(frozen,[frozen,'/usr/bin/python',SUPERVISOR],dict(os.environ,HOME=account.pw_dir))
   except Exception as exc:print('GPU child setup failed:',exc,file=sys.stderr,flush=True);os._exit(1)
  def stop(sig,frame):
   if child:
    try:os.killpg(child,signal.SIGTERM)
    except ProcessLookupError:pass
  signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)
  print(json.dumps({'status':'exclusive_partition_active','child_host_pid':child,'logical_cpus':[7,15],'execution_cpu':7}),flush=True)
  while True:
   waited,status=os.waitpid(child,0)
   if waited==child:return os.waitstatus_to_exitcode(status)
 finally:
  # Any remaining job descendants belong only to this dedicated cgroup.
  if GROUP.exists():
   members=(GROUP/'cgroup.procs').read_text().split()
   if members:
    if (GROUP/'cgroup.kill').exists():write(GROUP/'cgroup.kill','1')
    else:
     for pid in members:
      try:os.kill(int(pid),signal.SIGKILL)
      except ProcessLookupError:pass
    for unused in range(50):
     if not (GROUP/'cgroup.procs').read_text().strip():break
     time.sleep(.1)
   write(GROUP/'cpuset.cpus.partition','member')
   GROUP.rmdir()
   print('Exclusive CPU reservation released',flush=True)
if __name__=='__main__':
 try:sys.exit(main())
 except Exception as exc:print(str(exc),file=sys.stderr,flush=True);sys.exit(1)
