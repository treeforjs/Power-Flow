from pathlib import Path
import json,os,signal,time,subprocess,shutil,argparse,hashlib,fcntl
W=Path('/home/treefor/sim-workspace/powerflow')
OUT=W/'outputs/pic-profile';OUT.mkdir(exist_ok=True)
M=json.loads((W/'outputs/full60-gpu-launch.json').read_text())
PRODUCTION=Path(M['output'])
FROZEN='/home/treefor/projects/runtimes/arch-base-20261028/run-frozen'
LIBS=Path('/home/treefor/projects/runtimes/mfem-cuda-benchmark-20261029/cuda-libs/nvidia')
ENV=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',RAIKOU_PIC_FIELD_BACKEND='cudss',HWLOC_COMPONENTS='-gl,-opencl,-cuda,-nvml,-rsmi,-levelzero')
LD=':'.join(str(LIBS/v/'lib') for v in ['cu12','cublas','cuda_runtime'])
def inspect(pid):
 p=Path('/proc')/str(pid);stat=(p/'stat').read_text().split();return {'pid':pid,'exe':os.readlink(p/'exe'),'args':[(v.decode()) for v in (p/'cmdline').read_bytes().split(b'\0') if v],'start_ticks':stat[21],'state':stat[2]}
def manual_pause_requested():
 p=W/'outputs/full60-manual-pause.json'
 return p.exists() and json.loads(p.read_text()).get('status')=='paused_by_explicit_request'
def freeze_seed(seed=None):
 seed=Path(seed) if seed is not None else OUT/'seed-current';seed.mkdir()
 for unused in range(5):
  meta=(PRODUCTION/'restart.meta').read_bytes();name=meta.decode().split()[0]
  try:
   with (PRODUCTION/name).open('rb') as f:data=f.read()
   break
  except FileNotFoundError:continue
 else:raise RuntimeError('Could not capture stable checkpoint')
 (seed/'restart.meta').write_bytes(meta);(seed/name).write_bytes(data)
 with (PRODUCTION/'history.csv').open() as f:header=next(f)
 (seed/'history.csv').write_text(header)
 elapsed=float(meta.decode().split()[-2]);origin=M['start_time_s']
 (seed/'seed.json').write_text(json.dumps({'production':str(PRODUCTION),'start_ns':(origin+elapsed)*1e9,'checkpoint':name,'checkpoint_sha256':hashlib.sha256(data).hexdigest()},indent=2)+'\n')
 return seed
def _run(label,mode='plain',executable=None,window_ps=1,seed=None,nsys=None,cpu=6,extra_env=None):
 seed=Path(seed) if seed else OUT/'seed-current'
 if not seed.exists():seed=freeze_seed()
 report=json.loads((seed/'seed.json').read_text())
 target=OUT/label
 if target.exists():raise RuntimeError('Output already exists '+str(target))
 shutil.copytree(seed,target)
 binary=executable or M['executable'];end_ns=report['start_ns']+window_ps*.001
 run_env=dict(ENV,**(extra_env or {}))
 args=[binary,M['fields'],M['thermal'],str(target),format(end_ns,'.17g'),M['voltage']]
 if mode=='gprofng':args=['gprofng','collect','app','-p','1','-a','off','-F','off','-o',str(target/'cpu.er'),*args]
 if mode=='nsys':args=[nsys,'profile','--trace=cuda','--sample=none','--cpuctxsw=none','--force-overwrite=false','--output='+str(target/'gpu'),*args]
 if cpu not in (6,7):raise ValueError('Current testing allocation is CPU 6 or CPU 7 only')
 if cpu==7:
  allocation=json.loads((W/'outputs/pic-testing-allocation.json').read_text())
  group=Path('/sys/fs/cgroup/powerflow-pic-gpu-reserved')
  events=dict(line.split() for line in (group/'cgroup.events').read_text().splitlines())
  if allocation['status']!='paused_job_frozen_core7_released_for_testing' or events.get('frozen')!='1' or not manual_pause_requested():
   raise RuntimeError('CPU7 testing requires the preserved, frozen full-run pause')
 command=[FROZEN,'nice','-n','15','taskset','-c',str(cpu),'env','LD_LIBRARY_PATH='+LD,*args]
 if mode=='nsys':
  command=['nice','-n','15','taskset','-c',str(cpu),'/home/treefor/projects/runtimes/pic-profilers-20261029/bin/capture-sealed-pic',str(target/'gpu'),'60',M['fields'],M['thermal'],str(target),format(end_ns,'.17g'),M['voltage']]
 native=[]
 for p in Path('/proc').iterdir():
  if not p.name.isdigit():continue
  try:
   item=inspect(int(p.name))
   if item['exe']==M['executable'] and M['output'] in item['args']:native.append(item)
  except (OSError,UnicodeDecodeError):pass
 if len(native)!=1:raise RuntimeError('Expected exactly one verified production process')
 original=native[0];was_stopped=original['state'] in ('T','t');fd=os.pidfd_open(original['pid'])
 if inspect(original['pid'])['start_ticks']!=original['start_ticks']:raise RuntimeError('PID identity changed')
 before=json.loads((PRODUCTION/'progress.json').read_text())
 receipt={'label':label,'mode':mode,'binary':binary,'seed':report,'window_ps':window_ps,'command':command,'production_identity':original,'production_before':before,'profile_cpu':cpu}
 child=None;paused=False;wall=time.monotonic()
 try:
  signal.pidfd_send_signal(fd,signal.SIGSTOP);paused=True;time.sleep(.2)
  with (target/'profiler.log').open('w') as log:
   child=subprocess.Popen(command,env=run_env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
   try:code=child.wait(timeout=180)
   except BaseException:
    os.killpg(child.pid,signal.SIGKILL);child.wait();raise
  receipt['returncode']=code
  receipt['wall_s']=time.monotonic()-wall
  if code:raise RuntimeError('Replay/profiler failed; see '+str(target/'profiler.log'))
 finally:
  if paused:
   if was_stopped or manual_pause_requested():
    receipt['production_left_stopped']=True
    receipt['production_pause_reason']='Preserve pre-existing stopped state or explicit manual pause request'
   else:
    signal.pidfd_send_signal(fd,signal.SIGCONT);receipt['production_resumed_unix_s']=time.time()
  os.close(fd)
  (target/'run-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
 return target
def run(*args,**kwargs):
 cpu=kwargs.get('cpu',args[6] if len(args)>6 else 6)
 with (OUT/'testing-gpu.lock').open('a') as gpu_lease:
  fcntl.flock(gpu_lease,fcntl.LOCK_EX)
  with (OUT/('testing-core7.lock' if cpu==7 else 'testing-core.lock')).open('a') as lease:
   fcntl.flock(lease,fcntl.LOCK_EX)
   return _run(*args,**kwargs)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('label');p.add_argument('--mode',choices=['plain','gprofng','nsys'],default='plain');p.add_argument('--binary');p.add_argument('--window-ps',type=float,default=1);p.add_argument('--seed');p.add_argument('--nsys');p.add_argument('--cpu',type=int,default=6);a=p.parse_args()
 print(run(a.label,a.mode,a.binary,a.window_ps,a.seed,a.nsys,a.cpu))
