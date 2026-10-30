from pathlib import Path
import os,json,hashlib,time,subprocess,signal,sys
w=Path('<workspace>')
manifest_path=w/'outputs/full60-gpu-launch.json'
m=json.loads(manifest_path.read_text());out=Path(m['output']);status=w/'outputs/full60-gpu-supervisor.json'
def write(path,value):
 p=path.with_suffix(path.suffix+'.tmp');p.write_text(json.dumps(value,indent=2)+'\n');p.replace(path)
def record(phase,**values):write(status,dict(phase=phase,updated_unix_s=time.time(),output=str(out),**values))
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def cpus(value):
 result=set()
 for item in value.strip().split(','):
  if not item:continue
  a,*b=item.split('-');result.update(range(int(a),int(b[0])+1) if b else [int(a)])
 return result
def reservation():
 group=Path('/sys/fs/cgroup/powerflow-pic-gpu-reserved')
 if (group/'cpuset.cpus.partition').read_text().strip()!='isolated':raise RuntimeError('Exclusive CPU partition is not isolated')
 if cpus((group/'cpuset.cpus.exclusive.effective').read_text())!={7,15}:raise RuntimeError('Exclusive physical core not reserved')
 if set(os.sched_getaffinity(0))!={7}:raise RuntimeError('GPU job must be pinned toCPU7')
 own=Path('/proc/self/cgroup').read_text().strip()
 if '/powerflow-pic-gpu-reserved' not in own:raise RuntimeError('GPU job is outside the exclusive cgroup')
 return dict(partition='isolated',reserved_logical_cpus=[7,15],process_affinity=list(os.sched_getaffinity(0)),own_cgroup=own,checked_unix_s=time.time())
try:
 proof=reservation()
 for path,expected in m['input_sha256'].items():
  if sha(path)!=expected:raise RuntimeError('Frozen launch input changed: '+path)
 if sha(m['executable'])!=m['executable_sha256']:raise RuntimeError('Qualified executable changed')
 # The staged initial checkpoint is immutable evidence. New checkpoints may replace
 # this file only after the verified executable has accepted the state.
 if not (out/'progress.json').exists() and sha(out/m['predecessor_checkpoint'])!=m['checkpoint_sha256']:raise RuntimeError('Copied checkpoint changed')
 write(out/'core-reservation.json',proof)
 env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',HWLOC_COMPONENTS='-gl,-opencl,-cuda,-nvml,-rsmi,-levelzero',RAIKOU_PIC_FIELD_BACKEND='cudss')
 libs=Path('/home/treefor/projects/runtimes/mfem-cuda-benchmark-20261029/cuda-libs/nvidia')
 env['LD_LIBRARY_PATH']=':'.join(str(libs/v/'lib') for v in ['cu12','cublas','cuda_runtime'])
 command=[m['executable'],m['fields'],m['thermal'],str(out),'60',m['voltage']]
 child=None
 def stop(sig,frame):
  if child is not None and child.poll() is None:child.terminate()
 signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)
 with (out/'native.log').open('a') as log:
  child=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,env=env)
  record('running_full60_gpu_exclusive_core',launcher_pid_in_frozen_namespace=child.pid,mpi_ranks=1,backend='cudss',reserved_logical_cpus=[7,15],execution_affinity=[7],resume_time_s=m['resume_time_s'],requested_end_s=60e-9)
  last_check=time.monotonic()
  while child.poll() is None:
   time.sleep(2)
   if time.monotonic()-last_check>=30:
    try:proof=reservation();write(out/'core-reservation.json',proof)
    except Exception:
     child.terminate();raise
    last_check=time.monotonic()
  code=child.returncode
 if code:raise RuntimeError(f'GPU full run exited{code}; latest atomic checkpoint is retained')
 p=json.loads((out/'progress.json').read_text())
 if p['phase']!='completed' or abs(p['accepted_time_s']-60e-9)>1e-20:raise RuntimeError('No accepted60ns completion')
 record('completed_full60_gpu_exclusive_core',accepted_time_s=p['accepted_time_s'],needs_render=True)
except Exception as exc:
 record('failed',error=str(exc));raise
