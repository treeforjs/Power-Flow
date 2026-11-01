#!/usr/bin/env python3
from pathlib import Path
import fcntl,subprocess,json,time,os,hashlib,csv
from contextlib import contextmanager
ROOT=Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs/ss304-material-feedback-pilot-20261029')
RUNTIME=Path('/home/treefor/projects/runtimes/pic-electrode-heat-feedback-20261029')
BINARY=RUNTIME/'build-solid-material-pilot/solid-material-feedback-pilot'
FROZEN='/home/treefor/projects/runtimes/arch-base-20261028/run-frozen'
LOCK=Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-profile/testing-core7.lock')
@contextmanager
def lease():
 with LOCK.open('a')as f:
  fcntl.flock(f,fcntl.LOCK_EX)
  try:yield
  finally:fcntl.flock(f,fcntl.LOCK_UN)
records=[]
def execute(name,nx,ny,outside,gap,dt,steps=8):
 with lease():
  dest=ROOT/name;dest.mkdir(exist_ok=True);n=int((dest/'CURRENT').read_text())if(dest/'CURRENT').exists()else -1
  cmd=[FROZEN,str(BINARY),str(dest),str(nx),str(ny),str(outside),str(gap),str(dt),str(steps)]
  start=time.monotonic();p=subprocess.run(cmd,capture_output=True,text=True);elapsed=time.monotonic()-start
  log=dest/f'chunk-after-{n:06d}.log';log.write_text(p.stdout+p.stderr)
  record={'case':name,'command':cmd,'exit':p.returncode,'wall_s':elapsed,'prior_step':n,'log':str(log)};records.append(record)
  (ROOT/'chunk-commands.json').write_text(json.dumps(records,indent=2)+'\n')
  now=int((dest/'CURRENT').read_text())if(dest/'CURRENT').exists()else -1
  (ROOT/'progress.json').write_text(json.dumps({'case':name,'last_accepted_step':now,'last_exit':p.returncode,'last_log':str(log),'completed_chunks':len(records)},indent=2)+'\n')
  print(name,p.returncode,round(elapsed,3),p.stdout.splitlines()[-2:],flush=True)
  if p.returncode:raise RuntimeError('pilot failed: '+str(log))
 return dest
# Actual deterministic chunk/restart proof, same inputs/solver; not manufactured
# equivalent physics. Compare all native heat+field checkpoint bytes after8steps.
a=execute('mesh32-dt250ps',32,4,20,12,2.5e-10,4)
b=execute('restart-monolithic-eight',32,4,20,12,2.5e-10,8)
with lease():
 pa=a/'checkpoints/step_000008';pb=b/'checkpoints/step_000008';comparison={}
 for file in ['upper.bin','lower.bin','fields.bin','row.csv']:
  comparison[file]={'identical':(pa/file).read_bytes()==(pb/file).read_bytes(),'bytes':(pa/file).stat().st_size,'sha256':hashlib.sha256((pa/file).read_bytes()).hexdigest()}
 assert all(v['identical']for v in comparison.values())
 (ROOT/'restart-proof.json').write_text(json.dumps(comparison,indent=2)+'\n')
for name,nx,ny,no,ng,dt in [('mesh32-dt250ps',32,4,20,12,2.5e-10),('mesh64-dt250ps',64,8,40,24,2.5e-10),('mesh64-dt125ps',64,8,40,24,1.25e-10)]:
 while True:
  with lease():
   dest=ROOT/name;done=(dest/'STOP.json').exists()or(dest/'DONE.json').exists()
  if done:break
  execute(name,nx,ny,no,ng,dt)
  time.sleep(.05)
 with lease():
  step=int((dest/'CURRENT').read_text());rows=[]
  for i in range(step+1):rows.append((dest/'checkpoints'/f'step_{i:06d}'/'row.csv').read_text().strip())
  # Header is shared with native driver and kept verbatim in this script.
  header='step,time_s,current_A,upper_current_A,lower_current_A,predicted_loop_voltage_V,min_T_K,max_T_K,max_phase,max_liquid_fraction,thermal_energy_J,hot_headroom_J,min_node_hot_headroom_J,magnetic_energy_J,joule_increment_J,terminal_work_increment_J,BE_magnetic_damping_J,EM_residual_J,thermal_residual_J,joule_mapping_error_J,cumulative_work_J,cumulative_joule_J,cumulative_damping_J,sigma_min_S_m,sigma_max_S_m,linear_relative_residual,linear_iterations'
  (dest/'trajectory.csv').write_text(header+'\n'+'\n'.join(rows)+'\n')
  raw=list(csv.DictReader((dest/'trajectory.csv').open()));last=raw[-1]
  summary={'case':name,'mesh':{'nx':nx,'ny':ny,'outside':no,'gap':ng},'max_dt_s':dt,'accepted_steps':step,'last_accepted':last,'stop':json.loads((dest/'STOP.json').read_text())if(dest/'STOP.json').exists()else None,'full_plasma_run':False,'GPU_work':False}
  (dest/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('CASE_SUMMARY',json.dumps(summary),flush=True)
