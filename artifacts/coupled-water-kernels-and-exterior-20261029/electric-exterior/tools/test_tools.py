#!/usr/bin/env python3
"""Bounded CPU-only contract/corruption tests. Never invokes dispatch.main."""
import json,math,os,pathlib,shutil,signal,struct,subprocess,sys,tempfile,time
from unittest import mock
import compare_fields as compare
import dispatch
from static_contract import O,F,FIELDS,PARTICLES,read_json,write_json,uniform_audit,sha

def rejected(fn):
 try:fn()
 except (AssertionError,FileNotFoundError,ValueError):return True
 raise AssertionError('Invalid input accepted')

def main():
 assert os.sched_getaffinity(0)=={7} and os.getpriority(os.PRIO_PROCESS,0)>=15
 tests=[]
 with tempfile.TemporaryDirectory(prefix='cpu-contract-',dir=O) as td:
  t=pathlib.Path(td);ps=t/'particles.bin';a=t/'a.bin';b=t/'b.bin'
  positions=[(.004,.004,0),(.004,.0073,0),(.0019,.0071,0),(.008,.008,0)]
  old=[(3,4,12,0,0,2,1,.5),(0,2,0,0,1,0,1,.5),(1,0,0,1,0,0,1,.5),(1,0,0,0,0,0,1,.5)]
  new=[(6,8,12,0,0,2,1,.5),(0,3,0,0,1,0,1,.5),(2,0,0,1,0,0,1,.5),(1,0,0,0,0,0,1,.5)]
  ps.write_bytes(b''.join(PARTICLES.pack(i,*r,0,0,0,-1.602176634e-19,9.1093837e-31,i+1) for i,r in enumerate(positions)))
  a.write_bytes(b''.join(FIELDS.pack(*v) for v in old));b.write_bytes(b''.join(FIELDS.pack(*v) for v in new))
  result=compare.compare(ps,a,b,4);all_=result['regions']['all']
  assert math.isclose(all_['weighted_RMS_delta_E_V_m'],math.sqrt(3),rel_tol=1e-15)
  assert math.isclose(all_['weighted_RMS_delta_Eparallel_V_m'],math.sqrt(5/6),rel_tol=1e-15)
  assert math.isclose(all_['weighted_RMS_delta_Eperpendicular_V_m'],math.sqrt(25/6),rel_tol=1e-15)
  assert all_['zero_B_macroparticles']==1 and result['regions']['finite_width_gap']['macroparticles']==1 and result['regions']['above_upper_back']['macroparticles']==1 and result['regions']['foil_corner_vacuum_250um']['macroparticles']==1
  tests.append('Manufactured weighted vector/parallel/perpendicular statistics and geometric masks')
  original=b.read_bytes();b.write_bytes(original[:-1]);rejected(lambda:compare.compare(ps,a,b,4));b.write_bytes(original)
  values=list(new[0]);values[6]=float('nan');b.write_bytes(FIELDS.pack(*values)+original[FIELDS.size:]);rejected(lambda:compare.compare(ps,a,b,4))
  values=list(new[0]);values[3]=-0.;b.write_bytes(FIELDS.pack(*values)+original[FIELDS.size:]);rejected(lambda:compare.compare(ps,a,b,4));b.write_bytes(original)
  rejected(lambda:compare.compare(ps,a,t/'missing.bin',4));tests.append('Truncated/missing/nonfinite field records and signed-zero B changes fail closed')
  # Exercise the independent uniform CSR/nodal/work reader on pre-existing
  # accepted fresh nodes, using a clearly labeled CPU ledger surrogate only.
  surrogate=t/'uniform-reader-surrogate';surrogate.mkdir();accepted=O.parent/'cathode-emission-audit-20261029/static-domain/9mm'
  for name in ['result.json','control-gates.json','resources.json','phi_nodes.bin','psi_nodes.bin','particle_fields.bin']:shutil.copy2(accepted/name,surrogate/name)
  shutil.copy2(O/'uniform-control-sealed/work-and-layout.json',surrogate/'work-and-layout.json')
  checked=uniform_audit(surrogate,3);assert checked['passed_additive_control'] and not checked['original_strict_passed']
  rejected(lambda:uniform_audit(surrogate,0))
  tests.append('Existing accepted fresh-node surrogate: independent CSR residuals, source work, nodal checks; strict failure retained')
  w=read_json(surrogate/'work-and-layout.json');w['work_checks_passed']=False;write_json(surrogate/'work-and-layout.json',w);rejected(lambda:uniform_audit(surrogate,3));tests.append('Failed work flag rejected; exit3 alone cannot qualify uniform control')
  request=read_json(O/'gpu-request.json');mapping=dispatch.cases_from_request(request);assert len(mapping)==13
  rejected(lambda:dispatch.approve_followup(t/'missing-approval.json',['nested-105-40'],request));tests.append('Every requested case has explicit argv; missing follow-up approval rejected')
  # PIDFD test child only sleeps on CPU. It is a new owned child, never PIC.
  cmd=[sys.executable,'-c','import time;time.sleep(10)',str(t/'unique-owned-output')]
  child=subprocess.Popen(cmd,start_new_session=True);fd=os.pidfd_open(child.pid)
  try:
   deadline=time.monotonic()+2
   while time.monotonic()<deadline:
    try:
     info=dispatch.identity(child.pid,pathlib.Path(sys.executable).resolve(),cmd);break
    except (AssertionError,FileNotFoundError):time.sleep(.01)
   else:raise AssertionError('Owned CPU child identity not observed')
   assert dispatch.descendant(child.pid,os.getpid()) and info['start_ticks']
   rejected(lambda:dispatch.identity(child.pid,pathlib.Path(sys.executable).resolve(),cmd[:-1]+['wrong-output']))
   original_stat=dispatch.stat_info;calls=[0]
   def reused(pid):
    value=original_stat(pid);calls[0]+=1
    if calls[0]>1:value['start_ticks']=str(int(value['start_ticks'])+1)
    return value
   with mock.patch.object(dispatch,'stat_info',reused):rejected(lambda:dispatch.identity(child.pid,pathlib.Path(sys.executable).resolve(),cmd))
   tests.append('Exact executable/output/start identity, ancestry and PID reuse rejection on owned CPU sleep child')
  finally:
   signal.pidfd_send_signal(fd,signal.SIGTERM);child.wait(timeout=2);os.close(fd)
 write_json(O/'tools-tests.json',dict(status='passed',tests=tests,GPU_execution=False,PIC_or_Poisson_execution=False,production_actions=False,CPU=[7],nice=15,uniform_reader_test_scope='CPU-only surrogate with existing accepted fresh binaries and separately copied CPU ledger; not a newly qualified GPU receipt'))
 print(f'{len(tests)} bounded CPU tool contract tests PASS; no GPU/PIC/production action')

if __name__=='__main__':main()
