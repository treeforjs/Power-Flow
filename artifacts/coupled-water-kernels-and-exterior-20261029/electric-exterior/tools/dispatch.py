#!/usr/bin/env python3
"""Parent-invoked static GPU dispatcher. Importing/testing never launches GPU."""
import argparse,fcntl,importlib.util,json,os,pathlib,select,signal,subprocess,time
from static_contract import O,P,F,read_json,write_json,record,sha,verify_record,verify_pins,check_result,uniform_audit
W=O.parent;TOOLS=pathlib.Path('/home/treefor/projects/runtimes/pic-cuda-gather-investigation-20261029/tools')
HOST_GUARD=TOOLS/'host_guard.py';LAUNCHER=TOOLS/'run-frozen-with-host-proc'
PRODUCTION=278902
LIMITS=dict(host_RSS_bytes=4*1024**3,global_GPU_bytes=10*1024**3,GPU_increase_bytes=6*1024**3,output_bytes=64*1024**2,wall_s=120)

def stat_info(pid):
 p=pathlib.Path('/proc')/str(pid);parts=(p/'stat').read_text().rsplit(')',1)[1].split()
 return dict(pid=pid,state=parts[0],ppid=int(parts[1]),start_ticks=parts[19],rss_bytes=int(parts[21])*os.sysconf('SC_PAGE_SIZE'))
def identity(pid,exe,args):
 assert pid!=PRODUCTION,'Production may never be a candidate'
 p=pathlib.Path('/proc')/str(pid);a=stat_info(pid)
 assert os.readlink(p/'exe')==str(exe),'Candidate executable identity changed'
 assert [v.decode() for v in (p/'cmdline').read_bytes().split(b'\0') if v]==args,'Candidate exact argv/output identity changed'
 b=stat_info(pid);assert a['start_ticks']==b['start_ticks'],'PID reused during identity check'
 return b
def descendant(pid,ancestor):
 for _ in range(32):
  if pid==ancestor:return True
  if pid<=1:return False
  pid=stat_info(pid)['ppid']
 return False
def owned_processes(root_pid):
 result=[]
 for p in pathlib.Path('/proc').iterdir():
  if p.name.isdigit():
   try:
    pid=int(p.name)
    if descendant(pid,root_pid):result.append(stat_info(pid))
   except (FileNotFoundError,ProcessLookupError,PermissionError):pass
 return result
def native_info(exe,args,owner):
 result=[]
 for info in owned_processes(owner):
  try:
   if os.readlink(pathlib.Path('/proc')/str(info['pid'])/'exe')==str(exe):result.append(identity(info['pid'],exe,args))
  except (FileNotFoundError,ProcessLookupError,PermissionError):pass
 return result
def gpu_memory():
 s=subprocess.check_output(['/usr/bin/nvidia-smi','--id=0','--query-gpu=memory.used','--format=csv,noheader,nounits'],text=True,timeout=2).strip()
 assert s.isdigit(),'Unexpected GPU-memory response';return int(s)*1024**2
def output_size(paths):
 return sum(p.stat().st_size for root in paths if root.exists() for p in ([root] if root.is_file() else root.rglob('*')) if p.is_file())
def production_hold():
 spec=importlib.util.spec_from_file_location('static_production_hold',TOOLS/'production_hold.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
 saved=read_json(W/'full60-manual-pause.json')['native_identity'];assert int(saved['pid'])==PRODUCTION
 fd=os.pidfd_open(PRODUCTION);poller=select.poll();poller.register(fd,select.POLLIN)
 def check():
  assert not poller.poll(0),'Original production process exited'
  p=pathlib.Path('/proc')/str(PRODUCTION)
  result=module._check_identity(saved,(p/'stat').read_text(),(p/'status').read_text(),os.readlink(p/'exe'),(p/'cgroup').read_text(),dict(line.split() for line in (module.CGROUP/'cgroup.events').read_text().splitlines()))
  assert not poller.poll(0),'Production exited during identity check';return result
 return fd,check
def cases_from_request(request):
 result={'uniform9':request['uniform_control']['argv']}
 for stage in request['stages']:
  for entry in stage['cases']:
   argv=entry['argv'];name=pathlib.Path(argv[-3]).name;assert name not in result;result[name]=argv
 for name,argv in result.items():
  assert argv[-6:]==[str(P/'build/electric_exterior_gpu'),'gpu',str(F),str(O/'parent-static-runs'/name),argv[-2],argv[-1]]
 return result
def approve_followup(path,selected,request):
 approval=read_json(path)
 assert approval['approved_for_bounded_static_followup'] is True and approval['no_particle_evolution'] is True and approval['production_resume_authorized'] is False
 assert approval['gpu_request_sha256']==sha(O/'gpu-request.json')
 assert approval['dispatcher_bindings_sha256']==sha(O/'dispatcher-bindings.json')
 verify_record(approval['uniform_dispatch_receipt']);uniform=read_json(approval['uniform_dispatch_receipt']['path'])
 assert uniform['status']=='uniform_additive_match_original_strict_failure_preserved' and uniform['original_strict_passed'] is False
 assert uniform['cases'][0]['validation']['passed_additive_control'] is True
 assert selected and len(set(selected))==len(selected) and set(selected)<=set(approval['approved_case_ids']) and 'uniform9' not in selected
 assert set(approval['approved_case_ids'])<=set(cases_from_request(request))-{'uniform9'}
 return dict(file=record(path),approval=approval)

def run_case(name,argv,check_hold):
 output=pathlib.Path(argv[-3]);guard=O/'parent-static-guards'/name
 stdout=O/'parent-static-logs'/(name+'.stdout');stderr=stdout.with_suffix('.stderr')
 assert not output.exists() and not guard.exists() and not stdout.exists() and not stderr.exists(),'Fresh output/guard/log paths required'
 stdout.parent.mkdir(parents=True,exist_ok=True)
 case=dict(case=name,output=str(output),guard=str(guard),pins_before=verify_pins(),hold_before=check_hold(),candidate_pidfds=[],limits=LIMITS)
 base=gpu_memory();assert base<=LIMITS['global_GPU_bytes'];case.update(GPU_baseline_bytes=base,sampled_peak_GPU_bytes=base,sampled_peak_candidate_RSS_bytes=0,sampled_peak_owned_sum_RSS_bytes=0,peak_output_bytes=0)
 native_argv=argv[-6:];exe=pathlib.Path(native_argv[0]);assert native_argv[3]==str(output)
 child_argv=[str(LAUNCHER),*argv[1:]]
 cmd=['timeout','--signal=TERM','--kill-after=2s','120s','python3',str(HOST_GUARD),'--out',str(guard),'--',*child_argv]
 case['command']=cmd
 environment=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',OMP_DYNAMIC='FALSE',PYTHONNOUSERSITE='1')
 for key in ['LD_PRELOAD','PYTHONPATH','PYTHONHOME']:environment.pop(key,None)
 fds={};proc=None;owner_fd=None;started=time.monotonic();last_gpu=started;failure=None
 def stop_owned():
  for pid,(fd,tick) in fds.items():
   assert pid!=PRODUCTION
   try:signal.pidfd_send_signal(fd,signal.SIGKILL)
   except ProcessLookupError:pass
  if proc is not None and proc.poll() is None and owner_fd is not None:
   # This fd was opened on this Popen's private timeout/session leader only.
   try:signal.pidfd_send_signal(owner_fd,signal.SIGTERM)
   except ProcessLookupError:pass
   try:proc.wait(timeout=4)
   except subprocess.TimeoutExpired:
    signal.pidfd_send_signal(owner_fd,signal.SIGKILL);proc.wait(timeout=2)
 try:
  with stdout.open('x') as so,stderr.open('x') as se:
   proc=subprocess.Popen(cmd,stdout=so,stderr=se,env=environment,start_new_session=True);owner_fd=os.pidfd_open(proc.pid)
   case['owned_timeout_pid']=proc.pid;case['owned_timeout_start_ticks']=stat_info(proc.pid)['start_ticks']
   while proc.poll() is None:
    check_hold()
    for item in native_info(exe,native_argv,proc.pid):
     pid=item['pid']
     if pid not in fds:
      fd=os.pidfd_open(pid);repeated=identity(pid,exe,native_argv);assert repeated['start_ticks']==item['start_ticks'];fds[pid]=(fd,item['start_ticks']);case['candidate_pidfds'].append(dict(host_pid=pid,start_ticks=item['start_ticks'],exe=str(exe),exact_argv=native_argv))
     assert fds[pid][1]==item['start_ticks'];case['sampled_peak_candidate_RSS_bytes']=max(case['sampled_peak_candidate_RSS_bytes'],item['rss_bytes'])
    assert len(fds)<=1,'Multiple candidate identities in one case'
    tree=sum(p['rss_bytes'] for p in owned_processes(proc.pid));case['sampled_peak_owned_sum_RSS_bytes']=max(case['sampled_peak_owned_sum_RSS_bytes'],tree)
    size=output_size([output,guard,stdout,stderr]);case['peak_output_bytes']=max(case['peak_output_bytes'],size)
    elapsed=time.monotonic()-started
    if time.monotonic()-last_gpu>=.5:
     now=gpu_memory();last_gpu=time.monotonic();case['sampled_peak_GPU_bytes']=max(case['sampled_peak_GPU_bytes'],now)
    assert elapsed<LIMITS['wall_s'],'Static case120s wall cap'
    assert tree<=LIMITS['host_RSS_bytes'],'Owned candidate-tree4GiB sampled RSS cap'
    assert case['sampled_peak_GPU_bytes']<=LIMITS['global_GPU_bytes'],'Global10GiB sampled GPU cap'
    assert case['sampled_peak_GPU_bytes']-base<=LIMITS['GPU_increase_bytes'],'6GiB sampled GPU increase cap'
    assert size<=LIMITS['output_bytes'],'64MiB output cap'
    time.sleep(.02 if not fds else .1)
   code=proc.returncode;case['returncode']=code;assert len(fds)==1,'No exact candidate PIDFD/executable/output/start identity was observed; fail closed'
  case['wall_s']=time.monotonic()-started
  assert case['wall_s']<122,'Timeout cleanup exceeded bounded wall allowance'
  case['GPU_after_bytes']=gpu_memory();case['sampled_peak_GPU_bytes']=max(case['sampled_peak_GPU_bytes'],case['GPU_after_bytes'])
  assert case['sampled_peak_GPU_bytes']<=LIMITS['global_GPU_bytes'] and case['sampled_peak_GPU_bytes']-base<=LIMITS['GPU_increase_bytes']
  assert output_size([output,guard,stdout,stderr])<=LIMITS['output_bytes']
  g=read_json(guard/'host-guard.json');assert not g['production_signalled'] and g['guard_error'] is None and g['child_returncode']==code
  assert g['before']['cgroup_frozen'] and g['after']['cgroup_frozen'] and g['production_pidfd_retained_until_child_exit']
  case['host_guard']=g
  if name=='uniform9':case['validation']=uniform_audit(output,code)
  else:
   assert code==0;case['validation']=check_result(output,False)
   # Magnetic input/mapping must remain bitwise identical per particle.
   from static_contract import FIELDS
   prior=pathlib.Path(next(row['path'] for row in read_json(O/'gpu-request.json')['uniform_control']['accepted_fresh_control'] if pathlib.Path(row['path']).name=='particle_fields.bin'))
   with prior.open('rb') as a,(output/'particle_fields.bin').open('rb') as b:
    for _ in range(46726):assert a.read(FIELDS.size)[24:48]==b.read(FIELDS.size)[24:48],'Native B changed'
   case['B_bitwise_fixed']=True
  case['output_receipts']=[record(p) for p in sorted(output.iterdir()) if p.is_file()]
  return case
 except BaseException as exc:
  failure=repr(exc);case['error']=failure;stop_owned();raise
 finally:
  stop_owned()
  case['finished_unix_s']=time.time()
  try:case['pins_after']=verify_pins();case['hold_after']=check_hold()
  except BaseException as exc:
   case['final_verification_error']=repr(exc)
   if failure is None:raise
  finally:
   (O/'parent-static-case-receipts').mkdir(exist_ok=True)
   write_json(O/'parent-static-case-receipts'/(name+'.json'),case)
   for fd,_ in fds.values():os.close(fd)
   if owner_fd is not None:os.close(owner_fd)

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('mode',choices=['verify-only','uniform-only','approved-followup']);parser.add_argument('--approval',type=pathlib.Path);parser.add_argument('--case',action='append',default=[]);parser.add_argument('--receipt',type=pathlib.Path,required=True);args=parser.parse_args()
 assert os.sched_getaffinity(0)=={7} and os.getpriority(os.PRIO_PROCESS,0)>=15,'Invoke under taskset-c7/nice15'
 def interrupted(signum,frame):raise InterruptedError(f'Dispatch interrupted by signal{signum}; stopping only owned static child')
 signal.signal(signal.SIGTERM,interrupted)
 assert args.receipt.resolve().is_relative_to(O.resolve()) and not args.receipt.exists(),'Fresh receipt under owned output required'
 request=read_json(O/'gpu-request.json');mapping=cases_from_request(request)
 state=dict(mode=args.mode,status='preparing',cases=[],GPU_executed=False,production_signalled=False,production_resume_authorized=False,no_particle_evolution=True,started_unix_s=time.time(),CPU=[7],nice=os.getpriority(os.PRIO_PROCESS,0))
 pidfd=None
 try:
  if args.mode=='verify-only':
   assert not args.case and args.approval is None
   state['pins']=verify_pins();pidfd,hold=production_hold();state['hold']=hold();state['status']='read_only_pins_and_frozen_production_verified';return
  if args.mode=='uniform-only':assert not args.case and args.approval is None;selected=['uniform9']
  else:
   assert args.approval is not None;state['explicit_followup_approval']=approve_followup(args.approval,args.case,request);selected=args.case
  args.receipt.parent.mkdir(parents=True,exist_ok=True)
  with (W/'pic-profile/testing-gpu.lock').open('a') as gl,(W/'pic-profile/testing-core7.lock').open('a') as cl:
   fcntl.flock(gl,fcntl.LOCK_EX);fcntl.flock(cl,fcntl.LOCK_EX)
   state['pins_before']=verify_pins();pidfd,hold=production_hold();state['hold_before']=hold()
   for name in selected:
    assert name in mapping
    state['GPU_executed']=True;state['cases'].append(run_case(name,mapping[name],hold));write_json(args.receipt,state)
   state['hold_after']=hold();state['pins_after']=verify_pins()
  if args.mode=='uniform-only':state.update(status='uniform_additive_match_original_strict_failure_preserved',original_strict_passed=False)
  else:state['status']='explicitly_approved_static_cases_passed_pending_field_analysis'
 except BaseException as exc:state.update(status='failed_closed',error=repr(exc));raise
 finally:
  if pidfd is not None:os.close(pidfd)
  state['finished_unix_s']=time.time();args.receipt.parent.mkdir(parents=True,exist_ok=True);write_json(args.receipt,state)

if __name__=='__main__':main()
