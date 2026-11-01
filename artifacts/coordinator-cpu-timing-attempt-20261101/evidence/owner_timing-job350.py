#!/usr/bin/python -I
"""Single owner-core timing attempt, cooperating with the existing study SMT locks.

No physical-core lease is acquired here: SimQueue already owns that reservation.
A bounded sibling-lock wait may fail; it never signals another job or retries the
numerical benchmark. The comparison still checks actual SMT load and memory.
"""
import fcntl,hashlib,importlib.util,json,os,sqlite3,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
spec=importlib.util.spec_from_file_location('timing_runner',ROOT/'run_v8.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
LOCKS=Path('/home/treefor/projects/.workspace/outputs/cpu-isa-study-20261031/locks')

def bounded_lock(handle,kind,seconds=30):
    deadline=time.monotonic()+seconds
    while True:
        try:fcntl.flock(handle,kind|fcntl.LOCK_NB);return True
        except BlockingIOError:
            if time.monotonic()>=deadline:return False
            time.sleep(.2)

def owner_admission(kind):
    if not os.getenv('SIMQ_JOB_ID'):raise RuntimeError('scheduled owner job required')
    admission=m.admission(kind)
    cpu=admission['cpu_affinity'][0]
    if cpu not in (6,7) or os.getenv('SIMQ_BORROW')=='1':raise RuntimeError('owner Power-Flow core required')
    return admission,cpu

def check_job(job,cpu):
    if job['id']!=int(os.environ['SIMQ_JOB_ID']) or job['state']!='running' or job['partition']!='powerflow' or job['borrow'] or not job['exclusive_timing'] or job['unit_name']!=os.environ['SIMQ_UNIT'] or job['cpu_list']!=str(cpu) or job['mem_bytes']!=int(os.environ['SIMQ_MEM']):
        raise RuntimeError('live exclusive owner allocation not verified')

def host_main(label,tag):
    if m.inside_frozen():raise RuntimeError('host scheduler verification required before frozen entry')
    admission,cpu=owner_admission('exclusive-owner-timing-host')
    # Read-only scheduler evidence: an environment assertion alone cannot prove
    # the scheduler reserved this physical core against sibling borrowing.
    # This must precede run-frozen, which intentionally hides ~/.local.
    with sqlite3.connect('file:/home/treefor/.local/state/simqueue/queue.sqlite3?mode=ro',uri=True) as db:
        db.row_factory=sqlite3.Row
        row=db.execute('SELECT id,state,partition,cpu_list,mem_bytes,borrow,exclusive_timing,unit_name,request_sha256 FROM jobs WHERE id=?',(int(os.environ['SIMQ_JOB_ID']),)).fetchone()
        if row is None:raise RuntimeError('scheduled job record absent')
        job=dict(row)
    check_job(job,cpu)
    cgroup=Path('/proc/self/cgroup').read_text()
    unified=[x.split(':',2)[2] for x in cgroup.splitlines() if x.startswith('0::')]
    if len(unified)!=1:raise RuntimeError('unified job cgroup required')
    memory_file=Path('/sys/fs/cgroup')/unified[0].lstrip('/')/'memory.max'
    memory_max=memory_file.read_text().strip()
    if memory_max!=str(job['mem_bytes']):raise RuntimeError('actual cgroup memory cap differs from scheduler')
    proof={'schema':'coordinator-owner-host-admission-v1','scheduler_job':job,'host_admission':admission,
      'host_pid':os.getpid(),'host_cgroup':cgroup,'memory_max_path':str(memory_file),'memory_max':memory_max,
      'launcher_sha256':m.sha(Path(__file__)),'comparison_script_sha256':m.sha(ROOT/'tools/compare_candidate.py'),
      'label':label,'tag':tag,'verified_unix':time.time()}
    path=m.OUT/(tag+'-host-admission.json')
    with path.open('x') as f:json.dump(proof,f,indent=2);f.write('\n')
    command=[str(m.F),'/usr/bin/python','-I',str(Path(__file__).resolve()),label,tag,'--frozen-proof',str(path)]
    result=subprocess.run(command,cwd=ROOT)
    receipt={'schema':'coordinator-owner-host-exit-v1','proof_sha256':m.sha(path),
      'launcher_sha256':m.sha(Path(__file__)),'command':command,'returncode':result.returncode,'finished_unix':time.time()}
    (m.OUT/(tag+'-host-exit.json')).write_text(json.dumps(receipt,indent=2)+'\n')
    return result.returncode

def frozen_main(label,tag,proof_path):
    if not m.inside_frozen():raise RuntimeError('frozen scheduled owner job required')
    admission,cpu=owner_admission('exclusive-owner-timing')
    path=Path(proof_path)
    if path!=m.OUT/(tag+'-host-admission.json'):raise RuntimeError('unexpected host proof path')
    proof=json.loads(path.read_text());job=proof['scheduler_job'];check_job(job,cpu)
    if proof['label']!=label or proof['tag']!=tag or proof['host_cgroup']!=Path('/proc/self/cgroup').read_text():raise RuntimeError('host/frozen job continuity differs')
    if not 0<=time.time()-proof['verified_unix']<60:raise RuntimeError('host scheduler verification is stale')
    if proof['launcher_sha256']!=m.sha(Path(__file__)) or proof['comparison_script_sha256']!=m.sha(ROOT/'tools/compare_candidate.py'):raise RuntimeError('timing script changed after host verification')
    if Path(proof['memory_max_path']).read_text().strip()!=str(job['mem_bytes']):raise RuntimeError('frozen actual memory cap differs')
    siblings=Path('/sys/devices/system/cpu/cpu%d/topology/thread_siblings_list'%cpu).read_text().strip()
    if set(map(int,siblings.split(',')))!={cpu,cpu+8}:raise RuntimeError('unexpected SMT topology')
    wanted=LOCKS/'timing-wanted.lock';sibling=LOCKS/('cpu%d.lock'%(cpu+8))
    out=m.OUT/(tag+'-launcher.json')
    record={'schema':'coordinator-cooperative-owner-timing-v2','admission':admission,'scheduler_job':job,
      'host_proof_path':str(path),'host_proof_sha256':m.sha(path),'host_admission':proof,
      'thread_siblings_list':siblings,'wanted_lock':str(wanted),'sibling_lock':str(sibling),'lock_wait_limit_s':30,
      'launcher_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      'comparison_script_sha256':hashlib.sha256((ROOT/'tools/compare_candidate.py').read_bytes()).hexdigest(),
      'replaces_canceled_queued_jobs':[327,338],'started_unix':time.time(),'status':'waiting_for_sibling'}
    def save():out.write_text(json.dumps(record,indent=2)+'\n')
    save()
    with wanted.open('a') as w,sibling.open('a') as s:
        if not bounded_lock(w,fcntl.LOCK_SH):
            record['status']='refused_wanted_lock';save();return 75
        start=time.monotonic()
        if not bounded_lock(s,fcntl.LOCK_EX):
            record['status']='refused_sibling_lock';record['wait_s']=time.monotonic()-start;save();return 75
        record['wait_s']=time.monotonic()-start
        fcntl.flock(w,fcntl.LOCK_UN)
        record['status']='sibling_lock_held';save()
        # The actual load/memory gates remain inside this unchanged comparison.
        result=subprocess.run([sys.executable,'-I',str(ROOT/'tools/compare_candidate.py'),label,tag,'timing'],cwd=ROOT)
        record['returncode']=result.returncode
        record['status']='completed' if result.returncode==0 else 'comparison_refused_or_failed'
        record['finished_unix']=time.time();save()
        return result.returncode

if __name__=='__main__':
    args=sys.argv[1:]
    if len(args)==2:raise SystemExit(host_main(*args))
    if len(args)==4 and args[2]=='--frozen-proof':raise SystemExit(frozen_main(args[0],args[1],args[3]))
    raise SystemExit('usage: owner_timing.py LABEL TAG [--frozen-proof PATH]')
