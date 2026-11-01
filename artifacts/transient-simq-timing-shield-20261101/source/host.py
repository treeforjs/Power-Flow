#!/usr/bin/python -I
"""Host adapter. 'plan' reads only; 'canary' requires reviewed root staging.

No arbitrary command, CPU set, cgroup path, output path or database is accepted.
The independent child watchdog owns all writes and rollback; parent death
closes its liveness pipe. Neither process signals or moves any job.
"""
from __future__ import annotations
import argparse, fcntl, hashlib, importlib.util, json, os, re, signal, sqlite3, stat, sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ORIGINAL = Path('/home/treefor/projects/runtimes/transient-simq-timing-shield-v1-20261101')
CG = Path('/sys/fs/cgroup')
STATE = Path('/home/treefor/.local/state/simqueue')
LEASE = Path('/home/treefor/projects/.workspace/outputs/pic-profile/testing-core.lock')
STUDY = Path('/home/treefor/projects/.workspace/outputs/cpu-isa-study-20261031/locks')
MEM = 4 * 1024**3
UID = 1000


def require(value, message):
    if not value:
        raise RuntimeError(message)


def root_safe(p):
    # Privileged execution must be staged under a root-owned, non-writable
    # directory. Never execute a helper imported from a mutable project tree.
    for q in (p,) + tuple(p.parents):
        s = q.lstat()
        require(not stat.S_ISLNK(s.st_mode) and s.st_uid == 0 and not s.st_mode & 0o022,
                'privileged source staging is not root-owned and protected: ' + str(q))


if os.geteuid() == 0:
    root_safe(Path(__file__).absolute())
    root_safe(HERE/'shield.py')
spec = importlib.util.spec_from_file_location('shield', HERE/'shield.py')
core = importlib.util.module_from_spec(spec)
spec.loader.exec_module(core)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def source_identity():
    manifest = {n: sha(HERE/n) for n in ('host.py', 'shield.py')}
    return hashlib.sha256(json.dumps(manifest, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def redacted(value):
    if isinstance(value, dict):
        return {k: '[redacted]' if k == 'launch_token' else redacted(v) for k,v in value.items()}
    if isinstance(value, list):
        return [redacted(v) for v in value]
    return value


def report(value):
    try:
        print(json.dumps(redacted(value), sort_keys=True), flush=True)
    except BrokenPipeError:
        pass


def start_ticks(pid):
    return int((Path('/proc')/str(pid)/'stat').read_text().rsplit(')', 1)[1].split()[19])


def proc_cgroup(pid):
    rows = (Path('/proc')/str(pid)/'cgroup').read_text().splitlines()
    unified = [r[3:] for r in rows if r.startswith('0::')]
    require(len(unified) == 1, 'unified process cgroup missing')
    return unified[0]


def validate_row(row, job, unit):
    require(row is not None and row['id'] == job and row['state'] == 'running', 'exact job must be running')
    token = row['launch_token']
    require(isinstance(token, str) and re.fullmatch('[0-9a-f]{32}', token), 'invalid launch token')
    require(unit == 'simq-job-%d-%s.service' % (job, token[:12]) and row['unit_name'] == unit, 'unit/token mismatch')
    require(row['partition'] == 'powerflow' and row['cpu_list'] == '6' and row['borrow'] == 0
            and row['exclusive_timing'] == 1, 'CPU6 exclusive owner allocation required')
    require(row['mem_bytes'] == MEM and 0 <= row['swap_bytes'] <= 1024**3
            and 90 <= row['time_limit_s'] <= 180 and row['launch_backend'] == 'systemd-user', 'canary resource limits differ')
    require(row['gpus_requested'] == 0, 'canary cannot hold GPU allocation')
    expected = ['/usr/bin/python', '-I', str(ORIGINAL/'canary_wait.py'), '--seconds', '75']
    require(json.loads(row['argv_json']) == expected, 'only the exact waiting canary is supported')
    require(json.loads(row['lock_paths_json']) == [{'mode': 'exclusive', 'path': str(LEASE), 'resource': 'cpu:6'}],
            'job lease declaration differs')


class Host:
    def __init__(self, job, unit, canary_sha, journal=None, locks=()):
        self.job, self.unit, self.canary_sha = job, unit, canary_sha
        self.journal = journal
        self.locks = locks
        self.initial_wrapper = None

    def file(self, p):
        require(not p.startswith('/') and '..' not in p.split('/'), 'invalid cgroup path')
        return CG/p

    def read(self, p):
        return self.file(p).read_text().strip()

    def exists(self, p):
        return self.file(p).exists()

    def inode(self, p):
        s = self.file(p).stat()
        return [s.st_dev, s.st_ino]

    def write(self, p, value):
        require(os.geteuid() == 0 and self.journal is not None, 'no armed journal; writes refused')
        # Engine supplies only its fixed ancestry/job paths and these controls.
        group, _, name = p.rpartition('/')
        if not group:
            name = p
        require(group in core.ANCESTORS + (core.BASE+'/'+self.unit,), 'write outside exact ancestry')
        require(name in ('cgroup.subtree_control', 'cpuset.cpus.exclusive', 'cpuset.cpus', 'cpuset.cpus.partition'),
                'unsupported control write')
        fd = os.open(self.file(p), os.O_WRONLY | os.O_CLOEXEC | os.O_NOFOLLOW)
        try:
            data = (value+'\n').encode()
            require(os.write(fd, data) == len(data), 'short control write')
        finally:
            os.close(fd)

    def audit(self, value):
        value = {'monotonic_s': time.monotonic(), **value}
        line = (json.dumps(value, sort_keys=True)+'\n').encode()
        if self.journal is not None:
            require(os.write(self.journal, line) == len(line), 'short audit write')
            os.fsync(self.journal)
        report(value)  # Public output never exposes the full launch token.

    def row(self):
        columns = ('id,state,partition,cpu_list,borrow,exclusive_timing,mem_bytes,swap_bytes,time_limit_s,'
                   'gpus_requested,launch_backend,unit_name,launch_token,argv_json,lock_paths_json,pid,boot_id')
        with sqlite3.connect('file:'+str(STATE/'queue.sqlite3')+'?mode=ro', uri=True) as db:
            db.row_factory = sqlite3.Row
            row = db.execute('SELECT '+columns+' FROM jobs WHERE id=?', (self.job,)).fetchone()
            return dict(row) if row else None

    def allocation(self, group, active):
        for p,fd in self.locks:
            held, current = os.fstat(fd), p.stat()
            require((held.st_dev,held.st_ino) == (current.st_dev,current.st_ino), 'cooperative lock file replaced')
        row = self.row()
        validate_row(row, self.job, self.unit)
        require(sha(ORIGINAL/'canary_wait.py') == self.canary_sha, 'canary source changed')
        require(Path('/sys/devices/system/cpu/cpu6/topology/thread_siblings_list').read_text().strip() == '6,14',
                'SMT topology differs')
        require(row['boot_id'] == Path('/proc/sys/kernel/random/boot_id').read_text().strip(), 'job boot differs')
        pid = row['pid']
        require(isinstance(pid, int) and pid > 1, 'wrapper PID missing')
        ticks = start_ticks(pid)
        age = float(Path('/proc/uptime').read_text().split()[0]) - ticks/os.sysconf('SC_CLK_TCK')
        require(0 <= age < (40 if self.initial_wrapper is None else 65), 'canary activation/continuation window expired')
        require(proc_cgroup(pid) == '/'+group, 'wrapper escaped job cgroup')
        require(os.sched_getaffinity(pid) == {6}, 'wrapper affinity changed')
        status = (Path('/proc')/str(pid)/'status').read_text().splitlines()
        uids = next(x for x in status if x.startswith('Uid:')).split()[1:]
        require(all(int(x) == UID for x in uids), 'wrapper user differs')
        command = (Path('/proc')/str(pid)/'cmdline').read_bytes().split(b'\0')
        wrapper = b'/home/treefor/projects/runtimes/simqueue-0.6.2-20261031/source/simqueue/jobwrap.py'
        require(wrapper in command and b'run' in command, 'PID is not the queue wrapper')
        require(self.read(group+'/memory.max') == str(MEM), 'actual job memory cap differs')
        require(self.read(group+'/memory.swap.max') == str(row['swap_bytes']), 'actual swap cap differs')
        available = int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines()
                             if line.startswith('MemAvailable:'))) * 1024
        require(available >= 2*MEM, 'less than cap plus 4 GiB physical memory headroom')
        require(self.read(group+'/cgroup.freeze') == '0', 'job is frozen')
        require('populated 1' in self.read(group+'/cgroup.events'), 'job cgroup is empty')
        members = self.members(group)
        require(pid in members, 'wrapper absent from job')
        for member in members:
            for task in (Path('/proc')/str(member)/'task').iterdir():
                require(os.sched_getaffinity(int(task.name)) == {6}, 'job thread escaped CPU6')
        s = LEASE.stat()
        holders = []
        for line in Path('/proc/locks').read_text().splitlines():
            words = line.split()
            if '->' in words:
                continue
            if len(words) >= 8 and words[1:4] == ['FLOCK', 'ADVISORY', 'WRITE']:
                devino = words[5].split(':')
                if len(devino) == 3 and (int(devino[0],16),int(devino[1],16),int(devino[2])) == (os.major(s.st_dev),os.minor(s.st_dev),s.st_ino):
                    holders.append(int(words[4]))
        require(holders == [pid] or (len(holders) == 1 and holders[0] in members), 'CPU6 lease is not held by this job')
        self.initial_wrapper = (pid,ticks)
        return {'job': row, 'wrapper_start_ticks': ticks, 'cgroup_inode': self.inode(group),
                'lease_inode': [s.st_dev,s.st_ino], 'canary_sha256': self.canary_sha}

    def members(self, group):
        p = CG/group
        paths = [p/'cgroup.procs'] + list(p.glob('**/cgroup.procs'))
        result = set()
        for f in paths:
            result.update(map(int, f.read_text().split()))
        return result

    def foreign_workers(self, group):
        online = core.cpus(Path('/sys/devices/system/cpu/online').read_text())
        require(online == set(range(16)), 'online CPU topology differs')
        for p in Path('/proc').glob('[0-9]*'):
            try:
                command = (p/'cmdline').read_bytes()
                if not command:  # kernel per-CPU/IRQ activity is measured, never modified
                    continue
                cg = proc_cgroup(int(p.name))
                require(cg != '/', 'root-cgroup userspace cannot be shielded safely')
                if cg == '/'+group or cg.startswith('/'+group+'/'):
                    continue
                for t in (p/'task').iterdir():
                    affinity = os.sched_getaffinity(int(t.name))
                    require(not affinity & core.PAIR or affinity == online,
                            'foreign restricted-affinity worker overlaps CPU6/14: '+t.name)
            except (FileNotFoundError, ProcessLookupError):
                continue

    def outsiders(self, group):
        for f in CG.glob('**/cpuset.cpus.effective'):
            rel = str(f.parent.relative_to(CG))
            if rel == '.' or rel in core.ANCESTORS or rel == group or rel.startswith(group+'/'):
                continue
            try:
                require(not core.cpus(f.read_text()) & core.PAIR, 'outsider can use shield CPUs: '+rel)
            except FileNotFoundError:
                continue
        # A branch without a delegated cpuset file inherits the closest checked
        # ancestor. Every root branch has a file once root enables cpuset.


def lock_file(p, mode):
    fd = os.open(p, os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW)
    try:
        s = os.fstat(fd)
        require(stat.S_ISREG(s.st_mode) and s.st_uid == UID, 'invalid cooperative lock file')
        fcntl.flock(fd, mode | fcntl.LOCK_NB)
        return fd
    except BaseException:
        os.close(fd)
        raise


def monitor(transaction, alive, stopped, clock=time.monotonic, sleep=time.sleep):
    deadline = clock()+20
    while clock() < deadline and not stopped[0]:
        try:
            if os.read(alive,1) == b'':
                raise RuntimeError('parent control pipe closed')
        except BlockingIOError:
            pass
        transaction.verify()
        sleep(.2)
    require(not stopped[0], 'watchdog interrupted')


def stop_handlers(stopped):
    for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
        signal.signal(sig, lambda *_: stopped.__setitem__(0, True))


def watchdog(args, alive):
    # This process performs all writes and retains locks until rollback ends.
    stopped = [False]
    stop_handlers(stopped)
    os.set_blocking(alive, False)
    require(proc_cgroup(os.getpid()) != '/'+core.BASE+'/'+args.unit, 'watchdog cannot live inside job')
    os.sched_setaffinity(0, os.sched_getaffinity(0) - core.PAIR)
    require(not os.sched_getaffinity(0) & core.PAIR, 'watchdog affinity overlaps shield')
    locks, journal, transaction = [], None, None
    code = 1
    try:
        locks.append((STUDY/'timing-wanted.lock',lock_file(STUDY/'timing-wanted.lock', fcntl.LOCK_SH)))
        locks.append((STUDY/'cpu14.lock',lock_file(STUDY/'cpu14.lock', fcntl.LOCK_EX)))
        # Nonblocking acquisition fails promptly; no waiting while owning a core.
        directory = Path('/run')/('powerflow-timing-canary-%d'%args.job)
        directory.mkdir(mode=0o700)
        root_safe(directory)
        journal = os.open(directory/'journal.jsonl', os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC | os.O_NOFOLLOW, 0o600)
        backend = Host(args.job, args.unit, args.canary_sha256, journal, tuple(locks))
        backend.audit({'event':'watchdog_ready','source_identity':source_identity(),'pid':os.getpid(),
                       'cgroup':proc_cgroup(os.getpid()),'affinity':sorted(os.sched_getaffinity(0)),
                       'hard_hold_limit_s':20,'journal':str(directory/'journal.jsonl')})
        transaction = core.Shield(backend,args.unit)
        transaction.setup()
        monitor(transaction,alive,stopped)
        code = 0
    except BaseException as error:
        if transaction:
            try:
                transaction.b.audit({'event':'refused','error':str(error)})
            except Exception as audit_error:
                report({'event':'refused_audit_unavailable','error':str(error),'audit_error':str(audit_error)})
        else:
            report({'event':'refused_before_setup','error':str(error)})
    finally:
        if transaction:
            errors = transaction.rollback()
            try:
                transaction.b.audit({'event':'cleanup','errors':errors,'clean':not errors})
            except Exception as audit_error:
                errors.append('final cleanup audit unavailable: '+str(audit_error))
                report({'event':'cleanup_audit_unavailable','errors':errors,'clean':False})
            if errors:
                code = 2
        if journal is not None:
            os.close(journal)
        for _,fd in reversed(locks):
            os.close(fd)
        os.close(alive)
    return code


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=('plan','canary'))
    p.add_argument('--job',type=int,required=True)
    p.add_argument('--unit',required=True)
    p.add_argument('--source-sha256',required=True)
    p.add_argument('--canary-sha256',required=True)
    args = p.parse_args()
    require(__debug__, 'optimized Python guards are forbidden')
    require(args.job > 0 and re.fullmatch('simq-job-[0-9]+-[0-9a-f]{12}\\.service',args.unit), 'invalid job/unit')
    require(args.source_sha256 == source_identity(), 'helper source identity differs')
    require(re.fullmatch('[0-9a-f]{64}',args.canary_sha256), 'invalid canary pin')
    if args.mode == 'plan':
        t = core.Shield(Host(args.job,args.unit,args.canary_sha256), args.unit)
        t.preflight()
        report({'status':'read_only_preflight_passed','allocation':t.identity})
        return 0
    require(os.geteuid() == 0, 'live canary requires reviewed privileged staging')
    read, write = os.pipe2(os.O_CLOEXEC)
    pid = os.fork()
    if pid == 0:
        os.close(write)
        os._exit(watchdog(args,read))
    os.close(read)
    # Parent never writes controllers. Any parent exit closes the heartbeat
    # pipe; child independently rolls back. No signal is sent to any worker.
    try:
        _, status = os.waitpid(pid,0)
        return os.waitstatus_to_exitcode(status)
    finally:
        os.close(write)


if __name__ == '__main__':
    raise SystemExit(main())
