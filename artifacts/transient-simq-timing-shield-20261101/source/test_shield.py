#!/usr/bin/python -I
"""Metadata-only tests: all writes target the in-memory simulated filesystem."""
import contextlib, copy, importlib.util, io, json, os, signal, unittest
from unittest.mock import patch
from pathlib import Path

HERE = Path(__file__).resolve().parent
def module(name):
    spec = importlib.util.spec_from_file_location(name,HERE/(name+'.py'))
    mod = importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
core, host = module('shield'), module('host')
UNIT='simq-job-123-abcdef012345.service'
GROUP=core.BASE+'/'+UNIT


class Fake:
    def __init__(self):
        self.files={};self.ids={};self.next_id=1;self.events=[];self.writes=[]
        self.fail_before=None;self.fail_after=None;self.fail_audit=None;self.persistent_audit_failure=False;self.count=0
        self.bad_foreign=False;self.bad_outside=False;self.job_identity={'job':123,'unit':UNIT,'pid_start':5}
        self.groups=list(core.ANCESTORS)+[GROUP,'system.slice','init.scope',core.BASE+'/other.service']
        self.put('cgroup.controllers','cpuset cpu io memory pids')
        self.put('cpuset.cpus.isolated','')
        for g in self.groups:
            self.put(core.path(g,'cgroup.subtree_control'),'cpu memory pids')
            if g:
                self.put(core.path(g,'cgroup.type'),'domain')
                self.put(core.path(g,'cgroup.procs'),'101' if g==GROUP else '')
        self.baseline=copy.deepcopy(self.files)

    def put(self,p,v):
        self.files[p]=v
        if p not in self.ids:self.ids[p]=self.next_id;self.next_id+=1
    def read(self,p):return self.files[p]
    def exists(self,p):return p in self.files
    def inode(self,p):return [1,self.ids[p]]
    def audit(self,v):
        if self.persistent_audit_failure:raise OSError('persistent audit storage failure')
        if self.fail_audit==len(self.events):self.fail_audit=None;raise OSError('injected audit failure')
        self.events.append(copy.deepcopy(v))
    def allocation(self,g,active):
        core.require(g==GROUP,'group differs')
        core.require(core.path(g,'cgroup.procs') in self.files,'job disappeared')
        return copy.deepcopy(self.job_identity)
    def foreign_workers(self,g):core.require(not self.bad_foreign,'foreign restricted worker')
    def outsiders(self,g):
        core.require(not self.bad_outside,'outsider still overlaps')
        for other in self.groups:
            if not other or other in core.ANCESTORS or other==GROUP:continue
            p=core.path(other,'cpuset.cpus.effective')
            if p in self.files:core.require(not core.cpus(self.files[p])&core.PAIR,'effective outsider overlap')
    def children(self,g):
        return [c for c in self.groups if c and c.rpartition('/')[0]==g]
    def materialize(self,g):
        for c in self.children(g):
            for n,v in [('cpuset.cpus',''),('cpuset.cpus.exclusive',''),('cpuset.cpus.partition','member'),
                        ('cpuset.cpus.effective','0-15'),('cpuset.cpus.exclusive.effective','')]:
                self.put(core.path(c,n),v)
    def recompute(self):
        isolated=self.files.get(GROUP+'/cpuset.cpus.partition')=='isolated'
        self.files['cpuset.cpus.isolated']='6,14' if isolated else ''
        for g in self.groups:
            p=core.path(g,'cpuset.cpus.effective')
            if p in self.files:
                self.files[p]='6,14' if isolated and g==GROUP else ('0-5,7-13,15' if isolated and g not in core.ANCESTORS else '0-15')
            p=core.path(g,'cpuset.cpus.exclusive.effective')
            if p in self.files:self.files[p]=self.files[core.path(g,'cpuset.cpus.exclusive')]
    def write(self,p,v):
        self.count+=1
        if self.count==self.fail_before:raise OSError('injected pre-write failure')
        self.writes.append((p,v))
        if p.endswith('cgroup.subtree_control'):
            enabled=set(self.files[p].split())
            if v=='+cpuset':
                enabled.add('cpuset');self.materialize(p.rpartition('/')[0])
            elif v=='-cpuset':
                enabled.remove('cpuset')
                for child in self.children(p.rpartition('/')[0]):
                    core.require('cpuset' not in self.files[core.path(child,'cgroup.subtree_control')].split(),'child still uses controller')
                    for name in list(self.files):
                        if name.startswith(child+'/cpuset.') and '/' not in name[len(child)+1:]:
                            del self.files[name];del self.ids[name]
            else:raise RuntimeError('unexpected controller command')
            self.files[p]=' '.join(c for c in self.files['cgroup.controllers'].split() if c in enabled)
        else:self.files[p]=v
        self.recompute()
        if self.count==self.fail_after:raise OSError('injected post-write failure')


class TransactionTests(unittest.TestCase):
    def make(self):
        b=Fake();return b,core.Shield(b,UNIT)
    def test_read_only_preflight(self):
        b,t=self.make();t.preflight();self.assertEqual(b.files,b.baseline);self.assertFalse(b.writes)
    def test_full_setup_restore_exact(self):
        b,t=self.make();t.setup();self.assertTrue(t.active);self.assertEqual(t.rollback(),[])
        self.assertEqual(b.files,b.baseline)
    def test_every_prewrite_failure_rolls_back(self):
        b,t=self.make();t.setup();n=b.count
        self.assertEqual(n,12)
        for point in range(1,n+1):
            with self.subTest(point=point):
                b,t=self.make();b.fail_before=point
                with self.assertRaises(OSError):t.setup()
                self.assertEqual(t.rollback(),[]);self.assertEqual(b.files,b.baseline)
    def test_every_postwrite_failure_rolls_back(self):
        for point in range(1,13):
            with self.subTest(point=point):
                b,t=self.make();b.fail_after=point
                with self.assertRaises(OSError):t.setup()
                self.assertEqual(t.rollback(),[]);self.assertEqual(b.files,b.baseline)
    def test_audit_failure_before_and_after_each_write(self):
        for point in range(25):
            with self.subTest(point=point):
                b,t=self.make();b.fail_audit=point
                with self.assertRaises(OSError):t.setup()
                self.assertEqual(t.rollback(),[]);self.assertEqual(b.files,b.baseline)
    def test_persistent_audit_failure_during_cleanup_still_restores(self):
        b,t=self.make();t.setup();b.persistent_audit_failure=True
        errors=t.rollback();self.assertEqual(len(errors),12)
        self.assertTrue(all('audit unavailable' in e for e in errors));self.assertEqual(b.files,b.baseline)
    def test_persistent_audit_failure_and_conflict_preserves_conflict(self):
        b,t=self.make();t.setup();b.persistent_audit_failure=True
        b.files[GROUP+'/cpuset.cpus']='6'
        errors=t.rollback();self.assertEqual(len(errors),3)
        self.assertEqual(b.files[GROUP+'/cpuset.cpus.partition'],'member')
        self.assertEqual(b.files[GROUP+'/cpuset.cpus'],'6')
        self.assertIn('cpuset',b.files['cgroup.subtree_control'])
    def test_foreign_worker_refuses_without_writes(self):
        b,t=self.make();b.bad_foreign=True
        with self.assertRaises(core.Refusal):t.setup()
        self.assertFalse(b.writes)
    def test_other_partition_refuses_without_writes(self):
        b,t=self.make();b.files['cpuset.cpus.isolated']='7,15'
        with self.assertRaises(core.Refusal):t.setup()
        self.assertFalse(b.writes)
    def test_unavailable_controller_refuses_without_writes(self):
        b,t=self.make();b.files['cgroup.controllers']='cpu memory pids'
        with self.assertRaises(core.Refusal):t.setup()
        self.assertFalse(b.writes)
    def test_ancestor_direct_process_refuses(self):
        b,t=self.make();b.files[core.ANCESTORS[2]+'/cgroup.procs']='99'
        with self.assertRaises(core.Refusal):t.setup()
        self.assertFalse(b.writes)
    def test_ancestor_partition_or_exclusive_refuses(self):
        for name,v in [('cpuset.cpus.partition','root'),('cpuset.cpus.exclusive','6,14')]:
            b,t=self.make();b.put(core.ANCESTORS[1]+'/'+name,v)
            with self.assertRaises(core.Refusal):t.setup()
            self.assertFalse(b.writes)
    def test_existing_controller_preserved(self):
        b,t=self.make();b.write('cgroup.subtree_control','+cpuset');baseline=copy.deepcopy(b.files)
        t.setup();self.assertEqual(t.rollback(),[]);self.assertEqual(b.files,baseline)
    def test_verification_outsider_failure_rolls_back(self):
        b,t=self.make();b.bad_outside=True
        with self.assertRaises(core.Refusal):t.setup()
        self.assertEqual(t.rollback(),[]);self.assertEqual(b.files,b.baseline)
    def test_changed_live_identity_refuses(self):
        b,t=self.make();t.setup();b.job_identity['pid_start']=6
        with self.assertRaises(core.Refusal):t.verify()
        self.assertEqual(t.rollback(),[]);self.assertEqual(b.files,b.baseline)
    def test_partition_degradation_refuses_and_does_not_override(self):
        b,t=self.make();t.setup();p=GROUP+'/cpuset.cpus.partition';b.files[p]='isolated invalid (parent change)'
        with self.assertRaises(core.Refusal):t.verify()
        count=len(b.writes);errors=t.rollback()
        self.assertTrue(errors);self.assertEqual(len(b.writes),count)
        self.assertEqual(b.files[p],'isolated invalid (parent change)')
    def test_concurrent_same_path_change_stops_ancestor_cleanup(self):
        b,t=self.make();t.setup();p=GROUP+'/cpuset.cpus';b.files[p]='6'
        errors=t.rollback();self.assertTrue(errors)
        self.assertEqual(b.files[p],'6');self.assertIn('cpuset',b.files['cgroup.subtree_control'])
    def test_replaced_inode_is_not_touched(self):
        b,t=self.make();t.setup();p=GROUP+'/cpuset.cpus.partition';b.ids[p]+=10000
        count=len(b.writes);self.assertTrue(t.rollback());self.assertEqual(len(b.writes),count)
    def test_job_disappearance_allows_ancestor_cleanup(self):
        b,t=self.make();t.setup()
        for p in list(b.files):
            if p.startswith(GROUP+'/'):del b.files[p];del b.ids[p]
        b.recompute();b.groups.remove(GROUP)
        self.assertEqual(t.rollback(),[])
        expected={p:v for p,v in b.baseline.items() if not p.startswith(GROUP+'/')}
        self.assertEqual(b.files,expected)
    def test_cleanup_failure_does_not_change_parent(self):
        b,t=self.make();t.setup();b.fail_before=b.count+1
        self.assertTrue(t.rollback());self.assertEqual(b.files[GROUP+'/cpuset.cpus.partition'],'isolated')
        self.assertIn('cpuset',b.files['cgroup.subtree_control'])
    def test_parent_controller_concurrent_change_preserved(self):
        b,t=self.make();t.setup();b.files['cgroup.subtree_control']='cpuset cpu io memory pids'
        self.assertTrue(t.rollback());self.assertEqual(b.files['cgroup.subtree_control'],'cpuset cpu io memory pids')
    def test_cpu_ranges(self):
        self.assertEqual(core.cpus('0-2,6,14'),{0,1,2,6,14})
        for v in ('4-2','0-99999','1-2-3'):
            with self.assertRaises(core.Refusal):core.cpus(v)


class RowTests(unittest.TestCase):
    def row(self):
        return {'id':123,'state':'running','launch_token':'abcdef012345'+'a'*20,'unit_name':UNIT,
                'partition':'powerflow','cpu_list':'6','borrow':0,'exclusive_timing':1,
                'mem_bytes':host.MEM,'swap_bytes':1024**3,'time_limit_s':120,
                'launch_backend':'systemd-user','gpus_requested':0,
                'argv_json':json.dumps(['/usr/bin/python','-I',str(host.ORIGINAL/'canary_wait.py'),'--seconds','75']),
                'lock_paths_json':json.dumps([{'mode':'exclusive','path':str(host.LEASE),'resource':'cpu:6'}])}
    def test_exact_row(self):host.validate_row(self.row(),123,UNIT)
    def test_each_bad_allocation(self):
        changes={'id':124,'state':'starting','launch_token':'0'*32,'unit_name':'other.service',
                 'partition':'spheral','cpu_list':'7','borrow':1,'exclusive_timing':0,
                 'mem_bytes':host.MEM+1,'swap_bytes':2*1024**3,'time_limit_s':181,
                 'launch_backend':'other','gpus_requested':1,
                 'argv_json':json.dumps(['/usr/bin/python','simulation.py']), 'lock_paths_json':'[]'}
        for k,v in changes.items():
            with self.subTest(field=k):
                r=self.row();r[k]=v
                with self.assertRaises(RuntimeError):host.validate_row(r,123,UNIT)
    def test_unarmed_host_write_always_refused(self):
        b=host.Host(123,UNIT,'0'*64)
        with self.assertRaises(RuntimeError):b.write('cgroup.subtree_control','+cpuset')
    def test_public_audit_redacts_nested_launch_token(self):
        token='1234567890abcdef1234567890abcdef'
        value={'allocation':{'job':{'launch_token':token}},'list':[{'launch_token':token}]}
        out=io.StringIO()
        with contextlib.redirect_stdout(out):host.Host(123,UNIT,'0'*64).audit(value)
        self.assertNotIn(token,out.getvalue());self.assertEqual(out.getvalue().count('[redacted]'),2)
        self.assertEqual(value['allocation']['job']['launch_token'],token)


class MonitorTests(unittest.TestCase):
    def setUp(self):
        self.r,self.w=os.pipe();os.set_blocking(self.r,False)
        self.b=Fake();self.t=core.Shield(self.b,UNIT);self.t.setup()
    def tearDown(self):
        os.close(self.r)
        if self.w is not None:os.close(self.w)
    def test_parent_death_refuses_and_restores(self):
        os.close(self.w);self.w=None
        try:
            with self.assertRaisesRegex(RuntimeError,'parent control pipe closed'):
                host.monitor(self.t,self.r,[False],sleep=lambda _:None)
        finally:self.assertEqual(self.t.rollback(),[])
        self.assertEqual(self.b.files,self.b.baseline)
    def test_signal_request_refuses_and_restores(self):
        try:
            with self.assertRaisesRegex(RuntimeError,'watchdog interrupted'):
                host.monitor(self.t,self.r,[True],sleep=lambda _:None)
        finally:self.assertEqual(self.t.rollback(),[])
        self.assertEqual(self.b.files,self.b.baseline)
    def test_bounded_normal_monitor_and_restore(self):
        times=iter([0,0,5,10,15,20])
        host.monitor(self.t,self.r,[False],clock=lambda:next(times),sleep=lambda _:None)
        self.assertEqual(self.t.rollback(),[]);self.assertEqual(self.b.files,self.b.baseline)
    def test_lost_allocation_refuses_and_restores(self):
        self.b.job_identity['pid_start']=999
        try:
            with self.assertRaises(core.Refusal):host.monitor(self.t,self.r,[False],sleep=lambda _:None)
        finally:self.assertEqual(self.t.rollback(),[])
        self.assertEqual(self.b.files,self.b.baseline)
    def test_term_int_hup_handlers_request_cleanup(self):
        handlers={};stopped=[False]
        with patch.object(host.signal,'signal',lambda sig,fn:handlers.__setitem__(sig,fn)):
            host.stop_handlers(stopped)
        self.assertEqual(set(handlers),{signal.SIGTERM,signal.SIGINT,signal.SIGHUP})
        for handler in handlers.values():
            stopped[0]=False;handler(None,None);self.assertTrue(stopped[0])
    def test_actual_term_to_own_fake_watchdog_restores(self):
        # Only this test's newly forked child receives a signal. Its backend
        # is the in-memory fake; no cgroup or other worker is touched.
        reader,writer=os.pipe();pid=os.fork()
        if pid==0:
            os.close(reader);stopped=[False];host.stop_handlers(stopped)
            os.write(writer,b'R')
            refusal=None
            try:host.monitor(self.t,self.r,stopped)
            except RuntimeError as error:refusal=str(error)
            errors=self.t.rollback()
            result={'refusal':refusal,'errors':errors,'restored':self.b.files==self.b.baseline}
            os.write(writer,json.dumps(result).encode());os.close(writer);os._exit(0)
        os.close(writer)
        self.assertEqual(os.read(reader,1),b'R')
        os.kill(pid,signal.SIGTERM)
        waited,status=os.waitpid(pid,0)
        result=json.loads(os.read(reader,8192));os.close(reader)
        self.assertEqual(waited,pid);self.assertEqual(os.waitstatus_to_exitcode(status),0)
        self.assertEqual(result,{'refusal':'watchdog interrupted','errors':[],'restored':True})


if __name__=='__main__':unittest.main(verbosity=2)
