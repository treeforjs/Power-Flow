import ast
import importlib.util
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
import unittest

ROOT=Path(__file__).resolve().parents[1]


def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    obj=importlib.util.module_from_spec(spec);spec.loader.exec_module(obj)
    return obj


class RunnerTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(prefix='runner-audit-',dir=ROOT/'evidence')
        self.addCleanup(self.tmp.cleanup);self.dir=Path(self.tmp.name)
        self.runner=module(ROOT/'source/run_jobs.py','runner_test')
        self.runner.QUEUED=True
        self.runner.CORES=[min(os.sched_getaffinity(0))]
        self.runner.FROZEN='/usr/bin/env' # test already runs inside the allocated frozen process
        self.runner.BINDING='fixture-binding'
        self.binary=self.dir/'synthetic-swarm'
        self.blob=self.dir/'blob';self.blob.write_bytes(b'synthetic blob')
        self.output=self.dir/'output.json'
        self.job={'name':'fixture','output':str(self.output),'log':str(self.dir/'run.log'),'cwd':str(self.dir),
                  'argv':[str(self.binary),str(self.blob),str(self.output),'100,200','3.22e22','300','2000','123']}
        self.jobs=self.dir/'jobs.json';self.jobs.write_text(json.dumps([self.job]))
        self.receipt=self.dir/'receipt.json'

    def script(self,delay=0):
        numbers=self.runner._proof.V.POINT_NUMBERS
        counts=self.runner._proof.V.POINT_COUNTS
        self.binary.write_text(f'''#!/usr/bin/python3
import json,hashlib,sys,time,os
from pathlib import Path
Path({str(self.dir/'child.pid')!r}).write_text(str(os.getpid()))
with open({str(self.dir/'executions')!r},'a') as f:f.write('run\\n')
time.sleep({delay!r})
rows=[]
for en in sys.argv[3].split(','):
 r={{k:1.0 for k in {numbers!r}}};r.update({{k:0 for k in {counts!r}}})
 r.update(E_over_N_Td=float(en),process_events={{'synthetic':0}});rows.append(r)
Path(sys.argv[2]).write_text(json.dumps(dict(blob_sha256=hashlib.sha256(Path(sys.argv[1]).read_bytes()).hexdigest(),
 N_m3=float(sys.argv[4]),T_gas_K=float(sys.argv[5]),electrons=int(sys.argv[6]),seed=int(sys.argv[7]),points=rows)))
''')
        self.binary.chmod(0o755)

    def test_every_staged_python_parses(self):
        for path in (ROOT/'source').glob('*.py'):
            with self.subTest(path=path.name):ast.parse(path.read_text(),filename=str(path))

    def test_publish_then_verified_resume_does_not_rerun(self):
        self.script()
        self.assertEqual(self.runner.main(str(self.jobs),str(self.receipt)),0)
        self.assertTrue(self.runner._proof.completed(self.job,self.runner.BINDING))
        first=self.output.read_bytes()
        self.assertEqual(self.runner.main(str(self.jobs),str(self.receipt)),0)
        self.assertEqual((self.dir/'executions').read_text(),'run\n')
        self.assertEqual(self.output.read_bytes(),first)
        self.assertTrue(json.loads(self.receipt.read_text())['jobs'][0]['resumed'])

    def test_termination_reaps_own_child_and_keeps_public_output(self):
        self.script(delay=30)
        self.output.write_bytes(b'previous authoritative output\n')
        code=f'''import importlib.util,os
s=importlib.util.spec_from_file_location('r',{str(ROOT/'source/run_jobs.py')!r});r=importlib.util.module_from_spec(s);s.loader.exec_module(r)
r.QUEUED=True;r.CORES=[min(os.sched_getaffinity(0))];r.FROZEN='/usr/bin/env';r.BINDING='fixture-binding'
r.main({str(self.jobs)!r},{str(self.receipt)!r})
'''
        manager=subprocess.Popen([sys.executable,'-I','-c',code],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        try:
            deadline=time.monotonic()+5
            while not (self.dir/'child.pid').exists() and time.monotonic()<deadline:
                if manager.poll() is not None:self.fail(manager.communicate()[1].decode())
                time.sleep(.02)
            self.assertTrue((self.dir/'child.pid').exists())
            pid=int((self.dir/'child.pid').read_text())
            manager.send_signal(signal.SIGTERM);manager.communicate(timeout=8)
            with self.assertRaises(ProcessLookupError):os.kill(pid,0)
            self.assertEqual(self.output.read_bytes(),b'previous authoritative output\n')
            self.assertFalse(Path(str(self.output)+'.complete.json').exists())
            self.assertFalse(list(self.dir.glob('*.partial')))
        finally:
            if manager.poll() is None:manager.kill();manager.communicate(timeout=5)
            self.clean_toy_child()

    def clean_toy_child(self):
        path=self.dir/'child.pid'
        if path.exists():
            try:os.killpg(int(path.read_text()),signal.SIGKILL)
            except ProcessLookupError:pass

    def test_run_checks_timeout_reaps_manager_child_sessions(self):
        self.script(delay=30)
        self.output.write_bytes(b'previous authoritative output\n')
        code=f'''import importlib.util,os
s=importlib.util.spec_from_file_location('r',{str(ROOT/'source/run_jobs.py')!r});r=importlib.util.module_from_spec(s);s.loader.exec_module(r)
r.QUEUED=True;r.CORES=[min(os.sched_getaffinity(0))];r.FROZEN='/usr/bin/env';r.BINDING='fixture-binding'
r.main({str(self.jobs)!r},{str(self.receipt)!r})
'''
        tree=ast.parse((ROOT/'source/run_checks.py').read_text())
        actual_class=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='Run')
        scope={'Path':Path,'R':self.dir,'os':os,'signal':signal,'subprocess':subprocess,'time':time}
        exec(compile(ast.Module(body=[actual_class],type_ignores=[]),'actual_run_checks_class','exec'),scope)
        run=scope['Run'].__new__(scope['Run']);run.out=self.dir
        try:
            with self.assertRaises(subprocess.TimeoutExpired):
                run.step('manager-timeout',[sys.executable,'-I','-c',code],kind='host',timeout=1)
            self.assertTrue((self.dir/'child.pid').exists())
            with self.assertRaises(ProcessLookupError):os.kill(int((self.dir/'child.pid').read_text()),0)
            self.assertEqual(self.output.read_bytes(),b'previous authoritative output\n')
            self.assertFalse(Path(str(self.output)+'.complete.json').exists())
            self.assertFalse(list(self.dir.glob('*.partial')))
        finally:self.clean_toy_child()


if __name__=='__main__':unittest.main(verbosity=2)
