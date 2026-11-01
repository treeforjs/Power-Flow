"""Synthetic preparation guards; does not mutate numerical sources or inputs."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('actual_prep',ROOT/'tools/prep_stage.py')
P = importlib.util.module_from_spec(spec); spec.loader.exec_module(P)


class Guards(unittest.TestCase):
    def test_external_build_flags_refuse_before_publication(self):
        before = {str(p.relative_to(ROOT)):P.sha(p) for p in (ROOT/'evidence').glob('prep-*.json')}
        result = subprocess.run([sys.executable,'-I',str(ROOT/'tools/prep_stage.py'),'guard-probe'],
            env=dict(os.environ,PF_BUILD_ISA_FLAGS='-ffast-math'),capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('unqualified build flags refused',result.stderr)
        self.assertFalse((ROOT/'evidence/prep-guard-probe.json').exists())
        after = {str(p.relative_to(ROOT)):P.sha(p) for p in (ROOT/'evidence').glob('prep-*.json')}
        self.assertEqual(before,after)

    def test_added_changed_missing_sources_and_dependencies_refuse(self):
        original = P.R
        with tempfile.TemporaryDirectory(prefix='prep-guard-') as directory:
            r = Path(directory); P.R = r
            try:
                for d in ('include','device','tests','tools','data/sources','data/extracted','adoption','dep'):(r/d).mkdir(parents=True)
                for name in ('run_checks.py','run_prep.py','include/fixed.h','dep/dep.h','data/sources/input.txt'):(r/name).write_text('fixed\n')
                pins = {'source_files':{n:P.sha(r/n) for n in ('run_checks.py','run_prep.py','include/fixed.h')},
                        'input_files':{'data/sources/input.txt':P.sha(r/'data/sources/input.txt')},
                        'dependencies':{str(r/'dep/dep.h'):P.sha(r/'dep/dep.h')},'dependency_roots':[str(r/'dep')]}
                (r/'adoption/prep-inputs.json').write_text(json.dumps(pins))
                P.verify_sources()
                for path in ('include/shadow.h','data/sources/new.txt','dep/shadow.h'):
                    with self.subTest(path=path):
                        (r/path).write_text('extra\n')
                        with self.assertRaises(ValueError):P.verify_sources()
                        (r/path).unlink()
                p=r/'include/fixed.h';p.write_text('changed\n')
                with self.assertRaises(ValueError):P.verify_sources()
                p.unlink()
                with self.assertRaises(ValueError):P.verify_sources()
            finally:P.R=original


if __name__ == '__main__':
    if not os.environ.get('SIMQ_JOB_ID'):raise SystemExit('SimQueue allocation required')
    unittest.main(verbosity=2)
