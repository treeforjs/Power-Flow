#!/usr/bin/python -I
"""Fault-inject the actual run wrapper with temporary inputs and a zero-exit driver.

No scientific runtime file is modified. Execute inside an allocated frozen job.
"""
import importlib.util
import json
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location('coordinator_run', ROOT / 'run_v8.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def main(destination):
    records = []
    with tempfile.TemporaryDirectory(prefix='run-integrity-', dir=ROOT) as directory:
        temporary = Path(directory)
        for fault in ('none', 'dependency', 'source', 'binary', 'missing-dependency'):
            case = temporary / fault
            case.mkdir()
            source, dependency, binary = (case / n for n in ('source.txt', 'library.a', 'driver'))
            for f in (source, dependency, binary):
                f.write_bytes(b'original fixture bytes\n')
            m.OUT = case / 'outputs'
            build = m.OUT / 'build-fixture'
            build.mkdir(parents=True)
            m.select = lambda bindings: None
            m.admission = lambda *args: {'fixture_only': True}
            m.sealed_pins = lambda: {'dependency': m.sha(dependency)}
            m.own_sources = lambda: {'source': m.sha(source)}
            (build / 'build-receipt.json').write_text(json.dumps({
                'status': 'passed', 'binary': str(binary), 'binary_sha256': m.sha(binary),
                'own_sources_sha256': m.own_sources(), 'sealed': m.sealed_pins(), 'bindings': {},
            }))

            def execute(*args, **kwargs):
                if fault == 'missing-dependency':
                    dependency.unlink()
                elif fault != 'none':
                    {'source': source, 'dependency': dependency, 'binary': binary}[fault].write_bytes(b'changed during execution\n')
                return {'returncode': 0, 'fixture_only': True}

            m.execute = execute
            rc = m.run('fixture', 'run', 'result', [])
            receipt = json.loads((m.OUT / 'result/run-receipt.json').read_text())
            expected = 'passed' if fault == 'none' else 'failed'
            passed = (rc == 0) == (fault == 'none') and receipt['status'] == expected
            records.append({'fault': fault, 'passed': passed, 'wrapper_rc': rc, 'receipt': receipt})
    result = {'schema': 'coordinator-run-integrity-faults-v1', 'cases': records,
              'passed': sum(c['passed'] for c in records), 'total': len(records)}
    Path(destination).write_text(json.dumps(result, indent=2) + '\n')
    return 0 if all(c['passed'] for c in records) else 1


if __name__ == '__main__':
    raise SystemExit(main(*sys.argv[1:]))
