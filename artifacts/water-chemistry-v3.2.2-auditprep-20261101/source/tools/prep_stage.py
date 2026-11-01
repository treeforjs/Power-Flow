"""Bounded frozen preparation checks. Receipts never qualify a physical model."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

R = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('prep_validation', R/'tools/swarm_validation.py')
V = importlib.util.module_from_spec(spec); spec.loader.exec_module(V)


def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as stream:
        for b in iter(lambda: stream.read(1048576), b''): h.update(b)
    return h.hexdigest()


def verify_sources():
    pin = json.loads((R/'adoption/prep-inputs.json').read_text())
    current = {str(p.relative_to(R)) for d in ('include','device','tests','tools')
               for p in (R/d).rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    current.update(('run_checks.py','run_prep.py'))
    if current != set(pin['source_files']): raise ValueError('source file set changed')
    current_inputs = {str(p.relative_to(R)) for d in ('data/sources','data/extracted')
                      for p in (R/d).rglob('*') if p.is_file()}
    if current_inputs != {p for p in pin['input_files'] if p.startswith(('data/sources/','data/extracted/'))}:
        raise ValueError('input file set changed')
    for directory in pin['dependency_roots']:
        current_deps = {str(p) for p in Path(directory).rglob('*') if p.is_file()}
        if current_deps != {p for p in pin['dependencies'] if p.startswith(directory+'/')}:
            raise ValueError('dependency file set changed: '+directory)
    for section in ('source_files', 'input_files', 'dependencies'):
        for path, expected in pin[section].items():
            if sha(R/path) != expected: raise ValueError('pinned '+section+' changed: '+path)
    return sha(R/'adoption/prep-inputs.json')


def generated():
    return {str(p.relative_to(R)): sha(p) for p in sorted((R/'data/tables').glob('*')) if p.is_file()}


def prior(name, binding):
    q = json.loads((R/'evidence'/('prep-'+name+'.json')).read_text())
    if not q['passed'] or q['binding_sha256'] != binding: raise ValueError('missing/changed prerequisite: '+name)
    if name == 'regenerate' and q['generated_files'] != generated(): raise ValueError('regenerated tables changed')
    return q


def main(stage):
    if not os.environ.get('SIMQ_JOB_ID'): raise SystemExit('SimQueue allocation required')
    if len(os.sched_getaffinity(0)) != 1: raise SystemExit('exactly one allocated logical CPU required')
    if os.environ.get('PF_BUILD_ISA_FLAGS','').strip(): raise SystemExit('unqualified build flags refused; PF_BUILD_ISA_FLAGS must be empty')
    os.chdir(R)
    before = verify_sources()
    dest = R/'evidence'/('prep-'+stage+'.json')
    if dest.exists(): raise SystemExit('existing stage receipt preserved; use a new reviewed preparation version')
    receipt = {'schema':'chemistry-audit-preparation-v1', 'stage':stage, 'binding_sha256':before,
               'job_id':os.environ['SIMQ_JOB_ID'], 'cpus':sorted(os.sched_getaffinity(0)),
               'production_qualified':False, 'physics_model_validated':False, 'build_isa_flags':'', 'passed':False, 'entries':[]}
    def run(name, argv, extra_env=None):
        env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
        env['PF_BUILD_ISA_FLAGS'] = ''
        if stage != 'gpu': env['CUDA_VISIBLE_DEVICES'] = ''
        env.update(extra_env or {})
        log = R/'evidence'/('prep-'+stage+'-'+name+'.log')
        with log.open('w') as stream:
            result = subprocess.run(argv, stdout=stream, stderr=subprocess.STDOUT, env=env, check=False)
        receipt['entries'].append({'name':name,'argv':argv,'exit':result.returncode,'log_sha256':sha(log)})
        print(name, 'exit', result.returncode, flush=True)
        if result.returncode: raise RuntimeError(name+' failed; see '+str(log))
    tables = 'data/tables/chemistry-v3.2'
    try:
        if stage == 'regenerate':
            if generated(): raise ValueError('successor generated directory is not empty')
            run('spec', ['/usr/bin/python3','-I','tools/make_production_spec.py',tables+'-spec.json'])
            if sha(tables+'-spec.json') != sha('reference/v3.2/'+tables+'-spec.json'):
                raise ValueError('scientific specification changed; no automatic model adoption')
            run('tables', ['/usr/bin/python3','-I','tools/build_tables.py',tables+'-spec.json',tables])
            run('policy', ['/usr/bin/python3','-I','tools/make_coverage_policy.py',tables+'-manifest.json',tables+'-coverage-policy.json'])
            run('report', ['/usr/bin/python3','-I','tools/report_processes.py',tables+'-manifest.json','data/tables/PROCESSES.txt'])
            run('strict-policy', ['/usr/bin/python3','-I','tests/test_audit_policy.py'])
            receipt['generated_files'] = generated()
        elif stage in ('cpu', 'sanitized'):
            prior('regenerate', before)
            mode = 'normal' if stage == 'cpu' else 'sanitized'
            flags = {} if stage == 'cpu' else {'ASAN_OPTIONS':'detect_leaks=1:halt_on_error=1','UBSAN_OPTIONS':'halt_on_error=1:print_stacktrace=1'}
            outputs = {}
            for name, source in (('audit_regressions','tests/test_audit_regressions.cpp'),('test_core','tests/test_core.cpp')):
                binary = 'build/bin/'+name+'-'+stage
                run('build-'+name,['sh','tools/dev_build.sh',mode,source,binary]); outputs[binary]=sha(binary)
                argv = [binary,tables+'.blob']
                if name == 'test_core': argv += [tables+'-checks.json','evidence/core-'+stage,'10000']
                run(name,argv,flags)
                if sha(binary) != outputs[binary]: raise ValueError('tested executable changed')
            receipt['binaries'] = outputs
        elif stage == 'gpu-build':
            prior('regenerate', before)
            run('build',['sh','tools/dev_build_gpu.sh','build/bin/gpu_parity'])
            receipt['binaries'] = {'build/bin/gpu_parity':sha('build/bin/gpu_parity')}
        elif stage == 'gpu':
            if not os.environ.get('SIMQ_GPUS'): raise ValueError('GPU allocation required')
            prior('regenerate', before); build = prior('gpu-build', before)
            if build['binaries'] != {'build/bin/gpu_parity':sha('build/bin/gpu_parity')}: raise ValueError('GPU executable differs from build')
            receipt['binaries'] = build['binaries']
            for i in (1,2):
                run('process'+str(i),['build/bin/gpu_parity',tables+'.blob','evidence/gpu-process'+str(i),
                    '/home/treefor/projects/runtimes/pic-cuda-gather-investigation-20261029/tools/nvrtc/nvidia/cuda_nvrtc/lib/libnvrtc.so.12',
                    'GPU-a7b9b458-167f-cd57-fec9-3dea0995b6e9','1000','1'])
            a=json.loads((R/'evidence/gpu-process1/gpu-parity.json').read_text())
            b=json.loads((R/'evidence/gpu-process2/gpu-parity.json').read_text())
            if [c.get('record_digest_gpu') for c in a['cases']] != [c.get('record_digest_gpu') for c in b['cases']]:
                raise ValueError('fresh-process device digests differ')
            if receipt['binaries'] != {'build/bin/gpu_parity':sha('build/bin/gpu_parity')}: raise ValueError('GPU executable changed')
        else: raise ValueError('unsupported preparation stage')
        if before != verify_sources(): raise ValueError('source/input changed during preparation')
        if stage != 'regenerate': prior('regenerate',before)
        receipt['passed'] = True
    except BaseException as exc:
        receipt['failure'] = type(exc).__name__+': '+str(exc)
        V.atomic_json(dest,receipt)
        raise
    V.atomic_json(dest,receipt)
    print('PREPARATION PASSED; production and model qualification remain false',flush=True)


if __name__ == '__main__': main(sys.argv[1])
