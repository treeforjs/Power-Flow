#!/usr/bin/python -I
"""Build and run the Hall corrective v3 qualification.

  build:     core-lease-20261030/both-cores.sh /usr/bin/python -I run_hall_v3.py build LABEL [ubsan]
  cpu:       core-lease-20261030/any-core.sh --wait /usr/bin/python -I run_hall_v3.py cpu LABEL
  multistep: core-lease-20261030/any-core.sh --wait /usr/bin/python -I run_hall_v3.py multistep LABEL STEPS BIRTH_EVERY BIRTHS
  gpu:       flock <testing-gpu.lock> flock <testing-gpu-UUID.lock> core-lease-20261030/any-core.sh --wait \
                 /usr/bin/python -I run_hall_v3.py gpu LABEL CPU_RUN_DIR

Every compile and run step executes inside arch-base-20261028/run-frozen. Outputs:
/home/treefor/projects/.workspace/outputs/hall-corrective-v3-20261030/LABEL/. Sealed inputs (Hall 8dddec0,
hall-gpu-resident-v1 common/push, native candidate003, Hall corrective v2 sources) are hashed before and after
every stage and are only read.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import fcntl, hashlib, json, os, platform, subprocess, sys, time

P = Path(__file__).resolve().parent
R = Path('/home/treefor/projects/runtimes')
OUT = Path('/home/treefor/projects/.workspace/outputs/hall-corrective-v3-20261030')
N = R / 'pic-reaction-transactions-20261029/native-full-time-candidate003'
H = R / 'hall-coupled-transport-20261029/source'
HD = H / 'development/hall-canonical-driver-20261029'
I = R / 'pic-water-full-time-subflow-20261029/inputs'
M = R / 'mfem-4.10-pic-20261029'
G = R / 'hall-gpu-resident-v1-20261030'
HV2 = R / 'hall-corrective-v2-20261030'
F = R / 'arch-base-20261028/run-frozen'
LEASE_DIR = Path('/home/treefor/projects/.workspace/outputs/pic-profile')
LEASES = {6: LEASE_DIR / 'testing-core.lock', 7: LEASE_DIR / 'testing-core7.lock'}
GPU_UUID = 'GPU-a7b9b458-167f-cd57-fec9-3dea0995b6e9'
GPU_LOCKS = [LEASE_DIR / 'testing-gpu.lock', LEASE_DIR / ('testing-gpu-%s.lock' % GPU_UUID)]
NVRTC = R / 'pic-cuda-gather-investigation-20261029/tools/nvrtc/nvidia/cuda_nvrtc/lib/libnvrtc.so.12'
NVRTC_INC = R / 'pic-cuda-gather-investigation-20261029/tools/nvrtc/nvidia/cuda_nvrtc/include'
CUDA_INC = R / 'mfem-cuda-benchmark-20261029/cuda-libs/nvidia/cuda_runtime/include'
HALL_COMMIT = '8dddec0610982778bf02bdfe697b70a37129c012'
IDENTITY = 'c045c42d2f57082e3fb5e12172c440d2470e681a38df9ac406bd7e063609a08d'
NATIVE_LIB = N / 'build/libpic_native_candidate003.a'
NATIVE_LIB_SHA = 'd2568af9fac13f8aa056ae4a9a50bb2a07eee05d7218d1c4b79aef21317d8a20'
TRACE = H / 'development/hall-finite-foil-em-20261029/data/shot_17133_current.csv'
TRACE_SHA = '2c54d2dbaec6b6c8fcf2e45b5abd837aa72999f5fd14217f530277ffd5949076'
LIBS = [M / 'install-mpi/lib/libmfem.a', M / 'tpl/hypre/lib/libHYPRE.a', M / 'tpl/metis/lib/libmetis.a', M / 'tpl/GKlib/lib/libGKlib.a']
GPU_SEALED = [G / 'common/hall_gpu_common.cpp', G / 'common/hall_gpu_common.hpp', G / 'common/hall_gpu_primitives.cpp',
              G / 'common/hall_gpu_primitives.hpp', G / 'common/device/deterministic_primitives.cuh', G / 'push/hall_gpu_push.cpp',
              G / 'push/hall_gpu_push.hpp', G / 'push/push_core.hpp', G / 'push/push_abi.hpp', G / 'push/push_bodies.hpp',
              G / 'push/push_kernels.cu', G / 'push/push_fold.hpp', G / 'push/hall_adapter.hpp']
V2_EVIDENCE = Path('/home/treefor/projects/.workspace/outputs/hall-corrective-v2-20261030/b5/cpu-182242')
ENV = {'PATH': '/usr/bin:/bin', 'LC_ALL': 'C', 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1',
       'PYTHONNOUSERSITE': '1'}


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def fail(m):
    raise RuntimeError(m)


def held(lock):
    with open(lock, 'a') as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return True
        fcntl.flock(handle, fcntl.LOCK_UN)
        return False


def admission(stage):
    cores = os.sched_getaffinity(0)
    if not cores <= {6, 7} or os.getpriority(os.PRIO_PROCESS, 0) < 15:
        fail('only cores 6/7 at nice >= 15 are admitted; affinity %s' % sorted(cores))
    if stage != 'build' and len(cores) != 1:
        fail('run stages need exactly one leased core')
    for c in cores:
        if not held(LEASES[c]):
            fail('core %d lease is not held by the launching wrapper' % c)
    if stage == 'gpu' and not all(held(lk) for lk in GPU_LOCKS):
        fail('both GPU lock files must be held')
    if any(os.getenv(k) != '1' for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']):
        fail('one thread per numerical library required')
    if Path('/proc/sys/kernel/osrelease').read_text().strip() != (R / 'arch-base-20261028/kernel.release').read_text().strip():
        fail('kernel release differs from the frozen runtime')
    return {'stage': stage, 'cpu_affinity': sorted(cores), 'nice': os.getpriority(os.PRIO_PROCESS, 0),
            'core_leases_held': {str(c): held(LEASES[c]) for c in LEASES}, 'gpu_locks_held': [held(lk) for lk in GPU_LOCKS],
            'PF_CORE': os.getenv('PF_CORE'), 'python': sys.executable, 'python_version': platform.python_version(),
            'frozen_launcher_sha256': sha(F), 'time_unix': time.time()}


def sealed_pins():
    rec = {}
    manifest = json.loads((N / 'SOURCE-SHA256.json').read_text())
    if manifest['source_identity'] != IDENTITY:
        fail('native source identity differs')
    for name, digest in manifest['files'].items():
        if sha(N / name) != digest:
            fail('native sealed file changed: ' + name)
    if sha(NATIVE_LIB) != NATIVE_LIB_SHA:
        fail('native candidate003 library differs')
    rec['native'] = {'source_identity': IDENTITY, 'library_sha256': NATIVE_LIB_SHA, 'files_verified': len(manifest['files'])}
    head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=H, text=True).strip()
    dirty = subprocess.check_output(['git', 'status', '--porcelain'], cwd=H, text=True)
    if head != HALL_COMMIT or dirty:
        fail('Hall repository is not the clean corrective commit')
    rec['hall'] = {'commit': head, 'git_clean': True}
    rec['gpu_sealed_sources'] = {str(p): sha(p) for p in GPU_SEALED}
    rec['hall_v2_sources'] = {str(p): sha(p) for p in sorted((HV2 / 'source').rglob('*')) if p.is_file()}
    if sha(TRACE) != TRACE_SHA:
        fail('C19 trace differs')
    rec['trace_sha256'] = TRACE_SHA
    rec['nvrtc'] = {'library': str(NVRTC), 'sha256': sha(NVRTC)}
    return rec


def own_sources():
    files = sorted(p for d in (P / 'source', P / 'tests') for p in d.rglob('*') if p.is_file() and p.suffix in ('.hpp', '.cpp', '.cu', '.py'))
    return {str(p): sha(p) for p in files + [Path(__file__).resolve()]}


def execute(name, cmd, log_dir, timeout, cwd=None, env=None):
    log = log_dir / (name + '.log')
    start = time.monotonic()
    with log.open('w') as h:
        try:
            rc = subprocess.run([str(F), *map(str, cmd)], cwd=cwd or P, env=env or ENV, stdout=h, stderr=subprocess.STDOUT, timeout=timeout).returncode
        except subprocess.TimeoutExpired:
            rc = 'timeout'
    r = {'name': name, 'argv': [str(F), *map(str, cmd)], 'returncode': rc, 'elapsed_s': time.monotonic() - start,
         'log': str(log), 'log_sha256': sha(log)}
    print(name, rc, round(r['elapsed_s'], 1), flush=True)
    return r


def isystem():
    out = []
    for d in [N / 'source/include', I, I / 'include', N / 'inputs', H / 'development', M / 'install-mpi/include',
              M / 'tpl/hypre/include', M / 'tpl/metis/include']:
        out += ['-isystem', str(d)]
    return out


COMMON = ['-std=c++20', '-O2', '-DNDEBUG', '-ffp-contract=off', '-fno-fast-math', '-ffunction-sections', '-fdata-sections']
UBSAN = ['-std=c++20', '-O1', '-g', '-DNDEBUG', '-ffp-contract=off', '-fno-fast-math', '-fsanitize=undefined', '-fno-sanitize-recover=all',
         '-fno-omit-frame-pointer']


def build(label, variant='normal'):
    adm = admission('build')
    out = OUT / label / 'build'
    out.mkdir(parents=True, exist_ok=False)
    b = P / 'build' / label
    b.mkdir(parents=True, exist_ok=False)
    pins = sealed_pins()
    before = own_sources()
    base = UBSAN if variant == 'ubsan' else COMMON
    gpu_inc = ['-I' + str(G), '-I' + str(G / 'push'), '-isystem', str(CUDA_INC), '-isystem', str(NVRTC_INC), '-Wno-c++26-extensions']
    v3 = ['/usr/bin/c++', *base, '-I' + str(P / 'source'), '-I' + str(P / 'tests'), '-I' + str(G), *isystem()]
    v2 = ['/usr/bin/c++', *COMMON, '-I' + str(HV2 / 'source'), '-I' + str(P / 'tests'), *isystem()]
    sealed = ['/usr/bin/c++', *base, *isystem()]
    objs = {'motion_provider_v3': (v3, P / 'source/motion_provider.cpp'),
            'kinetic_boundary': (sealed, HD / 'kinetic_boundary.cpp'),
            'test_hall_v3': (v3, P / 'tests/test_hall_v3.cpp')}
    links = [('test_hall_v3', ['test_hall_v3', 'motion_provider_v3', 'kinetic_boundary'], [])]
    if variant == 'normal':
        objs.update({'motion_provider_v2': (v2, HV2 / 'source/motion_provider.cpp'),
                     'v2_physical_dump': (v2, P / 'tests/v2_physical_dump.cpp'),
                     'gpu_backends_v3': (v3 + gpu_inc, P / 'source/gpu/gpu_backends.cpp'),
                     'hall_gpu_push_v3': (['/usr/bin/c++', *COMMON, '-I' + str(P / 'source'), *gpu_inc], P / 'source/gpu/hall_gpu_push_v3.cpp'),
                     'hall_gpu_common': (['/usr/bin/c++', *COMMON, *gpu_inc], G / 'common/hall_gpu_common.cpp'),
                     'hall_gpu_primitives': (['/usr/bin/c++', *COMMON, *gpu_inc], G / 'common/hall_gpu_primitives.cpp'),
                     'test_hall_v3_gpu': (v3 + gpu_inc, P / 'tests/test_hall_v3_gpu.cpp')})
        links += [('v2_physical_dump', ['v2_physical_dump', 'motion_provider_v2', 'kinetic_boundary'], []),
                  ('test_hall_v3_gpu', ['test_hall_v3_gpu', 'gpu_backends_v3', 'hall_gpu_push_v3', 'hall_gpu_common', 'hall_gpu_primitives',
                                        'motion_provider_v3', 'kinetic_boundary'], ['-ldl'])]
    receipt = {'schema': 'hall-corrective-v3-build-receipt-v1', 'label': label, 'variant': variant, 'admission': adm, 'sealed': pins, 'commands': []}
    status = 'passed'
    try:
        jobs = [('compile-' + k, [*flags, '-MD', '-MF', b / (k + '.d'), '-c', src, '-o', b / (k + '.o')]) for k, (flags, src) in objs.items()]
        workers = 2 if len(adm['cpu_affinity']) == 2 else 1
        with ThreadPoolExecutor(max_workers=workers) as pool:
            results = list(pool.map(lambda j: execute(j[0], j[1], out, 3600), jobs))
        receipt['commands'] += results
        for r in results:
            if r['returncode'] != 0:
                fail(r['name'] + ' failed:\n' + Path(r['log']).read_text()[-8000:])
        link = ['/usr/bin/mpic++', '-Wl,--gc-sections'] + (['-fsanitize=undefined'] if variant == 'ubsan' else [])
        tail = [NATIVE_LIB, *LIBS, '-lcrypto']
        for name, parts, extra in links:
            r = execute('link-' + name, [*link, *[b / (p + '.o') for p in parts], *tail, *extra, '-o', b / name], out, 1800)
            receipt['commands'].append(r)
            if r['returncode'] != 0:
                fail(r['name'] + ' failed:\n' + Path(r['log']).read_text()[-8000:])
        deps = set()
        for k in objs:
            text = (b / (k + '.d')).read_text().replace('\\\n', ' ')
            for line in text.splitlines():
                if ':' in line:
                    deps |= {str(Path(x).resolve()) for x in line.split(':', 1)[1].split()}
        receipt['dependency_sha256'] = {d: sha(d) for d in sorted(deps)}
        v3deps = {d for k in objs if k.endswith('v3') or k == 'test_hall_v3' for d in (b / (k + '.d')).read_text().replace('\\\n', ' ').split()
                  if d.endswith('.hpp') or d.endswith('.h')}
        leaked = sorted(d for d in v3deps if 'hall-field-correction-20261030' in d or d.endswith('hall-canonical-driver-20261029/motion_receipt.hpp')
                        or 'hall-corrective-v2-20261030' in d)
        receipt['v3_objects_read_sealed_corrective_v2_or_sealed_motion_header'] = leaked
        if leaked:
            fail('a v3 object includes a sealed corrective, a v2 or the sealed motion receipt header: ' + ', '.join(leaked))
        receipt['binaries_sha256'] = {n: sha(b / n) for n, _, _ in links}
        receipt['build_dir'] = str(b)
        receipt['own_sources_sha256'] = before
        if own_sources() != before or sealed_pins() != pins:
            fail('sources or sealed inputs changed during the build')
    except Exception as e:
        status = 'failed'
        receipt['error'] = str(e)
    receipt['status'] = status
    (out / 'build-receipt.json').write_text(json.dumps(receipt, indent=1) + '\n')
    print(json.dumps({'build': str(out), 'status': status, 'error': receipt.get('error', '')[:4000]}))
    return 0 if status == 'passed' else 1


def load_build(label):
    br = json.loads((OUT / label / 'build/build-receipt.json').read_text())
    if br['status'] != 'passed':
        fail('build did not pass')
    for n, d in br['binaries_sha256'].items():
        if sha(Path(br['build_dir']) / n) != d:
            fail('binary differs from its build receipt: ' + n)
    if own_sources() != br['own_sources_sha256']:
        fail('v3 sources changed after the build')
    return br


def compare_rows(a, b, keys):
    ref = {r['fixture']: r for r in a}
    rows, ok = [], True
    for r in b:
        s = ref.get(r['fixture'])
        same = s is not None and all(s.get(k) == r.get(k) for k in keys)
        ok = ok and same
        rows.append({'fixture': r['fixture'], 'entry': r.get('entry'), 'identical': same})
    return ok and len(a) == len(b), rows


def run_cpu(label, ubsan_label=None):
    adm = admission('cpu')
    br = load_build(label)
    pins = sealed_pins()
    b = Path(br['build_dir'])
    out = OUT / label / ('cpu-' + time.strftime('%H%M%S'))
    out.mkdir(parents=True)
    cmds = [execute('components', [b / 'test_hall_v3', 'components', TRACE, out / 'components'], out, 7200),
            execute('parity-v3', [b / 'test_hall_v3', 'parity-dump', TRACE, out / 'parity-v3.json'], out, 3600),
            execute('parity-v2', [b / 'v2_physical_dump', TRACE, out / 'parity-v2.json'], out, 3600),
            execute('motion-dump-v3', [b / 'test_hall_v3', 'motion-dump', TRACE, out / 'motion-dump-v3.json'], out, 3600),
            execute('sink-dump-cpu', [b / 'test_hall_v3', 'sink-dump', TRACE, out / 'sink-dump-cpu.json'], out, 3600)]
    if ubsan_label:
        ub = load_build(ubsan_label)
        cmds.append(execute('components-ubsan', [Path(ub['build_dir']) / 'test_hall_v3', 'components', TRACE, out / 'components-ubsan'], out, 14400))
    summary = {'schema': 'hall-corrective-v3-cpu-receipt-v1', 'label': label, 'admission': adm, 'sealed': pins, 'commands': cmds,
               'build_receipt': str(OUT / label / 'build/build-receipt.json')}
    for key, sub in (('components', 'components'), ('components_ubsan', 'components-ubsan')):
        comp = out / sub / 'components.json'
        if comp.exists():
            c = json.loads(comp.read_text())
            summary[key] = {'failed': c['failed'], 'passed': sum(1 for x in c['cases'] if x['passed']), 'cases': len(c['cases']),
                            'names': [{'group': x['group'], 'name': x['name'], 'passed': x['passed']} for x in c['cases']]}
    if ubsan_label:
        a, u = out / 'components/components.json', out / 'components-ubsan/components.json'
        if a.exists() and u.exists():
            ca, cu = json.loads(a.read_text()), json.loads(u.read_text())
            summary['ubsan_same_outcomes'] = [x['passed'] for x in ca['cases']] == [x['passed'] for x in cu['cases']]
            log = (out / 'components-ubsan.log').read_text()
            summary['ubsan_runtime_errors'] = log.count('runtime error')
    try:
        p2, p3 = json.loads((out / 'parity-v2.json').read_text()), json.loads((out / 'parity-v3.json').read_text())
        keys = ('physical_sha256', 'motion_record_text_sha256', 'motion_seal_v2', 'motion_record_v3_sha256', 'motion_arrays_sha256', 'iterations',
                'combined_balance_J', 'field_heat_balance_J', 'field_receipt', 'refusal')
        rows = []
        for x, y in zip(p2, p3):
            rows.append({'scenario': x['scenario'], 'step': x['step'], 'identical': all(x.get(k) == y.get(k) for k in keys),
                         'differences': [k for k in keys if x.get(k) != y.get(k)]})
        summary['v2_physical_parity'] = {'identical': len(p2) == len(p3) and all(r['identical'] for r in rows), 'steps': rows}
    except Exception as e:
        summary['v2_physical_parity'] = {'identical': False, 'error': str(e)}
    try:
        keys = ('record_sha256', 'seal_sha256', 'record_v3_sha256', 'after_sha256', 'refusal', 'after_refusal')
        a = json.loads((V2_EVIDENCE / 'motion-dump-v2.json').read_text())
        bb = json.loads((out / 'motion-dump-v3.json').read_text())
        ok = len(a) == len(bb) and all(x['fixture'] == y['fixture'] and x['entry'] == y['entry'] and all(x.get(k) == y.get(k) for k in keys) for x, y in zip(a, bb))
        s = json.loads((V2_EVIDENCE / 'motion-dump-sealed.json').read_text())
        sk = ('record_sha256', 'seal_sha256', 'after_sha256', 'refusal', 'after_refusal')
        sealed_ok = len(s) == len([r for r in bb if r['entry'] == 'MotionAssembly::cpu']) and \
            all(x['fixture'] == y['fixture'] and x.get(k) == y.get(k) for x, y in zip(s, [r for r in bb if r['entry'] == 'MotionAssembly::cpu']) for k in sk)
        summary['v2_motion_records_identical'] = {'v2_dump': str(V2_EVIDENCE / 'motion-dump-v2.json'), 'identical_to_v2': ok,
                                                  'identical_to_sealed_8dddec0': sealed_ok, 'fixtures': len(bb)}
    except Exception as e:
        summary['v2_motion_records_identical'] = {'identical_to_v2': False, 'error': str(e)}
    summary['all_passed'] = all(c['returncode'] == 0 for c in cmds) and summary.get('v2_physical_parity', {}).get('identical', False) and \
        summary.get('v2_motion_records_identical', {}).get('identical_to_v2', False)
    if sealed_pins() != pins:
        fail('sealed inputs changed during the run')
    (out / 'cpu-receipt.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps({'run': str(out), 'all_passed': summary['all_passed']}))
    return 0 if summary['all_passed'] else 1


def run_multistep(label, steps, every, births):
    adm = admission('cpu')
    br = load_build(label)
    pins = sealed_pins()
    out = OUT / label / ('multistep-' + time.strftime('%H%M%S'))
    out.mkdir(parents=True)
    r = execute('multistep', [Path(br['build_dir']) / 'test_hall_v3', 'multistep', TRACE, out, steps, every, births], out, 6 * 3600)
    summary = {'schema': 'hall-corrective-v3-multistep-receipt-v1', 'label': label, 'admission': adm, 'sealed': pins, 'command': r}
    f = out / 'multistep-summary.json'
    if f.exists():
        summary['summary'] = json.loads(f.read_text())
        summary['ledger_sha256'] = sha(out / 'step-ledger.jsonl')
    (out / 'multistep-receipt.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps({'run': str(out), 'returncode': r['returncode']}))
    return 0 if r['returncode'] == 0 else 1


def run_gpu(label, cpu_dir):
    adm = admission('gpu')
    br = load_build(label)
    pins = sealed_pins()
    b = Path(br['build_dir'])
    out = OUT / label / ('gpu-' + time.strftime('%H%M%S'))
    out.mkdir(parents=True)
    env = dict(ENV)
    cmds = [execute('gpu-a', [b / 'test_hall_v3_gpu', TRACE, NVRTC, GPU_UUID, out / 'a'], out, 7200, env=env),
            execute('gpu-b', [b / 'test_hall_v3_gpu', TRACE, NVRTC, GPU_UUID, out / 'b'], out, 7200, env=env)]
    summary = {'schema': 'hall-corrective-v3-gpu-receipt-v1', 'label': label, 'admission': adm, 'sealed': pins, 'commands': cmds}
    for k in ('a', 'b'):
        f = out / k / 'gpu.json'
        if f.exists():
            g = json.loads(f.read_text())
            summary[k] = {'failed': g.get('failed'), 'error': g.get('error'), 'cases': [{'group': x['group'], 'name': x['name'], 'passed': x['passed']} for x in g.get('cases', [])],
                          'device': g.get('device'), 'teardown': g.get('teardown')}
    keys = ('record_text_sha256', 'record_v3_sha256', 'seal', 'arrays_sha256', 'after_sha256', 'sinks', 'sink_section_sha256', 'sink_kinds', 'motion_refusal',
            'candidate_field_sha256', 'candidate_heat_sha256', 'candidate_ledgers_sha256', 'iterations', 'combined_balance_J', 'energy_residual_J', 'refusal')
    try:
        cpu = json.loads((Path(cpu_dir) / 'sink-dump-cpu.json').read_text())
        for k in ('a', 'b'):
            gpu = json.loads((out / k / 'sink-dump-gpu.json').read_text())
            same, rows = compare_rows(cpu, gpu, keys)
            summary['sink_parity_' + k] = {'identical': same, 'fixtures': rows}
        v2 = json.loads((Path(cpu_dir) / 'motion-dump-v3.json').read_text())
        gv = json.loads((out / 'a' / 'gpu-motion-dump.json').read_text())
        ref = [r for r in v2 if r['entry'] == 'MotionAssembly::cpu']
        mk = ('record_sha256', 'seal_sha256', 'record_v3_sha256', 'after_sha256', 'refusal', 'after_refusal')
        summary['v2_fixture_parity'] = {'identical': len(ref) == len(gv) and all(all(x.get(q) == y.get(q) for q in mk) for x, y in zip(ref, gv)), 'fixtures': len(gv)}
    except Exception as e:
        summary['compare_error'] = str(e)
    summary['all_passed'] = all(c['returncode'] == 0 for c in cmds) and summary.get('sink_parity_a', {}).get('identical', False) and \
        summary.get('sink_parity_b', {}).get('identical', False) and summary.get('v2_fixture_parity', {}).get('identical', False)
    (out / 'gpu-receipt.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps({'run': str(out), 'all_passed': summary['all_passed']}))
    return 0 if summary['all_passed'] else 1


if __name__ == '__main__':
    a = sys.argv
    if len(a) < 3 or a[1] not in ('build', 'cpu', 'gpu', 'multistep') or '/' in a[2]:
        raise SystemExit(__doc__)
    if a[1] == 'build':
        raise SystemExit(build(a[2], a[3] if len(a) > 3 else 'normal'))
    if a[1] == 'cpu':
        raise SystemExit(run_cpu(a[2], a[3] if len(a) > 3 else None))
    if a[1] == 'multistep':
        raise SystemExit(run_multistep(a[2], a[3], a[4], a[5]))
    raise SystemExit(run_gpu(a[2], a[3]))
