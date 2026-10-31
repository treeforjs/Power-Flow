#!/usr/bin/python -I
"""Build and run the Hall corrective v2 component qualification.

  build: core-lease-20261030/both-cores.sh /usr/bin/python -I run_hall_v2.py build LABEL
         (a single leased core from any-core.sh is accepted when the second core is leased elsewhere)
  cpu:   core-lease-20261030/any-core.sh --wait /usr/bin/python -I run_hall_v2.py cpu LABEL
  gpu:   flock <testing-gpu.lock> flock <testing-gpu-UUID.lock> core-lease-20261030/any-core.sh --wait \
             /usr/bin/python -I run_hall_v2.py gpu LABEL

Every compile and run step executes inside arch-base-20261028/run-frozen. Outputs:
/home/treefor/projects/.workspace/outputs/hall-corrective-v2-20261030/LABEL/. Sealed inputs (Hall 8dddec0,
hall-gpu-resident-v1 common/push, native candidate003) are hashed before and after every stage.
"""
from pathlib import Path
import fcntl, hashlib, json, os, platform, subprocess, sys, time

P = Path(__file__).resolve().parent
R = Path('/home/treefor/projects/runtimes')
OUT = Path('/home/treefor/projects/.workspace/outputs/hall-corrective-v2-20261030')
N = R / 'pic-reaction-transactions-20261029/native-full-time-candidate003'
H = R / 'hall-coupled-transport-20261029/source'
HC = H / 'development/hall-field-correction-20261030'
HD = H / 'development/hall-canonical-driver-20261029'
I = R / 'pic-water-full-time-subflow-20261029/inputs'
M = R / 'mfem-4.10-pic-20261029'
G = R / 'hall-gpu-resident-v1-20261030'
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
              G / 'push/push_kernels.cu', G / 'push/push_fold.hpp', G / 'push/hall_adapter.hpp', G / 'field/src/kernels.hpp',
              G / 'field/src/gpu_cg.hpp']
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
    if stage == 'gpu' and not all(held(l) for l in GPU_LOCKS):
        fail('both GPU lock files must be held')
    if any(os.getenv(k) != '1' for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']):
        fail('one thread per numerical library required')
    if Path('/proc/sys/kernel/osrelease').read_text().strip() != (R / 'arch-base-20261028/kernel.release').read_text().strip():
        fail('kernel release differs from the frozen runtime')
    return {'stage': stage, 'cpu_affinity': sorted(cores), 'nice': os.getpriority(os.PRIO_PROCESS, 0),
            'core_leases_held': {str(c): held(LEASES[c]) for c in LEASES}, 'gpu_locks_held': [held(l) for l in GPU_LOCKS],
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
    if sha(TRACE) != TRACE_SHA:
        fail('C19 trace differs')
    rec['trace_sha256'] = TRACE_SHA
    rec['nvrtc'] = {'library': str(NVRTC), 'sha256': sha(NVRTC)}
    return rec


def own_sources():
    files = sorted(p for d in (P / 'source', P / 'tests') for p in d.rglob('*') if p.is_file() and p.suffix in ('.hpp', '.cpp'))
    return {str(p): sha(p) for p in files + [Path(__file__).resolve()]}


def execute(name, cmd, log_dir, timeout, cwd=None):
    log = log_dir / (name + '.log')
    start = time.monotonic()
    with log.open('w') as h:
        try:
            rc = subprocess.run([str(F), *map(str, cmd)], cwd=cwd or P, env=ENV, stdout=h, stderr=subprocess.STDOUT, timeout=timeout).returncode
        except subprocess.TimeoutExpired:
            rc = 'timeout'
    r = {'name': name, 'argv': [str(F), *map(str, cmd)], 'returncode': rc, 'elapsed_s': time.monotonic() - start,
         'log': str(log), 'log_sha256': sha(log)}
    print(name, rc, round(r['elapsed_s'], 1), flush=True)
    return r


def isystem(sealed_corrective):
    dirs = [N / 'source/include', I, I / 'include', N / 'inputs']
    if sealed_corrective:
        dirs.append(HC)
    dirs += [H / 'development', M / 'install-mpi/include', M / 'tpl/hypre/include', M / 'tpl/metis/include']
    out = []
    for d in dirs:
        out += ['-isystem', str(d)]
    return out


COMMON = ['-std=c++20', '-O2', '-DNDEBUG', '-ffp-contract=off', '-fno-fast-math', '-ffunction-sections', '-fdata-sections']


def build(label):
    adm = admission('build')
    out = OUT / label / 'build'
    out.mkdir(parents=True, exist_ok=False)
    b = P / 'build' / label
    b.mkdir(parents=True, exist_ok=False)
    pins = sealed_pins()
    before = own_sources()
    v2 = ['/usr/bin/c++', *COMMON, '-I' + str(P / 'source'), '-I' + str(P / 'tests'), *isystem(False)]
    sealed = ['/usr/bin/c++', *COMMON, '-I' + str(P / 'tests'), *isystem(True)]
    gpu_flags = ['-I' + str(G), '-isystem', str(CUDA_INC), '-isystem', str(NVRTC_INC), '-Wno-c++26-extensions']
    objs = {
        'motion_provider_v2': (v2, P / 'source/motion_provider.cpp'),
        'kinetic_boundary': (sealed, HD / 'kinetic_boundary.cpp'),
        'sealed_native_arena_motion': (sealed, HD / 'native_arena_motion.cpp'),
        'test_hall_v2': (v2, P / 'tests/test_hall_v2.cpp'),
        'repro_hall_node_plane_v2': (v2, P / 'tests/repro_hall_node_plane_v2.cpp'),
        'sealed_motion_reference': (sealed, P / 'tests/sealed_motion_reference.cpp'),
        'gpu_backends': (v2 + gpu_flags, P / 'source/gpu/gpu_backends.cpp'),
        'hall_gpu_common': (['/usr/bin/c++', *COMMON, *gpu_flags], G / 'common/hall_gpu_common.cpp'),
        'hall_gpu_primitives': (['/usr/bin/c++', *COMMON, *gpu_flags], G / 'common/hall_gpu_primitives.cpp'),
        'hall_gpu_push': (['/usr/bin/c++', *COMMON, *gpu_flags], G / 'push/hall_gpu_push.cpp'),
        'test_hall_v2_gpu': (v2 + gpu_flags, P / 'tests/test_hall_v2_gpu.cpp'),
    }
    receipt = {'schema': 'hall-corrective-v2-build-receipt-v1', 'label': label, 'admission': adm, 'sealed': pins, 'commands': []}
    status = 'passed'
    try:
        for k, (flags, src) in objs.items():
            r = execute('compile-' + k, [*flags, '-MD', '-MF', b / (k + '.d'), '-c', src, '-o', b / (k + '.o')], out, 1800)
            receipt['commands'].append(r)
            if r['returncode'] != 0:
                fail(r['name'] + ' failed:\n' + Path(r['log']).read_text()[-5000:])
        link = ['/usr/bin/mpic++', '-Wl,--gc-sections']
        tail = [NATIVE_LIB, *LIBS, '-lcrypto']
        for name, parts, extra in [
                ('test_hall_v2', ['test_hall_v2', 'motion_provider_v2', 'kinetic_boundary'], []),
                ('repro_hall_node_plane_v2', ['repro_hall_node_plane_v2', 'motion_provider_v2', 'kinetic_boundary'], []),
                ('sealed_motion_reference', ['sealed_motion_reference', 'sealed_native_arena_motion', 'kinetic_boundary'], []),
                ('test_hall_v2_gpu', ['test_hall_v2_gpu', 'gpu_backends', 'hall_gpu_common', 'hall_gpu_primitives', 'hall_gpu_push',
                                      'motion_provider_v2', 'kinetic_boundary'], ['-ldl'])]:
            r = execute('link-' + name, [*link, *[b / (p + '.o') for p in parts], *tail, *extra, '-o', b / name], out, 1800)
            receipt['commands'].append(r)
            if r['returncode'] != 0:
                fail(r['name'] + ' failed:\n' + Path(r['log']).read_text()[-5000:])
        deps = set()
        for k in objs:
            text = (b / (k + '.d')).read_text().replace('\\\n', ' ')
            for line in text.splitlines():
                if ':' in line:
                    deps |= {str(Path(x).resolve()) for x in line.split(':', 1)[1].split()}
        receipt['dependency_sha256'] = {d: sha(d) for d in sorted(deps)}
        leaked = [d for d in deps if 'hall-field-correction-20261030' in d and 'sealed' not in d]
        v2deps = {d for k in objs if k not in ('sealed_motion_reference', 'sealed_native_arena_motion', 'kinetic_boundary')
                  for d in (b / (k + '.d')).read_text().replace('\\\n', ' ').split() if d.endswith('.hpp') or d.endswith('.h')}
        receipt['v2_objects_read_sealed_corrective_or_sealed_motion_header'] = sorted(
            d for d in v2deps if 'hall-field-correction-20261030' in d or d.endswith('hall-canonical-driver-20261029/motion_receipt.hpp')
            or 'hall-motion-v2/motion_receipt.hpp' in d)
        if receipt['v2_objects_read_sealed_corrective_or_sealed_motion_header']:
            fail('a v2 object includes the sealed corrective or sealed motion receipt header')
        del leaked
        receipt['binaries_sha256'] = {n: sha(b / n) for n in ('test_hall_v2', 'repro_hall_node_plane_v2', 'sealed_motion_reference', 'test_hall_v2_gpu')}
        receipt['build_dir'] = str(b)
        receipt['own_sources_sha256'] = before
        if own_sources() != before or sealed_pins() != pins:
            fail('sources or sealed inputs changed during the build')
    except Exception as e:
        status = 'failed'
        receipt['error'] = str(e)
    receipt['status'] = status
    (out / 'build-receipt.json').write_text(json.dumps(receipt, indent=1) + '\n')
    print(json.dumps({'build': str(out), 'status': status, 'error': receipt.get('error')}))
    return 0 if status == 'passed' else 1


def load_build(label):
    br = json.loads((OUT / label / 'build/build-receipt.json').read_text())
    if br['status'] != 'passed':
        fail('build did not pass')
    for n, d in br['binaries_sha256'].items():
        if sha(Path(br['build_dir']) / n) != d:
            fail('binary differs from its build receipt: ' + n)
    if own_sources() != br['own_sources_sha256']:
        fail('v2 sources changed after the build')
    return br


def compare_dumps(a_path, b_path, key_entry):
    a = json.loads(Path(a_path).read_text())
    b = json.loads(Path(b_path).read_text())
    ref = {r['fixture']: r for r in a}
    rows, ok = [], True
    for r in b:
        if key_entry and r.get('entry') != key_entry:
            continue
        s = ref.get(r['fixture'])
        same = s is not None and all(s.get(k) == r.get(k) for k in ('record_sha256', 'seal_sha256', 'after_sha256', 'refusal', 'after_refusal'))
        ok = ok and same
        rows.append({'fixture': r['fixture'], 'entry': r.get('entry'), 'identical': same, 'refusal': r.get('refusal'),
                     'reference_wall_s': s.get('wall_s') if s else None, 'wall_s': r.get('wall_s')})
    return ok, rows


def run_cpu(label):
    adm = admission('cpu')
    br = load_build(label)
    pins = sealed_pins()
    b = Path(br['build_dir'])
    out = OUT / label / ('cpu-' + time.strftime('%H%M%S'))
    out.mkdir(parents=True)
    cmds = [execute('components', [b / 'test_hall_v2', 'components', TRACE, out / 'components'], out, 3600),
            execute('repro-node-plane-v2', [b / 'repro_hall_node_plane_v2', TRACE], out, 1800),
            execute('motion-dump-v2', [b / 'test_hall_v2', 'motion-dump', TRACE, out / 'motion-dump-v2.json'], out, 3600),
            execute('motion-dump-sealed', [b / 'sealed_motion_reference', TRACE, out / 'motion-dump-sealed.json'], out, 3600)]
    summary = {'schema': 'hall-corrective-v2-cpu-receipt-v1', 'label': label, 'admission': adm, 'sealed': pins, 'commands': cmds,
               'build_receipt': str(OUT / label / 'build/build-receipt.json')}
    comp = out / 'components/components.json'
    if comp.exists():
        c = json.loads(comp.read_text())
        summary['components'] = {'failed': c['failed'], 'cases': [{'group': x['group'], 'name': x['name'], 'passed': x['passed']} for x in c['cases']],
                                 'field_state_hash_timing_s': c.get('field_state_hash_timing_s'), 'configuration_sha256': c.get('configuration_sha256')}
    if (out / 'motion-dump-v2.json').exists() and (out / 'motion-dump-sealed.json').exists():
        rows = []
        ok = True
        for entry in ('MotionAssembly::cpu', 'prepare_canonical_motion'):
            same, r = compare_dumps(out / 'motion-dump-sealed.json', out / 'motion-dump-v2.json', entry)
            ok = ok and same
            rows += r
        summary['motion_v2_cpu_equals_sealed'] = {'identical': ok, 'fixtures': rows}
    summary['all_passed'] = all(c['returncode'] == 0 for c in cmds) and summary.get('motion_v2_cpu_equals_sealed', {}).get('identical', False)
    if sealed_pins() != pins:
        fail('sealed inputs changed during the run')
    (out / 'cpu-receipt.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps({'run': str(out), 'all_passed': summary['all_passed']}))
    return 0 if summary['all_passed'] else 1


def run_gpu(label, cpu_dir=None):
    adm = admission('gpu')
    br = load_build(label)
    pins = sealed_pins()
    b = Path(br['build_dir'])
    out = OUT / label / ('gpu-' + time.strftime('%H%M%S'))
    out.mkdir(parents=True)
    cmds = [execute('gpu-a', [b / 'test_hall_v2_gpu', TRACE, NVRTC, GPU_UUID, out / 'a'], out, 3600),
            execute('gpu-b', [b / 'test_hall_v2_gpu', TRACE, NVRTC, GPU_UUID, out / 'b'], out, 3600)]
    summary = {'schema': 'hall-corrective-v2-gpu-receipt-v1', 'label': label, 'admission': adm, 'sealed': pins, 'commands': cmds}
    for k in ('a', 'b'):
        f = out / k / 'gpu.json'
        if f.exists():
            g = json.loads(f.read_text())
            summary[k] = {'failed': g.get('failed'), 'error': g.get('error'), 'cases': [{'group': x['group'], 'name': x['name'], 'passed': x['passed']} for x in g.get('cases', [])],
                          'device': g.get('device'), 'teardown': g.get('teardown')}
    try:
        ga = json.loads((out / 'a/gpu.json').read_text())
        gb = json.loads((out / 'b/gpu.json').read_text())
        def hashes(g):
            return [c['evidence'] for c in g['cases'] if c['group'] == 'closed-loop' and 'three runs' in c['name']][0]
        summary['fresh_process_gpu_candidate_identical'] = [r['candidate_field_sha256'] for r in hashes(ga)] == [r['candidate_field_sha256'] for r in hashes(gb)]
    except Exception as e:
        summary['fresh_process_gpu_candidate_identical'] = False
        summary['compare_error'] = str(e)
    if cpu_dir:
        same, rows = compare_dumps(Path(cpu_dir) / 'motion-dump-sealed.json', out / 'a/gpu-motion-dump.json', None)
        summary['gpu_motion_equals_sealed_reference'] = {'identical': same, 'fixtures': rows}
    summary['all_passed'] = all(c['returncode'] == 0 for c in cmds) and summary['fresh_process_gpu_candidate_identical'] and \
        (not cpu_dir or summary['gpu_motion_equals_sealed_reference']['identical'])
    (out / 'gpu-receipt.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps({'run': str(out), 'all_passed': summary['all_passed']}))
    return 0 if summary['all_passed'] else 1


if __name__ == '__main__':
    if len(sys.argv) < 3 or sys.argv[1] not in ('build', 'cpu', 'gpu') or '/' in sys.argv[2]:
        raise SystemExit('usage: run_hall_v2.py build|cpu|gpu LABEL [CPU_RUN_DIR]')
    if sys.argv[1] == 'build':
        raise SystemExit(build(sys.argv[2]))
    if sys.argv[1] == 'cpu':
        raise SystemExit(run_cpu(sys.argv[2]))
    raise SystemExit(run_gpu(sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else None))
