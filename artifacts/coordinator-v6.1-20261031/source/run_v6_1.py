#!/usr/bin/python -I
"""Coordinator v6.1 build and run driver (receipted; every compile, link and run inside the frozen launcher).

  build:    core-lease-20261030/both-cores.sh /usr/bin/python -I run_v6_1.py build LABEL
  run:      <leases held by the launcher: any-core.sh, flock on both GPU lock files, or SimQueue>
            /usr/bin/python -I run_v6_1.py run LABEL MODE OUTDIR [driver options...]
            MODE = run | checks | vacuum | profile; OUTDIR absolute (a new directory) or a name under the v6.1 output root.
  baseline: the sealed coordinator v6 binary of record (build q3, sha256 verified against its build receipt) with the
            same inputs and admission, for before/after timings and identity comparisons:
            /usr/bin/python -I run_v6_1.py baseline MODE OUTDIR [driver options...]

Admission: CPU affinity within cores 6-7 with every core lease held, one thread per library, the frozen
kernel release; GPU options additionally need both GPU lock files held. Sealed inputs are hashed before
and after every build and run and are only read.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import fcntl, hashlib, json, os, platform, subprocess, sys, time

P = Path(__file__).resolve().parent
R = Path('/home/treefor/projects/runtimes')
OUT = Path('/home/treefor/projects/.workspace/outputs/pic-water-common-commit-v6.1-20261031')
V6 = R / 'pic-water-common-commit-v6-20261031'
V6_OUT = Path('/home/treefor/projects/.workspace/outputs/pic-water-common-commit-v6-20261031')
V6_LABEL = 'q3'
V6_BINARY_SHA = 'eb1936f1a92557dd1f37ee53b02c6d1b1ec269e340d4e588cfdb03891d1f4b05'
F = R / 'arch-base-20261028/run-frozen'
N6 = R / 'pic-reaction-transactions-20261029/native-full-time-candidate006'
N6_IDENTITY = 'aa4d8f03100e64bed759c739013d91cac1bd66c8a9b2ab6bc3105af7973452c9'
N6_LIB = N6 / 'build/libpic_native_v6.a'
N6_LIB_SHA = '1c2525958457dfe5ff16def07341cd000625e1eefaf38d44a054766efc3d86bc'
H = R / 'hall-coupled-transport-20261029/source'
HALL_DEV = H / 'development'
HD = HALL_DEV / 'hall-canonical-driver-20261029'
HB = HD / 'bindings/heat-v2'
HALL_COMMIT = '8dddec0610982778bf02bdfe697b70a37129c012'
HV3 = R / 'hall-corrective-v3-20261030'
G = R / 'hall-gpu-resident-v1-20261030'
S6 = R / 'pic-electrode-heat-source-v6-20261030'
S61 = R / 'pic-electrode-heat-source-v6.1-20261031'
S61N = S61 / 'native/mfem_electrode_heat_source_v6_1'
V3K = R / 'pic-electrode-heat-feedback-20261029/native/mfem_electrode_heat_source_v3'
A1 = R / 'pic-arrival-owner-v1-20261030'
A11 = R / 'pic-arrival-owner-v1.1-20261030'
CV1 = R / 'pic-water-parent-coverage-v1-20261030'
CV2 = R / 'pic-water-parent-coverage-v2-20261030'
V3 = R / 'pic-water-common-commit-v3-20261030'
I = R / 'pic-water-full-time-subflow-20261029/inputs'
M = R / 'mfem-4.10-pic-20261029'
LIBS = [M / 'install-mpi/lib/libmfem.a', M / 'tpl/hypre/lib/libHYPRE.a', M / 'tpl/metis/lib/libmetis.a', M / 'tpl/GKlib/lib/libGKlib.a']
NVRTC = R / 'pic-cuda-gather-investigation-20261029/tools/nvrtc/nvidia/cuda_nvrtc/lib/libnvrtc.so.12'
NVRTC_INC = R / 'pic-cuda-gather-investigation-20261029/tools/nvrtc/nvidia/cuda_nvrtc/include'
CUDA_INC = R / 'mfem-cuda-benchmark-20261029/cuda-libs/nvidia/cuda_runtime/include'
GPU_UUID = 'GPU-a7b9b458-167f-cd57-fec9-3dea0995b6e9'
LEASE_DIR = Path('/home/treefor/projects/.workspace/outputs/pic-profile')
LEASES = {6: LEASE_DIR / 'testing-core.lock', 7: LEASE_DIR / 'testing-core7.lock'}
GPU_LOCKS = [LEASE_DIR / 'testing-gpu.lock', LEASE_DIR / ('testing-gpu-%s.lock' % GPU_UUID)]
TRACE = HALL_DEV / 'hall-finite-foil-em-20261029/data/shot_17133_current.csv'
TRACE_SHA = '2c54d2dbaec6b6c8fcf2e45b5abd837aa72999f5fd14217f530277ffd5949076'
DATA = R / 'pic-water-exact-count-subflow-20261029/inputs/h2o-gas.json'
DATA_SHA = '840663cdbcf4c39f1ab6f5407f9ecacb232a07c3c48002a4a182c8187792e4bf'
BINDING = Path('/home/treefor/projects/.workspace/outputs/pic-water-exact-count-subflow-20261029/run-binding.json')
OUTER3 = R / 'h2o-collision-data-20261029/data/collisions/h2o-gas-20261029-outer3.json'
OUTER3_SHA = 'c436020d3f2d07902fb059cfc0a269b2e914a3cbb2789300093d3631df896551'
KERNEL_BINDING = CV2 / 'kernel-binding-v2.json'
ENV = {'PATH': '/usr/bin:/bin', 'LC_ALL': 'C', 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1',
       'OMP_THREAD_LIMIT': '1', 'CUDA_VISIBLE_DEVICES': '', 'PYTHONNOUSERSITE': '1'}
GPU_ENV = {k: v for k, v in ENV.items() if k != 'CUDA_VISIBLE_DEVICES'}
COMMON = ['-std=c++20', '-O2', '-ffp-contract=off', '-fno-fast-math', '-ffunction-sections', '-fdata-sections']
WRAPS = ['_ZN6raikou12mfem_backend6GapPic12solve_fieldsEv', '_ZN6raikou12mfem_backend6GapPic7advanceEd',
         '_ZN6raikou12mfem_backend6GapPic13advance_emptyEd',
         '_ZN6raikou12mfem_backend6GapPic16advance_observedEdRKSt8functionIFvRKNS0_11GapParticleEidEE']


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def fail(message):
    raise RuntimeError(message)


def held(lock):
    with open(lock, 'a') as probe:
        try:
            fcntl.flock(probe, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return True
        fcntl.flock(probe, fcntl.LOCK_UN)
        return False


def admission(kind, gpu=False):
    if not __debug__:
        fail('optimized Python guard execution refused')
    cores = os.sched_getaffinity(0)
    if not cores or not cores <= {6, 7}:
        fail('affinity %s not admitted (Power-Flow cores 6-7 only)' % sorted(cores))
    if any(os.getenv(k) != '1' for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']):
        fail('one thread per library required')
    for core in cores:
        if not held(LEASES[core]):
            fail('core %d lease must be held by the launcher; no work started' % core)
    if gpu and not all(held(x) for x in GPU_LOCKS):
        fail('both GPU lock files must be held by the launcher')
    if Path('/proc/sys/kernel/osrelease').read_text().strip() != (R / 'arch-base-20261028/kernel.release').read_text().strip():
        fail('kernel release differs from the frozen runtime')
    return {'kind': kind, 'cpu_affinity': sorted(cores), 'nice': os.getpriority(os.PRIO_PROCESS, 0), 'gpu': gpu,
            'gpu_uuid': GPU_UUID if gpu else None, 'leases_verified_held': [str(LEASES[c]) for c in sorted(cores)] + ([str(x) for x in GPU_LOCKS] if gpu else []),
            'simq_job': os.getenv('SIMQ_JOB_ID'), 'python': platform.python_version(), 'frozen_launcher_sha256': sha(F)}


def tree(paths, suffixes=('.hpp', '.h', '.cpp', '.cu', '.cuh', '.json', '.txt', '.py')):
    out = {}
    for root in paths:
        root = Path(root)
        files = [root] if root.is_file() else sorted(p for p in root.rglob('*') if p.is_file() and p.suffix in suffixes and 'build' not in p.parts
                                                      and '__pycache__' not in p.parts)
        for p in files:
            out[str(p)] = sha(p)
    return out


def sealed_pins():
    rec = {}
    manifest = json.loads((N6 / 'SOURCE-SHA256.json').read_text())
    if manifest['source_identity'] != N6_IDENTITY:
        fail('native v6 source identity differs')
    for name, digest in manifest['files'].items():
        if sha(N6 / name) != digest:
            fail('native v6 sealed file changed: ' + name)
    if sha(N6_LIB) != N6_LIB_SHA:
        fail('native v6 core archive differs')
    rec['native_v6'] = {'identity': N6_IDENTITY, 'files_verified': len(manifest['files']), 'library': str(N6_LIB), 'library_sha256': N6_LIB_SHA}
    head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=H, text=True).strip()
    dirty = subprocess.check_output(['git', 'status', '--porcelain'], cwd=H, text=True)
    if head != HALL_COMMIT or dirty:
        fail('Hall repository is not the clean corrective commit')
    rec['hall_commit'] = head
    for name, digest, path in [('trace', TRACE_SHA, TRACE), ('water_data', DATA_SHA, DATA), ('outer3', OUTER3_SHA, OUTER3)]:
        if sha(path) != digest:
            fail(name + ' differs from its pinned digest')
    binding = json.loads(BINDING.read_text())
    if binding['dataset_sha256'] != DATA_SHA:
        fail('water data binding differs')
    rec['inputs'] = {str(TRACE): TRACE_SHA, str(DATA): DATA_SHA, str(OUTER3): OUTER3_SHA, str(BINDING): sha(BINDING),
                     str(KERNEL_BINDING): sha(KERNEL_BINDING)}
    rec['components_sha256'] = tree([HV3 / 'source', HV3 / 'native-v6-link/hall_v3_flows_adapter.hpp', S6 / 'native', S6 / 'binding',
                                     S61 / 'native', S61 / 'binding',
                                     A11 / 'include', A1 / 'include', CV2 / 'source', CV1 / 'source', V3 / 'source/native_parent_subflow.hpp',
                                     V3 / 'source/common_commit_contract.hpp', V3K, G / 'common', G / 'push', N6 / 'source/adapters',
                                     HB / 'source-kernel'])
    rec['nvrtc'] = {'library': str(NVRTC), 'sha256': sha(NVRTC)}
    return rec


def own_sources():
    # post-processing tools (tools/) are not part of the binary and are not pinned
    return tree([P / 'source', P / 'tests', Path(__file__).resolve()])


def execute(name, cmd, log_dir, timeout, env=None, cwd=None):
    log = log_dir / (name + '.log')
    start = time.monotonic()
    with log.open('w') as handle:
        try:
            r = subprocess.run([str(F), *map(str, cmd)], cwd=cwd or P, env=env or ENV, stdout=handle, stderr=subprocess.STDOUT, timeout=timeout)
            rc = r.returncode
        except subprocess.TimeoutExpired:
            rc = 'timeout'
    out = {'name': name, 'argv': [str(F), *map(str, cmd)], 'returncode': rc, 'elapsed_s': time.monotonic() - start, 'log': str(log),
           'log_sha256': sha(log)}
    print(name, rc, round(out['elapsed_s'], 1), flush=True)
    return out


def includes():
    sys_dirs = [N6 / 'source/include', I, I / 'include', HB, HB / 'source-kernel/include', N6 / 'inputs', HALL_DEV, HV3 / 'source', G, V3K, R,
                S61 / 'binding', S6 / 'binding', A11 / 'include', A1 / 'include', V3 / 'source', CV1 / 'source', CV2 / 'source', N6 / 'source/adapters',
                HV3 / 'native-v6-link', M / 'install-mpi/include', M / 'tpl/hypre/include', M / 'tpl/metis/include']
    out = []
    for d in sys_dirs:
        out += ['-isystem', str(d)]
    return out


def build(label):
    adm = admission('build')
    out = OUT / ('build-' + label)
    out.mkdir(parents=True, exist_ok=False)
    b = P / 'build' / label
    b.mkdir(parents=True, exist_ok=False)
    pins, before = sealed_pins(), own_sources()
    defs = ['-DHEAT_SOURCE_EXACT_COUNT_HEADER="raikou/exact_particle_count.h"', '-DPF_ARRIVAL_EXACT_COUNT_HEADER="raikou/exact_particle_count.h"']
    gpu_inc = ['-I' + str(G), '-I' + str(G / 'push'), '-isystem', str(CUDA_INC), '-isystem', str(NVRTC_INC), '-Wno-c++26-extensions']
    native_inc = ['-isystem', str(N6 / 'source/include'), '-isystem', str(I / 'include'), '-isystem', str(N6 / 'inputs')]
    mfem_inc = [x for d in ('install-mpi/include', 'tpl/hypre/include', 'tpl/metis/include') for x in ('-isystem', str(M / d))]
    ours = ['/usr/bin/mpic++', *COMMON, '-Wall', '-Wextra', '-Wpedantic', '-Werror', '-I' + str(P / 'source'), *defs, *includes()]
    hall = ['/usr/bin/c++', *COMMON, '-DNDEBUG', '-I' + str(HV3 / 'source'), '-I' + str(G), *native_inc, '-isystem', str(HALL_DEV), *mfem_inc]
    gpu_plain = ['/usr/bin/c++', *COMMON, '-DNDEBUG', *gpu_inc]
    kernel = ['/usr/bin/c++', *COMMON, '-DNDEBUG', '-I' + str(HB / 'source-kernel/include')]
    fast = ['/usr/bin/c++', *COMMON, '-DNDEBUG', '-Wall', '-Wextra', '-Werror', '-I' + str(S61N), '-I' + str(HB / 'source-kernel/include')]
    objects = {
        'driver_v6_1': ([*ours, *gpu_inc], P / 'tests/driver_v6_1.cpp'),
        'fast_emission': (fast, S61N / 'fast_emission.cpp'),
        'alloc_hook': (['/usr/bin/c++', *COMMON], P / 'tests/alloc_hook.cpp'),
        'motion_provider_v3': (hall, HV3 / 'source/motion_provider.cpp'),
        'kinetic_boundary': (['/usr/bin/c++', *COMMON, '-DNDEBUG', *native_inc, *mfem_inc], HD / 'kinetic_boundary.cpp'),
        'gpu_backends_v3': ([*hall, *gpu_inc], HV3 / 'source/gpu/gpu_backends.cpp'),
        'hall_gpu_push_v3': ([*gpu_plain, '-I' + str(HV3 / 'source'), *native_inc, '-isystem', str(HALL_DEV), *mfem_inc], HV3 / 'source/gpu/hall_gpu_push_v3.cpp'),
        'hall_gpu_common': (gpu_plain, G / 'common/hall_gpu_common.cpp'),
        'hall_gpu_primitives': (gpu_plain, G / 'common/hall_gpu_primitives.cpp'),
        'surface_emission': (kernel, HB / 'source-kernel/src/surface_emission.cpp'),
        'micro_layer': (kernel, HB / 'source-kernel/src/micro_layer.cpp'),
    }
    jobs = [('compile-' + k, [*flags, '-MD', '-MF', b / (k + '.d'), '-c', src, '-o', b / (k + '.o')]) for k, (flags, src) in objects.items()]
    receipt = {'schema': 'coordinator-v6.1-build-receipt-v1', 'label': label, 'admission': adm, 'sealed': pins, 'environment': ENV,
               'started_unix': time.time(), 'commands': []}
    status = 'passed'
    try:
        with ThreadPoolExecutor(max_workers=len(adm['cpu_affinity'])) as pool:
            results = list(pool.map(lambda j: execute(j[0], j[1], out, 3600), jobs))
        receipt['commands'] += results
        for r in results:
            if r['returncode'] != 0:
                fail(r['name'] + ' failed:\n' + Path(r['log']).read_text()[-12000:])
        exe = b / 'driver_v6_1'
        link = ['/usr/bin/mpic++', *['-Wl,--wrap=' + s for s in WRAPS], '-Wl,--gc-sections', *[b / (k + '.o') for k in objects], N6_LIB, *LIBS,
                '-lcrypto', '-ldl', '-o', exe]
        r = execute('link', link, out, 1800)
        receipt['commands'].append(r)
        if r['returncode'] != 0:
            fail('link failed:\n' + Path(r['log']).read_text()[-12000:])
        deps = set()
        for k in objects:
            for line in (b / (k + '.d')).read_text().replace('\\\n', ' ').splitlines():
                if ':' in line:
                    deps |= {str(Path(x).resolve()) for x in line.split(':', 1)[1].split()}
        receipt['dependency_sha256'] = {d: sha(d) for d in sorted(deps)}
        json_hpp = {d: sha(d) for d in deps if d.endswith('/json.hpp')}
        receipt['single_nlohmann_json_hpp'] = len(set(json_hpp.values())) == 1
        heat = sorted(d for d in deps if d.endswith('electrode_heat_feedback.hpp'))
        receipt['heat_v2_engine_headers'] = heat
        if len(heat) != 1 or not receipt['single_nlohmann_json_hpp']:
            fail('duplicate heat-v2 engine or json headers in the link: %s %s' % (heat, json_hpp))
        rogue = sorted(d for d in deps if 'hall-field-correction-20261030' in d or 'hall-corrective-v2-20261030' in d or
                       d.endswith('hall-canonical-driver-20261029/motion_receipt.hpp') or 'hall-motion-v2/motion_receipt.hpp' in d or
                       'native-full-time-candidate003' in d or 'native-full-time-candidate005' in d)
        receipt['sealed_superseded_headers_in_link'] = rogue
        if rogue:
            fail('superseded component headers leaked into the link: ' + ', '.join(rogue))
        receipt['objects_sha256'] = {k: sha(b / (k + '.o')) for k in objects}
        receipt['binary'] = str(exe)
        receipt['binary_sha256'] = sha(exe)
        receipt['own_sources_sha256'] = before
        if own_sources() != before or sealed_pins() != pins:
            fail('sources or sealed inputs changed during the build')
        receipt['inputs_unchanged_after_build'] = True
    except Exception as exc:
        status = 'failed'
        receipt['error'] = str(exc)
    receipt['status'] = status
    receipt['finished_unix'] = time.time()
    (out / 'build-receipt.json').write_text(json.dumps(receipt, indent=1) + '\n')
    print(json.dumps({'build': str(out), 'status': status, 'error': receipt.get('error', '')[-6000:]}))
    return 0 if status == 'passed' else 1


def run(label, mode, outdir, args):
    gpu = any(a == 'gpu' for a in args)
    adm = admission('run', gpu)
    br = json.loads((OUT / ('build-' + label) / 'build-receipt.json').read_text())
    if br['status'] != 'passed' or sha(br['binary']) != br['binary_sha256']:
        fail('build not usable')
    if own_sources() != br['own_sources_sha256']:
        fail('coordinator sources changed after the build')
    pins = sealed_pins()
    if pins != br['sealed']:
        fail('sealed inputs differ from the build')
    out = Path(outdir) if outdir.startswith('/') else OUT / outdir
    out.mkdir(parents=True, exist_ok=False)
    extra = ['--uuid', GPU_UUID, '--nvrtc', str(NVRTC)] if gpu else []
    cmd = [br['binary'], mode, TRACE, DATA, BINDING, OUTER3, KERNEL_BINDING, out, N6_IDENTITY, *args, *extra]
    receipt = {'schema': 'coordinator-v6.1-run-receipt-v1', 'label': label, 'mode': mode, 'args': args, 'admission': adm, 'sealed': pins,
               'binary_sha256': br['binary_sha256'], 'build_receipt': str(OUT / ('build-' + label) / 'build-receipt.json'), 'started_unix': time.time()}
    (out / 'run-receipt.json').write_text(json.dumps(receipt, indent=1) + '\n')
    r = execute(mode, cmd, out, 7 * 24 * 3600, env=GPU_ENV if gpu else ENV)
    receipt['command'] = r
    receipt['finished_unix'] = time.time()
    receipt['sealed_unchanged_after_run'] = sealed_pins() == pins
    receipt['files_sha256'] = {f.name: sha(f) for f in sorted(out.iterdir()) if f.is_file() and f.name != 'run-receipt.json'}
    (out / 'run-receipt.json').write_text(json.dumps(receipt, indent=1) + '\n')
    print(json.dumps({'run': str(out), 'returncode': r['returncode']}))
    return 0 if r['returncode'] == 0 else 1


def baseline(mode, outdir, args):
    gpu = any(a == 'gpu' for a in args)
    adm = admission('baseline', gpu)
    receipt_path = V6_OUT / ('build-' + V6_LABEL) / 'build-receipt.json'
    br = json.loads(receipt_path.read_text())
    binary = V6 / 'build' / V6_LABEL / 'driver_v6'
    if br['status'] != 'passed' or br['binary_sha256'] != V6_BINARY_SHA or sha(binary) != V6_BINARY_SHA:
        fail('sealed coordinator v6 binary of record differs from its receipt')
    pins = sealed_pins()
    out = Path(outdir) if outdir.startswith('/') else OUT / outdir
    out.mkdir(parents=True, exist_ok=False)
    extra = ['--uuid', GPU_UUID, '--nvrtc', str(NVRTC)] if gpu else []
    cmd = [binary, mode, TRACE, DATA, BINDING, OUTER3, KERNEL_BINDING, out, N6_IDENTITY, *args, *extra]
    receipt = {'schema': 'coordinator-v6.1-baseline-receipt-v1', 'baseline': 'coordinator v6 build ' + V6_LABEL, 'mode': mode, 'args': args,
               'admission': adm, 'sealed': pins, 'binary': str(binary), 'binary_sha256': V6_BINARY_SHA, 'build_receipt': str(receipt_path),
               'started_unix': time.time()}
    (out / 'run-receipt.json').write_text(json.dumps(receipt, indent=1) + '\n')
    r = execute(mode, cmd, out, 7 * 24 * 3600, env=GPU_ENV if gpu else ENV, cwd=V6)
    receipt['command'] = r
    receipt['finished_unix'] = time.time()
    receipt['sealed_unchanged_after_run'] = sealed_pins() == pins
    receipt['files_sha256'] = {f.name: sha(f) for f in sorted(out.iterdir()) if f.is_file() and f.name != 'run-receipt.json'}
    (out / 'run-receipt.json').write_text(json.dumps(receipt, indent=1) + '\n')
    print(json.dumps({'baseline': str(out), 'returncode': r['returncode']}))
    return 0 if r['returncode'] == 0 else 1


if __name__ == '__main__':
    a = sys.argv
    if len(a) >= 3 and a[1] == 'build' and '/' not in a[2]:
        raise SystemExit(build(a[2]))
    if len(a) >= 5 and a[1] == 'run' and a[3] in ('run', 'checks', 'vacuum', 'profile'):
        raise SystemExit(run(a[2], a[3], a[4], a[5:]))
    if len(a) >= 4 and a[1] == 'baseline' and a[2] in ('run', 'checks', 'vacuum'):
        raise SystemExit(baseline(a[2], a[3], a[4:]))
    raise SystemExit(__doc__)
