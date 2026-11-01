#!/usr/bin/python -I
"""Source v7.3 surface-birth and inherited qualification. Only live one-core SimQueue jobs;
all compiler and numerical processes use the frozen runtime. No GPU work.
Commands: build LABEL [normal|ubsan|asan], run LABEL [SANITIZER_LABEL ...].
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import fcntl, hashlib, json, os, platform, subprocess, sys, time

P = Path(__file__).resolve().parent
R = Path('/home/treefor/projects/runtimes')
OUT = Path('/home/treefor/projects/.workspace/outputs/source-owner-v7.3-surface-birth-20261101')
S61 = R / 'pic-electrode-heat-source-v6.1-20261031'
S6 = R / 'pic-electrode-heat-source-v6-20261030'
N = R / 'pic-reaction-transactions-20261029/native-full-time-candidate003'
H = R / 'hall-coupled-transport-20261029/source'
HD = H / 'development/hall-canonical-driver-20261029'
HB = P / 'vendor/heat-v2'
I = R / 'pic-water-full-time-subflow-20261029/inputs'
M = R / 'mfem-4.10-pic-20261029'
G = R / 'hall-gpu-resident-v1-20261030'
HV3 = R / 'hall-corrective-v3-20261030'
S5 = R / 'pic-electrode-heat-source-v5-20261030/native/mfem_electrode_heat_source_v5'
V3K = P / 'vendor/v3'
A1 = R / 'pic-arrival-owner-v1-20261030'
A11 = R / 'pic-arrival-owner-v1.1-20261030'
F = R / 'arch-base-20261028/run-frozen'
N84 = R / 'pic-reaction-transactions-20261029/native-full-time-v8.4-water-scope-20261101'
CHEM = R / 'pic-water-chemistry-v3.1.1-20261031'
BLOB = CHEM / 'data/tables/chemistry-v3.1.blob'
LEGACY_HB = HD / 'bindings/heat-v2'
V7N = P / 'native/mfem_electrode_heat_source_v7'
V61N = P / 'vendor/fast'
LEASE_DIR = Path('/home/treefor/projects/.workspace/outputs/pic-profile')
LEASES = {0: LEASE_DIR / 'testing-core0.lock', 1: LEASE_DIR / 'testing-core1.lock', 3: LEASE_DIR / 'testing-core3.lock',
          6: LEASE_DIR / 'testing-core.lock', 7: LEASE_DIR / 'testing-core7.lock'}
RETURN_FLAG = LEASE_DIR / 'RETURN-CORES-0-3'
GPU_LOCKS = [LEASE_DIR / 'testing-gpu.lock', LEASE_DIR / 'testing-gpu-GPU-a7b9b458-167f-cd57-fec9-3dea0995b6e9.lock']
GPU_UUID = 'GPU-a7b9b458-167f-cd57-fec9-3dea0995b6e9'
NVRTC_DIR = R / 'pic-cuda-gather-investigation-20261029/tools/nvrtc/nvidia/cuda_nvrtc'
NVRTC = NVRTC_DIR / 'lib/libnvrtc.so.12'
CUDA_INC = R / 'mfem-cuda-benchmark-20261029/cuda-libs/nvidia/cuda_runtime/include'
NATIVE_LIB = N / 'build/libpic_native_candidate003.a'
LIBS = [M / 'install-mpi/lib/libmfem.a', M / 'tpl/hypre/lib/libHYPRE.a', M / 'tpl/metis/lib/libmetis.a', M / 'tpl/GKlib/lib/libGKlib.a']
TRACE = H / 'development/hall-finite-foil-em-20261029/data/shot_17133_current.csv'
ENV = {'PATH': '/usr/bin:/bin', 'LC_ALL': 'C', 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1', 'PYTHONNOUSERSITE': '1', 'PF_V8_CHEM': 'v3.1'}
VARIANTS = {
    'normal': ['-std=c++20', '-O2', '-DNDEBUG', '-ffp-contract=off', '-fno-fast-math'],
    'ubsan': ['-std=c++20', '-O1', '-g', '-DNDEBUG', '-ffp-contract=off', '-fno-fast-math', '-fsanitize=undefined', '-fno-sanitize-recover=all', '-fno-omit-frame-pointer'],
    'asan': ['-std=c++20', '-O1', '-g', '-DNDEBUG', '-ffp-contract=off', '-fno-fast-math', '-fsanitize=address', '-fno-omit-frame-pointer'],
}
LINK_SAN = {'normal': [], 'ubsan': ['-fsanitize=undefined'], 'asan': ['-fsanitize=address']}


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def held(lock):
    with open(lock, 'a') as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return True
        fcntl.flock(handle, fcntl.LOCK_UN)
        return False


def admission(stage):
    if not __debug__ or stage == 'gpu':
        raise RuntimeError('CPU-only guarded source mass qualification required')
    cores = os.sched_getaffinity(0)
    unit, assigned = os.getenv('SIMQ_UNIT', ''), os.getenv('SIMQ_CPU_LIST', '')
    if not os.getenv('SIMQ_JOB_ID', '').isdigit() or not unit or unit not in Path('/proc/self/cgroup').read_text():
        raise RuntimeError('live SimQueue job cgroup required')
    if len(cores) != 1 or assigned != str(next(iter(cores))) or int(os.getenv('SIMQ_MEM', '0')) <= 0:
        raise RuntimeError('one allocated CPU and explicit memory cap required')
    core, borrow = next(iter(cores)), os.getenv('SIMQ_BORROW') == '1'
    if core not in (4, 5, 6, 7) or (core in (4, 5) and not borrow):
        raise RuntimeError('Power-Flow owner core or explicitly borrowed core 4/5 required')
    lease = LEASE_DIR / ('borrow-core%d.lock' % core) if borrow else LEASES[core]
    if not held(lease):
        raise RuntimeError('scheduled CPU lease is not held')
    if Path('/proc/sys/kernel/osrelease').read_text().strip() != (R / 'arch-base-20261028/kernel.release').read_text().strip():
        raise RuntimeError('kernel differs from frozen runtime')
    mem = {x.split(':')[0]: int(x.split()[1])*1024 for x in Path('/proc/meminfo').read_text().splitlines() if ':' in x and len(x.split())>=2 and x.split()[1].isdigit()}
    if mem['MemAvailable'] < 5*1024**3:
        raise RuntimeError('less than 5 GiB physical memory available; no build/run started')
    return {'stage': stage, 'cpu_affinity': sorted(cores), 'simq_job': os.environ['SIMQ_JOB_ID'],
            'simq_unit': unit, 'borrowed': borrow, 'memory_ceiling_MiB': int(os.environ['SIMQ_MEM']),
            'mem_available_bytes': mem['MemAvailable'], 'lease_verified_held': str(lease),
            'nice': os.getpriority(os.PRIO_PROCESS, 0), 'frozen_launcher_sha256': sha(F), 'time_unix': time.time()}


def sealed():
    dirs = (S5, V3K, HB, A1 / 'include', A11 / 'include', S6 / 'native', S6 / 'binding', S6 / 'tests', S61 / 'native', S61 / 'binding', G / 'common')
    files = [p for d in dirs for p in sorted(d.rglob('*')) if p.is_file() and p.suffix in ('.hpp', '.h', '.cpp', '.cuh', '.inc', '.json', '.txt')]
    return {str(p): sha(p) for p in files + [NATIVE_LIB, N84 / 'build/libpic_native_v8.a', N84 / 'build/kinetic-boundary.o', BLOB, F, R/'arch-base-20261028/root-files.json', R/'arch-base-20261028/root-files.sha256', *LIBS]}


def own():
    files = sorted(p for p in P.rglob('*') if p.is_file() and p.suffix in ('.hpp', '.h', '.cpp', '.cu', '.inc', '.py', '.json', '.txt') and 'build' not in p.parts)
    return {str(p): sha(p) for p in files}


def execute(name, cmd, log_dir, timeout, env=None):
    log = log_dir / (name + '.log')
    t = time.monotonic()
    with log.open('w') as h:
        try:
            rc = subprocess.run([str(F), *map(str, cmd)], cwd=P, env=env or ENV, stdout=h, stderr=subprocess.STDOUT, timeout=timeout).returncode
        except subprocess.TimeoutExpired:
            rc = 'timeout'
    print(name, rc, round(time.monotonic() - t, 1), flush=True)
    return {'name': name, 'argv': [str(F), *map(str, cmd)], 'returncode': rc, 'elapsed_s': time.monotonic() - t, 'log': str(log), 'log_sha256': sha(log)}


def owner_flags(base, compat=False):
    inc = (['-I' + str(P / 'compat-v6')] if compat else []) + ['-I' + str(R), '-I' + str(P / 'binding'), '-I' + str(V7N), '-I' + str(V61N),
           '-I' + str(HB / 'source-kernel/include'), '-I' + str(HV3 / 'source'), '-I' + str(G), '-I' + str(HB), '-I' + str(V3K),
           '-I' + str(A11 / 'include'), '-I' + str(A1 / 'include'),
           '-DHEAT_SOURCE_EXACT_COUNT_HEADER="raikou/exact_particle_count.h"', '-DPF_ARRIVAL_EXACT_COUNT_HEADER="raikou/exact_particle_count.h"']
    for d in [N / 'source/include', I, I / 'include', N / 'inputs', H / 'development', M / 'install-mpi/include', M / 'tpl/hypre/include', M / 'tpl/metis/include']:
        inc += ['-isystem', str(d)]
    return ['/usr/bin/c++', *base, *inc]


def light_flags(base):
    return ['/usr/bin/c++', *base, '-I' + str(V7N), '-I' + str(V61N), '-I' + str(HB / 'source-kernel/include'), '-I' + str(A1 / 'include'),
            '-isystem', str(N / 'source/include'), '-isystem', str(N / 'inputs')]


def gpu_flags(base):
    return ['/usr/bin/c++', *base, '-Wno-c++26-extensions', '-I' + str(G), '-I' + str(V7N), '-I' + str(A1 / 'include'), '-isystem', str(N / 'source/include'),
            '-isystem', str(N / 'inputs'),
            '-isystem', str(CUDA_INC), '-isystem', str(NVRTC_DIR / 'include')]


EXES = {
    'test_source_v7': ('owner', ['test_source_v7']),
    'test_source_v6_on_v7': ('owner', ['test_source_v6_on_v7']),
    'test_surface_v7': ('light', ['test_surface_v7']),
    'test_wall_gpu_v7': ('gpu', ['test_wall_gpu_v7', 'hall_gpu_common']),
}


def dependency_path(logical):
    # Compiler .d paths refer to the frozen filesystem, not host /usr.
    path = Path(logical)
    if not path.is_absolute():
        path = P / path
    if path.parts[1] in ('usr', 'lib', 'lib64', 'bin', 'sbin', 'etc'):
        return R / 'arch-base-20261028/root' / path.relative_to('/')
    return path.resolve()


def build(label, variant):
    adm = admission('build')
    if variant not in VARIANTS:
        raise SystemExit('variant must be one of ' + ', '.join(VARIANTS))
    out = OUT / label / 'build'
    out.mkdir(parents=True, exist_ok=False)
    b = P / 'build' / label
    b.mkdir(parents=True, exist_ok=False)
    pins, before = sealed(), own()
    base = VARIANTS[variant]
    se_flags = ['/usr/bin/c++', *base, '-I' + str(HB / 'source-kernel/include')]
    native_flags = json.loads((P / 'native-registry-fixture-flags.json').read_text())['flags']
    native_flags += base[1:]
    native_flags += ['-I'+str(N84/'tests'), '-I'+str(P/'adapters'), '-I'+str(HB/'source-kernel/include'), '-I'+str(V61N), '-I'+str(V7N), '-I'+str(CHEM/'include/pf_chem_v3')]
    objs = {'test_registry_mass': (native_flags, P / 'tests/test_registry_mass.cpp'),
            'legacy_kernel_new': (se_flags, P / 'tests/legacy_kernel_bytes.cpp'),
            'legacy_kernel_old': (['/usr/bin/c++', *base, '-I'+str(LEGACY_HB/'source-kernel/include')], P / 'tests/legacy_kernel_bytes.cpp'),
            'surface_emission_old': (['/usr/bin/c++', *base, '-I'+str(LEGACY_HB/'source-kernel/include')], LEGACY_HB / 'source-kernel/src/surface_emission.cpp'),
            'micro_layer_old': (['/usr/bin/c++', *base, '-I'+str(LEGACY_HB/'source-kernel/include')], LEGACY_HB / 'source-kernel/src/micro_layer.cpp'),
            'test_source_v7': (owner_flags(base), P / 'tests/test_source_v7.cpp'),
            'test_source_v6_on_v7': (owner_flags(base, True), S6 / 'tests/test_source_v6.cpp'),
            'test_surface_v7': (light_flags(base), P / 'tests/test_surface_v7.cpp'),
            'fast_emission': ([*se_flags, '-I' + str(V61N)], V61N / 'fast_emission.cpp'),
            'emission_hot_v7': (se_flags, V7N / 'emission_hot_v7.cpp'),
            'motion_provider_v3': (owner_flags(base), HV3 / 'source/motion_provider.cpp'),
            'kinetic_boundary': (owner_flags(base), HD / 'kinetic_boundary.cpp'),
            'surface_emission': (se_flags, HB / 'source-kernel/src/surface_emission.cpp'),
            'micro_layer': (se_flags, HB / 'source-kernel/src/micro_layer.cpp')}
    exes = ['test_registry_mass','legacy_kernel_new','legacy_kernel_old','test_source_v7', 'test_source_v6_on_v7', 'test_surface_v7']
    if False:  # Wall GPU remains unmodified and outside this CPU-only patch.
        objs['test_wall_gpu_v7'] = (gpu_flags(base), P / 'tests/test_wall_gpu_v7.cpp')
        objs['hall_gpu_common'] = (gpu_flags(base), G / 'common/hall_gpu_common.cpp')
        exes.append('test_wall_gpu_v7')
    receipt = {'schema': 'source-v7-build-receipt-v1', 'label': label, 'variant': variant, 'admission': adm, 'sealed_sha256': pins, 'commands': []}
    status = 'passed'
    try:
        jobs = [('compile-' + k, [*f, '-MD', '-MF', b / (k + '.d'), '-c', src, '-o', b / (k + '.o')]) for k, (f, src) in objs.items()]
        with ThreadPoolExecutor(max_workers=max(1, len(adm['cpu_affinity']))) as pool:
            results = list(pool.map(lambda j: execute(j[0], j[1], out, 3600), jobs))
        receipt['commands'] += results
        for r in results:
            if r['returncode'] != 0:
                raise RuntimeError(r['name'] + ' failed:\n' + Path(r['log']).read_text()[-8000:])
        sealed_objs = [b / (k + '.o') for k in ('fast_emission', 'emission_hot_v7', 'motion_provider_v3', 'kinetic_boundary', 'surface_emission', 'micro_layer')]
        links = {'test_registry_mass': ['/usr/bin/mpic++', '-Wl,--gc-sections', *LINK_SAN[variant], b/'test_registry_mass.o', b/'fast_emission.o', b/'emission_hot_v7.o', b/'surface_emission.o', b/'micro_layer.o', N84/'build/kinetic-boundary.o', N84/'build/libpic_native_v8.a', *LIBS, '-lcrypto', '-ldl'],
                 'legacy_kernel_new': ['/usr/bin/c++', *LINK_SAN[variant], b/'legacy_kernel_new.o', b/'surface_emission.o', b/'micro_layer.o'],
                 'legacy_kernel_old': ['/usr/bin/c++', *LINK_SAN[variant], b/'legacy_kernel_old.o', b/'surface_emission_old.o', b/'micro_layer_old.o'],
                 'test_source_v7': ['/usr/bin/mpic++', '-Wl,--gc-sections', *LINK_SAN[variant], b / 'test_source_v7.o', *sealed_objs, NATIVE_LIB, *LIBS, '-lcrypto'],
                 'test_source_v6_on_v7': ['/usr/bin/mpic++', '-Wl,--gc-sections', *LINK_SAN[variant], b / 'test_source_v6_on_v7.o', *sealed_objs, NATIVE_LIB, *LIBS, '-lcrypto'],
                 'test_surface_v7': ['/usr/bin/c++', *LINK_SAN[variant], b / 'test_surface_v7.o', b / 'fast_emission.o', b / 'emission_hot_v7.o', b / 'surface_emission.o', '-lcrypto'],
                 'test_wall_gpu_v7': ['/usr/bin/c++', b / 'test_wall_gpu_v7.o', b / 'hall_gpu_common.o', '-ldl', '-lcrypto']}
        for exe in exes:
            r = execute('link-' + exe, [*links[exe], '-o', b / exe], out, 1800)
            receipt['commands'].append(r)
            if r['returncode'] != 0:
                raise RuntimeError('link failed:\n' + Path(r['log']).read_text()[-8000:])
        deps = set()
        for k in objs:
            for line in (b / (k + '.d')).read_text().replace('\\\n', ' ').splitlines():
                if ':' in line:
                    deps |= set(line.split(':', 1)[1].split())
        receipt['dependency_origin'] = {d: str(dependency_path(d)) for d in sorted(deps)}
        receipt['dependency_sha256'] = {actual: sha(actual) for actual in receipt['dependency_origin'].values()}
        receipt['sealed_v6_owner_headers_in_binaries'] = sorted(d for d in deps if '/pic-electrode-heat-source-v6-20261030/native/' in d)
        receipt['binary_sha256'] = {exe: sha(b / exe) for exe in exes}
        receipt['build_dir'] = str(b)
        receipt['own_sha256'] = before
        if sealed() != pins or own() != before:
            raise RuntimeError('sources changed during the build')
    except Exception as e:
        status = 'failed'
        receipt['error'] = str(e)
    receipt['status'] = status
    (out / 'build-receipt.json').write_text(json.dumps(receipt, indent=1) + '\n')
    print(json.dumps({'build': str(out), 'status': status, 'error': receipt.get('error', '')[:3000]}))
    return 0 if status == 'passed' else 1


def usable(lab):
    br = json.loads((OUT / lab / 'build/build-receipt.json').read_text())
    if br['status'] != 'passed' or br['own_sha256'] != own() or br['sealed_sha256'] != sealed():
        raise RuntimeError('source/dependency capture differs from qualified build')
    for path, digest in br['dependency_sha256'].items():
        if sha(Path(path)) != digest:
            raise RuntimeError('compiled dependency changed: '+path)
    for exe, digest in br['binary_sha256'].items():
        if sha(Path(br['build_dir']) / exe) != digest:
            raise RuntimeError('binary changed: '+exe)
    return br


def summarize(f):
    j = json.loads(f.read_text())
    return {'failed': j['failed'], 'passed': sum(1 for c in j['cases'] if c['passed']), 'cases': len(j['cases']), 'error': j.get('error'),
            'names': [{'group': c['group'], 'name': c['name'], 'passed': c['passed']} for c in j['cases']]}


def run(label, extra):
    adm = admission('run')
    out = OUT / label / ('run-' + time.strftime('%H%M%S'))
    out.mkdir(parents=True)
    summary = {'schema': 'source-v7-run-receipt-v1', 'label': label, 'admission': adm, 'sealed_sha256': sealed(), 'own_sha256': own(), 'commands': []}
    plan = [(label, 'normal')] + [(x, 'sanitizer-' + json.loads((OUT / x / 'build/build-receipt.json').read_text())['variant']) for x in extra]
    for lab, sub in plan:
        br = usable(lab)
        env = dict(ENV)
        if 'asan' in sub:
            env['ASAN_OPTIONS'] = 'detect_leaks=0:abort_on_error=1'   # MPI/MFEM static state is not freed at exit
        registry_out = out / (sub + '-test_registry_mass')
        registry_out.mkdir()
        rr = execute(sub+'-test_registry_mass', [Path(br['build_dir'])/'test_registry_mass', registry_out], out, 600, env)
        summary['commands'].append(rr)
        if rr['returncode'] != 0:
            raise RuntimeError('actual registry mass qualification failed: '+rr['log'])
        summary[sub+'/test_registry_mass'] = summarize(registry_out/'registry-mass.json')
        env['PF_SOURCE_MASS_CONFIG'] = str(registry_out/'mass-configuration.json')
        for legacy_exe in ('legacy_kernel_new','legacy_kernel_old'):
            r = execute(sub+'-'+legacy_exe,[Path(br['build_dir'])/legacy_exe,out/(sub+'-'+legacy_exe+'.bin')],out,600,env)
            summary['commands'].append(r)
        same = (out/(sub+'-legacy_kernel_new.bin')).read_bytes() == (out/(sub+'-legacy_kernel_old.bin')).read_bytes()
        summary[sub+'/legacy_bytes']={'equal':same,'new_sha256':sha(out/(sub+'-legacy_kernel_new.bin')),'old_sha256':sha(out/(sub+'-legacy_kernel_old.bin'))}
        if not same:
            raise RuntimeError('default emission/layer bytes differ from sealed legacy kernels')
        for exe, args, result in (('test_surface_v7', [], 'surface-v7.json'), ('test_source_v7', [], 'source-v7.json'),
                                  ('test_source_v6_on_v7', [TRACE], 'source-v6.json')):
            d = out / (sub + '-' + exe)
            d.mkdir()
            r = execute(sub + '-' + exe, [Path(br['build_dir']) / exe, *args, d], out, 7200, env)
            summary['commands'].append(r)
            if (d / result).exists():
                summary[sub + '/' + exe] = summarize(d / result)
    summary['all_passed'] = all(c['returncode'] == 0 for c in summary['commands']) and own() == summary['own_sha256'] and sealed() == summary['sealed_sha256']
    for lab, _ in plan:
        usable(lab)
    (out / 'run-receipt.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps({'run': str(out), 'all_passed': summary['all_passed']}))
    return 0 if summary['all_passed'] else 1


def gpu(label):
    adm = admission('gpu')
    br = usable(label)
    out = OUT / label / ('gpu-' + time.strftime('%H%M%S'))
    out.mkdir(parents=True)
    env = dict(ENV, CUDA_VISIBLE_DEVICES=GPU_UUID)
    summary = {'schema': 'source-v7-gpu-receipt-v1', 'label': label, 'admission': adm, 'sealed_sha256': sealed(), 'own_sha256': own(),
               'nvrtc_library': str(NVRTC), 'nvrtc_library_sha256': sha(NVRTC), 'device_uuid_required': GPU_UUID, 'commands': []}
    d = out / 'test_wall_gpu_v7'
    d.mkdir()
    r = execute('gpu-test_wall_gpu_v7', [Path(br['build_dir']) / 'test_wall_gpu_v7', NVRTC, GPU_UUID, V7N, A1 / 'include', d], out, 3600, env)
    summary['commands'].append(r)
    if (d / 'wall-gpu-v7.json').exists():
        summary['test_wall_gpu_v7'] = summarize(d / 'wall-gpu-v7.json')
    summary['all_passed'] = r['returncode'] == 0
    (out / 'gpu-receipt.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps({'gpu': str(out), 'all_passed': summary['all_passed']}))
    return 0 if summary['all_passed'] else 1


if __name__ == '__main__':
    a = sys.argv
    if len(a) < 3 or a[1] not in ('build', 'run', 'gpu') or '/' in a[2]:
        raise SystemExit(__doc__)
    if a[1] == 'build':
        raise SystemExit(build(a[2], a[3] if len(a) > 3 else 'normal'))
    if a[1] == 'run':
        raise SystemExit(run(a[2], a[3:]))
    raise SystemExit(gpu(a[2]))
