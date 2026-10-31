#!/usr/bin/python -I
"""Coordinator v5 (Stage 1 flow on Hall corrective v2, native candidate003) build and run driver.

Builds hold both core leases (6,7) when available, else one leased core; runs hold one core lease.
GPU stages additionally hold both GPU lock files:
  cd /home/treefor/projects/runtimes/pic-water-common-commit-v5-20261030
  core-lease-20261030/both-cores.sh /usr/bin/python -I run_v5.py build LABEL
  core-lease-20261030/any-core.sh --wait /usr/bin/python -I run_v5.py cpu LABEL [SWEEP_BUDGET_S] [WATER_WEIGHT]
  flock GPU.lock flock GPU-UUID.lock core-lease-20261030/any-core.sh --wait /usr/bin/python -I run_v5.py gpu LABEL
  ... run_v5.py multistep LABEL cpu|gpu CARRIERS MAX_STEPS [PARITY_EVERY]
  ... run_v5.py hall-loop LABEL cpu|gpu CARRIERS STEPS [PARITY_EVERY]

Every compile, link and test command runs inside the frozen launcher. Sealed inputs are hashed
before and after; nothing outside the v5 runtime directory and the v5 output directory is written.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import fcntl, hashlib, json, os, platform, resource, subprocess, sys, time

P = Path(__file__).resolve().parent
R = Path('/home/treefor/projects/runtimes')
OUT = Path('/home/treefor/projects/.workspace/outputs/pic-water-common-commit-v5-20261030')
N = R / 'pic-reaction-transactions-20261029/native-full-time-candidate003'
H = R / 'hall-coupled-transport-20261029/source'
HC = H / 'development/hall-field-correction-20261030'
HD = H / 'development/hall-canonical-driver-20261029'
I = R / 'pic-water-full-time-subflow-20261029/inputs'
V3 = R / 'pic-water-common-commit-v3-20261030'
M = R / 'mfem-4.10-pic-20261029'
F = R / 'arch-base-20261028/run-frozen'
LEASE_DIR = Path('/home/treefor/projects/.workspace/outputs/pic-profile')
LEASES = {6: LEASE_DIR / 'testing-core.lock', 7: LEASE_DIR / 'testing-core7.lock'}
HV2 = R / 'hall-corrective-v2-20261030'
G = R / 'hall-gpu-resident-v1-20261030'
GPU_UUID = 'GPU-a7b9b458-167f-cd57-fec9-3dea0995b6e9'
GPU_LOCKS = [LEASE_DIR / 'testing-gpu.lock', LEASE_DIR / ('testing-gpu-%s.lock' % GPU_UUID)]
NVRTC = R / 'pic-cuda-gather-investigation-20261029/tools/nvrtc/nvidia/cuda_nvrtc/lib/libnvrtc.so.12'
NVRTC_INC = R / 'pic-cuda-gather-investigation-20261029/tools/nvrtc/nvidia/cuda_nvrtc/include'
CUDA_INC = R / 'mfem-cuda-benchmark-20261029/cuda-libs/nvidia/cuda_runtime/include'
STAGE1 = Path('/home/treefor/sim-workspace/powerflow/outputs/pic-water-common-commit-v4-20261030/run-s1d-6-164057')
IDENTITY = 'c045c42d2f57082e3fb5e12172c440d2470e681a38df9ac406bd7e063609a08d'
# Native selection (build option): candidate003 (Stage 1 parity) or native v5 (candidate005).
NATIVES = {'003': {'N': R / 'pic-reaction-transactions-20261029/native-full-time-candidate003',
                   'IDENTITY': 'c045c42d2f57082e3fb5e12172c440d2470e681a38df9ac406bd7e063609a08d',
                   'LIB': 'build/libpic_native_candidate003.a',
                   'LIB_SHA': 'd2568af9fac13f8aa056ae4a9a50bb2a07eee05d7218d1c4b79aef21317d8a20', 'define': []},
           '005': {'N': R / 'pic-reaction-transactions-20261029/native-full-time-candidate005',
                   'IDENTITY': 'cca72822c9b61900aae7572f29c5d57314f45aef17c4b01c9fd2346c6fa3c883',
                   'LIB': 'build/libpic_native_v5.a',
                   'LIB_SHA': '715729114e5a7d3975bd7ee7db92c570e61f78c907f3e999812b34e7a74159f2', 'define': ['-DCOORDINATOR_NATIVE_V5=1']}}
NATIVE_WHICH = '003'


def select_native(which):
    global N, IDENTITY, NATIVE_LIB, NATIVE_LIB_SHA, NATIVE_WHICH
    c = NATIVES[which]
    N, IDENTITY, NATIVE_LIB, NATIVE_LIB_SHA, NATIVE_WHICH = c['N'], c['IDENTITY'], c['N'] / c['LIB'], c['LIB_SHA'], which
NATIVE_LIB = N / 'build/libpic_native_candidate003.a'
NATIVE_LIB_SHA = 'd2568af9fac13f8aa056ae4a9a50bb2a07eee05d7218d1c4b79aef21317d8a20'
HALL_COMMIT = '8dddec0610982778bf02bdfe697b70a37129c012'
TRACE = H / 'development/hall-finite-foil-em-20261029/data/shot_17133_current.csv'
TRACE_SHA = '2c54d2dbaec6b6c8fcf2e45b5abd837aa72999f5fd14217f530277ffd5949076'
DATA = R / 'pic-water-exact-count-subflow-20261029/inputs/h2o-gas.json'
BINDING = Path('/home/treefor/sim-workspace/powerflow/outputs/pic-water-exact-count-subflow-20261029/run-binding.json')
V3_SEALED = {'source/native_parent_subflow.hpp': '013a7e722fdbf8edf5f10c27d991de3a54482b444803d0dc776567c875ba22eb',
             'source/common_commit_contract.hpp': '7666c724ea4c6db5e8518322274aa998da8057cca8466114895ba70ab1578b20'}
LIBS = [M / 'install-mpi/lib/libmfem.a', M / 'tpl/hypre/lib/libHYPRE.a', M / 'tpl/metis/lib/libmetis.a', M / 'tpl/GKlib/lib/libGKlib.a']
NATIVE_WRAPS = ['_ZN6raikou12mfem_backend6GapPic12solve_fieldsEv', '_ZN6raikou12mfem_backend6GapPic7advanceEd',
                '_ZN6raikou12mfem_backend6GapPic13advance_emptyEd',
                '_ZN6raikou12mfem_backend6GapPic16advance_observedEdRKSt8functionIFvRKNS0_11GapParticleEidEE']
TIMING_WRAPS_UNUSED = ['_ZN11hall_driver24prepare_canonical_motionERKN5axial4MeshERKN6raikou12mfem_backend13GapPicOptionsESt4spanIKNS5_14ArenaRateEntryELm18446744073709551615EERKNS_13MotionBindingES9_IKmLm18446744073709551615EERKSt6vectorIdSaIdEESM_',
                '_ZN11hall_driver18total_field_sha256B5cxx11ERKSt6vectorIdSaIdEES4_',
                '_ZN11hall_driver23canonical_census_sha256B5cxx11ESt4spanIKN6raikou12mfem_backend14ArenaRateEntryELm18446744073709551615EE',
                '_ZN11hall_driver18motion_mesh_sha256B5cxx11ERKN5axial4MeshERKN6raikou12mfem_backend13GapPicOptionsE']
ENV = {'PATH': '/usr/bin:/bin', 'LC_ALL': 'C', 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1',
       'MKL_NUM_THREADS': '1', 'BLIS_NUM_THREADS': '1', 'NUMEXPR_NUM_THREADS': '1', 'VECLIB_MAXIMUM_THREADS': '1',
       'OMP_THREAD_LIMIT': '1', 'CUDA_VISIBLE_DEVICES': '', 'PYTHONNOUSERSITE': '1'}
GPU_ENV = {k: v for k, v in ENV.items() if k != 'CUDA_VISIBLE_DEVICES'}
RLIMIT_AS = 3 * 1024 ** 3
REJECT_EXPECT = {1: 'read-only object', 2: 'read-only object', 3: 'read-only object',
                 4: 'private within this context', 5: 'private within this context',
                 6: 'use of deleted function', 7: 'try_rollback_closed_parent_joint_v5',
                 8: 'incomplete type', 9: 'private within this context'}
GPU_SEALED = [G / 'common/hall_gpu_common.cpp', G / 'common/hall_gpu_common.hpp', G / 'common/hall_gpu_primitives.cpp',
              G / 'common/hall_gpu_primitives.hpp', G / 'common/device/deterministic_primitives.cuh', G / 'push/hall_gpu_push.cpp',
              G / 'push/hall_gpu_push.hpp', G / 'push/push_core.hpp', G / 'push/push_abi.hpp', G / 'push/push_bodies.hpp',
              G / 'push/push_kernels.cu', G / 'push/push_fold.hpp', G / 'push/hall_adapter.hpp']


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


def admission(kind):
    if not __debug__ or os.getenv('PYTHONOPTIMIZE'):
        fail('optimized Python guard execution refused')
    cores = os.sched_getaffinity(0)
    allowed = [{6, 7}, {6}, {7}] if kind == 'build' else [{6}, {7}]
    if cores not in allowed or os.getpriority(os.PRIO_PROCESS, 0) != 15:
        fail('affinity %s / nice %d not admitted for %s' % (sorted(cores), os.getpriority(os.PRIO_PROCESS, 0), kind))
    if any(os.getenv(k) != '1' for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']):
        fail('one thread per library required')
    for core in cores:
        if not held(LEASES[core]):
            fail('core %d lease must be held by the launching wrapper; no work started' % core)
    gpu = kind == 'gpu'
    if gpu and not all(held(l) for l in GPU_LOCKS):
        fail('both GPU lock files must be held by the launching wrapper')
    if Path('/proc/sys/kernel/osrelease').read_text().strip() != (R / 'arch-base-20261028/kernel.release').read_text().strip():
        fail('kernel release differs from the frozen runtime')
    if not gpu:
        resource.setrlimit(resource.RLIMIT_AS, (RLIMIT_AS, RLIMIT_AS))
    return {'kind': kind, 'cpu_affinity': sorted(cores), 'PF_CORE': os.getenv('PF_CORE'), 'nice': 15,
            'threads': 1, 'gpu': gpu, 'gpu_uuid': GPU_UUID if gpu else None,
            'leases_verified_held': [str(LEASES[c]) for c in sorted(cores)] + ([str(l) for l in GPU_LOCKS] if gpu else []),
            'rlimit_as_bytes': None if gpu else RLIMIT_AS, 'python': sys.executable, 'python_version': platform.python_version(),
            'frozen_launcher': str(F), 'frozen_launcher_sha256': sha(F)}


def sealed_pins():
    """Admit the sealed components; returns their identity record."""
    record = {}
    manifest = json.loads((N / 'SOURCE-SHA256.json').read_text())
    if manifest['source_identity'] != IDENTITY:
        fail('native source identity differs')
    for name, digest in manifest['files'].items():
        if sha(N / name) != digest:
            fail('native sealed file changed: ' + name)
    if sha(NATIVE_LIB) != NATIVE_LIB_SHA:
        fail('native candidate003 library differs')
    record['native'] = {'source_identity': IDENTITY, 'files_verified': len(manifest['files']),
                        'library': str(NATIVE_LIB), 'library_sha256': NATIVE_LIB_SHA,
                        'manifest_sha256': sha(N / 'SOURCE-SHA256.json')}
    head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=H, text=True).strip()
    dirty = subprocess.check_output(['git', 'status', '--porcelain'], cwd=H, text=True)
    if head != HALL_COMMIT or dirty:
        fail('Hall repository is not the clean corrective commit')
    record['hall'] = {'commit': head, 'git_clean': True}
    provider = {}
    for name in ('motion_receipt.hpp', 'native_arena_motion.cpp', 'kinetic_boundary.hpp', 'kinetic_boundary.cpp'):
        a, b = sha(N / 'inputs/hall-motion-v2' / name), sha(HD / name)
        if a != b:
            fail('native hall-motion-v2 input differs from the Hall provider: ' + name)
        provider[name] = a
    record['hall_motion_v2_byte_identity'] = provider
    for name, digest in V3_SEALED.items():
        if sha(V3 / name) != digest:
            fail('V3 sealed source changed: ' + name)
    record['v3_sources'] = V3_SEALED
    pins = json.loads((HC / 'runtime-pins.json').read_text())
    base = Path(pins['frozen_root']).parent
    if sha(base / 'root-files.json') != pins['frozen_root_manifest_sha256'] or sha(base / 'run-frozen') != pins['launcher_sha256']:
        fail('frozen runtime manifest or launcher changed')
    root = Path(pins['frozen_root'])
    for name, data in pins['selected_files'].items():
        p = root / name
        if 'sha256' in data:
            if sha(p) != data['sha256']:
                fail('frozen runtime file changed: ' + name)
        elif 'symlink' in data:
            if not p.is_symlink() or os.readlink(p) != data['symlink']:
                fail('frozen runtime link changed: ' + name)
        else:
            fail('unknown frozen runtime pin entry')
    for name, digest in pins['static_libraries'].items():
        if sha(name) != digest:
            fail('pinned static library changed: ' + name)
    record['frozen_runtime'] = {'runtime_pins': str(HC / 'runtime-pins.json'), 'runtime_pins_sha256': sha(HC / 'runtime-pins.json'),
                                'selected_files_verified': len(pins['selected_files']),
                                'static_libraries': pins['static_libraries']}
    if sha(TRACE) != TRACE_SHA:
        fail('C19 trace differs')
    binding = json.loads(BINDING.read_text())
    if sha(DATA) != binding['dataset_sha256']:
        fail('water data differs from its binding')
    record['gpu_sealed_sources'] = {str(p): sha(p) for p in GPU_SEALED}
    record['nvrtc'] = {'library': str(NVRTC), 'sha256': sha(NVRTC)}
    record['inputs'] = {'trace': str(TRACE), 'trace_sha256': TRACE_SHA, 'water_data': str(DATA),
                        'water_data_sha256': sha(DATA), 'binding': str(BINDING), 'binding_sha256': sha(BINDING)}
    return record


def tree_hashes(paths):
    out = {}
    for root in paths:
        root = Path(root)
        files = [root] if root.is_file() else sorted(p for p in root.rglob('*') if p.is_file())
        for p in files:
            out[str(p)] = sha(p)
    return out


def own_sources():
    return tree_hashes([P / 'source', P / 'tests', Path(__file__).resolve(), P / 'compare_stage1.py',
                        HV2 / 'source', HV2 / 'run_hall_v2.py'])


def frozen(cmd):
    return [str(F), *map(str, cmd)]


def execute(name, cmd, log_dir, timeout, expected=0, env=None):
    log = log_dir / (name + '.log')
    start = time.monotonic()
    with log.open('w') as handle:
        try:
            r = subprocess.run(frozen(cmd), cwd=P, env=env or ENV, stdout=handle, stderr=subprocess.STDOUT, timeout=timeout)
            rc = r.returncode
        except subprocess.TimeoutExpired:
            rc = 'timeout'
    return {'name': name, 'argv': frozen(cmd), 'returncode': rc, 'expected_returncode': expected,
            'ok': (rc == expected) if expected != 'nonzero' else (isinstance(rc, int) and rc != 0),
            'elapsed_s': time.monotonic() - start, 'log': str(log), 'log_sha256': sha(log)}


def dependencies(d_file):
    text = Path(d_file).read_text().replace('\\\n', ' ')
    deps = set()
    for line in text.splitlines():
        if ':' not in line:
            continue
        for item in line.split(':', 1)[1].split():
            deps.add(str(Path(item).resolve()))
    return sorted(deps)


def build(label, which='003'):
    select_native(which)
    adm = admission('build')
    OUT.mkdir(exist_ok=True)
    out = OUT / ('build-' + label)
    out.mkdir(exist_ok=False)
    bdir = P / 'build' / label
    bdir.mkdir(parents=True, exist_ok=False)
    pins = sealed_pins()
    before = own_sources()
    common = ['-std=c++20', '-O2', '-ffp-contract=off', '-fno-fast-math', '-ffunction-sections', '-fdata-sections']
    isystem = []
    for d in (N / 'source/include', I, I / 'include', N / 'inputs', HV2 / 'source', H / 'development',
              M / 'install-mpi/include', M / 'tpl/hypre/include', M / 'tpl/metis/include'):
        isystem += ['-isystem', str(d)]
    gpu_inc = ['-I' + str(G), '-isystem', str(CUDA_INC), '-isystem', str(NVRTC_INC)]
    ours = ['/usr/bin/mpic++', *common, '-Wall', '-Wextra', '-Wpedantic', '-Werror',
            '-I' + str(P / 'source'), '-I' + str(V3 / 'source'),
            '-DPIC_REACTION_NATIVE_SOURCE_SHA256="' + IDENTITY + '"', *NATIVES[which]['define'], *isystem]
    # Hall v2 and GPU translation units: Hall's option set (no warning flags, NDEBUG).
    hall = ['/usr/bin/c++', *common, '-DNDEBUG', *isystem]
    gpu = ['/usr/bin/c++', *common, '-DNDEBUG', *gpu_inc, '-Wno-c++26-extensions']
    objects = {'hall_motion_v2': (hall, HV2 / 'source/motion_provider.cpp'),
               'kinetic_boundary': (hall, HD / 'kinetic_boundary.cpp'),
               'hall_gpu_backends_v2': ([*hall, *gpu_inc], HV2 / 'source/gpu/gpu_backends.cpp'),
               'hall_gpu_common': (gpu, G / 'common/hall_gpu_common.cpp'),
               'hall_gpu_primitives': (gpu, G / 'common/hall_gpu_primitives.cpp'),
               'hall_gpu_push': (gpu, G / 'push/hall_gpu_push.cpp'),
               'hall_motion_adapter': (ours, N / 'source/src/mfem_pic_hall_motion_adapter.cpp'),
               'alloc_hook': (ours, P / 'tests/alloc_hook.cpp'),
               'test_stage1_v5': ([*ours, *gpu_inc], P / 'tests/test_stage1_v5.cpp')}
    jobs = [('compile-' + k, [*flags, '-MD', '-MF', str(bdir / (k + '.d')), '-c', str(src), '-o', str(bdir / (k + '.o'))])
            for k, (flags, src) in objects.items()]
    receipt = {'schema': 'stage1-v5-build-receipt-v1', 'label': label, 'admission': adm, 'sealed': pins,
               'native': {'which': which, 'source_identity': IDENTITY, 'library': str(NATIVE_LIB), 'library_sha256': NATIVE_LIB_SHA},
               'started_unix': time.time(), 'environment': ENV, 'commands': [],
               'decisions': {
                   'hall_v2': 'Hall corrective v2 (hall-corrective-v2-20261030/source) replaces the 8dddec0 corrective headers and the '
                       'sealed native_arena_motion.cpp; it defines hall_driver::MotionReceipt (same layout, extra v2 assembler friend) '
                       'and the unchanged hall_driver motion functions. Sealed kinetic_boundary.cpp, axial/foil field headers and '
                       'the v2 heat owner are used unchanged from 8dddec0.',
                   'native_adapter': 'native candidate003 mfem_pic_hall_motion_adapter.cpp is compiled unchanged; it only forward-declares '
                       'hall_driver::MotionReceipt and calls hall_driver::motion_record.',
                   'gpu': 'hall-gpu-resident-v1 common/ and push/ are compiled unchanged; the field CG runs on the same session.',
                   'mfem_gap_pic_not_linked_separately': 'pic_boris and pic_kinetic_energy_J come from libpic_native_candidate003.a.'}}

    def finish(status, error=None):
        receipt['status'] = status
        receipt['error'] = error
        receipt['finished_unix'] = time.time()
        (out / 'build-receipt.json').write_text(json.dumps(receipt, indent=1) + '\n')
        print(json.dumps({'build': str(out), 'status': status, 'error': error}))

    try:
        workers = 2 if len(adm['cpu_affinity']) == 2 else 1
        with ThreadPoolExecutor(max_workers=workers) as pool:
            results = list(pool.map(lambda j: execute(j[0], j[1], out, 1800), jobs))
        receipt['commands'] += results
        for r in results:
            if not r['ok']:
                fail(r['name'] + ' failed:\n' + Path(r['log']).read_text()[-6000:])
        exe = bdir / 'test_stage1_v5'
        link = ['/usr/bin/mpic++', '-rdynamic', *['-Wl,--wrap=' + s for s in NATIVE_WRAPS], '-Wl,--gc-sections',
                *[bdir / (k + '.o') for k in objects], NATIVE_LIB, *LIBS, '-lcrypto', '-ldl', '-o', exe]
        r = execute('link', link, out, 1800)
        receipt['commands'].append(r)
        if not r['ok']:
            fail('link failed:\n' + Path(r['log']).read_text()[-6000:])
        decl = [('declaration', [*ours, *gpu_inc, '-fsyntax-only', P / 'tests/declaration_v5.cpp'], 0)]
        decl += [('reject-%d' % k, [*ours, *gpu_inc, '-DREJECT_CASE=%d' % k, '-fsyntax-only', P / 'tests/declaration_v5.cpp'], 'nonzero')
                 for k in sorted(REJECT_EXPECT)]
        with ThreadPoolExecutor(max_workers=workers) as pool:
            results = list(pool.map(lambda j: execute(j[0], j[1], out, 900, j[2]), decl))
        rejects = {}
        for r in results:
            if r['name'].startswith('reject-'):
                k = int(r['name'].split('-')[1])
                text = Path(r['log']).read_text()
                r['expected_diagnostic'] = REJECT_EXPECT[k]
                r['diagnostic_found'] = REJECT_EXPECT[k] in text
                r['ok'] = r['ok'] and r['diagnostic_found']
                rejects[k] = r['ok']
        receipt['commands'] += results
        receipt['declaration'] = {'compiles': results[0]['ok'], 'reject_cases': rejects,
                                  'all_rejects_fail_for_the_stated_reason': all(rejects.values())}
        deps = {k: dependencies(bdir / (k + '.d')) for k in objects}
        allfiles = sorted({d for v in deps.values() for d in v})
        receipt['dependency_sha256'] = {d: sha(d) for d in allfiles}
        json_headers = {k: [d for d in v if Path(d).name in ('json.hpp', 'json.h')] for k, v in deps.items()}
        identities = {d: sha(d) for v in json_headers.values() for d in v}
        hpp = {d for d in identities if d.endswith('json.hpp')}
        receipt['json_header_identity'] = {'per_tu': json_headers, 'sha256': identities,
                                           'single_nlohmann_json_hpp': len(hpp) == 1 and len({identities[d] for d in hpp}) == 1}
        if not receipt['json_header_identity']['single_nlohmann_json_hpp']:
            fail('more than one nlohmann json.hpp in the link: ' + ', '.join(sorted(hpp)))
        rogue = [d for d in allfiles if 'native-v6.2' in d or 'hall-field-correction-20261030' in d or
                 d.endswith('hall-canonical-driver-20261029/motion_receipt.hpp') or 'hall-motion-v2/motion_receipt.hpp' in d]
        receipt['sealed_corrective_or_sealed_motion_receipt_headers_in_link'] = rogue
        if rogue:
            fail('sealed 8dddec0 corrective / motion receipt or native-v6.2 headers leaked into the link: ' + ', '.join(rogue))
        receipt['libraries_sha256'] = {str(x): sha(x) for x in [NATIVE_LIB, *LIBS]}
        receipt['objects_sha256'] = {k: sha(bdir / (k + '.o')) for k in objects}
        receipt['binary'] = str(exe)
        receipt['binary_sha256'] = sha(exe)
        receipt['own_sources_sha256'] = before
        if own_sources() != before or sealed_pins() != pins:
            fail('sources or sealed inputs changed during the build')
        receipt['inputs_unchanged_after_build'] = True
        if not receipt['declaration']['compiles'] or not receipt['declaration']['all_rejects_fail_for_the_stated_reason']:
            finish('declaration-failed')
            return 1
        finish('passed')
        return 0
    except Exception as exc:
        finish('failed', str(exc))
        return 1


def checked_build(label):
    build_receipt = json.loads((OUT / ('build-' + label) / 'build-receipt.json').read_text())
    select_native(build_receipt.get('native', {}).get('which', '003'))
    if build_receipt['status'] != 'passed':
        fail('build did not pass')
    exe = Path(build_receipt['binary'])
    if sha(exe) != build_receipt['binary_sha256']:
        fail('test binary differs from its build receipt')
    if own_sources() != build_receipt['own_sources_sha256']:
        fail('v5 or Hall v2 sources changed after the build')
    return build_receipt, exe


def run_record(name, r, d, adm, build_receipt, pins, result_name, extras):
    rec = {'schema': 'stage1-v5-run-receipt-v1', 'command': r, 'admission': adm,
           'binary_sha256': build_receipt['binary_sha256'], 'inputs_sha256': {
               str(TRACE): TRACE_SHA, str(DATA): pins['inputs']['water_data_sha256'], str(BINDING): pins['inputs']['binding_sha256']}}
    result = d / result_name
    if result.exists():
        res = json.loads(result.read_text())
        rec['result_file'] = str(result)
        rec['result_sha256'] = sha(result)
        rec['summary'] = res.get('summary')
        rec['cases'] = res.get('cases')
    for extra in extras:
        if (d / extra).exists():
            rec.setdefault('files', {})[extra] = {'path': str(d / extra), 'sha256': sha(d / extra)}
    return rec


def run(label, sweep_budget, water_weight):
    adm = admission('run')
    build_receipt, exe = checked_build(label)
    pins = sealed_pins()
    if pins != build_receipt['sealed']:
        fail('sealed inputs differ from the build')
    out = OUT / ('run-' + label + '-cpu-' + str(adm['cpu_affinity'][0]) + '-' + time.strftime('%H%M%S'))
    summary_native = build_receipt.get('native')
    out.mkdir(exist_ok=False)
    runs = {}
    summary = {'schema': 'stage1-v5-run-summary-v1', 'label': label, 'admission': adm, 'sealed': pins,
               'backends': {'field': 'cpu-pcg', 'push': 'cpu'}, 'native': summary_native,
               'build_receipt': str(OUT / ('build-' + label) / 'build-receipt.json'),
               'binary_sha256': build_receipt['binary_sha256'], 'started_unix': time.time()}
    status = 'passed'
    try:
        for name, mode, arg, timeout in [('t1-a', 't1', water_weight, 3600), ('t1-b', 't1', water_weight, 3600),
                                         ('matrix', 'matrix', sweep_budget, 4 * 3600)]:
            d = out / name
            d.mkdir()
            r = execute(name, [exe, mode, TRACE, DATA, BINDING, d, IDENTITY, repr(arg)], out, timeout, 0)
            rec = run_record(name, r, d, adm, build_receipt, pins, 't1-result.json' if mode == 't1' else 'matrix-result.json',
                             ('t1-receipt.json', 't1-timings.json', 'allocation-sweep.json'))
            (out / (name + '-receipt.json')).write_text(json.dumps(rec, indent=1) + '\n')
            runs[name] = rec
            print(name, r['returncode'], round(r['elapsed_s'], 1), flush=True)
        a, b = out / 't1-a/t1-receipt.json', out / 't1-b/t1-receipt.json'
        same = a.exists() and b.exists() and a.read_bytes() == b.read_bytes()
        summary['t1_two_fresh_processes_byte_identical'] = same
        summary['t1_receipt_sha256'] = [sha(a) if a.exists() else None, sha(b) if b.exists() else None]
        if a.exists() and (out / 'matrix/matrix-result.json').exists():
            t1 = json.loads(a.read_text())
            mx = json.loads((out / 'matrix/matrix-result.json').read_text())
            ref = (mx.get('reference_after_step') or {}).get('digest')
            summary['matrix_reference_equals_t1_receipt'] = ref == t1.get('P1', {}).get('after', {}).get('digest')
        per_test = {}
        for name, rec in runs.items():
            for case in rec.get('cases') or []:
                t = per_test.setdefault(case['test'], {'passed': 0, 'failed': 0, 'failed_cases': []})
                if case['passed']:
                    t['passed'] += 1
                else:
                    t['failed'] += 1
                    t['failed_cases'].append(name + ': ' + case['name'])
        summary['per_test'] = per_test
        summary['declaration_from_build'] = build_receipt.get('declaration')
        summary['runs'] = {k: {'returncode': v['command']['returncode'], 'elapsed_s': v['command']['elapsed_s'],
                               'receipt': str(out / (k + '-receipt.json'))} for k, v in runs.items()}
        if (out / 't1-a').exists():
            cmp = subprocess.run(['/usr/bin/python', '-I', str(P / 'compare_stage1.py'), str(STAGE1 / 't1-a'), str(out / 't1-a'),
                                  str(out / 'stage1-identity-mapping.json'), *(['--native-v5'] if NATIVE_WHICH == '005' else [])],
                                 capture_output=True, text=True, env=ENV)
            summary['stage1_comparison'] = {'returncode': cmp.returncode, 'output': cmp.stdout[-4000:] + cmp.stderr[-2000:],
                                            'mapping': str(out / 'stage1-identity-mapping.json')}
        if own_sources() != build_receipt['own_sources_sha256'] or sealed_pins() != pins:
            fail('sources or sealed inputs changed during the run')
        summary['inputs_unchanged_after_run'] = True
        if any(v['command']['returncode'] != 0 for v in runs.values()) or not same or \
                summary.get('stage1_comparison', {}).get('returncode') != 0:
            status = 'completed-with-failures'
    except Exception as exc:
        status = 'failed'
        summary['error'] = str(exc)
    summary['status'] = status
    summary['finished_unix'] = time.time()
    (out / 'stage1-v5-receipt.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps({'run': str(out), 'status': status}))
    return 0 if status == 'passed' else 1


def run_points(label, k0, k1):
    adm = admission('run')
    build_receipt, exe = checked_build(label)
    pins = sealed_pins()
    out = OUT / ('sweep-points-' + label + '-' + time.strftime('%H%M%S'))
    out.mkdir(parents=True, exist_ok=False)
    r = execute('sweep-points', [exe, 'sweep-point', TRACE, DATA, BINDING, out, IDENTITY, str(k0), str(k1)], out, 4 * 3600, 0)
    rec = {'schema': 'stage1-v5-sweep-points-receipt-v1', 'command': r, 'admission': adm, 'sealed': pins,
           'binary_sha256': build_receipt['binary_sha256'], 'native': build_receipt.get('native')}
    (out / 'sweep-points-receipt.json').write_text(json.dumps(rec, indent=1) + '\n')
    print(json.dumps({'run': str(out), 'returncode': r['returncode']}))
    return 0


def gpu_flags():
    return ['--field', 'gpu', '--push', 'gpu', '--uuid', GPU_UUID, '--nvrtc', str(NVRTC)]


def run_gpu(label, water_weight):
    adm = admission('gpu')
    build_receipt, exe = checked_build(label)
    pins = sealed_pins()
    out = OUT / ('run-' + label + '-gpu-' + str(adm['cpu_affinity'][0]) + '-' + time.strftime('%H%M%S'))
    out.mkdir(exist_ok=False)
    summary = {'schema': 'stage1-v5-gpu-run-summary-v1', 'label': label, 'admission': adm, 'sealed': pins,
               'backends': {'field': 'gpu1-factored-cg', 'push': 'gpu-push-v1'}, 'binary_sha256': build_receipt['binary_sha256']}
    runs = {}
    for name in ('gpu-step-a', 'gpu-step-b'):
        d = out / name
        d.mkdir()
        r = execute(name, [exe, 'gpu-step', TRACE, DATA, BINDING, d, IDENTITY, repr(water_weight), *gpu_flags()], out, 3600, 0, GPU_ENV)
        rec = run_record(name, r, d, adm, build_receipt, pins, 'gpu-step-result.json',
                         ('gpu-step-receipt.json', 'gpu-step-deterministic.json', 'gpu-teardown.json'))
        (out / (name + '-receipt.json')).write_text(json.dumps(rec, indent=1) + '\n')
        runs[name] = rec
        print(name, r['returncode'], round(r['elapsed_s'], 1), flush=True)
    a, b = out / 'gpu-step-a/gpu-step-deterministic.json', out / 'gpu-step-b/gpu-step-deterministic.json'
    summary['two_fresh_processes_byte_identical'] = a.exists() and b.exists() and a.read_bytes() == b.read_bytes()
    summary['runs'] = {k: {'returncode': v['command']['returncode'], 'summary': v.get('summary')} for k, v in runs.items()}
    summary['status'] = 'passed' if summary['two_fresh_processes_byte_identical'] and all(
        v['command']['returncode'] == 0 for v in runs.values()) else 'completed-with-failures'
    (out / 'gpu-step-summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps({'run': str(out), 'status': summary['status']}))
    return 0 if summary['status'] == 'passed' else 1


def run_series(kind, label, backend, carriers, steps, parity, extra=()):
    adm = admission('gpu' if backend == 'gpu' else 'run')
    build_receipt, exe = checked_build(label)
    pins = sealed_pins()
    name = '%s-%s-%s-%d' % (kind, backend, carriers, int(time.time()) % 100000)
    out = OUT / ('series-' + label) / name
    out.mkdir(parents=True, exist_ok=False)
    flags = gpu_flags() if backend == 'gpu' else ['--field', 'cpu', '--push', 'cpu']
    r = execute(kind, [exe, kind, TRACE, DATA, BINDING, out, IDENTITY, str(carriers), str(steps), str(parity), *flags, *extra], out,
                6 * 3600, 0, GPU_ENV if backend == 'gpu' else ENV)
    files = {}
    for f in sorted(out.iterdir()):
        if f.is_file():
            files[f.name] = sha(f)
    rec = {'schema': 'stage1-v5-series-receipt-v1', 'kind': kind, 'label': label, 'backend': backend, 'carriers': carriers,
           'native': build_receipt.get('native'), 'extra_flags': list(extra),
           'max_steps': steps, 'parity_every': parity, 'command': r, 'admission': adm, 'sealed': pins,
           'binary_sha256': build_receipt['binary_sha256'], 'files_sha256': files}
    (out / 'series-receipt.json').write_text(json.dumps(rec, indent=1) + '\n')
    print(json.dumps({'run': str(out), 'returncode': r['returncode']}))
    return 0 if r['returncode'] == 0 else 1


if __name__ == '__main__':
    if len(sys.argv) < 3 or sys.argv[1] not in ('build', 'cpu', 'gpu', 'multistep', 'hall-loop', 'points') or not sys.argv[2] or '/' in sys.argv[2]:
        raise SystemExit('usage: run_v5.py build|cpu|gpu|multistep|hall-loop LABEL [...]')
    if sys.argv[1] == 'build':
        raise SystemExit(build(sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else '003'))
    if sys.argv[1] == 'cpu':
        raise SystemExit(run(sys.argv[2], float(sys.argv[3]) if len(sys.argv) > 3 else 900.0,
                             float(sys.argv[4]) if len(sys.argv) > 4 else 1e18))
    if sys.argv[1] == 'points':
        raise SystemExit(run_points(sys.argv[2], int(sys.argv[3]), int(sys.argv[4])))
    if sys.argv[1] == 'gpu':
        raise SystemExit(run_gpu(sys.argv[2], float(sys.argv[3]) if len(sys.argv) > 3 else 1e18))
    if sys.argv[3] not in ('cpu', 'gpu'):
        raise SystemExit('backend must be cpu or gpu')
    raise SystemExit(run_series(sys.argv[1], sys.argv[2], sys.argv[3], int(float(sys.argv[4])), int(sys.argv[5]),
                                int(sys.argv[6]) if len(sys.argv) > 6 else 0, tuple(sys.argv[7:])))
