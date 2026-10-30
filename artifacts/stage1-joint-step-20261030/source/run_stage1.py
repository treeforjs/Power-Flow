#!/usr/bin/python -I
"""Stage 1 (first real all-owner step) build and run driver for coordinator v4.

Builds hold both core leases (6,7); runs hold one core lease (6 or 7):
  cd /home/treefor/projects/runtimes/pic-water-common-commit-v4-20261030
  /home/treefor/projects/runtimes/core-lease-20261030/both-cores.sh env CUDA_VISIBLE_DEVICES= \
      /usr/bin/python -I run_stage1.py build LABEL
  /home/treefor/projects/runtimes/core-lease-20261030/any-core.sh --wait env CUDA_VISIBLE_DEVICES= \
      /usr/bin/python -I run_stage1.py run LABEL [SWEEP_BUDGET_S] [WATER_WEIGHT]

Every compile, link and test command runs inside the frozen launcher with
RLIMIT_AS 3 GiB, one thread per numerical library and no GPU. Sealed inputs are
hashed before and after; nothing outside the v4 runtime directory and the v4
output directory is written.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import fcntl, hashlib, json, os, platform, resource, subprocess, sys, time

P = Path(__file__).resolve().parent
R = Path('/home/treefor/projects/runtimes')
OUT = Path('/home/treefor/sim-workspace/powerflow/outputs/pic-water-common-commit-v4-20261030')
N = R / 'pic-reaction-transactions-20261029/native-full-time-candidate003'
H = R / 'hall-coupled-transport-20261029/source'
HC = H / 'development/hall-field-correction-20261030'
HD = H / 'development/hall-canonical-driver-20261029'
I = R / 'pic-water-full-time-subflow-20261029/inputs'
V3 = R / 'pic-water-common-commit-v3-20261030'
M = R / 'mfem-4.10-pic-20261029'
F = R / 'arch-base-20261028/run-frozen'
LEASES = {6: Path('/home/treefor/sim-workspace/powerflow/outputs/pic-profile/testing-core.lock'),
          7: Path('/home/treefor/sim-workspace/powerflow/outputs/pic-profile/testing-core7.lock')}
IDENTITY = 'c045c42d2f57082e3fb5e12172c440d2470e681a38df9ac406bd7e063609a08d'
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
TIMING_WRAPS = ['_ZN11hall_driver24prepare_canonical_motionERKN5axial4MeshERKN6raikou12mfem_backend13GapPicOptionsESt4spanIKNS5_14ArenaRateEntryELm18446744073709551615EERKNS_13MotionBindingES9_IKmLm18446744073709551615EERKSt6vectorIdSaIdEESM_',
                '_ZN11hall_driver18total_field_sha256B5cxx11ERKSt6vectorIdSaIdEES4_',
                '_ZN11hall_driver23canonical_census_sha256B5cxx11ESt4spanIKN6raikou12mfem_backend14ArenaRateEntryELm18446744073709551615EE',
                '_ZN11hall_driver18motion_mesh_sha256B5cxx11ERKN5axial4MeshERKN6raikou12mfem_backend13GapPicOptionsE']
ENV = {'PATH': '/usr/bin:/bin', 'LC_ALL': 'C', 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1',
       'MKL_NUM_THREADS': '1', 'BLIS_NUM_THREADS': '1', 'NUMEXPR_NUM_THREADS': '1', 'VECLIB_MAXIMUM_THREADS': '1',
       'OMP_THREAD_LIMIT': '1', 'CUDA_VISIBLE_DEVICES': '', 'PYTHONNOUSERSITE': '1'}
RLIMIT_AS = 3 * 1024 ** 3
REJECT_EXPECT = {1: 'read-only object', 2: 'read-only object', 3: 'read-only object',
                 4: 'private within this context', 5: 'private within this context',
                 6: 'use of deleted function', 7: 'try_rollback_closed_parent_joint_v4',
                 8: 'incomplete type', 9: 'private within this context'}


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def fail(message):
    raise RuntimeError(message)


def admission(kind):
    if not __debug__ or os.getenv('PYTHONOPTIMIZE'):
        fail('optimized Python guard execution refused')
    cores = os.sched_getaffinity(0)
    allowed = [{6, 7}] if kind == 'build' else [{6}, {7}]
    if cores not in allowed or os.getpriority(os.PRIO_PROCESS, 0) != 15:
        fail('affinity %s / nice %d not admitted for %s' % (sorted(cores), os.getpriority(os.PRIO_PROCESS, 0), kind))
    if any(os.getenv(k) != '1' for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']) \
            or os.getenv('CUDA_VISIBLE_DEVICES') != '':
        fail('one thread per library and no GPU required')
    for core in cores:
        with open(LEASES[core], 'a') as probe:
            try:
                fcntl.flock(probe, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                continue
            fail('core %d lease must be held by the launching wrapper; no work started' % core)
    if Path('/proc/sys/kernel/osrelease').read_text().strip() != (R / 'arch-base-20261028/kernel.release').read_text().strip():
        fail('kernel release differs from the frozen runtime')
    resource.setrlimit(resource.RLIMIT_AS, (RLIMIT_AS, RLIMIT_AS))
    return {'kind': kind, 'cpu_affinity': sorted(cores), 'PF_CORE': os.getenv('PF_CORE'), 'nice': 15,
            'threads': 1, 'gpu': False, 'leases_verified_held': [str(LEASES[c]) for c in sorted(cores)],
            'rlimit_as_bytes': RLIMIT_AS, 'python': sys.executable, 'python_version': platform.python_version(),
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
    return tree_hashes([P / 'source', P / 'tests', Path(__file__).resolve()])


def frozen(cmd):
    return [str(F), *map(str, cmd)]


def execute(name, cmd, log_dir, timeout, expected=0):
    log = log_dir / (name + '.log')
    start = time.monotonic()
    with log.open('w') as handle:
        try:
            r = subprocess.run(frozen(cmd), cwd=P, env=ENV, stdout=handle, stderr=subprocess.STDOUT, timeout=timeout)
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


def build(label):
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
    for d in (N / 'source/include', I, I / 'include', N / 'inputs', HC, H / 'development',
              M / 'install-mpi/include', M / 'tpl/hypre/include', M / 'tpl/metis/include'):
        isystem += ['-isystem', str(d)]
    ours = ['/usr/bin/mpic++', *common, '-Wall', '-Wextra', '-Wpedantic', '-Werror',
            '-I' + str(P / 'source'), '-I' + str(V3 / 'source'),
            '-DPIC_REACTION_NATIVE_SOURCE_SHA256="' + IDENTITY + '"', *isystem]
    # Hall translation units: Hall's own option set (no warning flags, NDEBUG),
    # with the same native003 header resolution as every other TU in the link.
    hall = ['/usr/bin/c++', *common, '-DNDEBUG', *isystem]
    objects = {'native_arena_motion': (hall, HD / 'native_arena_motion.cpp'),
               'kinetic_boundary': (hall, HD / 'kinetic_boundary.cpp'),
               'hall_motion_adapter': (ours, N / 'source/src/mfem_pic_hall_motion_adapter.cpp'),
               'alloc_hook': (ours, P / 'tests/alloc_hook.cpp'),
               'test_stage1_joint': (ours, P / 'tests/test_stage1_joint.cpp')}
    jobs = [('compile-' + k, [*flags, '-MD', '-MF', str(bdir / (k + '.d')), '-c', str(src), '-o', str(bdir / (k + '.o'))])
            for k, (flags, src) in objects.items()]
    receipt = {'schema': 'stage1-v4-build-receipt-v1', 'label': label, 'admission': adm, 'sealed': pins,
               'started_unix': time.time(), 'environment': ENV, 'commands': [],
               'decisions': {
                   'hall_tu_include_resolution': 'Hall native_arena_motion.cpp and kinetic_boundary.cpp use Hall compile options '
                       '(-O2 -DNDEBUG, no warning flags) but resolve raikou/* to native candidate003 headers, so the whole link '
                       'has one definition of ArenaRateEntry/ArenaCount (native-v6.2 headers declare ArenaCount=electrode_heat::Count).',
                   'mfem_gap_pic_not_linked_separately': 'pic_boris and pic_kinetic_energy_J are already defined by native.o inside '
                       'libpic_native_candidate003.a, compiled from mfem_gap_pic.cpp whose bytes (fa1511c2...) equal native-v6.2 '
                       'mfem_gap_pic.cpp; linking that source again would duplicate every GapPic definition.',
                   'timing_wraps': 'four Hall motion-TU functions are link-wrapped only to accumulate wall time and call counts; '
                       'each wrapper calls the real function with unchanged arguments.'}}

    def finish(status, error=None):
        receipt['status'] = status
        receipt['error'] = error
        receipt['finished_unix'] = time.time()
        (out / 'build-receipt.json').write_text(json.dumps(receipt, indent=1) + '\n')
        print(json.dumps({'build': str(out), 'status': status, 'error': error}))

    try:
        with ThreadPoolExecutor(max_workers=2) as pool:
            results = list(pool.map(lambda j: execute(j[0], j[1], out, 1800), jobs))
        receipt['commands'] += results
        for r in results:
            if not r['ok']:
                fail(r['name'] + ' failed:\n' + Path(r['log']).read_text()[-6000:])
        nm = subprocess.run(frozen(['/usr/bin/nm', '--defined-only', bdir / 'native_arena_motion.o']), cwd=P, env=ENV,
                            capture_output=True, text=True, check=True).stdout
        defined = {line.split()[-1] for line in nm.splitlines() if line.strip()}
        missing = [s for s in TIMING_WRAPS if s not in defined]
        if missing:
            fail('timing wrap symbols are not defined by the Hall motion object: ' + ', '.join(missing))
        exe = bdir / 'test_stage1_joint'
        link = ['/usr/bin/mpic++', '-rdynamic', *['-Wl,--wrap=' + s for s in NATIVE_WRAPS + TIMING_WRAPS], '-Wl,--gc-sections',
                *[bdir / (k + '.o') for k in objects], NATIVE_LIB, *LIBS, '-lcrypto', '-o', exe]
        r = execute('link', link, out, 1800)
        receipt['commands'].append(r)
        if not r['ok']:
            fail('link failed:\n' + Path(r['log']).read_text()[-6000:])
        decl = [('declaration', [*ours, '-fsyntax-only', P / 'tests/declaration_v4.cpp'], 0)]
        decl += [('reject-%d' % k, [*ours, '-DREJECT_CASE=%d' % k, '-fsyntax-only', P / 'tests/declaration_v4.cpp'], 'nonzero')
                 for k in sorted(REJECT_EXPECT)]
        with ThreadPoolExecutor(max_workers=2) as pool:
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
        # Dependency closure from the -MD files: every header actually read.
        deps = {}
        for k in objects:
            deps[k] = dependencies(bdir / (k + '.d'))
        allfiles = sorted({d for v in deps.values() for d in v})
        receipt['dependency_sha256'] = {d: sha(d) for d in allfiles}
        json_headers = {k: [d for d in v if Path(d).name in ('json.hpp', 'json.h')] for k, v in deps.items()}
        identities = {d: sha(d) for v in json_headers.values() for d in v}
        hpp = {d for d in identities if d.endswith('json.hpp')}
        receipt['json_header_identity'] = {'per_tu': json_headers, 'sha256': identities,
                                           'single_nlohmann_json_hpp': len(hpp) == 1 and len({identities[d] for d in hpp}) == 1}
        if not receipt['json_header_identity']['single_nlohmann_json_hpp']:
            fail('more than one nlohmann json.hpp in the link: ' + ', '.join(sorted(hpp)))
        rogue = [d for d in allfiles if 'native-v6.2' in d]
        receipt['native_v6_2_headers_in_link'] = rogue
        if rogue:
            fail('native-v6.2 headers leaked into the link: ' + ', '.join(rogue))
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


def run(label, sweep_budget, water_weight):
    adm = admission('run')
    build_receipt = json.loads((OUT / ('build-' + label) / 'build-receipt.json').read_text())
    if build_receipt['status'] != 'passed':
        fail('build did not pass')
    exe = Path(build_receipt['binary'])
    if sha(exe) != build_receipt['binary_sha256']:
        fail('test binary differs from its build receipt')
    pins = sealed_pins()
    if pins != build_receipt['sealed']:
        fail('sealed inputs differ from the build')
    if own_sources() != build_receipt['own_sources_sha256']:
        fail('v4 sources changed after the build')
    out = OUT / ('run-' + label + '-' + str(adm['cpu_affinity'][0]) + '-' + time.strftime('%H%M%S'))
    out.mkdir(exist_ok=False)
    base = [exe, '', TRACE, DATA, BINDING, '', IDENTITY]
    runs = {}
    summary = {'schema': 'stage1-v4-run-summary-v1', 'label': label, 'admission': adm, 'sealed': pins,
               'build_receipt': str(OUT / ('build-' + label) / 'build-receipt.json'),
               'binary_sha256': build_receipt['binary_sha256'], 'started_unix': time.time()}
    status = 'passed'
    try:
        for name, mode, arg, timeout in [('t1-a', 't1', water_weight, 3600), ('t1-b', 't1', water_weight, 3600),
                                         ('matrix', 'matrix', sweep_budget, 4 * 3600)]:
            d = out / name
            d.mkdir()
            cmd = [exe, mode, TRACE, DATA, BINDING, d, IDENTITY, repr(arg)]
            r = execute(name, cmd, out, timeout, 0)
            rec = {'schema': 'stage1-v4-run-receipt-v1', 'command': r, 'cpu_affinity': adm['cpu_affinity'],
                   'nice': 15, 'threads': 1, 'gpu': False, 'rlimit_as_bytes': RLIMIT_AS,
                   'binary_sha256': build_receipt['binary_sha256'], 'inputs_sha256': {
                       str(TRACE): TRACE_SHA, str(DATA): pins['inputs']['water_data_sha256'],
                       str(BINDING): pins['inputs']['binding_sha256']}}
            result = d / ('t1-result.json' if mode == 't1' else 'matrix-result.json')
            if result.exists():
                res = json.loads(result.read_text())
                rec['result_file'] = str(result)
                rec['result_sha256'] = sha(result)
                rec['summary'] = res.get('summary')
                rec['cases'] = res.get('cases')
            for extra in ('t1-receipt.json', 't1-timings.json', 'allocation-sweep.json'):
                if (d / extra).exists():
                    rec.setdefault('files', {})[extra] = {'path': str(d / extra), 'sha256': sha(d / extra)}
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
        if own_sources() != build_receipt['own_sources_sha256'] or sealed_pins() != pins:
            fail('sources or sealed inputs changed during the run')
        summary['inputs_unchanged_after_run'] = True
        if any(v['command']['returncode'] != 0 for v in runs.values()) or not same:
            status = 'completed-with-failures'
    except Exception as exc:
        status = 'failed'
        summary['error'] = str(exc)
    summary['status'] = status
    summary['finished_unix'] = time.time()
    (out / 'stage1-receipt.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps({'run': str(out), 'status': status}))
    return 0 if status == 'passed' else 1


if __name__ == '__main__':
    if len(sys.argv) < 3 or sys.argv[1] not in ('build', 'run') or not sys.argv[2] or '/' in sys.argv[2]:
        raise SystemExit('usage: run_stage1.py build|run LABEL [SWEEP_BUDGET_S] [WATER_WEIGHT]')
    if sys.argv[1] == 'build':
        raise SystemExit(build(sys.argv[2]))
    budget = float(sys.argv[3]) if len(sys.argv) > 3 else 900.0
    weight = float(sys.argv[4]) if len(sys.argv) > 4 else 1e18
    raise SystemExit(run(sys.argv[2], budget, weight))
