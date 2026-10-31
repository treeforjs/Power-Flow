#!/usr/bin/python -I
"""Stage 1/v5 acceptance (T1-T7) of the unchanged coordinator v5 sources linked against Hall corrective v3.1
(CPU backends), and the comparison with the coordinator v5 run on Hall corrective v3 (run-p3-221714): accepted
physical state bit for bit, every digest and every non-digest value.

  build:   core-lease-20261030/both-cores.sh /usr/bin/python -I run_coordinator_parity.py build LABEL
  run:     core-lease-20261030/any-core.sh --wait /usr/bin/python -I run_coordinator_parity.py run LABEL [SWEEP_BUDGET_S]

Coordinator v5 (runtimes/pic-water-common-commit-v5-20261030) is read in place and compiled unchanged; the only
difference to its own build-c2 is the Hall owner: -isystem hall-corrective-v3.1/source, -include v2_compat.hpp
(namespace hall_corrective_v2 = hall_corrective_v3) and the v3 motion provider / GPU objects. Native candidate003.
Reference: Hall v3 coordinator parity run-p3-221714 (coordinator v5 on Hall v3, candidate003).
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import fcntl, hashlib, json, os, re, resource, subprocess, sys, time

P = Path(__file__).resolve().parent
R = Path('/home/treefor/projects/runtimes')
C5 = R / 'pic-water-common-commit-v5-20261030'
OUT = Path('/home/treefor/projects/.workspace/outputs/hall-corrective-v3.1-20261031/coordinator-parity')
REF = Path('/home/treefor/projects/.workspace/outputs/hall-corrective-v3-20261030/coordinator-parity/run-p3-221714')
N = R / 'pic-reaction-transactions-20261029/native-full-time-candidate003'
H = R / 'hall-coupled-transport-20261029/source'
HD = H / 'development/hall-canonical-driver-20261029'
I = R / 'pic-water-full-time-subflow-20261029/inputs'
V3 = R / 'pic-water-common-commit-v3-20261030'
M = R / 'mfem-4.10-pic-20261029'
G = R / 'hall-gpu-resident-v1-20261030'
F = R / 'arch-base-20261028/run-frozen'
LEASE_DIR = Path('/home/treefor/projects/.workspace/outputs/pic-profile')
LEASES = {6: LEASE_DIR / 'testing-core.lock', 7: LEASE_DIR / 'testing-core7.lock', 3: LEASE_DIR / 'testing-core3.lock'}   # core 3: temporary loan, CPU-only
NVRTC_INC = R / 'pic-cuda-gather-investigation-20261029/tools/nvrtc/nvidia/cuda_nvrtc/include'
CUDA_INC = R / 'mfem-cuda-benchmark-20261029/cuda-libs/nvidia/cuda_runtime/include'
IDENTITY = 'c045c42d2f57082e3fb5e12172c440d2470e681a38df9ac406bd7e063609a08d'
NATIVE_LIB = N / 'build/libpic_native_candidate003.a'
TRACE = H / 'development/hall-finite-foil-em-20261029/data/shot_17133_current.csv'
DATA = R / 'pic-water-exact-count-subflow-20261029/inputs/h2o-gas.json'
BINDING = Path('/home/treefor/projects/.workspace/outputs/pic-water-exact-count-subflow-20261029/run-binding.json')
LIBS = [M / 'install-mpi/lib/libmfem.a', M / 'tpl/hypre/lib/libHYPRE.a', M / 'tpl/metis/lib/libmetis.a', M / 'tpl/GKlib/lib/libGKlib.a']
NATIVE_WRAPS = ['_ZN6raikou12mfem_backend6GapPic12solve_fieldsEv', '_ZN6raikou12mfem_backend6GapPic7advanceEd',
                '_ZN6raikou12mfem_backend6GapPic13advance_emptyEd',
                '_ZN6raikou12mfem_backend6GapPic16advance_observedEdRKSt8functionIFvRKNS0_11GapParticleEidEE']
REJECT_EXPECT = {1: 'read-only object', 2: 'read-only object', 3: 'read-only object', 4: 'private within this context',
                 5: 'private within this context', 6: 'use of deleted function', 7: 'try_rollback_closed_parent_joint_v5',
                 8: 'incomplete type', 9: 'private within this context'}
ENV = {'PATH': '/usr/bin:/bin', 'LC_ALL': 'C', 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1',
       'BLIS_NUM_THREADS': '1', 'NUMEXPR_NUM_THREADS': '1', 'VECLIB_MAXIMUM_THREADS': '1', 'OMP_THREAD_LIMIT': '1',
       'CUDA_VISIBLE_DEVICES': '', 'PYTHONNOUSERSITE': '1'}
HEX = re.compile(r'^[0-9a-f]{64}$')


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def held(lock):
    with open(lock, 'a') as probe:
        try:
            fcntl.flock(probe, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return True
        fcntl.flock(probe, fcntl.LOCK_UN)
        return False


def admission(kind):
    cores = os.sched_getaffinity(0)
    if not (cores <= {6, 7} or cores == {3}) or os.getpriority(os.PRIO_PROCESS, 0) < 15:
        raise RuntimeError('affinity/nice not admitted')
    for c in cores:
        if not held(LEASES[c]):
            raise RuntimeError('core lease not held')
    if kind == 'run':
        resource.setrlimit(resource.RLIMIT_AS, (3 * 1024 ** 3, 3 * 1024 ** 3))
    return {'kind': kind, 'cpu_affinity': sorted(cores), 'nice': os.getpriority(os.PRIO_PROCESS, 0), 'frozen_launcher_sha256': sha(F)}


def tree(paths):
    out = {}
    for root in paths:
        root = Path(root)
        for p in ([root] if root.is_file() else sorted(x for x in root.rglob('*') if x.is_file() and 'build' not in x.relative_to(root).parts)):
            out[str(p)] = sha(p)
    return out


def execute(name, cmd, log_dir, timeout, expected=0):
    log = log_dir / (name + '.log')
    t = time.monotonic()
    with log.open('w') as h:
        try:
            rc = subprocess.run([str(F), *map(str, cmd)], cwd=C5, env=ENV, stdout=h, stderr=subprocess.STDOUT, timeout=timeout).returncode
        except subprocess.TimeoutExpired:
            rc = 'timeout'
    ok = (rc == expected) if expected != 'nonzero' else (isinstance(rc, int) and rc != 0)
    print(name, rc, round(time.monotonic() - t, 1), flush=True)
    return {'name': name, 'argv': [str(F), *map(str, cmd)], 'returncode': rc, 'ok': ok, 'elapsed_s': time.monotonic() - t,
            'log': str(log), 'log_sha256': sha(log)}


def build(label):
    adm = admission('build')
    out = OUT / ('build-' + label)
    out.mkdir(parents=True, exist_ok=False)
    b = P / 'build' / ('coordinator-' + label)
    b.mkdir(parents=True, exist_ok=False)
    sealed_before = tree([C5 / 'source', C5 / 'tests', V3 / 'source'])
    own_before = tree([P / 'source', Path(__file__).resolve()])
    common = ['-std=c++20', '-O2', '-ffp-contract=off', '-fno-fast-math', '-ffunction-sections', '-fdata-sections']
    isystem = []
    for d in (N / 'source/include', I, I / 'include', N / 'inputs', P / 'source', H / 'development', G, G / 'push',
              M / 'install-mpi/include', M / 'tpl/hypre/include', M / 'tpl/metis/include'):
        isystem += ['-isystem', str(d)]
    gpu_inc = ['-isystem', str(CUDA_INC), '-isystem', str(NVRTC_INC)]
    ours = ['/usr/bin/mpic++', *common, '-Wall', '-Wextra', '-Wpedantic', '-Werror', '-I' + str(C5 / 'source'), '-I' + str(V3 / 'source'),
            '-DPIC_REACTION_NATIVE_SOURCE_SHA256="' + IDENTITY + '"', *isystem, *gpu_inc, '-include', str(P / 'source/v2_compat.hpp')]
    hall = ['/usr/bin/c++', *common, '-DNDEBUG', *isystem, *gpu_inc, '-Wno-c++26-extensions']
    objects = {'hall_motion_v3': (hall, P / 'source/motion_provider.cpp'),
               'kinetic_boundary': (hall, HD / 'kinetic_boundary.cpp'),
               'hall_gpu_backends_v3': (hall, P / 'source/gpu/gpu_backends.cpp'),
               'hall_gpu_push_v3': (hall, P / 'source/gpu/hall_gpu_push_v3.cpp'),
               'hall_gpu_common': (hall, G / 'common/hall_gpu_common.cpp'),
               'hall_gpu_primitives': (hall, G / 'common/hall_gpu_primitives.cpp'),
               'hall_motion_adapter': (['/usr/bin/mpic++', *common, '-Wall', '-Wextra', '-Wpedantic', '-Werror', '-I' + str(C5 / 'source'),
                                        '-I' + str(V3 / 'source'), '-DPIC_REACTION_NATIVE_SOURCE_SHA256="' + IDENTITY + '"', *isystem],
                                       N / 'source/src/mfem_pic_hall_motion_adapter.cpp'),
               'alloc_hook': (ours, C5 / 'tests/alloc_hook.cpp'),
               'test_stage1_v5': (ours, C5 / 'tests/test_stage1_v5.cpp')}
    receipt = {'schema': 'hall-v3.1-coordinator-parity-build-v1', 'label': label, 'admission': adm, 'commands': [],
               'decision': 'coordinator v5 sources compiled unchanged; Hall owner = hall-corrective-v3.1 via -include v2_compat.hpp; native candidate003'}
    status = 'passed'
    try:
        jobs = [('compile-' + k, [*f, '-MD', '-MF', b / (k + '.d'), '-c', src, '-o', b / (k + '.o')]) for k, (f, src) in objects.items()]
        with ThreadPoolExecutor(max_workers=2 if len(adm['cpu_affinity']) == 2 else 1) as pool:
            res = list(pool.map(lambda j: execute(j[0], j[1], out, 3600), jobs))
        receipt['commands'] += res
        for r in res:
            if not r['ok']:
                raise RuntimeError(r['name'] + ' failed:\n' + Path(r['log']).read_text()[-8000:])
        exe = b / 'test_stage1_v5'
        link = ['/usr/bin/mpic++', '-rdynamic', *['-Wl,--wrap=' + s for s in NATIVE_WRAPS], '-Wl,--gc-sections',
                *[b / (k + '.o') for k in objects], NATIVE_LIB, *LIBS, '-lcrypto', '-ldl', '-o', exe]
        r = execute('link', link, out, 1800)
        receipt['commands'].append(r)
        if not r['ok']:
            raise RuntimeError('link failed:\n' + Path(r['log']).read_text()[-8000:])
        decl = [('declaration', [*ours, '-fsyntax-only', C5 / 'tests/declaration_v5.cpp'], 0)]
        decl += [('reject-%d' % k, [*ours, '-DREJECT_CASE=%d' % k, '-fsyntax-only', C5 / 'tests/declaration_v5.cpp'], 'nonzero') for k in sorted(REJECT_EXPECT)]
        with ThreadPoolExecutor(max_workers=2 if len(adm['cpu_affinity']) == 2 else 1) as pool:
            res = list(pool.map(lambda j: execute(j[0], j[1], out, 900, j[2]), decl))
        rejects = {}
        for r in res:
            if r['name'].startswith('reject-'):
                k = int(r['name'].split('-')[1])
                r['diagnostic_found'] = REJECT_EXPECT[k] in Path(r['log']).read_text()
                r['ok'] = r['ok'] and r['diagnostic_found']
                rejects[k] = r['ok']
        receipt['commands'] += res
        receipt['declaration'] = {'compiles': res[0]['ok'], 'reject_cases': rejects}
        deps = set()
        for k in objects:
            for line in (b / (k + '.d')).read_text().replace('\\\n', ' ').splitlines():
                if ':' in line:
                    deps |= {str(Path(x).resolve()) for x in line.split(':', 1)[1].split()}
        rogue = sorted(d for d in deps if 'hall-corrective-v2-20261030' in d or 'hall-corrective-v3-20261030' in d or 'hall-field-correction-20261030' in d or
                       d.endswith('hall-canonical-driver-20261029/motion_receipt.hpp'))
        receipt['hall_v2_v3_or_sealed_motion_headers_in_link'] = rogue
        if rogue:
            raise RuntimeError('Hall v2 / sealed headers in the link: ' + ', '.join(rogue))
        receipt['binary'] = str(exe)
        receipt['binary_sha256'] = sha(exe)
        receipt['coordinator_v5_sources_sha256'] = sealed_before
        receipt['hall_v3_1_sources_sha256'] = own_before
        if tree([C5 / 'source', C5 / 'tests', V3 / 'source']) != sealed_before or tree([P / 'source', Path(__file__).resolve()]) != own_before:
            raise RuntimeError('sources changed during the build')
        if not receipt['declaration']['compiles'] or not all(rejects.values()):
            status = 'declaration-failed'
    except Exception as e:
        status = 'failed'
        receipt['error'] = str(e)
    receipt['status'] = status
    (out / 'build-receipt.json').write_text(json.dumps(receipt, indent=1) + '\n')
    print(json.dumps({'build': str(out), 'status': status, 'error': receipt.get('error', '')[:3000]}))
    return 0 if status == 'passed' else 1


# ---------------------------------------------------------------- comparison (Hall v2 run -> Hall v3 run)
INTENTIONAL_TEXT = set()   # Hall v3 -> v3.1 CPU backends: no intentional text change
EMBEDDED_JSON = {'chemistry_checkpoint'}


class Diff:
    def __init__(self):
        self.mapping, self.roles, self.issues, self.texts, self.added = {}, {}, [], set(), set()

    def walk(self, a, b, path):
        if isinstance(a, dict) and isinstance(b, dict):
            for k in a:
                if k not in b:
                    self.issues.append(['missing key', path + '/' + k])
                elif k in EMBEDDED_JSON and isinstance(a[k], str) and isinstance(b[k], str):
                    self.walk(json.loads(a[k]), json.loads(b[k]), path + '/' + k)
                else:
                    self.walk(a[k], b[k], path + '/' + k)
            for k in b:
                if k not in a:
                    self.added.add(path.split('/')[-1] + '/' + k)
        elif isinstance(a, list) and isinstance(b, list):
            if len(a) != len(b):
                self.issues.append(['length', path, len(a), len(b)])
                return
            for i, (x, y) in enumerate(zip(a, b)):
                self.walk(x, y, '%s[%d]' % (path, i))
        elif isinstance(a, str) and isinstance(b, str) and HEX.match(a) and HEX.match(b):
            if a != b:
                self.mapping.setdefault(a, set()).add(b)
                self.roles.setdefault(a, set()).add(re.sub(r'\[\d+\]', '[]', path))
        elif type(a) is not type(b) or a != b:
            if isinstance(a, str) and (a, b) in INTENTIONAL_TEXT:
                self.texts.add((a, b))
            else:
                self.issues.append(['value', path, repr(a)[:200], repr(b)[:200]])


def load(p):
    return json.loads(Path(p).read_text())


def compare(old, new, out):
    d = Diff()
    report = {'schema': 'hall-v3.1-coordinator-parity-mapping-v1', 'hall_v3_run': str(old), 'hall_v3_1_run': str(new), 'populations': {}}
    ok = True
    for label in ('P1', 'P2', 'P1-stage1-z1mm', 'P2-stage1-z1mm'):
        pop = {}
        ho, hn = old / (label + '-after-owner-bytes/hall-save'), new / (label + '-after-owner-bytes/hall-save')
        co, cn = load(ho / 'checkpoint.json'), load(hn / 'checkpoint.json')
        pop['hall_field_words_identical'] = co['field'] == cn['field']
        pop['hall_charge_ledger_identical'] = co['charge'] == cn['charge']
        heat = sorted(p.name for p in (ho / 'heat').iterdir())
        pop['heat_slices_byte_identical'] = heat == sorted(p.name for p in (hn / 'heat').iterdir()) and all(
            (ho / 'heat' / f).read_bytes() == (hn / 'heat' / f).read_bytes() for f in heat)
        pop['heat_slices'] = len(heat)
        pop['v3_kinetic_ledgers_empty'] = all(t['count'] == 0 for t in cn.get('kinetic_ledgers', {}).get('sinks', [])) and cn.get('kinetic_ledgers', {}).get('loads') == 0
        for key in ('configuration_sha256', 'field_sha256', 'heat_sha256', 'motion_mesh_sha256', 'schema'):
            d.walk(co[key], cn[key], label + '/hall-checkpoint/' + key)
        d.walk(load(old / (label + '-after-owner-bytes/native-save.json')), load(new / (label + '-after-owner-bytes/native-save.json')), label + '/native-save')
        d.walk(load(old / (label + '-trial-audit.json')), load(new / (label + '-trial-audit.json')), label + '/trial-audit')
        ro, rn = load(old / 't1-receipt.json')[label], load(new / 't1-receipt.json')[label]
        for k in ro:
            if k in ('before', 'after'):
                a, b = dict(ro[k]), dict(rn[k])
                for x in (a, b):
                    x.pop('hall_save_sha256', None)
                d.walk(a, b, label + '/t1-receipt/' + k)
            elif k == 'publish':
                a = {x: y for x, y in ro[k].items() if not x.startswith('allocations_') or x == 'allocations_between_first_and_last_swap'}
                b = {x: y for x, y in rn[k].items() if not x.startswith('allocations_') or x == 'allocations_between_first_and_last_swap'}
                d.walk(a, b, label + '/t1-receipt/' + k)
            elif k in ('timings', 'phase_times', 'hall_timings') or k.endswith('_s') or k.endswith('_ms'):
                continue
            else:
                d.walk(ro[k], rn[k], label + '/t1-receipt/' + k)
        ok = ok and pop['hall_field_words_identical'] and pop['hall_charge_ledger_identical'] and pop['heat_slices_byte_identical'] and pop['v3_kinetic_ledgers_empty']
        report['populations'][label] = pop
    inconsistent = {k: sorted(v) for k, v in d.mapping.items() if len(v) != 1}
    inverse = {}
    for k, v in d.mapping.items():
        for x in v:
            inverse.setdefault(x, set()).add(k)
    collisions = {k: sorted(v) for k, v in inverse.items() if len(v) != 1}
    report['identity_mapping'] = sorted(({'old': k, 'new': sorted(v)[0], 'roles': sorted(d.roles[k])[:12], 'role_count': len(d.roles[k])}
                                         for k, v in d.mapping.items()), key=lambda r: r['roles'][0])
    report['mapping_size'] = len(d.mapping)
    report['mapping_inconsistent'] = inconsistent
    report['mapping_collisions'] = collisions
    report['intentional_text_changes'] = sorted(d.texts)
    report['added_keys'] = sorted(d.added)
    report['non_digest_differences'] = d.issues[:200]
    report['non_digest_difference_count'] = len(d.issues)
    ok = ok and not d.issues and not inconsistent and not collisions
    report['physical_state_bit_identical_and_mapping_consistent'] = ok
    out.write_text(json.dumps(report, indent=1) + '\n')
    return ok, report


def run(label, budget):
    adm = admission('run')
    br = load(OUT / ('build-' + label) / 'build-receipt.json')
    if br['status'] != 'passed' or sha(br['binary']) != br['binary_sha256']:
        raise RuntimeError('build not usable')
    out = OUT / ('run-' + label + '-' + time.strftime('%H%M%S'))
    out.mkdir(parents=True)
    exe = br['binary']
    summary = {'schema': 'hall-v3.1-coordinator-parity-run-v1', 'label': label, 'admission': adm, 'binary_sha256': br['binary_sha256'],
               'reference_run': str(REF), 'runs': {}}
    for name, mode, arg in (('t1-a', 't1', 1e18), ('t1-b', 't1', 1e18), ('matrix', 'matrix', float(budget))):
        d = out / name
        d.mkdir()
        r = execute(name, [exe, mode, TRACE, DATA, BINDING, d, IDENTITY, repr(arg)], out, 4 * 3600)
        res = d / ('t1-result.json' if mode == 't1' else 'matrix-result.json')
        summary['runs'][name] = {'command': r, 'summary': load(res).get('summary') if res.exists() else None}
    a, b = out / 't1-a/t1-receipt.json', out / 't1-b/t1-receipt.json'
    summary['t1_two_fresh_processes_byte_identical'] = a.exists() and b.exists() and a.read_bytes() == b.read_bytes()
    if (out / 't1-a').exists():
        ok, rep = compare(REF / 't1-a', out / 't1-a', out / 'v3-to-v3.1-identity-mapping.json')
        summary['v3_parity'] = {'physical_state_bit_identical_and_mapping_consistent': ok, 'populations': rep['populations'],
                                'mapping_size': rep['mapping_size'], 'non_digest_difference_count': rep['non_digest_difference_count'],
                                'first_issues': rep['non_digest_differences'][:10], 'intentional_text_changes': rep['intentional_text_changes'],
                                'added_keys': rep['added_keys'], 'mapping': str(out / 'v3-to-v3.1-identity-mapping.json')}
    per = {}
    for name, v in summary['runs'].items():
        for t, s in (v.get('summary') or {}).items():
            x = per.setdefault(t, {'passed': 0, 'failed': 0})
            x['passed'] += s.get('passed', 0)
            x['failed'] += s.get('failed', 0)
    summary['per_test'] = per
    summary['declaration_from_build'] = br.get('declaration')
    summary['status'] = 'passed' if all(v['command']['returncode'] == 0 for v in summary['runs'].values()) and \
        summary['t1_two_fresh_processes_byte_identical'] and summary.get('v3_parity', {}).get('physical_state_bit_identical_and_mapping_consistent') else 'completed-with-failures'
    (out / 'parity-receipt.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps({'run': str(out), 'status': summary['status'], 'per_test': per}))
    return 0 if summary['status'] == 'passed' else 1


if __name__ == '__main__':
    a = sys.argv
    if len(a) < 3 or a[1] not in ('build', 'run') or '/' in a[2]:
        raise SystemExit(__doc__)
    if a[1] == 'build':
        raise SystemExit(build(a[2]))
    raise SystemExit(run(a[2], a[3] if len(a) > 3 else 300))
