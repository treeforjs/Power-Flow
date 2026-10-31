#!/usr/bin/python -I
"""Build and run the heat/source owner v6 binding qualification (H5, with Hall corrective v3).

  build: core-lease-20261030/both-cores.sh /usr/bin/python -I run_source_v6.py build LABEL [ubsan]
  run:   core-lease-20261030/any-core.sh --wait /usr/bin/python -I run_source_v6.py run LABEL [UBSAN_LABEL]
  growth: core-lease-20261030/any-core.sh --wait /usr/bin/python -I run_source_v6.py growth LABEL T1:STEPS1,T2:STEPS2 W0

Every compile and run step executes inside arch-base-20261028/run-frozen. Outputs:
/home/treefor/projects/.workspace/outputs/hall-corrective-v3-20261030/source-v6/LABEL/. Sealed inputs (source owner v5,
heat/source v3 kernel, Hall 8dddec0 heat-v2 bindings, native candidate003) are hashed before and after and only read.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import fcntl, hashlib, json, os, platform, subprocess, sys, time

P = Path(__file__).resolve().parent
R = Path('/home/treefor/projects/runtimes')
OUT = Path('/home/treefor/projects/.workspace/outputs/hall-corrective-v3-20261030/source-v6')
N = R / 'pic-reaction-transactions-20261029/native-full-time-candidate003'
H = R / 'hall-coupled-transport-20261029/source'
HD = H / 'development/hall-canonical-driver-20261029'
HB = HD / 'bindings/heat-v2'
I = R / 'pic-water-full-time-subflow-20261029/inputs'
M = R / 'mfem-4.10-pic-20261029'
G = R / 'hall-gpu-resident-v1-20261030'
HV3 = R / 'hall-corrective-v3-20261030'
S5 = R / 'pic-electrode-heat-source-v5-20261030/native/mfem_electrode_heat_source_v5'
V3K = R / 'pic-electrode-heat-feedback-20261029/native/mfem_electrode_heat_source_v3'
A1 = R / 'pic-arrival-owner-v1-20261030'
A11 = R / 'pic-arrival-owner-v1.1-20261030'
F = R / 'arch-base-20261028/run-frozen'
LEASE_DIR = Path('/home/treefor/projects/.workspace/outputs/pic-profile')
LEASES = {6: LEASE_DIR / 'testing-core.lock', 7: LEASE_DIR / 'testing-core7.lock'}
NATIVE_LIB = N / 'build/libpic_native_candidate003.a'
LIBS = [M / 'install-mpi/lib/libmfem.a', M / 'tpl/hypre/lib/libHYPRE.a', M / 'tpl/metis/lib/libmetis.a', M / 'tpl/GKlib/lib/libGKlib.a']
TRACE = H / 'development/hall-finite-foil-em-20261029/data/shot_17133_current.csv'
ENV = {'PATH': '/usr/bin:/bin', 'LC_ALL': 'C', 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1', 'PYTHONNOUSERSITE': '1'}
COMMON = ['-std=c++20', '-O2', '-DNDEBUG', '-ffp-contract=off', '-fno-fast-math']
UBSAN = ['-std=c++20', '-O1', '-g', '-DNDEBUG', '-ffp-contract=off', '-fno-fast-math', '-fsanitize=undefined', '-fno-sanitize-recover=all', '-fno-omit-frame-pointer']


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
    cores = os.sched_getaffinity(0)
    if not cores <= {6, 7} or os.getpriority(os.PRIO_PROCESS, 0) < 15:
        raise RuntimeError('only cores 6/7 at nice >= 15 are admitted')
    for c in cores:
        if not held(LEASES[c]):
            raise RuntimeError('core %d lease is not held' % c)
    return {'stage': stage, 'cpu_affinity': sorted(cores), 'nice': os.getpriority(os.PRIO_PROCESS, 0), 'python': platform.python_version(),
            'frozen_launcher_sha256': sha(F), 'time_unix': time.time()}


def sealed():
    files = [p for d in (S5, V3K, HB, A1 / 'include') for p in sorted(d.rglob('*')) if p.is_file() and p.suffix in ('.hpp', '.h', '.cpp', '.json', '.txt')]
    return {str(p): sha(p) for p in files + [NATIVE_LIB]}


def own():
    files = sorted(p for d in (P, A11) for p in d.rglob('*') if p.is_file() and p.suffix in ('.hpp', '.cpp', '.py', '.txt') and 'build' not in p.parts)
    return {str(p): sha(p) for p in files}


def execute(name, cmd, log_dir, timeout):
    log = log_dir / (name + '.log')
    t = time.monotonic()
    with log.open('w') as h:
        try:
            rc = subprocess.run([str(F), *map(str, cmd)], cwd=P, env=ENV, stdout=h, stderr=subprocess.STDOUT, timeout=timeout).returncode
        except subprocess.TimeoutExpired:
            rc = 'timeout'
    print(name, rc, round(time.monotonic() - t, 1), flush=True)
    return {'name': name, 'argv': [str(F), *map(str, cmd)], 'returncode': rc, 'elapsed_s': time.monotonic() - t, 'log': str(log), 'log_sha256': sha(log)}


def flags(base):
    inc = ['-I' + str(R), '-I' + str(P / 'binding'), '-I' + str(HV3 / 'source'), '-I' + str(G), '-I' + str(HB), '-I' + str(V3K),
           '-I' + str(A11 / 'include'), '-I' + str(A1 / 'include'),
           '-DHEAT_SOURCE_EXACT_COUNT_HEADER="raikou/exact_particle_count.h"', '-DPF_ARRIVAL_EXACT_COUNT_HEADER="raikou/exact_particle_count.h"']
    for d in [N / 'source/include', I, I / 'include', N / 'inputs', H / 'development', M / 'install-mpi/include', M / 'tpl/hypre/include', M / 'tpl/metis/include']:
        inc += ['-isystem', str(d)]
    return ['/usr/bin/c++', *base, *inc]


def build(label, variant):
    adm = admission('build')
    out = OUT / label / 'build'
    out.mkdir(parents=True, exist_ok=False)
    b = P / 'build' / label
    b.mkdir(parents=True, exist_ok=False)
    pins, before = sealed(), own()
    base = UBSAN if variant == 'ubsan' else COMMON
    objs = {'test_source_v6': (flags(base), P / 'tests/test_source_v6.cpp'),
            'motion_provider_v3': (flags(base), HV3 / 'source/motion_provider.cpp'),
            'kinetic_boundary': (flags(base), HD / 'kinetic_boundary.cpp'),
            'surface_emission': (['/usr/bin/c++', *base, '-I' + str(HB / 'source-kernel/include')], HB / 'source-kernel/src/surface_emission.cpp'),
            'micro_layer': (['/usr/bin/c++', *base, '-I' + str(HB / 'source-kernel/include')], HB / 'source-kernel/src/micro_layer.cpp')}
    receipt = {'schema': 'source-v6-build-receipt-v1', 'label': label, 'variant': variant, 'admission': adm, 'sealed_sha256': pins, 'commands': []}
    status = 'passed'
    try:
        jobs = [('compile-' + k, [*f, '-MD', '-MF', b / (k + '.d'), '-c', src, '-o', b / (k + '.o')]) for k, (f, src) in objs.items()]
        with ThreadPoolExecutor(max_workers=2 if len(adm['cpu_affinity']) == 2 else 1) as pool:
            results = list(pool.map(lambda j: execute(j[0], j[1], out, 3600), jobs))
        receipt['commands'] += results
        for r in results:
            if r['returncode'] != 0:
                raise RuntimeError(r['name'] + ' failed:\n' + Path(r['log']).read_text()[-8000:])
        link = ['/usr/bin/mpic++', '-Wl,--gc-sections', *(['-fsanitize=undefined'] if variant == 'ubsan' else []),
                *[b / (k + '.o') for k in objs], NATIVE_LIB, *LIBS, '-lcrypto', '-o', b / 'test_source_v6']
        r = execute('link', link, out, 1800)
        receipt['commands'].append(r)
        if r['returncode'] != 0:
            raise RuntimeError('link failed:\n' + Path(r['log']).read_text()[-8000:])
        deps = set()
        for k in objs:
            for line in (b / (k + '.d')).read_text().replace('\\\n', ' ').splitlines():
                if ':' in line:
                    deps |= {str(Path(x).resolve()) for x in line.split(':', 1)[1].split()}
        receipt['dependency_sha256'] = {d: sha(d) for d in sorted(deps)}
        heat = sorted(d for d in deps if d.endswith('electrode_heat_feedback.hpp'))
        receipt['heat_v2_engine_headers'] = heat
        if len(heat) != 1:
            raise RuntimeError('more than one heat-v2 engine header in the link: ' + ', '.join(heat))
        receipt['binary_sha256'] = sha(b / 'test_source_v6')
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


def run(label, ubsan_label=None):
    adm = admission('run')
    out = OUT / label / ('run-' + time.strftime('%H%M%S'))
    out.mkdir(parents=True)
    summary = {'schema': 'source-v6-run-receipt-v1', 'label': label, 'admission': adm, 'sealed_sha256': sealed(), 'commands': []}
    for lab, sub in ((label, 'normal'), (ubsan_label, 'ubsan')):
        if not lab:
            continue
        br = json.loads((OUT / lab / 'build/build-receipt.json').read_text())
        if br['status'] != 'passed' or sha(Path(br['build_dir']) / 'test_source_v6') != br['binary_sha256']:
            raise RuntimeError('build %s not usable' % lab)
        d = out / sub
        d.mkdir()
        r = execute(sub, [Path(br['build_dir']) / 'test_source_v6', TRACE, d], out, 7200)
        summary['commands'].append(r)
        f = d / 'source-v6.json'
        if f.exists():
            j = json.loads(f.read_text())
            summary[sub] = {'failed': j['failed'], 'passed': sum(1 for c in j['cases'] if c['passed']), 'cases': len(j['cases']), 'error': j.get('error'),
                            'names': [{'group': c['group'], 'name': c['name'], 'passed': c['passed']} for c in j['cases']]}
    if ubsan_label and 'normal' in summary and 'ubsan' in summary:
        summary['ubsan_runtime_errors'] = (out / 'ubsan.log').read_text().count('runtime error')
        summary['ubsan_same_outcomes'] = [c['passed'] for c in summary['normal']['names']] == [c['passed'] for c in summary['ubsan']['names']]
    summary['all_passed'] = all(c['returncode'] == 0 for c in summary['commands'])
    (out / 'run-receipt.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps({'run': str(out), 'all_passed': summary['all_passed']}))
    return 0 if summary['all_passed'] else 1


def growth(label, specs, w0):
    adm = admission('run')
    br = json.loads((OUT / label / 'build/build-receipt.json').read_text())
    if br['status'] != 'passed' or sha(Path(br['build_dir']) / 'test_source_v6') != br['binary_sha256']:
        raise RuntimeError('build not usable')
    out = OUT / label / ('growth-' + time.strftime('%H%M%S'))
    out.mkdir(parents=True)
    rec = {'schema': 'source-v6-growth-receipt-v1', 'admission': adm, 'w0': w0, 'runs': []}
    for spec in specs.split(','):
        T, steps = spec.split(':')
        r = execute('growth-%sK' % T, [Path(br['build_dir']) / 'test_source_v6', TRACE, out, 'growth', T, steps, w0], out, 6 * 3600)
        f = out / ('growth-%sK.json' % T)
        rec['runs'].append({'command': r, 'result': json.loads(f.read_text()) if f.exists() else None})
    (out / 'growth-receipt.json').write_text(json.dumps(rec, indent=1) + '\n')
    print(json.dumps({'run': str(out), 'ok': all(x['command']['returncode'] == 0 for x in rec['runs'])}))
    return 0


if __name__ == '__main__':
    a = sys.argv
    if len(a) < 3 or a[1] not in ('build', 'run', 'growth') or '/' in a[2]:
        raise SystemExit(__doc__)
    if a[1] == 'build':
        raise SystemExit(build(a[2], a[3] if len(a) > 3 else 'normal'))
    if a[1] == 'growth':
        raise SystemExit(growth(a[2], a[3], a[4]))
    raise SystemExit(run(a[2], a[3] if len(a) > 3 else None))
