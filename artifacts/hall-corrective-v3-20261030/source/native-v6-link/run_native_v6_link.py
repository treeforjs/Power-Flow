#!/usr/bin/python -I
"""Native v6 admission of Hall corrective v3 motion (flows record).

  build: core-lease-20261030/both-cores.sh /usr/bin/python -I run_native_v6_link.py build LABEL
  run:   core-lease-20261030/any-core.sh --wait /usr/bin/python -I run_native_v6_link.py run LABEL [STEPS]

Native candidate006 (identity aa4d8f03...) is read in place: its sealed sources (core arena and canonical
TUs, test fixture v6_fixture.hpp, arrival v1 births adapter) are hashed against SOURCE-SHA256.json and compiled
here with the flags of its build-components.py; the sealed exact-count object native-exact-count/build/native.o
is linked as native does. The Hall corrective v3 CPU motion provider (source/motion_provider.cpp) is compiled
against the native v6 headers (as native compiles the Hall corrective v2 provider). Every compile and run step
executes inside arch-base-20261028/run-frozen. Outputs: /home/treefor/projects/.workspace/outputs/
hall-corrective-v3-20261030/native-v6-link/LABEL/.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import fcntl, hashlib, json, os, subprocess, sys, time

P = Path(__file__).resolve().parent
HV3 = P.parent
R = Path('/home/treefor/projects/runtimes')
OUT = Path('/home/treefor/projects/.workspace/outputs/hall-corrective-v3-20261030/native-v6-link')
N6 = R / 'pic-reaction-transactions-20261029/native-full-time-candidate006'
N6_IDENTITY = 'aa4d8f03100e64bed759c739013d91cac1bd66c8a9b2ab6bc3105af7973452c9'
EXACT_O = R / 'pic-reaction-transactions-20261029/native-exact-count/build/native.o'
VENDOR = R / 'pic-reaction-transactions-20261029/source/vendor'
H = R / 'hall-coupled-transport-20261029/source'
HALL = H / 'development'
HALL_COMMIT = '8dddec0610982778bf02bdfe697b70a37129c012'
ARRIVAL = R / 'pic-arrival-owner-v1-20261030'
ARRIVAL11 = R / 'pic-arrival-owner-v1.1-20261030'
V3_COMMIT = R / 'pic-water-common-commit-v3-20261030'
G = R / 'hall-gpu-resident-v1-20261030'
M = R / 'mfem-4.10-pic-20261029'
F = R / 'arch-base-20261028/run-frozen'
LEASE_DIR = Path('/home/treefor/projects/.workspace/outputs/pic-profile')
LEASES = {6: LEASE_DIR / 'testing-core.lock', 7: LEASE_DIR / 'testing-core7.lock'}
LIBS = [M / 'install-mpi/lib/libmfem.a', M / 'tpl/hypre/lib/libHYPRE.a', M / 'tpl/metis/lib/libmetis.a', M / 'tpl/GKlib/lib/libGKlib.a', '-lcrypto']
ENV = {'PATH': '/usr/bin:/bin', 'LC_ALL': 'C', 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1', 'PYTHONNOUSERSITE': '1'}
SYSTEM = [x for d in ('install-mpi/include', 'tpl/hypre/include', 'tpl/metis/include') for x in ('-isystem', str(M / d))]
INC = ['-I' + str(N6 / 'source/include'), '-isystem', str(N6 / 'inputs'), '-I' + str(VENDOR)] + SYSTEM
FP = ['-ffp-contract=off', '-fno-fast-math', '-ffunction-sections', '-fdata-sections']
STRICT = ['-Wall', '-Wextra', '-Wpedantic', '-Werror']
DEF = ['-DPIC_REACTION_NATIVE_SOURCE_SHA256="' + N6_IDENTITY + '"']


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
    return {'stage': stage, 'cpu_affinity': sorted(cores), 'nice': os.getpriority(os.PRIO_PROCESS, 0), 'frozen_launcher_sha256': sha(F),
            'time_unix': time.time()}


def sealed():
    manifest = json.loads((N6 / 'SOURCE-SHA256.json').read_text())
    if manifest['source_identity'] != N6_IDENTITY:
        fail('native v6 source identity differs')
    for name, digest in manifest['files'].items():
        if sha(N6 / name) != digest:
            fail('native v6 sealed file changed: ' + name)
    tests = {str(p): sha(p) for p in [N6 / 'tests/v6_fixture.hpp']}
    head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=H, text=True).strip()
    dirty = subprocess.check_output(['git', 'status', '--porcelain'], cwd=H, text=True)
    if head != HALL_COMMIT or dirty:
        fail('Hall repository is not the clean corrective commit')
    return {'native_v6_identity': N6_IDENTITY, 'native_v6_files_verified': len(manifest['files']), 'native_v6_test_fixture': tests,
            'exact_count_object_sha256': sha(EXACT_O), 'arrival_v1_owner_sha256': sha(ARRIVAL / 'include/pf_arrival/arrival_owner.hpp'),
            'hall_commit': head}


def own():
    files = sorted(p for d in (HV3 / 'source', P, ARRIVAL11 / 'include') for p in d.rglob('*')
                   if p.is_file() and p.suffix in ('.hpp', '.cpp', '.py') and 'build' not in p.parts)
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


def build(label):
    adm = admission('build')
    out = OUT / label / 'build'
    out.mkdir(parents=True, exist_ok=False)
    b = P / 'build' / label
    b.mkdir(parents=True, exist_ok=False)
    pins, before = sealed(), own()
    cxx = ['/usr/bin/mpic++', '-std=c++20']
    hall_cc = ['/usr/bin/c++', '-std=c++20', '-O2', '-DNDEBUG'] + FP
    test_inc = ['-I' + str(N6 / 'tests'), '-I' + str(ARRIVAL / 'include'), '-I' + str(ARRIVAL11 / 'include'), '-I' + str(N6 / 'source/adapters'),
                '-I' + str(V3_COMMIT / 'source'),
                '-DPF_ARRIVAL_EXACT_COUNT_HEADER="raikou/exact_particle_count.h"', '-isystem', str(HALL), '-isystem', str(HV3 / 'source'),
                '-isystem', str(G), '-I' + str(P)]
    objs = {'arena': cxx + ['-O2', '-g0'] + FP + INC + DEF + STRICT + ['-c', N6 / 'source/src/mfem_pic_arena.cpp'],
            'canonical': cxx + ['-O2', '-g0'] + FP + INC + DEF + STRICT + ['-c', N6 / 'source/src/mfem_pic_canonical.cpp'],
            'kinetic_boundary': hall_cc + INC + ['-c', HALL / 'hall-canonical-driver-20261029/kinetic_boundary.cpp'],
            'motion_provider_v3': hall_cc + ['-I' + str(HV3 / 'source'), '-I' + str(G)] + INC + ['-isystem', str(HALL), '-c', HV3 / 'source/motion_provider.cpp'],
            'test_native_v6_admission': cxx + ['-O1', '-g0'] + FP + INC + DEF + STRICT + test_inc + ['-c', P / 'test_native_v6_admission.cpp']}
    receipt = {'schema': 'hall-v3-native-v6-link-build-v1', 'label': label, 'admission': adm, 'sealed': pins, 'commands': []}
    status = 'passed'
    try:
        jobs = [('compile-' + k, [*c, '-o', b / (k + '.o')]) for k, c in objs.items()]
        with ThreadPoolExecutor(max_workers=2 if len(adm['cpu_affinity']) == 2 else 1) as pool:
            results = list(pool.map(lambda j: execute(j[0], j[1], out, 3600), jobs))
        receipt['commands'] += results
        for r in results:
            if r['returncode'] != 0:
                raise RuntimeError(r['name'] + ' failed:\n' + Path(r['log']).read_text()[-8000:])
        link = ['/usr/bin/mpic++', '-Wl,--gc-sections', *[b / (k + '.o') for k in ('test_native_v6_admission', 'motion_provider_v3', 'kinetic_boundary',
                                                                                 'arena', 'canonical')], EXACT_O, *LIBS, '-o', b / 'test_native_v6_admission']
        r = execute('link', link, out, 1800)
        receipt['commands'].append(r)
        if r['returncode'] != 0:
            raise RuntimeError('link failed:\n' + Path(r['log']).read_text()[-8000:])
        receipt['binary_sha256'] = sha(b / 'test_native_v6_admission')
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


def run(label, steps):
    adm = admission('run')
    br = json.loads((OUT / label / 'build/build-receipt.json').read_text())
    if br['status'] != 'passed' or sha(Path(br['build_dir']) / 'test_native_v6_admission') != br['binary_sha256']:
        fail('build not usable')
    if own() != br['own_sha256']:
        fail('sources changed after the build')
    pins = sealed()
    out = OUT / label / ('run-' + time.strftime('%H%M%S'))
    out.mkdir(parents=True)
    r = execute('admission', [Path(br['build_dir']) / 'test_native_v6_admission', out, *([steps] if steps else [])], out, 6 * 3600)
    rec = {'schema': 'hall-v3-native-v6-link-run-v1', 'label': label, 'admission': adm, 'sealed': pins, 'command': r}
    f = out / 'native-v6-admission.json'
    if f.exists():
        j = json.loads(f.read_text())
        rec['status'] = j.get('status')
        rec['check_count'] = j.get('check_count')
        rec['failure'] = j.get('failure')
        rec['measurements'] = j.get('measurements')
    (out / 'run-receipt.json').write_text(json.dumps(rec, indent=1) + '\n')
    print(json.dumps({'run': str(out), 'status': rec.get('status'), 'checks': rec.get('check_count'), 'failure': rec.get('failure')}))
    return 0 if r['returncode'] == 0 else 1


if __name__ == '__main__':
    a = sys.argv
    if len(a) < 3 or a[1] not in ('build', 'run') or '/' in a[2]:
        raise SystemExit(__doc__)
    if a[1] == 'build':
        raise SystemExit(build(a[2]))
    raise SystemExit(run(a[2], a[3] if len(a) > 3 else None))
