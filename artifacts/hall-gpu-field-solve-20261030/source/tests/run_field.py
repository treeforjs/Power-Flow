"""Qualification runner for the Hall GPU field backend (gpu1) against the sealed
Hall corrective field model (commit 8dddec0). Receipts go to
/home/treefor/sim-workspace/powerflow/outputs/hall-gpu-resident-v1-20261030/field/LABEL.

Phases and launch lines (leases are taken by the launcher, verified here):
  build LABEL
      both-cores.sh /usr/bin/python -I run_field.py build LABEL
  prep LABEL MESH | guard LABEL MESH | report LABEL
      any-core.sh --wait /usr/bin/python -I run_field.py prep LABEL MESH
  gpu LABEL MESH
      flock .../pic-profile/testing-gpu.lock any-core.sh --wait /usr/bin/python -I run_field.py gpu LABEL MESH

Every compile and run step executes inside the pinned frozen launcher
(arch-base-20261028/run-frozen) with a cleaned environment. GPU access inside
that launcher works because it binds /dev and the frozen root carries the same
user-space driver (libcuda 615.78.08) as the host kernel module.
"""
from pathlib import Path
import fcntl, hashlib, json, os, platform, shutil, statistics, subprocess, sys, time

FIELD = Path(__file__).resolve().parents[1]
REPO = Path('/home/treefor/projects/runtimes/hall-coupled-transport-20261029/source')
COMMIT = '8dddec0610982778bf02bdfe697b70a37129c012'
MODULES = ['hall-field-correction-20261030', 'hall-canonical-driver-20261029',
           'hall-finite-foil-em-20261029', 'hall-coupled-transport-20261029']
OUT_ROOT = Path('/home/treefor/sim-workspace/powerflow/outputs/hall-gpu-resident-v1-20261030/field')
LEASES = Path('/home/treefor/sim-workspace/powerflow/outputs/pic-profile')
CORE_LOCK = {6: LEASES / 'testing-core.lock', 7: LEASES / 'testing-core7.lock'}
GPU_LOCK = LEASES / 'testing-gpu.lock'
FROZEN = '/home/treefor/projects/runtimes/arch-base-20261028/run-frozen'
DEP = Path('/home/treefor/projects/runtimes/mfem-4.10-pic-20261029')
LIBS = [DEP / 'install-mpi/lib/libmfem.a', DEP / 'tpl/hypre/lib/libHYPRE.a',
        DEP / 'tpl/metis/lib/libmetis.a', DEP / 'tpl/GKlib/lib/libGKlib.a']
CUDA_INCLUDE = Path('/home/treefor/projects/runtimes/mfem-cuda-benchmark-20261029/cuda-libs/nvidia/cuda_runtime/include')
NVRTC = Path('/home/treefor/projects/runtimes/pic-cuda-gather-investigation-20261029/tools/nvrtc/nvidia/cuda_nvrtc')
FROZEN_ROOT = Path('/home/treefor/projects/runtimes/arch-base-20261028/root')
CLEAN_ENV = {'PATH': '/usr/bin:/bin', 'LANG': 'C.UTF-8', 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1',
             'MKL_NUM_THREADS': '1', 'PYTHONNOUSERSITE': '1'}
# name: (transverse, axial cells, air cells, width m, CPU repeats, GPU repeats, address-space GiB, role)
MESHES = {
    'm29a2': (29, 2, 1, .0058, 20, 20, 4, 'smallest admissible box (5.8 mm); timing point for the crossover only'),
    'm45a2': (45, 2, 1, .009, 20, 20, 4, 'Hall box, 2 axial foil cells; timing point for the crossover'),
    'm45': (45, 4, 1, .009, 20, 20, 4, 'Hall default MeshSpec{45,4,1,.009}'),
    'm90': (90, 4, 1, .009, 10, 20, 8, 'refined transverse x2'),
    'm135a6': (135, 6, 1, .009, 5, 20, 12, 'refined, ~5e5 edges'),
    'm270a6': (270, 6, 1, .009, 3, 10, 16, 'refined, ~2e6 edges'),
}
GPU_PROCESSES = 3   # fresh processes with default scheduling (bit-identity across processes)
GPU_BATCH, GPU_GRAPH = 4, 1


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def git_blob(p):
    data = Path(p).read_bytes()
    return hashlib.sha1(b'blob %d\0' % len(data) + data).hexdigest()


def fail(message):
    raise RuntimeError(message)


def lease_held(path):
    """True when another open file description holds the lock (the launcher's)."""
    with open(path, 'a') as probe:
        try:
            fcntl.flock(probe, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return True
        fcntl.flock(probe, fcntl.LOCK_UN)
        return False


def admission(phase):
    aff = os.sched_getaffinity(0)
    if os.getpriority(os.PRIO_PROCESS, 0) < 15:
        fail('nice 15 required')
    if any(os.getenv(k) != '1' for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']):
        fail('one thread per library required')
    leases = {}
    if phase == 'build':
        if aff != {6, 7}:
            fail('build requires both Power-Flow cores {6,7}')
        for c in (6, 7):
            if not lease_held(CORE_LOCK[c]):
                fail('core %d lease must be held by the launcher' % c)
            leases['core%d' % c] = str(CORE_LOCK[c])
    else:
        if aff not in ({6}, {7}):
            fail('one Power-Flow core (6 or 7) required')
        core = next(iter(aff))
        if os.getenv('PF_CORE') not in (None, str(core)):
            fail('PF_CORE disagrees with affinity')
        if not lease_held(CORE_LOCK[core]):
            fail('core %d lease must be held by the launcher' % core)
        leases['core%d' % core] = str(CORE_LOCK[core])
        if phase == 'gpu':
            if not lease_held(GPU_LOCK):
                fail('GPU lease must be held by the launcher')
            leases['gpu'] = str(GPU_LOCK)
    if Path('/proc/sys/kernel/osrelease').read_text().strip() != '7.2.11-arch1-1':
        fail('kernel release differs from the frozen runtime')
    head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True).strip()
    dirty = subprocess.check_output(['git', 'status', '--porcelain'], cwd=REPO, text=True)
    if head != COMMIT or dirty:
        fail('Hall repository is not the clean corrective commit')
    pins_path = REPO / 'development/hall-field-correction-20261030/runtime-pins.json'
    pins = json.loads(pins_path.read_text())
    frozen_root = Path(pins['frozen_root'])
    base = frozen_root.parent
    if sha(base / 'root-files.json') != pins['frozen_root_manifest_sha256'] or sha(base / 'run-frozen') != pins['launcher_sha256']:
        fail('frozen runtime manifest or launcher changed')
    if phase == 'build':
        for name, data in pins['selected_files'].items():
            p = frozen_root / name
            if 'sha256' in data:
                if sha(p) != data['sha256']:
                    fail('frozen runtime file changed: ' + name)
            elif 'symlink' in data:
                if not p.is_symlink() or os.readlink(p) != data['symlink']:
                    fail('frozen runtime link changed: ' + name)
        for name, expected in pins['static_libraries'].items():
            if sha(name) != expected:
                fail('pinned static library changed: ' + name)
    return {'phase': phase, 'cpu_affinity': sorted(aff), 'nice': os.getpriority(os.PRIO_PROCESS, 0),
            'PF_CORE': os.getenv('PF_CORE'), 'leases_verified_held_by_launcher': leases, 'kernel_release': '7.2.11-arch1-1',
            'hall_commit': head, 'hall_git_clean': True, 'runtime_pins': str(pins_path), 'runtime_pins_sha256': sha(pins_path),
            'frozen_launcher': FROZEN, 'frozen_launcher_sha256': sha(FROZEN),
            'frozen_runtime_manifest_checked': phase == 'build'}


def run(name, argv, log, timeout, env, receipts, cwd):
    start = time.monotonic()
    with open(log, 'w') as f:
        r = subprocess.run(argv, cwd=cwd, env=env, stdout=f, stderr=subprocess.STDOUT, timeout=timeout)
    rec = {'name': name, 'argv': argv, 'returncode': r.returncode, 'elapsed_s': time.monotonic() - start,
           'log': str(log), 'log_sha256': sha(log), 'environment': env, 'cwd': str(cwd)}
    receipts.append(rec)
    return rec


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2) + '\n')


def device_identity():
    info = {}
    try:
        info['proc_driver_version'] = Path('/proc/driver/nvidia/version').read_text()
    except OSError as e:
        info['proc_driver_version_error'] = str(e)
    try:
        info['nvidia_smi'] = subprocess.check_output(
            ['nvidia-smi', '--query-gpu=index,name,uuid,driver_version,compute_cap,pci.bus_id,memory.total', '--format=csv'],
            text=True, timeout=30)
    except Exception as e:  # noqa: BLE001
        info['nvidia_smi_error'] = str(e)
    return info


def phase_build(label, out, receipt):
    out.mkdir(parents=True, exist_ok=False)
    source = out / 'source'
    tree = {}
    for line in subprocess.check_output(['git', 'ls-tree', '-r', COMMIT, '--', *['development/' + m for m in MODULES]],
                                        cwd=REPO, text=True).splitlines():
        meta, path = line.split('\t', 1)
        tree[path] = meta.split()[2]
    for name in MODULES:
        shutil.copytree(REPO / 'development' / name, source / 'development' / name, symlinks=True)
    files = sorted(p for p in source.rglob('*') if p.is_file() and not p.is_symlink())
    mismatch = [str(p) for p in files if tree.get(str(p.relative_to(source))) != git_blob(p)]
    if mismatch or len(files) != len(tree):
        fail('snapshot differs from the commit tree')
    fsrc = out / 'field-source'
    shutil.copytree(FIELD / 'src', fsrc / 'src')
    shutil.copytree(FIELD / 'tests', fsrc / 'tests', ignore=shutil.ignore_patterns('__pycache__'))
    if (FIELD / 'README.txt').exists():
        shutil.copy2(FIELD / 'README.txt', fsrc / 'README.txt')
    native = source / 'development/hall-canonical-driver-20261029/bindings/native-v6.2'
    common = ['-std=c++20', '-O2', '-DNDEBUG', '-ffp-contract=off', '-fno-fast-math', '-Wall', '-Wextra']
    build = out / 'build'
    build.mkdir()
    cpu = ['/usr/bin/c++', *common, '-isystem', str(source), '-I' + str(native / 'source/include'), '-isystem', str(native / 'inputs'),
           '-I' + str(native / 'vendor')]
    for d in ['install-mpi/include', 'tpl/hypre/include', 'tpl/metis/include']:
        cpu += ['-isystem', str(DEP / d)]
    cpu += [str(fsrc / 'src/cpu_main.cpp'), *map(str, LIBS), '/usr/lib/libmpi.so', '-lcrypto', '-o', str(build / 'hall_field_cpu')]
    gpu = ['/usr/bin/c++', *common, '-isystem', str(CUDA_INCLUDE), '-isystem', str(NVRTC / 'include'), '-isystem', str(native / 'vendor'),
           '-isystem', str(native / 'inputs'), str(fsrc / 'src/gpu_main.cpp'), '/usr/lib/libcuda.so.1', str(NVRTC / 'lib/libnvrtc.so.12'),
           '-Wl,-rpath,' + str(NVRTC / 'lib'), '-lcrypto', '-o', str(build / 'hall_field_gpu')]
    env = dict(CLEAN_ENV, CUDA_VISIBLE_DEVICES='')
    receipts = []
    procs = []
    for name, cmd in (('build-cpu', cpu), ('build-gpu', gpu)):
        argv = [FROZEN, 'prlimit', '--as=8589934592', '--cpu=900', '--', *cmd]
        log = out / (name + '.log')
        procs.append((name, argv, log, time.monotonic(), subprocess.Popen(argv, cwd=out, env=env, stdout=open(log, 'w'), stderr=subprocess.STDOUT)))
    for name, argv, log, start, p in procs:
        rc = p.wait(timeout=1200)
        receipts.append({'name': name, 'argv': argv, 'returncode': rc, 'elapsed_s': time.monotonic() - start, 'log': str(log),
                         'log_sha256': sha(log), 'environment': env, 'cwd': str(out)})
    receipt['commands'] = receipts
    if any(r['returncode'] for r in receipts):
        fail('build failed: ' + ' / '.join((out / (r['name'] + '.log')).read_text()[-3000:] for r in receipts if r['returncode']))
    for r in receipts:
        if Path(r['log']).stat().st_size:
            fail('build emitted warnings: ' + r['log'])
    inputs = {str(p): sha(p) for p in sorted(source.rglob('*')) if p.is_file() and not p.is_symlink()}
    inputs.update({str(p): sha(p) for p in sorted(fsrc.rglob('*')) if p.is_file()})
    for d in ['install-mpi/include', 'tpl/hypre/include', 'tpl/metis/include']:
        inputs.update({str(p): sha(p) for p in sorted((DEP / d).rglob('*')) if p.is_file()})
    inputs.update({str(p): sha(p) for p in LIBS})
    inputs.update({str(p): sha(p) for p in sorted(CUDA_INCLUDE.rglob('*')) if p.is_file()})
    inputs.update({str(p): sha(p) for p in sorted(NVRTC.rglob('*')) if p.is_file()})
    for lib in ['usr/lib/libcuda.so.615.78.08', 'usr/lib/libmpi.so', 'usr/bin/c++']:
        p = FROZEN_ROOT / lib
        inputs[str(p.resolve())] = sha(p)
    write_json(out / 'input-sha256.json', inputs)
    receipt['input_sha256_manifest'] = str(out / 'input-sha256.json')
    receipt['input_sha256_manifest_sha256'] = sha(out / 'input-sha256.json')
    receipt['input_files'] = len(inputs)
    receipt['binaries'] = {n: sha(build / n) for n in ['hall_field_cpu', 'hall_field_gpu']}
    receipt['trace'] = str(source / 'development/hall-finite-foil-em-20261029/data/shot_17133_current.csv')
    receipt['trace_sha256'] = sha(receipt['trace'])


def mesh_args(mesh):
    t, a, air, w = MESHES[mesh][:4]
    return [str(t), str(a), str(air), repr(w)]


def phase_prep_or_guard(phase, out, mesh, receipt):
    build = json.loads((out / 'build.json').read_text())
    binary = out / 'build/hall_field_cpu'
    if sha(binary) != build['binaries']['hall_field_cpu']:
        fail('CPU binary changed since build')
    mdir = out / 'meshes' / mesh
    if phase == 'prep':
        mdir.mkdir(parents=True, exist_ok=False)
        extra = str(MESHES[mesh][4])
    else:
        if (mdir / 'guard.json').exists():
            fail('guard already recorded for this mesh')
        gpu = json.loads((mdir / 'gpu-phase.json').read_text())
        if gpu.get('status') != 'passed':
            fail('GPU phase did not pass')
        extra = str(mdir / 'gpu-run-1')
    env = dict(CLEAN_ENV, CUDA_VISIBLE_DEVICES='')
    argv = [FROZEN, 'prlimit', '--as=%d' % (MESHES[mesh][6] << 30), '--cpu=3000', '--',
            str(binary), phase, build['trace'], str(mdir), *mesh_args(mesh), extra]
    receipts = []
    rec = run(phase, argv, mdir / (phase + '.log'), 3600, env, receipts, mdir)
    receipt['commands'] = receipts
    receipt['binary_sha256'] = build['binaries']['hall_field_cpu']
    receipt['trace_sha256'] = sha(build['trace'])
    if rec['returncode'] != 0:
        fail('%s failed: %s' % (phase, Path(rec['log']).read_text()[-3000:]))
    result = json.loads((mdir / (phase + '.json')).read_text())
    receipt['result'] = str(mdir / (phase + '.json'))
    receipt['result_sha256'] = sha(mdir / (phase + '.json'))
    if phase == 'prep':
        receipt['fixtures'] = {n: {'file': f['fixture_file'], 'sha256': f['fixture_sha256'], 'verified_sha256': sha(f['fixture_file'])}
                               for n, f in result['fixtures'].items()}
    else:
        if gpu['inputs_sha256'] != {n: sha(f['fixture_file']) for n, f in json.loads((mdir / 'prep.json').read_text())['fixtures'].items()}:
            fail('fixtures changed since the GPU phase')


def phase_gpu(out, mesh, receipt):
    build = json.loads((out / 'build.json').read_text())
    binary = out / 'build/hall_field_gpu'
    if sha(binary) != build['binaries']['hall_field_gpu']:
        fail('GPU binary changed since build')
    mdir = out / 'meshes' / mesh
    prep = json.loads((mdir / 'prep.json').read_text())
    inputs = {n: sha(f['fixture_file']) for n, f in prep['fixtures'].items()}
    for n, f in prep['fixtures'].items():
        if inputs[n] != f['fixture_sha256']:
            fail('fixture changed since prep: ' + n)
    receipt['binary_sha256'] = build['binaries']['hall_field_gpu']
    receipt['inputs_sha256'] = inputs
    receipt['device_identity_before'] = device_identity()
    env = dict(CLEAN_ENV)
    receipts = []
    gpu_repeats = MESHES[mesh][5]
    runs = [('gpu-run-%d' % (k + 1), 0, 1 if k == 0 else 0) for k in range(GPU_PROCESSES)] + [('gpu-run-blocking', 1, 0)]
    for name, blocking, sweep in runs:
        d = mdir / name
        if d.exists():
            fail('output exists: ' + str(d))
        argv = [FROZEN, 'prlimit', '--cpu=900', '--', str(binary), str(mdir), str(d), '3', str(gpu_repeats),
                str(GPU_BATCH), str(GPU_GRAPH), str(blocking), str(sweep)]
        d.mkdir(parents=True)
        rec = run(name, argv, mdir / (name + '.log'), 1800, env, receipts, mdir)
        receipt['commands'] = receipts
        if rec['returncode'] != 0:
            fail('%s failed (rc %d): %s' % (name, rec['returncode'], Path(rec['log']).read_text()[-3000:]))
    receipt['device_identity_after'] = device_identity()
    # Cross-process bit identity of every solution.
    digests = {}
    for name, _, _ in runs:
        g = json.loads((mdir / name / 'gpu.json').read_text())
        if g['cleanup']['issues'] or g['cleanup']['quarantined'] or not g['cleanup']['primary_released']:
            fail('cleanup issue in ' + name)
        for fx, f in g['fixtures'].items():
            for mode, m in f['modes'].items():
                digests.setdefault('%s/%s' % (fx, mode), {})[name] = {
                    'x_sha256': m['x_sha256'], 'repeat_identical': m['repeat_solves_bitwise_identical'],
                    'emulation_bitwise': m['x_vs_host_emulation']['bitwise_equal'], 'x_file_sha256': sha(m['x_file'])}
    identity = {k: {'processes': v, 'identical_across_processes': len({e['x_sha256'] for e in v.values()}) == 1,
                    'all_repeats_identical': all(e['repeat_identical'] for e in v.values()),
                    'equal_to_host_emulation': all(e['emulation_bitwise'] for e in v.values())} for k, v in digests.items()}
    receipt['bit_identity'] = identity
    if not all(v['identical_across_processes'] and v['all_repeats_identical'] for v in identity.values()):
        fail('GPU solutions are not bit-identical across solves/processes')


def phase_report(out, receipt):
    rows, fixtures = [], {}
    for mesh in MESHES:
        mdir = out / 'meshes' / mesh
        if not (mdir / 'guard.json').exists():
            continue
        prep = json.loads((mdir / 'prep.json').read_text())
        guard = json.loads((mdir / 'guard.json').read_text())
        gphase = json.loads((mdir / 'gpu-phase.json').read_text())
        g1 = json.loads((mdir / 'gpu-run-1' / 'gpu.json').read_text())
        gb = json.loads((mdir / 'gpu-run-blocking' / 'gpu.json').read_text())
        for fx, p in prep['fixtures'].items():
            gr = guard['fixtures'][fx]
            gf = g1['fixtures'][fx]
            entry = {'mesh': mesh, 'role': MESHES[mesh][7], 'edges': prep['mesh']['edges'], 'nnz_assembled': p['operator']['nnz_assembled'],
                     'fixture': fx, 'dt_s': p['dt_s'],
                     'csr_vs_matrix_free': [{'vector': c['vector'], 'assembled': c['assembled_csr_vs_matrix_free'],
                                             'factored_bitwise': c['factored_vs_matrix_free']['bitwise_equal']} for c in p['spmv_checks']],
                     'cpu_iterations': p['cpu_sealed_pcg']['record']['iterations_per_correction'],
                     'cpu_kappa': p['cpu_sealed_pcg']['record']['lanczos']['kappa_estimate'],
                     'cpu_sealed_ms': p['cpu_sealed_pcg']['timing']['median_ms'],
                     'cpu_plain_assembled_ms': p['cpu_plain_pcg_assembled']['timing']['median_ms'],
                     'cpu_plain_factored_ms': p['cpu_plain_pcg_factored']['timing']['median_ms'],
                     'sigma_rows': p['operator']['sigma_dependent_rows'], 'negative_controls': gr['negative_controls'], 'gpu': {}}
            for mode, m in gf['modes'].items():
                v = gr['variants']['gpu_' + mode]
                gd = v['guards']
                ritz = m['lanczos']
                ratio = v['independent_residual']['residual'] / v['independent_residual']['rhs_norm']
                entry['gpu'][mode] = {
                    'iterations': m['iterations_per_correction'], 'kappa': ritz['kappa_estimate'],
                    'forward_error_estimate': ritz['kappa_estimate'] * ratio,
                    'residual': v['independent_residual'],
                    'gauss': gd['absolute_gauss'], 'divB': gd['absolute_divB'], 'gauss_increment': gd['gauss_increment'],
                    'energy': gd['field_energy'], 'field_heat': gd['field_heat_energy'], 'joule': gd['joule_transfer'],
                    'x_vs_cpu_sealed': m['x_vs_cpu_sealed'],
                    'bit_identity': gphase['bit_identity']['%s/%s' % (fx, mode)],
                    'resident_ms': m['timing_resident']['median_ms'],
                    'transfers_ms': m['timing_with_rhs_solution_transfers']['median_ms'],
                    'sigma_refresh_transfers_ms': m['timing_with_sigma_refresh_and_transfers']['median_ms'],
                    'full_upload_transfers_ms': m['timing_with_full_operator_upload_and_transfers']['median_ms'],
                    'blocking_sync_resident_ms': gb['fixtures'][fx]['modes'][mode]['timing_resident']['median_ms']}
            entry['cpu_guards'] = {k: gr['variants']['cpu_sealed']['guards'][k] for k in ['absolute_gauss', 'absolute_divB', 'gauss_increment', 'field_energy', 'field_heat_energy']}
            entry['cpu_candidate_reproduced'] = gr['variants']['cpu_sealed']['guards']['after_sha256'] == p['cpu_path_guards']['after_sha256']
            fixtures['%s/%s' % (mesh, fx)] = entry
            rows.append(entry)
    receipt['fixtures'] = fixtures
    # Crossover (per fixture type): smallest edge count at which GPU beats CPU.
    cross = {}
    for fx in ['dc-10fs', 'dc-1fs', 'heated-10fs', 'heated-1fs']:
        pts = sorted((r for r in rows if r['fixture'] == fx), key=lambda r: r['edges'])
        cross[fx] = [{'edges': r['edges'], 'mesh': r['mesh'], 'cpu_sealed_ms': r['cpu_sealed_ms'], 'cpu_plain_assembled_ms': r['cpu_plain_assembled_ms'],
                      'gpu_factored_resident_ms': r['gpu']['factored']['resident_ms'], 'gpu_factored_sigma_refresh_transfers_ms': r['gpu']['factored']['sigma_refresh_transfers_ms'],
                      'gpu_assembled_resident_ms': r['gpu']['assembled']['resident_ms']} for r in pts]
    receipt['timing_by_fixture'] = cross


def main():
    if len(sys.argv) < 3:
        fail('usage: run_field.py build|prep|gpu|guard|report LABEL [MESH]')
    phase, label = sys.argv[1], sys.argv[2]
    if phase not in ('build', 'prep', 'gpu', 'guard', 'report') or not label or '/' in label or label in ('.', '..'):
        fail('bad phase or label')
    mesh = None
    if phase in ('prep', 'gpu', 'guard'):
        if len(sys.argv) != 4 or sys.argv[3] not in MESHES:
            fail('mesh required: ' + ', '.join(MESHES))
        mesh = sys.argv[3]
    elif len(sys.argv) != 3:
        fail('unexpected arguments')
    receipt = {'schema': 'hall-gpu-field-qualification-v1', 'label': label, 'phase': phase, 'mesh': mesh,
               'runner_argv': sys.argv, 'runner': __file__, 'runner_sha256': sha(__file__), 'python': sys.executable,
               'python_version': platform.python_version(), 'started_unix': time.time()}
    receipt['admission'] = admission(phase)
    out = OUT_ROOT / label
    if phase != 'build' and not (out / 'build.json').exists():
        fail('build receipt missing')
    path = out / ('build.json' if phase == 'build' else ('report.json' if phase == 'report' else 'meshes/%s/%s-phase.json' % (mesh, phase)))
    try:
        if phase == 'build':
            phase_build(label, out, receipt)
        elif phase in ('prep', 'guard'):
            phase_prep_or_guard(phase, out, mesh, receipt)
        elif phase == 'gpu':
            phase_gpu(out, mesh, receipt)
        else:
            phase_report(out, receipt)
        receipt['status'] = 'passed'
    except Exception as e:  # noqa: BLE001
        receipt['status'] = 'failed'
        receipt['error'] = str(e)
        if out.exists():
            path.parent.mkdir(parents=True, exist_ok=True)
            receipt['finished_unix'] = time.time()
            write_json(path, receipt)
        raise
    receipt['finished_unix'] = time.time()
    path.parent.mkdir(parents=True, exist_ok=True)
    write_json(path, receipt)
    print(phase, mesh or '', 'passed', path)


if __name__ == '__main__':
    main()
