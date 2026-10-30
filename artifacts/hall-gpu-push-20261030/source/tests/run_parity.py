"""Build and run the GPU push parity harness; write receipts.

Stages (each launched under the leases it needs; the runner verifies them):
  build   both-cores.sh python -I run_parity.py build LABEL
  cpu     any-core.sh --wait python -I run_parity.py cpu LABEL [fixture ...]
  gpu     flock <gpu lock> any-core.sh --wait python -I run_parity.py gpu LABEL [fixture ...]
  report  any-core.sh --wait python -I run_parity.py report LABEL
Every compile and run step executes inside the frozen launcher
(arch-base-20261028/run-frozen). The frozen root carries libcuda.so.615.78.08,
byte-identical to the host driver library (sealed in root-files.json); /dev is
bound, so GPU steps run inside the sandbox with the same NVRTC 12.8 library.
Outputs: OUT/LABEL/{build,cpu,gpu}/ plus receipt, parity and timing JSON.
"""
from pathlib import Path
import fcntl
import hashlib
import json
import os
import platform
import statistics
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import compare  # noqa: E402

OUT_ROOT = Path('/home/treefor/sim-workspace/powerflow/outputs/hall-gpu-resident-v1-20261030/push')
LEASES = Path('/home/treefor/sim-workspace/powerflow/outputs/pic-profile')
CORE_LOCK = {6: LEASES / 'testing-core.lock', 7: LEASES / 'testing-core7.lock'}
GPU_LOCK = LEASES / 'testing-gpu.lock'
HALL = Path('/home/treefor/projects/runtimes/hall-coupled-transport-20261029/source')
HALL_COMMIT = '8dddec0610982778bf02bdfe697b70a37129c012'
OLD = HALL / 'development/hall-canonical-driver-20261029'
NATIVE = OLD / 'bindings/native-v6.2'
FIELD = HALL / 'development/hall-field-correction-20261030'
TRACE = HALL / 'development/hall-finite-foil-em-20261029/data/shot_17133_current.csv'
DEP = Path('/home/treefor/projects/runtimes/mfem-4.10-pic-20261029')
LIBS = [DEP / 'install-mpi/lib/libmfem.a', DEP / 'tpl/hypre/lib/libHYPRE.a', DEP / 'tpl/metis/lib/libmetis.a', DEP / 'tpl/GKlib/lib/libGKlib.a']
NVRTC_DIR = Path('/home/treefor/projects/runtimes/pic-cuda-gather-investigation-20261029/tools/nvrtc/nvidia/cuda_nvrtc')
NVRTC = NVRTC_DIR / 'lib/libnvrtc.so.12'
CUDA_INC = Path('/home/treefor/projects/runtimes/mfem-cuda-benchmark-20261029/cuda-libs/nvidia/cuda_runtime/include')
FROZEN = Path('/home/treefor/projects/runtimes/arch-base-20261028/run-frozen')
GPU_UUID = 'GPU-a7b9b458-167f-cd57-fec9-3dea0995b6e9'
GPU_REPEATS_IN_PROCESS = 5
FRESH_PROCESSES = 3
FIXTURES = ['port_pair_w1e+09_seed', 'port_pair_w1e+09_final', 'port_pair_w1e+11_seed', 'port_pair_w1e+11_final',
            'port_pair_w1e+12_seed', 'port_pair_w1e+12_final',
            'syn_record_1e3', 'syn_record_1e4', 'syn_record_1e5', 'syn_record_1e6',
            'syn_boundary_1e3', 'syn_boundary_1e4', 'syn_boundary_1e5', 'syn_boundary_1e6', 'syn_metal_1e3']
ENV = {'PATH': '/usr/bin:/bin', 'LANG': 'C.UTF-8', 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1',
       'PYTHONNOUSERSITE': '1'}


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def held(lock):
    """True if another holder (the launching lease wrapper) holds the lock."""
    with open(lock, 'a') as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return True
        fcntl.flock(handle, fcntl.LOCK_UN)
        return False


def admission(stage):
    cpus = os.sched_getaffinity(0)
    nice = os.getpriority(os.PRIO_PROCESS, 0)
    if any(c not in (6, 7) for c in cpus):
        raise RuntimeError('only cores 6 and 7 are permitted; affinity is %s' % sorted(cpus))
    if stage == 'build':
        if cpus != {6, 7} or not held(CORE_LOCK[6]) or not held(CORE_LOCK[7]):
            raise RuntimeError('build requires both-cores.sh (cores 6,7 with both leases)')
    else:
        if len(cpus) != 1:
            raise RuntimeError('single-core stage requires any-core.sh')
        core = next(iter(cpus))
        if not held(CORE_LOCK[core]):
            raise RuntimeError('core %d lease is not held' % core)
        if stage == 'gpu' and not held(GPU_LOCK):
            raise RuntimeError('GPU lease is not held')
    if nice < 15:
        raise RuntimeError('nice 15 required')
    if any(os.getenv(k) != '1' for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']):
        raise RuntimeError('one thread per library required')
    if Path('/proc/sys/kernel/osrelease').read_text().strip() != '7.2.11-arch1-1':
        raise RuntimeError('kernel release differs from the frozen runtime')
    head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=HALL, text=True).strip()
    dirty = subprocess.check_output(['git', 'status', '--porcelain'], cwd=HALL, text=True)
    if head != HALL_COMMIT or dirty:
        raise RuntimeError('Hall repository is not the clean commit ' + HALL_COMMIT)
    pins = json.loads((FIELD / 'runtime-pins.json').read_text())
    base = Path(pins['frozen_root']).parent
    if sha(base / 'root-files.json') != pins['frozen_root_manifest_sha256'] or sha(base / 'run-frozen') != pins['launcher_sha256']:
        raise RuntimeError('frozen runtime manifest or launcher changed')
    return {'stage': stage, 'cpu_affinity': sorted(cpus), 'nice': nice, 'core_leases_held': {str(c): held(l) for c, l in CORE_LOCK.items()},
            'gpu_lease_held': held(GPU_LOCK), 'PF_CORE': os.getenv('PF_CORE'), 'hall_commit': head, 'hall_clean': True,
            'kernel_release': '7.2.11-arch1-1', 'frozen_launcher_sha256': sha(FROZEN), 'runtime_pins_sha256': sha(FIELD / 'runtime-pins.json'),
            'python': sys.executable, 'python_version': platform.python_version(), 'runner_sha256': sha(__file__), 'argv': sys.argv,
            'time_unix': time.time()}


def run_step(name, cmd, log_dir, env_extra=None, as_limit=None, timeout=3600, expected=0):
    pre = [str(FROZEN), 'prlimit', '--cpu=%d' % timeout]
    if as_limit:
        pre.append('--as=%d' % as_limit)
    argv = pre + ['--', *map(str, cmd)]
    env = dict(ENV, **(env_extra or {}))
    log = log_dir / (name + '.log')
    start = time.monotonic()
    with log.open('w') as f:
        try:
            rc = subprocess.run(argv, cwd=str(log_dir), env=env, stdout=f, stderr=subprocess.STDOUT, timeout=timeout).returncode
        except subprocess.TimeoutExpired:
            rc = 'timeout'
    entry = {'name': name, 'argv': argv, 'environment': env, 'returncode': rc, 'expected_returncode': expected,
             'elapsed_s': time.monotonic() - start, 'log': str(log), 'log_sha256': sha(log)}
    print(name, rc, round(entry['elapsed_s'], 2), flush=True)
    return entry


def stage_build(label):
    out = OUT_ROOT / label / 'build'
    out.mkdir(parents=True, exist_ok=False)
    receipt = {'admission': admission('build')}
    hall_flags = ['-std=c++20', '-O2', '-DNDEBUG', '-ffp-contract=off', '-fno-fast-math', '-ffunction-sections', '-fdata-sections',
                  '-I' + str(NATIVE / 'source/include'), '-isystem', str(NATIVE / 'inputs'), '-I' + str(NATIVE / 'vendor')]
    for d in ['install-mpi/include', 'tpl/hypre/include', 'tpl/metis/include']:
        hall_flags += ['-isystem', str(DEP / d)]
    cpu_cmd = ['/usr/bin/c++', *hall_flags, '-I' + str(FIELD), '-I' + str(ROOT), '-MD', '-MF', str(out / 'cpu_reference.d'),
               str(ROOT / 'tests/cpu_reference.cpp'), str(OLD / 'native_arena_motion.cpp'), str(OLD / 'kinetic_boundary.cpp'),
               str(NATIVE / 'source/src/mfem_gap_pic.cpp'), *map(str, LIBS), '/usr/lib/libmpi.so', '-lcrypto', '-Wl,--gc-sections',
               '-o', str(out / 'cpu_reference')]
    gpu_cmd = ['/usr/bin/c++', '-std=c++20', '-O2', '-DNDEBUG', '-ffp-contract=off', '-fno-fast-math', '-Wall', '-Wextra',
               '-Wno-c++26-extensions', '-I' + str(ROOT), '-isystem', str(CUDA_INC), '-isystem', str(NVRTC_DIR / 'include'),
               '-MD', '-MF', str(out / 'gpu_parity.d'), str(ROOT / 'tests/gpu_parity.cpp'), str(ROOT / 'common/hall_gpu_common.cpp'),
               str(ROOT / 'common/hall_gpu_primitives.cpp'), str(ROOT / 'push/hall_gpu_push.cpp'), '-ldl', '-lcrypto', '-o', str(out / 'gpu_parity')]
    loop_cmd = ['/usr/bin/c++', *hall_flags, '-Wno-c++26-extensions', '-I' + str(FIELD), '-I' + str(ROOT), '-isystem', str(CUDA_INC),
                '-isystem', str(NVRTC_DIR / 'include'), '-MD', '-MF', str(out / 'gpu_closed_loop.d'), str(ROOT / 'tests/gpu_closed_loop.cpp'),
                str(ROOT / 'common/hall_gpu_common.cpp'), str(ROOT / 'common/hall_gpu_primitives.cpp'), str(ROOT / 'push/hall_gpu_push.cpp'),
                str(OLD / 'native_arena_motion.cpp'), str(OLD / 'kinetic_boundary.cpp'), str(NATIVE / 'source/src/mfem_gap_pic.cpp'),
                *map(str, LIBS), '/usr/lib/libmpi.so', '-lcrypto', '-ldl', '-Wl,--gc-sections', '-o', str(out / 'gpu_closed_loop')]
    steps = []
    for name, cmd in [('build-cpu-reference', cpu_cmd), ('build-gpu-parity', gpu_cmd), ('build-gpu-closed-loop', loop_cmd)]:
        steps.append(run_step(name, cmd, out, {'CUDA_VISIBLE_DEVICES': ''}, as_limit=4 << 30, timeout=900))
        if steps[-1]['returncode'] != 0:
            receipt['steps'] = steps
            (out / 'build-receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
            raise RuntimeError('build failed: ' + name + '\n' + (out / (name + '.log')).read_text()[-5000:])
    inputs = set()
    for dep in ['cpu_reference.d', 'gpu_parity.d', 'gpu_closed_loop.d']:
        text = (out / dep).read_text().replace('\\\n', ' ')
        for token in text.split():
            if token.endswith(':') or not Path(token).is_file():
                continue
            inputs.add(str(Path(token).resolve()))
    inputs |= {str(p) for p in LIBS} | {str(NVRTC), '/usr/lib/libmpi.so', str(TRACE), str(FIELD / 'runtime-pins.json')}
    inputs |= {str(p) for p in (ROOT / 'common').rglob('*') if p.is_file()} | {str(p) for p in (ROOT / 'push').rglob('*') if p.is_file()}
    inputs |= {str(p) for p in (ROOT / 'tests').glob('*') if p.is_file() and p.suffix in ('.py', '.cpp', '.hpp')}
    manifest = {p: sha(p) for p in sorted(inputs)}
    (out / 'input-sha256.json').write_text(json.dumps(manifest, indent=2) + '\n')
    for s in steps:
        b = Path(s['argv'][s['argv'].index('-o') + 1])
        s['binary'] = str(b)
        s['binary_sha256'] = sha(b)
    receipt.update({'steps': steps, 'input_manifest': str(out / 'input-sha256.json'), 'input_manifest_sha256': sha(out / 'input-sha256.json'),
                    'input_files': len(manifest)})
    (out / 'build-receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')


def verify_inputs(label):
    manifest = json.loads((OUT_ROOT / label / 'build/input-sha256.json').read_text())
    changed = [p for p, h in manifest.items() if sha(p) != h]
    if changed:
        raise RuntimeError('inputs changed since the build: ' + ', '.join(changed[:5]))
    return len(manifest)


def stage_cpu(label, fixtures):
    base = OUT_ROOT / label
    out = base / 'cpu'
    out.mkdir(parents=True, exist_ok=True)
    receipt = {'admission': admission('cpu'), 'inputs_verified_before': verify_inputs(label)}
    exe = base / 'build/cpu_reference'
    steps = [run_step('cpu-selftest', [exe, 'selftest', out], out, {'CUDA_VISIBLE_DEVICES': ''}, as_limit=8 << 30, timeout=1800)]
    steps.append(run_step('cpu-reference', [exe, 'run', TRACE, out, *fixtures], out, {'CUDA_VISIBLE_DEVICES': ''}, as_limit=24 << 30, timeout=7200))
    receipt['steps'] = steps
    receipt['inputs_verified_after'] = verify_inputs(label)
    (out / ('cpu-receipt-%d.json' % int(time.time()))).write_text(json.dumps(receipt, indent=2) + '\n')
    if any(s['returncode'] != 0 for s in steps):
        raise RuntimeError('CPU stage failed; see logs in ' + str(out))


def stage_gpu(label, fixtures, closed_loop=True):
    base = OUT_ROOT / label
    out = base / 'gpu'
    out.mkdir(parents=True, exist_ok=True)
    receipt = {'admission': admission('gpu'), 'inputs_verified_before': verify_inputs(label), 'nvrtc_library': str(NVRTC),
               'nvrtc_library_sha256': sha(NVRTC), 'device_uuid_required': GPU_UUID, 'fresh_processes': FRESH_PROCESSES,
               'repeats_in_process': GPU_REPEATS_IN_PROCESS}
    exe = base / 'build/gpu_parity'
    steps = []
    for fx in fixtures:
        fdir = base / 'cpu' / fx
        if not (fdir / 'input.bin').is_file():
            raise RuntimeError('fixture missing: ' + fx)
        (out / fx).mkdir(exist_ok=True)
        for k in range(1, FRESH_PROCESSES + 1):
            stem = out / fx / ('run%d' % k)
            steps.append(run_step('gpu-%s-run%d' % (fx, k), [exe, fdir, stem, NVRTC, GPU_UUID, GPU_REPEATS_IN_PROCESS], out / fx,
                                  {'CUDA_VISIBLE_DEVICES': GPU_UUID}, timeout=3600))
            steps[-1]['fixture_input_sha256'] = sha(fdir / 'input.bin')
            if (stem.with_suffix('.bin')).is_file():
                steps[-1]['output_sha256'] = sha(stem.with_suffix('.bin'))
    if closed_loop:
        steps.append(run_step('gpu-closed-loop', [base / 'build/gpu_closed_loop', TRACE, NVRTC, GPU_UUID, out / 'closed-loop.json'], out,
                              {'CUDA_VISIBLE_DEVICES': GPU_UUID}, timeout=1800))
    receipt['steps'] = steps
    receipt['inputs_verified_after'] = verify_inputs(label)
    (out / ('gpu-receipt-%d.json' % int(time.time()))).write_text(json.dumps(receipt, indent=2) + '\n')


def stage_report(label):
    base = OUT_ROOT / label
    admission('report')
    cpu_summary = json.loads((base / 'cpu/cpu-summary.json').read_text())
    parity, timing = {}, {}
    seal_exe = base / 'build/cpu_reference'
    for fx in FIXTURES:
        fdir = base / 'cpu' / fx
        gdir = base / 'gpu' / fx
        if not fdir.is_dir():
            continue
        cj = json.loads((fdir / 'cpu.json').read_text())
        entry = {'carriers': cj['carriers'], 'mode': cj['mode'], 'sealed_ran': cj['sealed']['ran'], 'sealed_refusal': cj['sealed'].get('refusal', ''),
                 'replica_refusal': cj['replica']['refusal'], 'replica_refusal_rank': cj['replica']['refusal_rank']}
        if cj['sealed']['ran']:
            entry['replica_vs_sealed'] = compare.compare(str(fdir / 'sealed.bin'), str(fdir / 'replica.bin'), subset=True)
            entry['hostport_vs_sealed'] = compare.compare(str(fdir / 'sealed.bin'), str(fdir / 'hostport.bin'), subset=True)
        entry['hostport_vs_replica'] = compare.compare(str(fdir / 'replica.bin'), str(fdir / 'hostport.bin'))
        runs = sorted(gdir.glob('run*.bin')) if gdir.is_dir() else []
        entry['gpu_runs'] = []
        for r in runs:
            rj = json.loads(r.with_suffix('.json').read_text())
            item = {'run': r.name, 'sha256': sha(r), 'repeat_outputs_identical': rj.get('repeat_outputs_identical'),
                    'teardown_released': rj.get('teardown', {}).get('released'), 'quarantined': rj.get('teardown', {}).get('quarantined'),
                    'vs_replica': compare.compare(str(fdir / 'replica.bin'), str(r))}
            if cj['sealed']['ran']:
                item['vs_sealed'] = compare.compare(str(fdir / 'sealed.bin'), str(r), subset=True)
            entry['gpu_runs'].append(item)
        shas = {g['sha256'] for g in entry['gpu_runs']}
        entry['gpu_fresh_process_outputs_identical'] = len(entry['gpu_runs']) == FRESH_PROCESSES and len(shas) == 1
        if runs and cj['mode'] == 'record' and not cj['replica']['refusal']:
            res = subprocess.run([str(FROZEN), str(seal_exe), 'seal', str(fdir), str(runs[0])], capture_output=True, text=True,
                                 env=dict(ENV, CUDA_VISIBLE_DEVICES=''), cwd=str(gdir))
            try:
                entry['gpu_record_seal_sha256'] = json.loads(res.stdout)['seal_sha256']
            except Exception:
                entry['gpu_record_seal_error'] = res.stdout[-2000:] + res.stderr[-2000:]
            entry['replica_record_seal_sha256'] = cj['replica'].get('record_seal_sha256')
            entry['sealed_record_seal_sha256'] = cj['sealed'].get('seal_sha256')
            entry['gpu_seal_equals_replica'] = entry.get('gpu_record_seal_sha256') == entry['replica_record_seal_sha256']
            if entry['sealed_record_seal_sha256']:
                entry['gpu_seal_equals_sealed'] = entry.get('gpu_record_seal_sha256') == entry['sealed_record_seal_sha256']
        ok = all(g['vs_replica']['identical'] and g.get('vs_sealed', {'identical': True})['identical'] for g in entry['gpu_runs'])
        entry['gpu_bit_identical_to_cpu'] = bool(runs) and ok
        parity[fx] = entry
        # timing
        if runs:
            rj = [json.loads(r.with_suffix('.json').read_text()) for r in runs]
            per_process = []
            for j in rj:
                rs = j['runs']
                steady = rs[1:] if len(rs) > 1 else rs
                per_process.append({'first_call_total_wall_ms': rs[0]['total_wall_ms'],
                                    'median_total_wall_ms': statistics.median(x['total_wall_ms'] for x in steady),
                                    'median_field_and_run_wall_ms': statistics.median(x['field_and_run_wall_ms'] for x in steady),
                                    'median_gpu_phases_ms': {k: statistics.median(x['gpu_phases_ms'][k] for x in steady) for k in steady[0]['gpu_phases_ms']},
                                    'median_host_fold_ms': statistics.median(x['host_fold_ms'] for x in steady),
                                    'setup_ms': j['setup_ms'], 'compile_ms': j['compile']['compile_ms'],
                                    'own_peak_device_bytes': j['memory']['own_peak_device_bytes'],
                                    'device_peak_used_by_all_processes': j['memory']['device_peak_used_by_all_processes'],
                                    'own_pinned_host_bytes': j['memory']['own_pinned_host_bytes']})
            timing[fx] = {'carriers': cj['carriers'], 'mode': cj['mode'], 'cpu': {'sealed': cj['sealed'].get('timing'), 'replica_phases_ms': cj['replica']['phases_ms'],
                          'replica_total_ms': cj['replica']['total_ms'], 'census_sha256_ms': cj.get('census_sha256_ms'),
                          'record_seal_ms': cj['replica'].get('record_seal_ms'), 'host_order_census_ms': cj.get('host_order_census_ms')},
                          'gpu_processes': per_process,
                          'device': rj[0]['device'], 'compile': {k: rj[0]['compile'][k] for k in ['nvrtc_version', 'arch_option', 'image_kind', 'options', 'image_sha256', 'nvrtc_library_sha256', 'nvrtc_builtins']}}
    (base / 'parity.json').write_text(json.dumps(parity, indent=1) + '\n')
    lines = ['%-24s %9s %12s %12s %12s %12s %12s %10s %9s' % ('fixture', 'carriers', 'cpu_sealed', 'cpu_numeric', 'gpu_total', 'gpu_steady',
                                                           'gpu_events', 'speedup', 'peak_MiB')]
    for fx, t in timing.items():
        if t['mode'] != 'record':
            continue
        sealed = (t['cpu']['sealed'] or {}).get('prepare_canonical_motion_ms')
        free = (t['cpu']['sealed'] or {}).get('charge_of_endpoints_ms')
        ph = t['cpu']['replica_phases_ms']
        numeric = sum(ph.get(k, 0) for k in ['census_loop_ms', 'motion_loop_ms', 'deposit_loop_ms', 'budget_ms', 'free_charge_ms'])
        gt = statistics.median(pp['median_total_wall_ms'] for pp in t['gpu_processes'])
        gs = statistics.median(pp['median_field_and_run_wall_ms'] for pp in t['gpu_processes'])
        ge = statistics.median(pp['median_gpu_phases_ms']['gpu_total_run'] for pp in t['gpu_processes'])
        peak = max(pp['own_peak_device_bytes'] for pp in t['gpu_processes']) / 2 ** 20
        ref = (sealed + free) if sealed is not None and free is not None else numeric
        lines.append('%-24s %9d %12s %12.3f %12.3f %12.3f %12.3f %10.1f %9.1f' % (fx, t['carriers'], '%.3f' % (sealed + free) if sealed is not None and free is not None else ('cap' if t['carriers'] > 100000 else 'n/a'), numeric, gt, gs, ge, ref / gt, peak))
    lines.append('ms; cpu_sealed = prepare_canonical_motion + charge_of_endpoints (sealed, includes census hash and record JSON seal);')
    lines.append('cpu_numeric = replica phases (census loop, motion, deposits, budget, free charge); gpu_total = median wall of')
    lines.append('set_carriers + set_field + run incl. all transfers and host folds; gpu_steady = set_field + run (carriers resident);')
    lines.append('gpu_events = stream time of run(); speedup = cpu_sealed (or cpu_numeric above the cap) / gpu_total.')
    (base / 'timing-table.txt').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))
    (base / 'timing.json').write_text(json.dumps(timing, indent=1) + '\n')
    selftest = json.loads((base / 'cpu/selftest.json').read_text())
    summary = {'selftest': selftest, 'fixtures': {fx: {'gpu_bit_identical_to_cpu': p['gpu_bit_identical_to_cpu'],
                                                       'fresh_process_identical': p['gpu_fresh_process_outputs_identical'],
                                                       'replica_vs_sealed': p.get('replica_vs_sealed', {}).get('identical'),
                                                       'hostport_vs_replica': p['hostport_vs_replica']['identical'],
                                                       'refusal': p['replica_refusal'], 'sealed_refusal': p['sealed_refusal'],
                                                       'gpu_seal_equals_sealed': p.get('gpu_seal_equals_sealed'),
                                                       'gpu_seal_equals_replica': p.get('gpu_seal_equals_replica')} for fx, p in parity.items()}}
    loop = base / 'gpu/closed-loop.json'
    if loop.is_file():
        lj = json.loads(loop.read_text())
        summary['closed_loop'] = {'all_equal': lj.get('all_equal'), 'error': lj.get('error'),
                                  'cases': [{k: v for k, v in c.items()} for c in lj.get('cases', [])]}
    (base / 'parity-summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    stage, label = sys.argv[1], sys.argv[2]
    if not label or '/' in label or label in ('.', '..'):
        raise SystemExit('label required')
    rest = sys.argv[3:] or FIXTURES
    if stage == 'build':
        stage_build(label)
    elif stage == 'cpu':
        stage_cpu(label, sys.argv[3:])
    elif stage == 'gpu':
        stage_gpu(label, [f for f in rest if f != 'closed-loop'] if sys.argv[3:] else FIXTURES, closed_loop=not sys.argv[3:] or 'closed-loop' in sys.argv[3:])
    elif stage == 'report':
        stage_report(label)
    else:
        raise SystemExit('unknown stage')
