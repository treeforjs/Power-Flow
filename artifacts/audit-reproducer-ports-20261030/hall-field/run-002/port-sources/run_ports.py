"""Build and run the field audit reproducer port against Hall commit 8dddec0.

Launch (lease held externally, CPU6 only):
  flock /home/treefor/sim-workspace/powerflow/outputs/pic-profile/testing-core.lock \
    taskset -c 6 nice -n 15 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    CUDA_VISIBLE_DEVICES= PYTHONNOUSERSITE=1 /usr/bin/python -I \
    /home/treefor/projects/runtimes/audit-reproducer-ports-20261030/hall-field/run_ports.py NEW_LABEL

Every compile/run step is executed as
  run-frozen prlimit --as=2147483648 --cpu=180 -- CMD   (240 s wall timeout)
with the flags and include paths of development/hall-field-correction-20261030/run_components.py,
once without and once with -fsanitize=undefined -fno-sanitize-recover=undefined.
Nothing in the Hall repository, the audit directory or earlier evidence is written.
"""
from pathlib import Path
import fcntl, hashlib, json, os, platform, shutil, subprocess, sys, time

HERE = Path(__file__).resolve().parent
REPO = Path('/home/treefor/projects/runtimes/hall-coupled-transport-20261029/source')
COMMIT = '8dddec0610982778bf02bdfe697b70a37129c012'
MODULES = ['hall-field-correction-20261030', 'hall-canonical-driver-20261029',
           'hall-finite-foil-em-20261029', 'hall-coupled-transport-20261029']
AUDIT = Path('/home/treefor/projects/Power-Flow/artifacts/pre-integration-audit-20261030')
PROBE = AUDIT / 'reviews/field'
OUT_ROOT = Path('/home/treefor/sim-workspace/powerflow/outputs/audit-reproducer-ports-20261030/hall-field')
LEASE = '/home/treefor/sim-workspace/powerflow/outputs/pic-profile/testing-core.lock'
FROZEN = '/home/treefor/projects/runtimes/arch-base-20261028/run-frozen'
DEP = Path('/home/treefor/projects/runtimes/mfem-4.10-pic-20261029')
LIBS = [DEP / 'install-mpi/lib/libmfem.a', DEP / 'tpl/hypre/lib/libHYPRE.a',
        DEP / 'tpl/metis/lib/libmetis.a', DEP / 'tpl/GKlib/lib/libGKlib.a']
CLEAN_ENV = {'PATH': '/usr/bin:/bin', 'LANG': 'C.UTF-8', 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1',
             'MKL_NUM_THREADS': '1', 'CUDA_VISIBLE_DEVICES': '', 'PYTHONNOUSERSITE': '1'}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def git_blob(p):
    data = Path(p).read_bytes()
    return hashlib.sha1(b'blob %d\0' % len(data) + data).hexdigest()


def fail(message):
    raise RuntimeError(message)


# ------------------------------------------------------------------ admission of the execution policy
if os.sched_getaffinity(0) != {6} or os.getpriority(os.PRIO_PROCESS, 0) != 15:
    fail('CPU6 affinity and nice 15 required')
if any(os.getenv(k) != '1' for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']) \
        or os.getenv('CUDA_VISIBLE_DEVICES') != '':
    fail('one thread per library and no GPU required')
with open(LEASE, 'a') as probe_handle:
    try:
        fcntl.flock(probe_handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        pass
    else:
        fail('CPU6 lease must be held by the launching flock; no work started')
if len(sys.argv) != 2 or not sys.argv[1] or '/' in sys.argv[1] or sys.argv[1] in ('.', '..'):
    fail('single new run label required')
label = sys.argv[1]

head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True).strip()
dirty = subprocess.check_output(['git', 'status', '--porcelain'], cwd=REPO, text=True)
if head != COMMIT or dirty:
    fail('Hall repository is not the clean corrective commit')
if Path('/proc/sys/kernel/osrelease').read_text().strip() != '7.2.11-arch1-1':
    fail('kernel release differs from the frozen runtime')

out = OUT_ROOT / label
out.mkdir(parents=True, exist_ok=False)
receipt = {'schema': 'hall-field-audit-port-receipt-v1', 'label': label, 'runner_argv': sys.argv,
           'runner_sha256': sha(__file__), 'python': sys.executable, 'python_version': platform.python_version(),
           'started_unix': time.time(), 'hall_commit': head, 'hall_git_clean': True,
           'kernel_release': '7.2.11-arch1-1', 'cpu_affinity': sorted(os.sched_getaffinity(0)),
           'nice': os.getpriority(os.PRIO_PROCESS, 0), 'lease': LEASE, 'lease_held_by': 'launching flock (verified busy)',
           'environment': CLEAN_ENV, 'frozen_launcher': FROZEN, 'frozen_launcher_sha256': sha(FROZEN)}

# ------------------------------------------------------------------ frozen runtime admission (same pins as run_components.py)
pins_path = REPO / 'development/hall-field-correction-20261030/runtime-pins.json'
pins = json.loads(pins_path.read_text())
frozen_root = Path(pins['frozen_root'])
base = frozen_root.parent
if sha(base / 'root-files.json') != pins['frozen_root_manifest_sha256'] or sha(base / 'run-frozen') != pins['launcher_sha256']:
    fail('frozen runtime manifest or launcher changed')
if Path('/proc/sys/kernel/osrelease').read_text().strip() != pins['kernel_release']:
    fail('kernel release differs from runtime pins')
for name, data in pins['selected_files'].items():
    p = frozen_root / name
    if 'sha256' in data:
        if sha(p) != data['sha256']:
            fail('frozen runtime file changed: ' + name)
    elif 'symlink' in data:
        if not p.is_symlink() or os.readlink(p) != data['symlink']:
            fail('frozen runtime link changed: ' + name)
    else:
        fail('unknown runtime pin entry')
for name, expected in pins['static_libraries'].items():
    if sha(name) != expected:
        fail('pinned static library changed: ' + name)
receipt['runtime_admission'] = {'status': 'passed', 'runtime_pins': str(pins_path), 'runtime_pins_sha256': sha(pins_path),
                                'runtime_entries': len(pins['selected_files']), 'static_libraries': pins['static_libraries']}

# ------------------------------------------------------------------ snapshot of the exact commit plus the port translation unit
source = out / 'source'
tree = {}
for line in subprocess.check_output(['git', 'ls-tree', '-r', COMMIT, '--', *['development/' + m for m in MODULES]],
                                    cwd=REPO, text=True).splitlines():
    meta, path = line.split('\t', 1)
    tree[path] = meta.split()[2]
for name in MODULES:
    shutil.copytree(REPO / 'development' / name, source / 'development' / name, symlinks=True)
snapshot_files = sorted(p for p in source.rglob('*') if p.is_file() and not p.is_symlink())
mismatch = [str(p) for p in snapshot_files if tree.get(str(p.relative_to(source))) != git_blob(p)]
if mismatch or len(snapshot_files) != len(tree):
    fail('snapshot differs from commit tree: %d mismatched, %d files vs %d blobs' % (len(mismatch), len(snapshot_files), len(tree)))
module = source / 'development/hall-field-correction-20261030'
old = source / 'development/hall-canonical-driver-20261029'
native = old / 'bindings/native-v6.2'
trace = source / 'development/hall-finite-foil-em-20261029/data/shot_17133_current.csv'
# Placed beside test_owner.cpp so that quoted-include resolution is identical to the corrective owner test.
port_tu = module / 'port_field_audit.cpp'
shutil.copy2(HERE / 'port_field_audit.cpp', port_tu)
provenance = out / 'port-sources'
provenance.mkdir()
for name in ['port_field_audit.cpp', 'run_ports.py']:
    shutil.copy2(HERE / name, provenance / name)

audit_copies = PROBE / 'inputs/hall-source/development/hall-finite-foil-em-20261029'
audit_identity = {name: {'audit_capture_sha256': sha(audit_copies / name),
                         'snapshot_sha256': sha(source / 'development/hall-finite-foil-em-20261029' / name)}
                  for name in ['axial_fields.hpp', 'axial_mesh.hpp', 'axial_particles.hpp', 'foil_em.hpp']}
for v in audit_identity.values():
    v['identical'] = v['audit_capture_sha256'] == v['snapshot_sha256']
if not all(v['identical'] for v in audit_identity.values()):
    fail('unchanged operator headers differ from the audit capture')

manifest = {str(p): sha(p) for p in source.rglob('*') if p.is_file() and not p.is_symlink()}
for d in ['install-mpi/include', 'tpl/hypre/include', 'tpl/metis/include']:
    manifest.update({str(p): sha(p) for p in (DEP / d).rglob('*') if p.is_file()})
manifest.update({str(p): sha(p) for p in LIBS})
manifest['/usr/lib/libmpi.so'] = sha('/usr/lib/libmpi.so')
(out / 'input-sha256.json').write_text(json.dumps(manifest, indent=2) + '\n')
headers = {str(p.relative_to(source)): sha(p) for p in sorted(source.rglob('*'))
           if p.is_file() and not p.is_symlink() and p.suffix in ('.hpp', '.h', '.cpp', '.inc')}
receipt['snapshot'] = {'directory': str(source), 'files': len(snapshot_files), 'verified_against_git_tree': COMMIT,
                       'port_translation_unit': str(port_tu), 'port_translation_unit_sha256': sha(port_tu),
                       'trace': str(trace), 'trace_sha256': sha(trace)}
receipt['source_header_sha256'] = headers
receipt['audit_operator_header_identity'] = audit_identity
receipt['audit_probe_reference'] = {str(PROBE / n): sha(PROBE / n) for n in
                                    ['probe_operators.cpp', 'operator-observations.json', 'AUDIT.txt', 'run_probe.py']}
receipt['audit_text_sha256'] = sha(AUDIT / 'AUDIT.txt')
receipt['input_manifest'] = {'path': str(out / 'input-sha256.json'), 'entries': len(manifest),
                             'sha256': sha(out / 'input-sha256.json')}

# ------------------------------------------------------------------ build and run, plain and UBSan
base_flags = ['-std=c++20', '-O2', '-DNDEBUG', '-ffp-contract=off', '-fno-fast-math', '-ffunction-sections', '-fdata-sections',
              '-I' + str(native / 'source/include'), '-isystem', str(native / 'inputs'), '-I' + str(native / 'vendor')]
build = out / 'build'
build.mkdir()
commands = []
for variant in ['plain', 'ubsan']:
    flags = list(base_flags)
    if variant == 'ubsan':
        flags += ['-fsanitize=undefined', '-fno-sanitize-recover=undefined']
    for d in ['install-mpi/include', 'tpl/hypre/include', 'tpl/metis/include']:
        flags += ['-isystem', str(DEP / d)]
    binary = build / ('port-' + variant)
    run_dir = out / ('run-' + variant)
    run_dir.mkdir()
    # mfem_gap_pic.cpp supplies pic_boris/pic_kinetic_energy_J for prepare_canonical_motion; the audit
    # probe compiled the same native-v6.2 source. The corrective owner test never references that path.
    commands.append(('build-' + variant, ['/usr/bin/c++', *flags, str(port_tu), str(old / 'native_arena_motion.cpp'),
                     str(old / 'kinetic_boundary.cpp'), str(native / 'source/src/mfem_gap_pic.cpp'),
                     *map(str, LIBS), '/usr/lib/libmpi.so', '-lcrypto',
                     '-Wl,--gc-sections', '-o', str(binary)], 0))
    commands.append(('run-' + variant, [str(binary), str(trace), str(run_dir)], 0))

receipt['commands'] = []
status = 'passed'
error = None
try:
    for name, cmd, expected in commands:
        actual = [FROZEN, 'prlimit', '--as=2147483648', '--cpu=180', '--', *cmd]
        start = time.monotonic()
        log_path = out / (name + '.log')
        with log_path.open('w') as log:
            try:
                r = subprocess.run(actual, cwd=out, env=CLEAN_ENV, stdout=log, stderr=subprocess.STDOUT, timeout=240)
                rc = r.returncode
            except subprocess.TimeoutExpired:
                rc = 'timeout'
        entry = {'name': name, 'argv': actual, 'returncode': rc, 'expected_returncode': expected,
                 'elapsed_s': time.monotonic() - start, 'log': str(log_path), 'log_sha256': sha(log_path),
                 'CPU': 6, 'nice': 15, 'threads': 1, 'lease': LEASE}
        if name.startswith('build-'):
            binary = Path(cmd[cmd.index('-o') + 1])
            entry['binary_sha256'] = sha(binary) if binary.exists() else None
        receipt['commands'].append(entry)
        (out / 'commands.json').write_text(json.dumps(receipt['commands'], indent=2) + '\n')
        print(name, rc, round(entry['elapsed_s'], 1), flush=True)
        if rc != expected:
            fail('step failed: ' + name + '\n' + log_path.read_text()[-6000:])
    changed = [p for p, h in manifest.items() if sha(p) != h]
    if changed:
        fail('inputs changed during the run: ' + ', '.join(changed[:5]))
    receipt['inputs_unchanged_after_run'] = True
except Exception as e:  # receipt is written on every exit path
    status, error = 'failed', str(e)

# ------------------------------------------------------------------ per-counterexample results
observations = {}
for variant in ['plain', 'ubsan']:
    path = out / ('run-' + variant) / 'port-observations.json'
    if path.exists():
        observations[variant] = json.loads(path.read_text())
receipt['observations'] = {v: {'path': str(out / ('run-' + v) / 'port-observations.json'),
                               'sha256': sha(out / ('run-' + v) / 'port-observations.json')} for v in observations}
if len(observations) == 2:
    receipt['plain_and_ubsan_observations_identical'] = observations['plain'] == observations['ubsan']


def outcome(entry, expect_accept):
    accepted = bool(entry.get('accepted'))
    return {'accepted': accepted, 'expected_accepted': expect_accept, 'matches_expected': accepted == expect_accept,
            'rejection': entry.get('rejection')}


if 'plain' in observations:
    o = observations['plain']
    g = o['A09_gauss']
    d = o['A09_divB']
    f = o['A10_field_step']
    orig = o['original_operator']
    results = {
        'fixture_identity': o['fixture_identity'],
        'A09_inconsistent_initial_charge': {
            'original': {'operator_accepts_particle_charge_operand': orig['A09_gauss']['api_operand_for_particle_charge'],
                         'field_vs_manufactured_charge_max_C': orig['initial_field_vs_manufactured_authored_charge_max_C'],
                         'mismatch_after_accepted_Fields_prepare_C': orig['A09_gauss']['field_vs_manufactured_charge_max_after_one_step_C'],
                         'gauss_increment_C': orig['A09_gauss']['gauss_increment_C']},
            'corrected': {k: outcome(g[k], k.startswith('control')) for k in
                          ['control_admit_field_without_free_charge', 'control_owner_empty_census',
                           'Model_admit_census_charge_dc_declaration', 'Model_admit_census_charge_absorbing_declaration',
                           'Model_prepare_from_inconsistent_state', 'Owner_actual_census_with_matching_ledger',
                           'Owner_actual_census_with_absorbing_ledger', 'Owner_actual_census_with_field_only_ledger',
                           'Owner_empty_census_with_charged_ledger']},
            'defect_vs_bound': {'max_defect_C': g['gauss_scan_dc_declaration']['max_abs_defect_C'],
                                'bound_at_max_C': g['gauss_scan_dc_declaration']['bound_at_max_defect_C'],
                                'solver_bound_C': g['solver_bound_C'],
                                'audit_mismatch_over_solver_bound': g['audit_mismatch_over_solver_bound'],
                                'nodes_violating_dc_declaration': g['gauss_scan_dc_declaration']['nodes_violating_bound'],
                                'nodes_violating_absorbing_declaration': g['gauss_scan_absorbing_declaration']['nodes_violating_bound']},
            'census_rows_constructed_by_port_not_native_arena': True},
        'A09_divB_50_T_per_m': {
            'original': orig['A09_divB'],
            'corrected': {k: outcome(d[k], k.endswith('control')) for k in
                          ['Model_admit', 'Model_prepare_from_perturbed_state', 'Owner_bootstrap_empty_census',
                           'restore_unmodified_checkpoint_control', 'restore_resealed_perturbed_checkpoint']},
            'perturbed_divB_scan': d['divB_scan_perturbed'],
            'smallest_detected_delta_T': min([s['delta_T'] for s in d['supplementary_detection_sweep'] if not s['Model_admit']['accepted']], default=None)},
        'A10_midpoint': {
            'original': orig['A10_midpoint'],
            'corrected_field_step': {k: f[k] for k in ['current_A', 'current_minus_probe_current_A', 'E_high_word_vs_authoritative',
                                                       'B_high_word_vs_authoritative', 'E_authoritative_equals_compensated_average',
                                                       'B_authoritative_equals_compensated_average', 'Maxwell_receipt_mid_E_equals_authoritative',
                                                       'joule_cell_power_from_authoritative',
                                                       'original_operator_same_current_high_word_vs_receipt_mid_E']},
            'closed_consumers': o['A10_closed_consumers'],
            'source_consumer': o['A10_source_consumer']},
        'gauss_bound_reconstruction': o['gauss_bound_reconstruction']}
    receipt['results'] = results

receipt['status'] = status
receipt['error'] = error
receipt['finished_unix'] = time.time()
receipt['scope'] = {'native_arena_constructed': False, 'actual_SSA': False, 'source_population': False,
                    'joint_native_publication': False, 'physical_run': False, 'hall_repo_written': False,
                    'audit_directory_written': False}
(out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
(out / 'run-status.json').write_text(json.dumps({'status': status, 'error': error}, indent=2) + '\n')
print('receipt', out / 'receipt.json', status, flush=True)
if status != 'passed':
    sys.exit(1)
