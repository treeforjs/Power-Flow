#!/usr/bin/python -I
"""Build and run the Hall node-plane reproducer with the sealed Hall recipe.

Compiled exactly like the Hall corrective owner test (development/
hall-field-correction-20261030/run_components.py flags, native-v6.2 headers)
plus native-v6.2 mfem_gap_pic.cpp for pic_boris/pic_kinetic_energy_J, as in the
audit field reproducer port. Launch under both core leases:
  /home/treefor/projects/runtimes/core-lease-20261030/both-cores.sh env CUDA_VISIBLE_DEVICES= \
      /usr/bin/python -I repro/run_repro.py LABEL
"""
from pathlib import Path
import fcntl, hashlib, json, os, resource, subprocess, sys, time

HERE = Path(__file__).resolve().parent
R = Path('/home/treefor/projects/runtimes')
H = R / 'hall-coupled-transport-20261029/source'
HD = H / 'development/hall-canonical-driver-20261029'
HC = H / 'development/hall-field-correction-20261030'
NV = HD / 'bindings/native-v6.2'
M = R / 'mfem-4.10-pic-20261029'
F = R / 'arch-base-20261028/run-frozen'
OUT = Path('/home/treefor/sim-workspace/powerflow/outputs/pic-water-common-commit-v4-20261030')
TRACE = H / 'development/hall-finite-foil-em-20261029/data/shot_17133_current.csv'
LEASES = {6: '/home/treefor/sim-workspace/powerflow/outputs/pic-profile/testing-core.lock',
          7: '/home/treefor/sim-workspace/powerflow/outputs/pic-profile/testing-core7.lock'}
ENV = {'PATH': '/usr/bin:/bin', 'LC_ALL': 'C', 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1',
       'MKL_NUM_THREADS': '1', 'CUDA_VISIBLE_DEVICES': '', 'PYTHONNOUSERSITE': '1'}
LIBS = [M / 'install-mpi/lib/libmfem.a', M / 'tpl/hypre/lib/libHYPRE.a', M / 'tpl/metis/lib/libmetis.a', M / 'tpl/GKlib/lib/libGKlib.a']


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


cores = os.sched_getaffinity(0)
if cores not in ({6, 7}, {6}, {7}) or os.getpriority(os.PRIO_PROCESS, 0) != 15 or os.getenv('CUDA_VISIBLE_DEVICES') != '':
    raise SystemExit('core 6/7 lease, nice 15 and no GPU required')
for core in cores:
    with open(LEASES[core], 'a') as probe:
        try:
            fcntl.flock(probe, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            continue
        raise SystemExit('core lease must be held by the launching wrapper')
resource.setrlimit(resource.RLIMIT_AS, (3 * 1024 ** 3, 3 * 1024 ** 3))
label = sys.argv[1]
head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=H, text=True).strip()
dirty = subprocess.check_output(['git', 'status', '--porcelain'], cwd=H, text=True)
if head != '8dddec0610982778bf02bdfe697b70a37129c012' or dirty:
    raise SystemExit('Hall repository is not the clean corrective commit')
out = OUT / ('repro-hall-node-plane-' + label)
out.mkdir(parents=True, exist_ok=False)
exe = out / 'repro_hall_node_plane'
flags = ['-std=c++20', '-O2', '-DNDEBUG', '-ffp-contract=off', '-fno-fast-math', '-ffunction-sections', '-fdata-sections',
         '-I' + str(NV / 'source/include'), '-isystem', str(NV / 'inputs'), '-I' + str(NV / 'vendor'), '-I' + str(HC)]
for d in ['install-mpi/include', 'tpl/hypre/include', 'tpl/metis/include']:
    flags += ['-isystem', str(M / d)]
sources = [HERE / 'repro_hall_node_plane.cpp', HD / 'native_arena_motion.cpp', HD / 'kinetic_boundary.cpp',
           NV / 'source/src/mfem_gap_pic.cpp']
commands = [('build', ['/usr/bin/c++', *flags, *map(str, sources), *map(str, LIBS), '/usr/lib/libmpi.so', '-lcrypto',
                       '-Wl,--gc-sections', '-o', str(exe)]),
            ('run', [str(exe), str(TRACE)])]
inputs = {str(p): sha(p) for p in [*sources, TRACE, *LIBS]}
inputs.update({str(p): sha(p) for p in HC.glob('*.hpp')})
inputs.update({str(p): sha(p) for p in [HD / 'motion_receipt.hpp', HD / 'kinetic_boundary.hpp', HD / 'axial_joint_v2.hpp',
                                         HD / 'axial_thermal_v2.hpp']})
receipt = {'schema': 'hall-node-plane-reproducer-receipt-v1', 'hall_commit': head, 'hall_git_clean': True,
           'cpu_affinity': sorted(cores), 'nice': 15, 'gpu': False, 'rlimit_as_bytes': 3 * 1024 ** 3,
           'environment': ENV, 'frozen_launcher_sha256': sha(F), 'inputs_sha256': inputs, 'commands': []}
status = 'passed'
for name, cmd in commands:
    log = out / (name + '.log')
    start = time.monotonic()
    with log.open('w') as handle:
        rc = subprocess.run([str(F), *cmd], cwd=HERE.parent, env=ENV, stdout=handle, stderr=subprocess.STDOUT, timeout=900).returncode
    receipt['commands'].append({'name': name, 'argv': [str(F), *cmd], 'returncode': rc, 'elapsed_s': time.monotonic() - start,
                                'log': str(log), 'log_sha256': sha(log)})
    if rc != 0:
        status = 'failed'
        break
if status == 'passed':
    text = (out / 'run.log').read_text()
    receipt['observations'] = json.loads(text[text.index('\n['):])
    receipt['binary_sha256'] = sha(exe)
receipt['inputs_unchanged_after_run'] = all(sha(p) == h for p, h in inputs.items())
receipt['status'] = status
(out / 'receipt.json').write_text(json.dumps(receipt, indent=1) + '\n')
print(json.dumps({'reproducer': str(out), 'status': status}))
