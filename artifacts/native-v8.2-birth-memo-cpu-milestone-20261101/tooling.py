"""Shared build/run helpers for native v8 (candidate008). Every compiler and
test process runs through the frozen launcher with cwd under
/home/treefor/projects, inside a caller-held core lease (core-lease-20261030
both-cores.sh / any-core.sh --wait for cores 6-7; the lent core 3 through its
testing-core3.lock with taskset and nice 15, never while RETURN-CORES-0-3
exists), at nice 15 with one thread. CPU processes see no GPU; the GPU forms
test runs under the two GPU locks. Normal processes run under RLIMIT_AS
(3 GiB unless a command states more); ASan/UBSan processes use
hard_rss_limit_mb instead."""
from pathlib import Path
import json, os, subprocess, hashlib, time, datetime

P = Path(__file__).resolve().parent
RT = P.parent
V5 = RT / 'native-full-time-candidate005'
V5_IDENTITY = 'cca72822c9b61900aae7572f29c5d57314f45aef17c4b01c9fd2346c6fa3c883'
V5_LIB_SHA = '715729114e5a7d3975bd7ee7db92c570e61f78c907f3e999812b34e7a74159f2'
V5_ADAPTER_SHA = '072e0f5fcdac870cc91a197e7667abc8d72aee88da41621ca6408c458040f518'
V6 = RT / 'native-full-time-candidate006'
V6_IDENTITY = 'aa4d8f03100e64bed759c739013d91cac1bd66c8a9b2ab6bc3105af7973452c9'
V6_LIB_SHA = '1c2525958457dfe5ff16def07341cd000625e1eefaf38d44a054766efc3d86bc'
V6_ADAPTER_SHA = '7fcac96a019414e23905804e181b330b73bef84918397df2923c65f1cf7551e2'
V6_V3_ADAPTER_SHA = '597c7e1eaa2992dd0e520b9793424cc41cd5936939ec2997e3b3d2834d2e648f'
EXACT = RT / 'native-exact-count'
OUT = Path('/home/treefor/projects/.workspace/outputs/native-v8.2-birth-memo-20261101')
V71 = RT / 'native-full-time-candidate007.1'
V71_IDENTITY = '86b4500195013def19e5c4e45e398dacf75c9c44b23007b02130e7b40178fad8'
V71_LIB_SHA = 'cce8d900a0106b3eb1882c84d4f95f1f3a96f28e031b94a26970d01ca539a741'
V71_OUT = Path('/home/treefor/projects/.workspace/outputs/native-v7.1-20261031')
CHEM = Path('/home/treefor/projects/runtimes/pic-water-chemistry-v3.1.1-20261031')
CHEM_V3 = Path('/home/treefor/projects/runtimes/pic-water-chemistry-v3-20261031')
CHEM_BLOB = CHEM_V3 / 'data/tables/chemistry-v3.blob'
CHEM_BLOB_SHA = '3ca3ea495630e6b6e125c05ebf45488fd89a899544a2aa138cefb9b1dd415b0a'
CHEM31_BLOB = CHEM / 'data/tables/chemistry-v3.1.blob'
CHEM31_BLOB_SHA = '70fcc2adc2fdaa62abed7ed02d68ffc09b1b18554eda8fa99c5451120ae454b8'
CHEM31_POLICY = CHEM / 'data/tables/chemistry-v3.1-coverage-policy.json'
CHEM31_POLICY_SHA = '15d44a91b8af93e0b14beb89434c55dc402e7196f10ca4bdaa2357564474b2d3'
CHEM_HEADER_SHA = {'chem_host.hpp': '87715119d2837b52b57162bd40a30687bc8cc109f6d170de5f8dd0f2ae11691e',
                   'chem_event.h': '95f460ec3111d93dbe9ca38e7737cbb95e80a3d3fcfee912e7745be3596d8b3b'}
MESHES = Path('/home/treefor/projects/runtimes/resolution-series-prep-20261031/meshes')
MERGE_SPLIT = Path('/home/treefor/projects/runtimes/resolution-series-prep-20261031/merge-split')
V7 = RT / 'native-full-time-candidate007'
V7_IDENTITY = '3f88fab5e285a13448ab709394dabc3cb65aab4820c76dec27c8121c0b042dc5'
V7_LIB_SHA = 'af2bd9e7dc4932421d7928af7716cc076af726c61834d1d2dc2d690582348dec'
V7_OUT = Path('/home/treefor/projects/.workspace/outputs/native-v7-20261031')
V6_OUT = Path('/home/treefor/projects/.workspace/outputs/pic-profile/'
              'reaction-transactions-20261029/native-full-time-candidate006')
GPU_COMMON = Path('/home/treefor/projects/runtimes/hall-gpu-resident-v1-20261030/common')
CUDA_INC = Path('/home/treefor/projects/runtimes/mfem-cuda-benchmark-20261029/cuda-libs/nvidia/cuda_runtime/include')
NVRTC_DIR = Path('/home/treefor/projects/runtimes/pic-cuda-gather-investigation-20261029/tools/nvrtc/nvidia/cuda_nvrtc')
NVRTC = NVRTC_DIR / 'lib/libnvrtc.so.12'
GPU_UUID = 'GPU-a7b9b458-167f-cd57-fec9-3dea0995b6e9'
LEASES = Path('/home/treefor/projects/.workspace/outputs/pic-profile')
# Hall v3.2 reference of motion-work-gross-v1 (read only; tolerance bit cross-check).
HALL_V32_WORK_GROSS = Path('/home/treefor/projects/runtimes/hall-corrective-v3.2-20261031/source/motion_work_gross_v1.hpp')
HALL_V32_SPEC = Path('/home/treefor/projects/runtimes/hall-corrective-v3.2-20261031/MOTION-WORK-GROSS-V1.txt')
FROZEN = '/home/treefor/projects/runtimes/arch-base-20261028/run-frozen'
MFEM = Path('/home/treefor/projects/runtimes/mfem-4.10-pic-20261029')
VENDOR = RT / 'source/vendor'
HALL = Path('/home/treefor/projects/runtimes/hall-coupled-transport-20261029/source/development')
HALL_DRIVER = HALL / 'hall-canonical-driver-20261029'
HALL_V2 = Path('/home/treefor/projects/runtimes/hall-corrective-v2-20261030')
ARRIVAL = Path('/home/treefor/projects/runtimes/pic-arrival-owner-v1-20261030')
# v8.1: merge-split-transaction-v1 component (CPU reference compiled in place; pinned).
MERGE_SPLIT = Path('/home/treefor/projects/runtimes/resolution-series-prep-20261031/merge-split')
MERGE_SPLIT_SHA = {'include/ms_core.h': 'e832545f159eded2afaeea771ae6b8f4d3926b55768c68c124791272ff801997',
                   'include/merge_split.hpp': 'e9f593837a4ec1e213754092ffa207ac3b49f9e12b64bd9fc797ba399d1e3663',
                   'src/merge_split.cpp': '906d33cc9e23455a91ef303a04fcf180bfb235188f26fc02a9d9fd98af623609'}
V3_COMMIT = Path('/home/treefor/projects/runtimes/pic-water-common-commit-v3-20261030')
S5 = Path('/home/treefor/projects/runtimes/pic-electrode-heat-source-v5-20261030')
S5_MOD = S5 / 'native/mfem_electrode_heat_source_v5'
S5_SEAL = Path('/home/treefor/projects/Power-Flow/artifacts/source-owner-v5-20261030/source-seal-v5.json')
SV3 = Path('/home/treefor/projects/runtimes/pic-electrode-heat-feedback-20261029/native/mfem_electrode_heat_source_v3')
SV2 = Path('/home/treefor/projects/runtimes/pic-electrode-heat-feedback-20261029/native/mfem_electrode_heat_feedback_v2')
COUNT = Path('/home/treefor/projects/runtimes/pic-water-exact-count-subflow-20261029/native/source/include/raikou/exact_particle_count.h')
COUNT_SHA = '46e01fe22e67cc9deb0da9320edef4041f4b3ef433de8ce1eca484463b405501'
HALL_COMMIT = '8dddec0610982778bf02bdfe697b70a37129c012'
RLIMIT_AS = 3 * 1024 ** 3
LIBS = [str(MFEM / 'install-mpi/lib/libmfem.a'), str(MFEM / 'tpl/hypre/lib/libHYPRE.a'),
        str(MFEM / 'tpl/metis/lib/libmetis.a'), str(MFEM / 'tpl/GKlib/lib/libGKlib.a'), '-lcrypto', '-ldl']
SYSTEM = []
for d in ['install-mpi/include', 'tpl/hypre/include', 'tpl/metis/include']:
    SYSTEM += ['-isystem', str(MFEM / d)]
TRAPS = ['_ZN6raikou12mfem_backend6GapPic12solve_fieldsEv', '_ZN6raikou12mfem_backend6GapPic7advanceEd',
         '_ZN6raikou12mfem_backend6GapPic13advance_emptyEd',
         '_ZN6raikou12mfem_backend6GapPic16advance_observedEdRKSt8functionIFvRKNS0_11GapParticleEidEE']
ASAN = ['-fsanitize=address,undefined', '-fno-omit-frame-pointer']
ASAN_ENV = {'ASAN_OPTIONS': 'detect_leaks=1:abort_on_error=1:hard_rss_limit_mb=3072',
            'UBSAN_OPTIONS': 'halt_on_error=1:print_stacktrace=1'}


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def resource_guard(build=False, gpu=False):
    cpus = os.sched_getaffinity(0)
    pairs = [{a, b} for a in (6, 7, 14, 15, 0, 1) for b in (6, 7, 14, 15, 0, 1) if a < b]
    allowed = [{6}, {7}, {0}, {1}, {15}] + (pairs if build else [])
    if os.environ.get('SIMQ_JOB_ID'):
        spec = os.environ.get('SIMQ_CPU_LIST', '')
        scheduled = set()
        for item in spec.split(','):
            if '-' in item:
                lo, hi = map(int, item.split('-')); scheduled.update(range(lo, hi + 1))
            elif item:
                scheduled.add(int(item))
        need(bool(os.environ.get('SIMQ_LAUNCH_TOKEN')) and bool(os.environ.get('SIMQ_UNIT')) and
             os.environ.get('SIMQ_PARTITION') == 'powerflow' and cpus == scheduled and len(cpus) <= 2 and
             os.environ['SIMQ_UNIT'] in Path('/proc/self/cgroup').read_text(),
             'SimQueue allocation must match the process affinity')
    else:
        need(cpus in allowed, f'core lease required (affinity {sorted(cpus)})')
    need(not (cpus & {0, 1, 2, 3} and (LEASES / 'RETURN-CORES-0-3').exists()), 'cores 0-3 must be returned')
    need(os.getpriority(os.PRIO_PROCESS, 0) >= 15, 'nice >= 15 required')
    for name in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']:
        need(os.environ.get(name) == '1', 'one thread required: ' + name)
    if not gpu:
        need(os.environ.get('CUDA_VISIBLE_DEVICES', None) == '', 'CUDA_VISIBLE_DEVICES must be empty')
    return sorted(cpus)


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def tree_hashes(root, suffixes=None):
    root = Path(root)
    return {str(p.relative_to(root)): sha256(p) for p in sorted(root.rglob('*'))
            if p.is_file() and 'build' not in p.relative_to(root).parts and '__pycache__' not in p.parts
            and 'snapshots' not in p.relative_to(root).parts and (suffixes is None or p.suffix in suffixes)}


def sealed_inputs():
    """Hashes of every sealed input read by v7; refuses if a pinned one changed."""
    v5 = json.loads((V5 / 'SOURCE-SHA256.json').read_text())
    need(v5['source_identity'] == V5_IDENTITY, 'native v5 identity changed')
    for name, digest in v5['files'].items():
        need(sha256(V5 / name) == digest, 'sealed native v5 file changed: ' + name)
    v6 = json.loads((V6 / 'SOURCE-SHA256.json').read_text())
    need(v6['source_identity'] == V6_IDENTITY, 'native v6 identity changed')
    for name, digest in v6['files'].items():
        need(sha256(V6 / name) == digest, 'sealed native v6 file changed: ' + name)
    need(sha256(V6 / 'build/libpic_native_v6.a') == V6_LIB_SHA, 'native v6 archive changed')
    need(sha256(V6 / 'build/libpic_hall_motion_adapter_v6.a') == V6_ADAPTER_SHA, 'native v6 adapter changed')
    need(sha256(V6 / 'build/libpic_hall_motion_v3_adapter_v6.a') == V6_V3_ADAPTER_SHA, 'native v6 v3 adapter changed')
    need(sha256(COUNT) == COUNT_SHA, 'ExactParticleCount header changed')
    head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=HALL.parent, text=True).strip()
    need(head == HALL_COMMIT, 'Hall checkout is not 8dddec0')
    for name in ['kinetic_boundary.cpp', 'kinetic_boundary.hpp', 'motion_receipt.hpp', 'native_arena_motion.cpp']:
        need(sha256(HALL_DRIVER / name) == sha256(P / 'inputs/hall-motion-v2' / name), 'Hall provider differs: ' + name)
    seal = json.loads(S5_SEAL.read_text())
    changed = [f['path'] for f in seal['files'] if sha256(S5 / f['path']) != f['sha256']]
    need(not changed, 'source v5 sealed files changed: %s' % changed[:4])
    v7 = json.loads((V7 / 'SOURCE-SHA256.json').read_text())
    need(v7['source_identity'] == V7_IDENTITY, 'native v7 identity changed')
    for name, digest in v7['files'].items():
        need(sha256(V7 / name) == digest, 'native v7 file changed: ' + name)
    need(sha256(V7 / 'build/libpic_native_v7.a') == V7_LIB_SHA, 'native v7 archive changed')
    v71 = json.loads((V71 / 'SOURCE-SHA256.json').read_text())
    need(v71['source_identity'] == V71_IDENTITY, 'native v7.1 identity changed')
    for name, digest in v71['files'].items():
        need(sha256(V71 / name) == digest, 'native v7.1 file changed: ' + name)
    need(sha256(V71 / 'build/libpic_native_v71.a') == V71_LIB_SHA, 'native v7.1 archive changed')
    need(sha256(CHEM_BLOB) == CHEM_BLOB_SHA, 'chemistry v3 blob changed')
    need(sha256(CHEM31_BLOB) == CHEM31_BLOB_SHA, 'chemistry v3.1 blob changed')
    need(sha256(CHEM31_POLICY) == CHEM31_POLICY_SHA, 'chemistry v3.1 coverage policy changed')
    for name, digest in CHEM_HEADER_SHA.items():
        need(sha256(CHEM / 'include/pf_chem_v3' / name) == digest, 'chemistry v3.1.1 header changed: ' + name)
    for name, digest in MERGE_SPLIT_SHA.items():
        need(sha256(MERGE_SPLIT / name) == digest, 'merge/split component file changed: ' + name)
    return {'merge_split_v1': {'path': str(MERGE_SPLIT), 'files': MERGE_SPLIT_SHA},
            'native_v71': {'identity': V71_IDENTITY, 'files_verified': len(v71['files']), 'archive_sha256': V71_LIB_SHA},
            'chemistry_v3': {'blob_sha256': CHEM_BLOB_SHA, 'blob_path': str(CHEM_BLOB)},
            'chemistry_v3_1': {'blob_sha256': CHEM31_BLOB_SHA, 'blob_path': str(CHEM31_BLOB),
                               'policy_sha256': CHEM31_POLICY_SHA, 'policy_path': str(CHEM31_POLICY),
                               'headers': tree_hashes(CHEM / 'include'), 'kernel_sha256': sha256(CHEM / 'device/chem_kernels.cu')},
'native_v5': {'identity': V5_IDENTITY, 'files_verified': len(v5['files'])},
            'native_v7': {'identity': V7_IDENTITY, 'files_verified': len(v7['files']), 'archive_sha256': V7_LIB_SHA},
            'native_v6': {'identity': V6_IDENTITY, 'files_verified': len(v6['files']),
                          'archive_sha256': V6_LIB_SHA, 'adapter_sha256': V6_ADAPTER_SHA,
                          'v3_adapter_sha256': V6_V3_ADAPTER_SHA},
            'exact_count_header_sha256': COUNT_SHA,
            'hall_commit': head,
            'hall_corrective_v2': tree_hashes(HALL_V2 / 'source', {'.hpp', '.cpp'}),
            'hall_gpu_resident_v1_common': tree_hashes(GPU_COMMON),
            'nvrtc_library_sha256': sha256(NVRTC),
            'hall_v32_motion_work_gross_v1_hpp_sha256': sha256(HALL_V32_WORK_GROSS),
            'hall_v32_motion_work_gross_v1_spec_sha256': sha256(HALL_V32_SPEC),
            'arrival_owner_v1': tree_hashes(ARRIVAL / 'include'),
            'common_commit_v3_contract_sha256': sha256(V3_COMMIT / 'source/common_commit_contract.hpp'),
            'source_v5': {'identity': seal['source_identity'], 'files_verified': len(seal['files'])}}


class Runner:
    def __init__(self, out, label):
        self.out = Path(out)
        self.out.mkdir(parents=True, exist_ok=True)
        self.label = label
        self.records = []

    def run(self, command, log, timeout=1800, env=None, sanitized=False, limit=True, check=True, rlimit=RLIMIT_AS):
        prefix = [FROZEN]
        if limit and not sanitized:
            prefix += ['prlimit', f'--as={rlimit}', '--']
        full = prefix + [str(x) for x in command]
        environment = dict(os.environ)
        if sanitized:
            environment.update(ASAN_ENV)
        if env:
            environment.update(env)
        start = time.monotonic()
        with (self.out / log).open('w') as handle:
            p = subprocess.run(full, stdout=handle, stderr=subprocess.STDOUT, timeout=timeout,
                               env=environment, cwd=str(P))
        record = {'log': log, 'command': full, 'returncode': p.returncode,
                  'elapsed_s': round(time.monotonic() - start, 3),
                  'log_sha256': sha256(self.out / log),
                  'cpus': sorted(os.sched_getaffinity(0)), 'nice': os.getpriority(os.PRIO_PROCESS, 0),
                  'rlimit_as': rlimit if (limit and not sanitized) else None,
                  'sanitizer_rss_limit_mb': 3072 if sanitized else None, 'utc': utc()}
        self.records.append(record)
        (self.out / (self.label + '-commands.json')).write_text(json.dumps(self.records, indent=2) + '\n')
        if check and p.returncode:
            print((self.out / log).read_text()[-12000:])
            raise RuntimeError(f'failed {log} exit {p.returncode}')
        return p.returncode


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()
