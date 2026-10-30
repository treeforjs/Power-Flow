"""Release the manually paused PIC partition without terminating its processes."""
from pathlib import Path
import hashlib
import json
import os
import time

W = Path('/home/treefor/sim-workspace/powerflow')
G = Path('/sys/fs/cgroup/powerflow-pic-gpu-reserved')

def cpus(value):
    result = set()
    for part in value.strip().split(','):
        if part:
            a, *b = part.split('-')
            result.update(range(int(a), int(b[0])+1) if b else [int(a)])
    return result

def write_record(path, value):
    tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(value, indent=2) + '\n')
    os.chown(tmp, 1000, 1000)
    tmp.replace(path)

def main():
    if os.geteuid() != 0:
        raise RuntimeError('Administrative authentication required for the root-owned cgroup')
    pause = json.loads((W/'outputs/full60-manual-pause.json').read_text())
    manifest = json.loads((W/'outputs/full60-gpu-launch.json').read_text())
    if pause['status'] != 'paused_by_explicit_request':
        raise RuntimeError('No explicit full-run pause is in force')
    identity = pause['native_identity']; proc = Path('/proc') / str(identity['pid'])
    stat = (proc/'stat').read_text().split()
    if os.readlink(proc/'exe') != manifest['executable'] or stat[21] != identity['start_ticks']:
        raise RuntimeError('Native process identity changed')
    args = [v.decode() for v in (proc/'cmdline').read_bytes().split(b'\0') if v]
    if args != identity['args'] or stat[2] not in ('T', 't'):
        raise RuntimeError('Expected the same stopped full simulation')
    if (proc/'cgroup').read_text().strip() != '0::/powerflow-pic-gpu-reserved':
        raise RuntimeError('Native process is outside the owned reservation')
    before = {n:(G/n).read_text().strip() for n in
              ['cpuset.cpus.partition','cpuset.cpus.exclusive','cpuset.cpus.effective','cgroup.freeze']}
    members = [int(v) for v in (G/'cgroup.procs').read_text().split()]
    if identity['pid'] not in members:
        raise RuntimeError('Native process is not a reservation member')
    # Freeze the guard as well as the already-SIGSTOP native program before
    # changing the partition, so no guard can terminate a deliberately paused job.
    (G/'cgroup.freeze').write_text('1\n')
    for unused in range(100):
        events = dict(line.split() for line in (G/'cgroup.events').read_text().splitlines())
        if events.get('frozen') == '1':
            break
        time.sleep(.05)
    else:
        raise RuntimeError('Cgroup did not reach the frozen state; partition unchanged')
    (G/'cpuset.cpus.partition').write_text('member\n')
    (G/'cpuset.cpus.exclusive').write_text('\n')
    available = cpus((Path('/sys/fs/cgroup/user.slice')/'cpuset.cpus.effective').read_text())
    if not {6,7} <= available:
        raise RuntimeError('Testing CPUs are not available; all production processes remain frozen')
    if (proc/'stat').read_text().split()[21] != identity['start_ticks']:
        raise RuntimeError('Native identity changed after transition')
    result = dict(status='paused_job_frozen_core7_released_for_testing', updated_unix_s=time.time(),
                  production_identity=identity, native_stopped=True, production_cgroup_frozen=True,
                  cgroup=str(G), partition_before=before, partition_after=(G/'cpuset.cpus.partition').read_text().strip(),
                  preserved_members=members, testing_logical_cpus=[6,7], gpu_benchmarks_serial=True,
                  cpu_builds='CPU6 with testing-core.lock; GPU benchmark CPU7 with testing-gpu.lock',
                  full_run_resume_policy='New explicit restart decision required. Re-establish and verify isolated7,15 before unfreezing guards and before SIGCONT of the native process.',
                  helper_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    write_record(W/'outputs/pic-testing-allocation.json', result)
    pause['core_reservation']='All job processes frozen; partition is member so CPU7 is available for testing. Full-run isolation must be restored before any future resume.'
    pause['testing_allocation_receipt']=str(W/'outputs/pic-testing-allocation.json')
    write_record(W/'outputs/full60-manual-pause.json', pause)
    print(json.dumps({'status':result['status'],'testing_cpus':[6,7],'production_native_pid':identity['pid']}))

if __name__ == '__main__':
    main()
