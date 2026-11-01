#!/usr/bin/python -I
"""Bounded, read-only SimQueue canary. No simulation, allocation or busy loop."""
import json, os, sys, time
from pathlib import Path

def require(value, message):
    if not value:
        raise RuntimeError(message)

def main():
    if sys.argv[1:] != ['--seconds', '75']:
        raise RuntimeError('only the bounded 75-second canary is supported')
    unit = os.environ['SIMQ_UNIT']
    group = Path('/sys/fs/cgroup') / ('user.slice/user-1000.slice/user@1000.service/simq.slice/' + unit)
    require(os.environ['SIMQ_CPU_LIST'] == '6', 'CPU6 allocation required; do not steal or force a core')
    require(os.environ['SIMQ_MEM'] == str(4 * 1024**3), '4 GiB allocation required')
    require(os.environ.get('SIMQ_BORROW', '0') == '0', 'borrowed allocation refused')
    require(os.sched_getaffinity(0) == {6}, 'CPU6 affinity required')
    require(Path('/proc/self/cgroup').read_text().strip() == '0::/' + str(group.relative_to('/sys/fs/cgroup')), 'cgroup differs')
    started = time.monotonic()
    print(json.dumps({'event': 'waiting', 'pid': os.getpid(), 'unit': unit}), flush=True)
    seen = []
    while time.monotonic() - started < 75:
        require(os.sched_getaffinity(0) == {6}, 'affinity changed')
        require((group/'memory.max').read_text().strip() == str(4 * 1024**3), 'memory cap changed')
        p = group/'cpuset.cpus.partition'
        value = p.read_text().strip() if p.exists() else 'controller_absent'
        if not seen or seen[-1] != value:
            seen.append(value)
            print(json.dumps({'event': 'state', 'elapsed_s': time.monotonic()-started, 'partition': value}), flush=True)
        time.sleep(.25)
    if 'isolated' not in seen or seen[-1] not in ('member', 'controller_absent'):
        raise RuntimeError('canary did not observe isolation followed by restoration')
    print(json.dumps({'event': 'canary_complete', 'states': seen}), flush=True)

if __name__ == '__main__':
    main()
