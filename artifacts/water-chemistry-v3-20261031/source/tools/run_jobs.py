#!/usr/bin/env python3
"""Run single-threaded CPU jobs on leased cores (Power-Flow cores 6-7 and the lent cores 0, 1, 3).

usage: run_jobs.py JOBS.json RECEIPT.json
JOBS.json: [{"name": ..., "argv": [...], "log": PATH, "cwd": DIR}, ...]
Each job runs under the frozen launcher, pinned to one core (taskset), nice 15, one thread, no GPU,
holding that core's lease file (flock) for its whole duration. No new job starts on cores 0-3 once
RETURN-CORES-0-3 exists. Cores 2, 4 and 5 are never used.
"""
import json
import os
import subprocess
import sys
import time

LOCKS = '/home/treefor/projects/.workspace/outputs/pic-profile'
FROZEN = '/home/treefor/projects/runtimes/arch-base-20261028/run-frozen'
LOCK_OF = {0: 'testing-core0.lock', 1: 'testing-core1.lock', 3: 'testing-core3.lock', 6: 'testing-core.lock', 7: 'testing-core7.lock'}
RETURN = os.path.join(LOCKS, 'RETURN-CORES-0-3')


def main(jobs_path, receipt_path):
    jobs = json.load(open(jobs_path))
    pending = list(range(len(jobs)))
    running = {}   # core -> (job index, Popen, start)
    results = {}
    while pending or running:
        for core, (j, p, t0) in list(running.items()):
            if p.poll() is not None:
                if p.returncode == 75:   # lease busy: requeue
                    pending.insert(0, j)
                else:
                    results[j] = {'name': jobs[j]['name'], 'core': core, 'exit': p.returncode, 'seconds': time.time() - t0}
                del running[core]
        for core in (6, 7, 0, 1, 3):
            if not pending:
                break
            if core in running:
                continue
            if core < 4 and os.path.exists(RETURN):
                continue
            j = pending.pop(0)
            job = jobs[j]
            cmd = ['flock', '-n', '-E', '75', os.path.join(LOCKS, LOCK_OF[core]), 'taskset', '-c', str(core), 'nice', '-n', '15',
                   'env', 'OMP_NUM_THREADS=1', 'CUDA_VISIBLE_DEVICES=', FROZEN] + job['argv']
            log = open(job['log'], 'w')
            p = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, cwd=job.get('cwd'))
            running[core] = (j, p, time.time())
        time.sleep(1.0)
    json.dump({'jobs': [results[i] for i in sorted(results)]}, open(receipt_path, 'w'), indent=1)
    bad = [r for r in results.values() if r['exit'] != 0]
    print('jobs', len(results), 'failed', len(bad))
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main(*sys.argv[1:3]))
