#!/usr/bin/env python3
"""Run single-threaded CPU jobs on leased cores (PFC_LEASE, default Power-Flow cores 6-7 and the lent cores 0, 1; 15 during a loan).

usage: run_jobs.py JOBS.json RECEIPT.json
JOBS.json: [{"name": ..., "argv": [...], "log": PATH, "cwd": DIR}, ...]
Each job runs under the frozen launcher, pinned to one core (taskset), nice 15, one thread, no GPU,
holding that core's lease file (flock) for its whole duration. No new job starts on cores 0-3 once
RETURN-CORES-0-3 exists. Cores 2-5 are never used.
"""
import json
import fcntl
import importlib.util
from pathlib import Path
import tempfile
import os
import signal
import subprocess
import sys
import time

LOCKS = '/home/treefor/projects/.workspace/outputs/pic-profile'
FROZEN = '/home/treefor/projects/runtimes/arch-base-20261028/run-frozen'
LOCK_OF = {0: 'testing-core0.lock', 1: 'testing-core1.lock', 6: 'testing-core.lock', 7: 'testing-core7.lock', 15: 'testing-cpu15.lock'}
NICE_OF = {15: '19'}   # logical CPU 15 during a core loan (shared by the Power-Flow agents)
QUEUED = os.environ.get('PFC_LEASE') == 'none'   # inside a SimQueue job: the queue holds the job's CPUs; no lease files
CORES = sorted(os.sched_getaffinity(0)) if QUEUED else [int(c) for c in os.environ.get('PFC_LEASE', '6 7 0 1').split()]
RETURN = os.path.join(LOCKS, 'RETURN-CORES-0-3')
_spec = importlib.util.spec_from_file_location('swarm_job_integrity', Path(__file__).with_name('swarm_job_integrity.py'))
_proof = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_proof)
BINDING = os.environ.get('PF_QUAL_BINDING_SHA256', 'unbound-diagnostic')


def main(jobs_path, receipt_path):
    jobs = json.load(open(jobs_path))
    pending = list(range(len(jobs)))
    results = {}
    running = {}   # core -> (job index, process, start, claim, temporary output, input identity)
    def terminated(signum, frame):
        raise InterruptedError('swarm manager received signal ' + str(signum))
    previous = {sig: signal.signal(sig, terminated) for sig in (signal.SIGTERM, signal.SIGINT)}
    try:
        while pending or running:
            for core, (j, p, t0, claim, temporary, before) in list(running.items()):
                if p.poll() is not None:
                    code = p.returncode
                    if code == 0:
                        try:
                            _proof.publish(temporary, jobs[j], before, BINDING)
                        except (OSError, ValueError, KeyError, TypeError) as exc:
                            print(jobs[j]['name'], 'completion refused:', str(exc), flush=True)
                            code = 1
                    if code == 75:   # core lease busy: requeue
                        pending.insert(0, j)
                    else:
                        results[j] = {'name': jobs[j]['name'], 'core': core, 'exit': code,
                                      'seconds': time.time() - t0, 'completion_verified': code == 0}
                    Path(temporary).unlink(missing_ok=True)
                    claim.close()
                    del running[core]
            for core in CORES:
                if not pending:
                    break
                if core in running:
                    continue
                if core < 4 and os.path.exists(RETURN) and not QUEUED:
                    continue
                j = pending.pop(0)
                job = jobs[j]
                _proof.V.need(job.get('output'), 'job lacks explicit output')
                claim = open(job['output'] + '.claim', 'a')
                try:
                    fcntl.flock(claim, fcntl.LOCK_EX | fcntl.LOCK_NB)
                except BlockingIOError:
                    claim.close()
                    pending.append(j)
                    continue
                # Recheck under the claim: another runner may have completed it.
                if _proof.completed(job, BINDING):
                    results[j] = {'name': job['name'], 'core': None, 'exit': 0, 'seconds': 0.0,
                                  'resumed': True, 'completion_verified': True}
                    claim.close()
                    continue
                before = _proof.identity(job, BINDING)
                fd, temporary = tempfile.mkstemp(prefix='.' + Path(job['output']).name + '.',
                                                 suffix='.partial', dir=Path(job['output']).parent)
                os.close(fd)
                argv = list(job['argv'])
                argv[2] = temporary
                cmd = ([] if QUEUED else ['flock', '-n', '-E', '75', os.path.join(LOCKS, LOCK_OF[core])]) + \
                    ['taskset', '-c', str(core), 'nice', '-n', NICE_OF.get(core, '15'), 'env', 'OMP_NUM_THREADS=1', 'CUDA_VISIBLE_DEVICES=', FROZEN] + argv
                log = open(job['log'], 'w')
                p = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, cwd=job.get('cwd'), start_new_session=True)
                log.close()
                running[core] = (j, p, time.time(), claim, temporary, before)
            time.sleep(1.0)
        _proof.V.atomic_json(receipt_path, {'jobs': [results[i] for i in sorted(results)]})
        bad = [r for r in results.values() if r['exit'] != 0]
        print('jobs', len(results), 'failed', len(bad))
        return 1 if bad else 0
    finally:
        # Only process groups created by this manager are signaled.
        for _, child, _, _, _, _ in running.values():
            if child.poll() is None:
                try: os.killpg(child.pid, signal.SIGTERM)
                except ProcessLookupError: pass
        for _, child, _, claim, temporary, _ in running.values():
            try:
                child.wait(timeout=5)
            except subprocess.TimeoutExpired:
                try: os.killpg(child.pid, signal.SIGKILL)
                except ProcessLookupError: pass
                child.wait(timeout=5)
            Path(temporary).unlink(missing_ok=True)
            claim.close()
        for sig, handler in previous.items(): signal.signal(sig, handler)



if __name__ == '__main__':
    sys.exit(main(*sys.argv[1:3]))
