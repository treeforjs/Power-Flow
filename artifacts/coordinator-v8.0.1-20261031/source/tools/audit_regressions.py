#!/usr/bin/python -I
"""End-to-end continuation regressions; execute inside an allocated frozen-runtime job.

Uses short physical coupled runs, including checkpoint replay and a second restart
after extending a time-limited run. No process/rank-count claim is made.
"""
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = Path('/home/treefor/projects/.workspace/outputs/pic-water-common-commit-v8.0.1-20261031')


def main(label, tag, expected_restore_mode):
    if expected_restore_mode not in ('direct', 'replay'):
        raise ValueError('expected restore mode must be direct or replay')
    records = []

    def checkpoint_hashes(directory):
        return {str(p.relative_to(directory)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in sorted(directory.rglob('*')) if p.is_file()}

    def run(name, *args, success=True):
        dest = OUT / (tag + '-' + name)
        restore = Path(args[args.index('--restore-checkpoint') + 1]) if '--restore-checkpoint' in args else None
        checkpoint_before = checkpoint_hashes(restore) if restore else None
        command = [sys.executable, '-I', str(ROOT / 'run_v8.py'), 'run', label, 'run', str(dest),
                   '--fixture', 'qual', '--dt-mode', 'fixed', *map(str, args)]
        result = subprocess.run(command, capture_output=True, text=True)
        log = (dest / 'run.log').read_text() if (dest / 'run.log').exists() else result.stdout + result.stderr
        receipt_path = dest / 'run-receipt.json'
        receipt = json.loads(receipt_path.read_text()) if receipt_path.exists() else {}
        if not all(receipt.get(key) is True for key in
                   ('sealed_unchanged_after_run', 'own_sources_unchanged_after_run', 'binary_unchanged_after_run')):
            raise RuntimeError(name + ': run integrity was not verified: ' + log[-1000:])
        if success and result.returncode:
            raise RuntimeError(name + ': ' + log[-2500:] + result.stderr[-1000:])
        if restore and checkpoint_hashes(restore) != checkpoint_before:
            raise RuntimeError(name + ': input checkpoint changed during restore')
        if not success:
            for filename in ('step-receipts.jsonl', 'step-digests.txt', 'ledger.bin', 'deterministic.bin'):
                p = dest / filename
                if p.exists() and p.stat().st_size:
                    raise RuntimeError(name + ': refused continuation published accepted data: ' + filename)
            if any((dest / 'checkpoints').glob('ckpt-*')):
                raise RuntimeError(name + ': refused continuation published a checkpoint')
        summary = json.loads((dest / 'run-summary.json').read_text()) if (dest / 'run-summary.json').exists() else {}
        return dest, result.returncode, summary, log

    def check(name, condition, **evidence):
        records.append(dict(name=name, passed=bool(condition), **evidence))
        print(('PASS ' if condition else 'FAIL ') + name, flush=True)

    def close(a, b):
        return math.isclose(a, b, rel_tol=0, abs_tol=1e-20)

    base, _, bs, _ = run('window-base', '--window', 4e-14, '--checkpoint-every', 2)
    check('initial time-limited run accepts four steps', bs['stop'] == 'end_time' and bs['steps_completed'] == 4)
    ck = lambda directory, step: directory / 'checkpoints' / ('ckpt-%012d' % step)
    accepted = bs['final_time_s']
    direct, _, ds, _ = run('direct-window', '--restore-checkpoint', ck(base, 0), '--window', 2e-14)
    check('direct genesis restore honors relative window', ds['steps_this_process'] == 2 and close(ds['final_time_s'], accepted - 2e-14)
          and ds['restored_from']['mode'] == 'direct')
    _, _, direct_end, _ = run('direct-end', '--restore-checkpoint', ck(base, 0), '--end-time', accepted - 2e-14)
    check('direct genesis restore honors absolute end', direct_end['steps_this_process'] == 2 and
          close(direct_end['final_time_s'], accepted - 2e-14) and direct_end['restored_from']['mode'] == 'direct')
    _, _, direct_keep, _ = run('direct-keep-target', '--restore-checkpoint', ck(base, 0))
    check('direct genesis restore retains saved endpoint without override', direct_keep['steps_completed'] == 4 and
          direct_keep['steps_this_process'] == 4 and close(direct_keep['final_time_s'], accepted) and
          direct_keep['restored_from']['mode'] == 'direct')
    _, _, crossing, _ = run('off-grid-window', '--restore-checkpoint', ck(base, 0), '--window', 2.5e-14)
    check('off-grid target stops at first accepted step crossing the threshold', crossing['steps_this_process'] == 3 and
          close(crossing['final_time_s'], accepted - 1e-14) and crossing['stop'] == 'end_time' and
          crossing['restored_from']['mode'] == 'direct')
    extended, _, es, _ = run('extend-window', '--restore-checkpoint', ck(base, 4), '--window', 2e-14, '--checkpoint-every', 2)
    check('completed window extends after verified restore', es['steps_this_process'] == 2 and close(es['final_time_s'], accepted + 2e-14)
          and es['restored_from']['mode'] == expected_restore_mode,
          restored_mode=es['restored_from']['mode'])
    _, _, chained, _ = run('extend-chained', '--restore-checkpoint', ck(extended, 6), '--end-time', accepted + 4e-14)
    check('second restart honors absolute end after prior extension', chained['steps_this_process'] == 2 and
          close(chained['final_time_s'], accepted + 4e-14) and chained['restored_from']['mode'] == expected_restore_mode,
          restored_mode=chained['restored_from']['mode'])
    _, _, keep, _ = run('keep-target', '--restore-checkpoint', ck(base, 2))
    check('omitted target retains saved endpoint', keep['steps_completed'] == 4 and close(keep['final_time_s'], accepted)
          and keep['restored_from']['mode'] == expected_restore_mode)
    _, _, steps, _ = run('step-target', '--restore-checkpoint', ck(base, 4), '--steps', 6)
    check('explicit steps replaces saved time limit', steps['steps_completed'] == 6 and steps['stop'] == 'steps'
          and steps['restored_from']['mode'] == expected_restore_mode)
    step_base, _, _, _ = run('step-base', '--steps', 4, '--checkpoint-every', 2)
    _, _, saved_steps, _ = run('keep-step-target', '--restore-checkpoint', ck(step_base, 2))
    check('omitted target retains saved step limit', saved_steps['steps_completed'] == 4 and saved_steps['steps_this_process'] == 2
          and saved_steps['restored_from']['mode'] == expected_restore_mode)
    _, rc, _, log = run('past-target', '--restore-checkpoint', ck(base, 4), '--end-time', accepted, success=False)
    check('non-advancing explicit end is refused', rc != 0 and 'later than the accepted checkpoint time' in log)
    _, rc, _, log = run('conflicting-target', '--window', 2e-14, '--steps', 2, success=False)
    check('ambiguous targets are refused', rc != 0 and 'choose only one' in log)
    result = dict(schema='coordinator-v8.0.1-audit-regressions-v1', cases=records,
                  passed=sum(r['passed'] for r in records), total=len(records), executed_processes_per_run=1,
                  expected_nonempty_restore_mode=expected_restore_mode,
                  restored_inputs_unchanged=True, refusal_published_no_accepted_data=True,
                  time_target_semantics='first accepted step at or beyond threshold; exact landing not guaranteed')
    (OUT / (tag + '-results.json')).write_text(json.dumps(result, indent=2) + '\n')
    return 0 if all(r['passed'] for r in records) else 1


if __name__ == '__main__':
    raise SystemExit(main(*sys.argv[1:]))
