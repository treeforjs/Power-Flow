#!/usr/bin/python -I
"""Coordinator v8 restart tests (A05/A06): whole-run checkpoint, restore into a fresh process, crash during a checkpoint
write, tamper and corruption refusals, rank-metadata independence. Every driver run goes through run_v8.py (receipted).

  restart_tests.py LABEL TAG [--resident] [--lane CMD...]
      LABEL  run_v8.py build label;  TAG  prefix of the output directories (OUT/TAG-*)
      --lane the command prefix that holds the core (and GPU) leases, e.g.
             flock .../testing-core1.lock taskset -c 1 nice -n 15 env OMP_NUM_THREADS=1 ... run-frozen /usr/bin/python -I
      --resident: the scenarios on the device-resident Hall owner (the lane must hold the GPU locks)

Scenarios (each: an uninterrupted reference run with checkpoints, then fresh-process restores and comparisons):
  qual    qualification fixture (births every step, sinks, ionization), fixed 10 fs, 30 steps, checkpoints every 10
  dense   plasma blob 8e10 (omega_p dt ~0.03-0.05), adaptive controller, 30 steps, checkpoints every 10 (mid-plasma)
  pilot   pilot fixture, adaptive controller with every owner cap at 1 ps and the production source interval preset
          (LTE history, source rate samples and controller state in the checkpoint), 60 steps, checkpoints every 20
Restores: every checkpoint of each scenario (step 0: native restores directly; later steps: native v7.1 refuses a nonempty
disk restore, so the verified replay from the genesis checkpoint runs) continued to the reference end; the continuation
must be bit-identical (deterministic ledger, step digests, step receipts, control ledger without timings, probes).
Crash: the qual reference rerun with a crash (exit 86) injected while writing the checkpoint of step 20 at each stage
(0: owner files written, 1: before the manifest, 2: manifest temporary written, 3: committed files before the directory
rename); recovery (--resume) must pick checkpoint 10 and continue bit-identically.
Tamper: copies of the qual checkpoint root with: truncated native.json, one flipped byte in hall/checkpoint.json, a
missing owner file, an owner removed from the manifest (resealed), a mismatched native owner version (resealed), one
flipped byte in the manifest body (not resealed), a flipped byte in arrival.bin with the manifest digest updated
(resealed: the owner's own checksum refuses), a flipped byte in the genesis checkpoint (replay base). All must be refused.
Rank metadata: --ranks 3 changes metadata only; every run executes one process.
This does not qualify a distributed MPI restart or scaling.
"""
import hashlib, json, os, shutil, subprocess, sys, time
from pathlib import Path

P = Path(__file__).resolve().parent.parent
OUT = Path('/home/treefor/projects/.workspace/outputs/pic-water-common-commit-v8.0.1-20261031')
sys.path.insert(0, str(P / 'tools'))
import compare_v8  # noqa: E402

CAPS = ['--hall-dt-max', '1e-12', '--source-dt-max', '1e-12', '--neutral-dt-max', '1e-12', '--charged-only-particle-cap', '--dt-cap', '1e-12']
SCENARIOS = {
    'qual': (['--fixture', 'qual', '--steps', '30', '--checkpoint-every', '10'], 30, [0, 10, 20]),
    'dense': (['--fixture', 'dense', '--blob-weight', '8e10', '--dt-mode', 'adaptive', '--steps', '30', '--checkpoint-every', '10'], 30, [10, 20]),
    'pilot': (['--fixture', 'pilot', '--dt-mode', 'adaptive', *CAPS, '--source-interval-preset', 'production', '--steps', '60', '--checkpoint-every', '20'],
              60, [20, 40]),
}


def main():
    args = sys.argv[1:]
    label, tag = args[0], args[1]
    resident = '--resident' in args
    lane = args[args.index('--lane') + 1:] if '--lane' in args else []
    extra = ['--hall-owner', 'resident'] if resident else []
    results = {'label': label, 'tag': tag, 'resident': resident, 'cases': []}

    def run(name, *driver_args, expect_ok=True):
        out = OUT / ('%s-%s' % (tag, name))
        cmd = [*(lane or [sys.executable, '-I']), str(P / 'run_v8.py'), 'run', label, 'run', str(out), *extra, *driver_args]
        t = time.time()
        r = subprocess.run(cmd, cwd=P, capture_output=True, text=True)
        log = (out / 'run.log').read_text() if (out / 'run.log').exists() else r.stdout + r.stderr
        rc = r.returncode
        if (out / 'run-receipt.json').exists():   # the driver's own exit code (run_v8.py maps every failure to 1)
            receipt = json.loads((out / 'run-receipt.json').read_text())
            if not all(receipt.get(key) is True for key in
                       ('sealed_unchanged_after_run', 'own_sources_unchanged_after_run', 'binary_unchanged_after_run')):
                raise RuntimeError(name + ': run integrity was not verified')
            rc = receipt.get('command', {}).get('returncode', rc)
        return out, rc, log, time.time() - t

    def record(case, ok, **ev):
        results['cases'].append({'case': case, 'passed': bool(ok), **ev})
        print(('pass ' if ok else 'FAIL ') + case, json.dumps(ev)[:400], flush=True)

    refs = {}
    for sc, (sargs, end, points) in SCENARIOS.items():
        ref, rc, log, dt = run(sc + '-ref', *sargs)
        refs[sc] = ref
        summ = json.loads((ref / 'run-summary.json').read_text()) if (ref / 'run-summary.json').exists() else {}
        record(sc + ': reference run with checkpoints', rc == 0 and summ.get('stop') == 'steps', rc=rc, wall_s=round(dt, 1),
               checkpoints=[c['step'] for c in summ.get('checkpoints', [])], totals=summ.get('totals'))
        for k in points:
            ck = ref / 'checkpoints' / ('ckpt-%012d' % k)
            out, rc, log, dt = run('%s-restore-%d' % (sc, k), '--restore-checkpoint', str(ck), '--steps', str(end))
            rest = json.loads((out / 'restore.json').read_text()) if (out / 'restore.json').exists() else {}
            cmp = compare_v8.identity(ref, out, None, k + 1) if rc == 0 else {'bit_identical': False}
            record('%s: restore at step %d in a fresh process continues bit-identically to step %d' % (sc, k, end), rc == 0 and cmp['bit_identical'],
                   rc=rc, mode=rest.get('mode'), replayed=rest.get('replayed_steps'), cross_checks=rest.get('cross_checks'),
                   restore_s=rest.get('restore_s'), rows=cmp.get('ledger', {}).get('rows_compared'), tail=log[-300:] if rc else '')
    # crash during a checkpoint write
    qargs = SCENARIOS['qual'][0]
    for stage in (0, 1, 2, 3):
        out, rc, log, dt = run('qual-crash-s%d' % stage, *qargs, '--crash-in-checkpoint', '20:%d' % stage)
        root = out / 'checkpoints'
        tmp = sorted(p.name for p in root.iterdir() if p.name.endswith('.tmp')) if root.exists() else []
        done = sorted(p.name for p in root.iterdir() if not p.name.endswith('.tmp')) if root.exists() else []
        res, rc2, log2, dt2 = run('qual-crash-s%d-resume' % stage, '--resume', str(root), '--steps', '30')
        rest = json.loads((res / 'restore.json').read_text()) if (res / 'restore.json').exists() else {}
        cmp = compare_v8.identity(refs['qual'], res, None, 11) if rc2 == 0 else {'bit_identical': False}
        record('crash at stage %d of the step-20 checkpoint write: recovery takes checkpoint 10 and continues bit-identically' % stage,
               rc == 86 and rc2 == 0 and rest.get('step') == 10 and cmp['bit_identical'] and tmp == ['ckpt-%012d.tmp' % 20],
               crash_rc=rc, temporary=tmp, committed=done, resumed_from_step=rest.get('step'), skipped=rest.get('skipped'), mode=rest.get('mode'))
    # tamper and corruption
    src = refs['qual'] / 'checkpoints'

    def tampered(name, edit):
        d = OUT / ('%s-tamper-%s' % (tag, name))
        if d.exists():
            shutil.rmtree(d)
        shutil.copytree(src, d)
        edit(d)
        return d

    def flip(p, at=None):
        b = bytearray(p.read_bytes())
        i = len(b) // 2 if at is None else at
        b[i] ^= 0x01
        p.write_bytes(bytes(b))

    def reseal(ck, change):
        m = json.loads((ck / 'manifest.json').read_text())
        change(m['body'])
        m['body_sha256'] = hashlib.sha256(json.dumps(m['body'], separators=(',', ':'), ensure_ascii=False, sort_keys=True).encode()).hexdigest()
        (ck / 'manifest.json').write_text(json.dumps(m, indent=1) + '\n')

    def set_file_digest(ck, rel):
        def ch(body):
            data = (ck / rel).read_bytes()
            body['files'][rel] = {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
        reseal(ck, ch)

    c10 = 'ckpt-%012d' % 10
    cases = {
        'truncated-native': lambda d: (d / c10 / 'native.json').write_bytes((d / c10 / 'native.json').read_bytes()[:1000]),
        'flipped-hall-byte': lambda d: flip(d / c10 / 'hall' / 'checkpoint.json'),
        'missing-owner-file': lambda d: (d / c10 / 'source-3.bin').unlink(),
        'owner-removed-from-manifest': lambda d: reseal(d / c10, lambda b: b['owners'].pop('arrival')),
        'mismatched-native-version': lambda d: reseal(d / c10, lambda b: b['owners']['native'].__setitem__('version', '3f88fab5e285a13448ab709394dabc3cb65aab4820c76dec27c8121c0b042dc5')),
        'mismatched-hall-binding': lambda d: reseal(d / c10, lambda b: b['owners']['hall']['version'].__setitem__('binding', 'v3.2')),
        'flipped-manifest-byte': lambda d: flip(d / c10 / 'manifest.json', 200),
        'flipped-arrival-byte-resealed': lambda d: (flip(d / c10 / 'arrival.bin'), set_file_digest(d / c10, 'arrival.bin')),
        'flipped-coordinator-state-resealed': lambda d: (flip(d / c10 / 'controller-ring.bin'), set_file_digest(d / c10, 'controller-ring.bin')),
        'flipped-genesis-checkpoint-byte': lambda d: flip(d / ('ckpt-%012d' % 0) / 'chemistry.json'),
    }
    for name, edit in cases.items():
        d = tampered(name, edit)
        out, rc, log, dt = run('tamper-%s-restore' % name, '--restore-checkpoint', str(d / c10), '--steps', '12')
        fatal = [l for l in log.splitlines() if 'fatal' in l or 'refus' in l.lower() or 'differs' in l or 'mismatch' in l]
        record('tamper %s: restore refused' % name, rc != 0 and not (out / 'run-summary.json').exists(), rc=rc, refusal=(fatal[-1] if fatal else log[-300:]))
    # Rank metadata only; this command does not launch MPI ranks.
    out, rc, log, dt = run('qual-restore-10-ranks3', '--restore-checkpoint', str(src / c10), '--steps', '30', '--ranks', '3')
    rest = json.loads((out / 'restore.json').read_text()) if (out / 'restore.json').exists() else {}
    cmp = compare_v8.identity(refs['qual'], out, None, 11) if rc == 0 else {'bit_identical': False}
    record('one-process restore with rank metadata changed from 1 to 3 continues bit-identically', rc == 0 and cmp['bit_identical'],
           ranks_recorded=rest.get('ranks_recorded'), ranks_now=rest.get('ranks_now'), mpi_ranks_executed=1,
           distributed_restart_qualified=False)
    results['passed'] = sum(c['passed'] for c in results['cases'])
    results['total'] = len(results['cases'])
    (OUT / ('%s-restart-results.json' % tag)).write_text(json.dumps(results, indent=1) + '\n')
    print('%d/%d passed' % (results['passed'], results['total']))
    return 0 if results['passed'] == results['total'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
