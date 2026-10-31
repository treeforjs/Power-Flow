"""Summary of a native v7.1 verification run (OUT/LABEL/verification-result.json and logs) as text.
usage: python3 scripts/results.py LABEL"""
import json
import sys
from pathlib import Path

OUT = Path('/home/treefor/projects/.workspace/outputs/native-v7.1-20261031')
label = sys.argv[1] if len(sys.argv) > 1 else 'final'
F = OUT / label
r = json.loads((F / 'verification-result.json').read_text())
art = json.loads((F / 'build-artifact-hashes.json').read_text())

out = []
w = out.append
w('Native v7.1 (candidate007.1) results -- verification run %r' % label)
w('Receipts: %s (receipt-commands.json, verification-result.json, build-artifact-hashes.json, logs)' % F)
w('Source identity %s; base native v7 %s (committed, unchanged: hashed before/after).' % (r['source_identity'], r['base_v7_identity']))
w('Core archive build/libpic_native_v71.a sha256 %s.' % art['build/libpic_native_v71.a'])
w('')
w('Suites (last line of each log)')
for k, v in sorted(r['suite_summaries'].items()):
    if not k.startswith(('timing-', 'gpu-forms', 'verify-exact')):
        w(' %-34s %s' % (k, v))
w('')
w('Host mode vs native v7 evidence: %s' % ('identical (suite outputs; checkpoints except source_sha256, anchors, file digest)'
                                          if not any(r['host_vs_v7_differences'].values()) else r['host_vs_v7_differences']))
w('Device modes vs host mode: %s' % ('identical (suite outputs; every checkpoint byte for byte)'
                                     if not any(r['device_vs_host_differences'].values()) else r['device_vs_host_differences']))
w('Parity trace (60 steps, forced events, K = 4): v7 vs v7.1 host %s; device_verified vs host %s; device vs host %s' % (
    r['parity']['v7_vs_host'] or 'identical except anchors/source identity', r['parity']['device_verified_vs_host'] or 'identical',
    r['parity']['device_vs_host'] or 'identical'))
w('')
w('Execution counters (device modes)')
for k, v in sorted(r['executions'].items()):
    w(' %-34s %s' % (k, json.dumps({x: y for x, y in v.items() if x != 'execution'}, sort_keys=True)))
w('')
w('GPU forms (RTX 3080 Ti, NVRTC, two fresh processes, identical non-timing output)')
for line in r['gpu_forms']:
    if any(t in line for t in ['gpu_fused_block_lists', 'gpu_v5_id_roots', 'gpu_field_a1', 'gpu_checks']):
        w(' ' + line.replace('RESULT ', ''))
w('')
w('Native per-step host cost (thread CPU, median of the steps after the first; ms/step [us/carrier])')
w(' mesh    carriers   v7 (before)          v7.1 host            v7.1 device')
t = r['timings']
for mesh, name in [(45, '200um'), (180, '50um')]:
    for n in [10000, 100000, 1000000]:
        cells = []
        for tag in ['v7-host', 'v71-host', 'v71-device']:
            x = t.get('%s-%d-%d' % (tag, mesh, n))
            cells.append('%9.3f [%.5f]' % (x['native_cpu_ms_per_step_median'], x['native_cpu_us_per_carrier']) if x else '-')
        w(' %-6s %9d   %-20s %-20s %s' % (name, n, *cells))
for mesh, name in [(45, '200um'), (180, '50um')]:
    x = t.get('v71-device_stream-%d-10000' % mesh)
    if x:
        w(' %-6s %9d   device lane with the v7 streaming field/current digests (hashed on the host): %.3f ms/step' % (
            name, 10000, x['native_cpu_ms_per_step_median']))
w('')
w('Phases (ms/step, thread CPU): begin / first exposure + fork / census / admission + accept / second exposure / seal + publish')
order = ['begin', 'first_exposure_fork', 'census', 'admission_accept', 'second_exposure', 'seal_publish']
for key in sorted(t):
    x = t[key]
    w(' %-22s %s%s' % (key, ' / '.join('%.3f' % x['phases'][p]['cpu_ms_per_step'] for p in order),
                       ('  (device admissions %s, host materialisations %s)' % (x.get('device_admissions'), x.get('host_row_materialisations'))
                        if 'device' in key else '')))
w('')
w('ASan/UBSan (host mode)')
for k, v in sorted(r['sanitized'].items()):
    w(' %-22s %s' % (k, v))
lr = r['long_run_200000']
w('')
w('Long run 200000 steps (births and sinks every step): %s' % json.dumps({k: v for k, v in lr.items() if k != 'samples'}, sort_keys=True))
print('\n'.join(out))
