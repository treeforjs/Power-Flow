"""Summary of a native v7 verification run (OUT/LABEL/verification-result.json and logs) as text.
usage: python3 scripts/results.py LABEL"""
import json
import re
import statistics
import sys
from pathlib import Path

OUT = Path('/home/treefor/projects/.workspace/outputs/native-v7-20261031')
label = sys.argv[1] if len(sys.argv) > 1 else 'final'
F = OUT / label
r = json.loads((F / 'verification-result.json').read_text())
P = Path(__file__).resolve().parent.parent
art = json.loads((F / 'build-artifact-hashes.json').read_text())


def field(line, key):
    m = re.search(r'\b%s=([^ ]+)' % re.escape(key), line)
    return m.group(1) if m else ''


out = []
w = out.append
w('Native v7 (candidate007) results -- verification run %r' % label)
w('Receipts: %s (receipt-commands.json, verification-result.json, build-artifact-hashes.json, logs)' % F)
w('Source identity %s; base native v6 %s (sealed, unchanged: hashed before/after).' % (r['source_identity'], r['base_v6_identity']))
w('Core archive build/libpic_native_v7.a sha256 %s.' % art['build/libpic_native_v7.a'])
w('')
w('v6 parity')
w(' legacy: %d checks; arena-tests.json identical to native v6 evidence and to ASan/UBSan; other outputs differ only in' % r['legacy_checks'])
w('   %s.' % sorted({k for v in r['legacy_output_differences_vs_v6'].values() for k in v}))
pc = r['ported_candidate003_controls']
w(' ported candidate003 controls: %d checks, names identical to v6, allocation-failure sweep %s (identical to v6);' % (pc['check_count'], pc['allocation_failure_sweep']))
w('   genesis checkpoint differs from v6 only in %s.' % r['genesis_checkpoint_differences_vs_v6'])
h = r['v5_history_controls']
w(' v5 history controls: %d normal, %d ASan/UBSan reduced; v6 checks replaced by the v7 motion form: %s;' % (h['normal'], h['sanitized_reduced'], h['replaced_v6_checks']))
w('   new: %s.' % h['new_v7_checks'])
f = r['v6_flows_ported']
fn = r['v6_flows_name_changes']
w(' v6 flows controls (ported): %d normal, %d ASan/UBSan reduced; %d v6 check names renamed/replaced by %d v7 names' % (
    f['normal'], f['sanitized_reduced'], len(fn['removed_v6_names']), len(fn['new_v7_names'])))
w('   (v3 record/JSON route checks -> v7 record checks; every v6 control kept in v7 form).')
s = r['v6_source_arrival']
w(' source v5 + arrival v1 sweep: %d checks, %d cases, %d admitted steps, %d exact stream equalities, %d births (= v6);' % (
    s['check_count'], s['cases'], s['admitted_steps'], s['exact_stream_equalities'], s['births']))
ss = r['v6_source_arrival_sanitized_reduced']
w('   ASan/UBSan reduced %d checks, %d equalities, %d births.' % (ss['check_count'], ss['exact_stream_equalities'], ss['births']))
p = r['physical_state_parity_trace']
w(' physical-state parity trace (sealed v6 JSON route vs v7 adapter, 60 steps, forced events, K = 4): physical state')
w('   (every native carrier with packet ID and state) identical at every step; genesis identical; differing fields %s.' % p['differing_identity_fields'])
lr = r['long_run_200000']
w(' long run 200000 steps (births and sinks every step): %s' % json.dumps({k: v for k, v in lr.items() if k != 'samples'}, sort_keys=True))
w('')
w('Exact sums (API-007.txt 1.1)')
for line in r['forms_cpu']:
    if 'adversarial' in line or 'partition' in line or 'kinetic n=' in line or 'field_work' in line:
        w(' ' + line.replace('RESULT ', ''))
for line in r['exact_rational_verification']:
    w(' ' + line.replace('RESULT ', ''))
w('')
w('Native forms vs sealed functions: %s' % json.dumps(r['forms_native']['measurements'], sort_keys=True))
w('')
w('motion-work-gross-v1 (amendment A1): %d checks (ASan/UBSan %d)' % (r['work_gross_v1']['check_count'], r['work_gross_v1']['sanitized']))
for c in r['work_gross_v1']['cases']:
    if 'control' not in c:
        w(' %s: sealed provider refused (%s)' % (c['case'], c.get('sealed_provider_refusal', '')))
        continue
    sc = c['control']['native_scales']
    w(' %s: admitted, tolerance %.4g J (G_K %.4g, C_K %.4g, G_W %.4g, Q %.4g, E^ %.4g), defect/tol %.3g, defect/v6-net-tol %.3g,' % (
        c['case'], sc['tolerance_J'], sc['G_K_J'], sc['C_K_J'], sc['G_W_J'], sc['Q_C'], sc['E_hat_V_m'],
        c['defect_over_gross_tolerance'], c['defect_over_v6_net_tolerance']))
    w('   bits equal to Hall reference: %s; tampers (margin |defect|/tol, refused):' % c['control']['tolerance_bits_equal_hall'])
    for t in c['tampers']:
        w('     %-58s %9.3g  %s  bits-equal %s' % (t['case'], t['defect_over_tolerance'], 'refused' if t['refusal'] else 'ADMITTED',
                                                  t['tolerance_bits_equal_hall']))
w('')
w('CPU vs GPU (RTX 3080 Ti, NVRTC, two fresh processes, identical non-timing output)')
for line in r['gpu']:
    if any(k in line for k in ['gpu_dot random', 'gpu_field_work', 'gpu_kinetic', 'gpu_block_list', 'gpu_census', 'resident_host', 'gpu_checks', 'kernel_source']):
        w(' ' + line.replace('RESULT ', ''))
w('')
w('Benchmark (median, ms; us/carrier in brackets)')
b = r['benchmark']
for n in [10000, 100000, 1000000]:
    for tag in ['v6', 'v7']:
        med = b['%s-%d' % (tag, n)]['median_s']
        w(' %s n=%d: %s' % (tag, n, ', '.join('%s %.3f [%.4f]' % (k[:-2], 1e3 * v, 1e6 * v / n) for k, v in sorted(med.items()))))
print('\n'.join(out))
