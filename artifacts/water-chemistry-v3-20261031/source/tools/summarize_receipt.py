#!/usr/bin/env python3
"""Summarize a qualification receipt into RESULTS.txt (and copy SWARM-VALIDATION.txt) in the component directory.

usage: summarize_receipt.py OUTPUT_LABEL_DIR
"""
import json
import os
import shutil
import sys

R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main(d):
    rc = json.load(open(os.path.join(d, 'receipt.json')))
    L = ['Qualification results of water chemistry v3 (receipt %s)' % os.path.join(d, 'receipt.json'), '']
    for st in rc['stages']:
        L.append('[%s] %s  %s' % ('PASS' if st['passed'] else 'FAIL', st['stage'], st['utc']))
        for e in st['entries']:
            L.append('    %-44s exit %s%s' % (e['name'], e.get('exit'), ('  %.1f s' % e['seconds']) if 'seconds' in e else ''))
        if st['stage'] == 'tables':
            L.append('    blob sha256 %s' % st.get('blob_sha256'))
    for sub in ('core-normal', 'core-sanitized'):
        p = os.path.join(d, sub, 'core-results.json')
        if os.path.exists(p):
            j = json.load(open(p))
            L.append('')
            L.append('%s: %d failures' % (sub, j['failures']))
            for r in j['results']:
                if r['check'].startswith('release/stock never negative'):
                    continue
                L.append('    %s %-62s %s' % ('PASS' if r['pass'] else 'FAIL', r['check'], r['detail'][:150]))
            lp = os.path.join(d, sub, 'ledger-1e6.json')
            if os.path.exists(lp):
                lg = json.load(open(lp))
                L.append('    ledger: %d events; in %.6g + %.6g eV; out survivor %.6g, birth %.6g, reservoir %.6g, gas %.6g; q %.6g; internal %s; defect %.3g eV (exact hex %s...); max event defect %.3g relative'
                         % (lg['events'], lg['ke_in_projectile_eV'], lg['ke_in_target_eV'], lg['ke_out_survivor_eV'], lg['ke_out_birth_eV'],
                            lg['ke_out_reservoir_eV'], lg['ke_out_gas_eV'], lg['q_reference_eV'], ['%.4g' % x for x in lg['internal_eV']],
                            lg['defect_eV'], lg['defect_exact_hex'][-24:], lg['max_relative_event_defect']))
    for k in (1, 2):
        p = os.path.join(d, 'gpu-process%d' % k, 'gpu-parity.json')
        if os.path.exists(p):
            j = json.load(open(p))
            L.append('')
            dev = j.get('device', {})
            comp = j.get('compile', {})
            L.append('gpu process %d: %s cc %s.%s, arch option %s (%s), nvrtc %s' % (k, dev.get('name'), dev.get('cc_major'), dev.get('cc_minor'),
                                                                                   comp.get('arch_option'), comp.get('image_kind'), comp.get('nvrtc_version', comp.get('nvrtc_major'))))
            for c in j.get('cases', []):
                if c['case'] == 'math':
                    L.append('    math parity: %d values, %d mismatches' % (c['values'], c['mismatches']))
                else:
                    L.append('    step seed %s step %s block %s%s: %d events, carrier/event mismatches %d/%d, record digest cpu==gpu %s, ledger cpu==gpu %s, kernel %.2f ms, cpu %.1f ms'
                             % (c['seed'], c['step'], c['block'], ' split' if c['split_launch'] else '', c['events'], c['carrier_mismatches'], c['event_mismatches'],
                                c['record_digest_cpu'] == c['record_digest_gpu'], c['ledger_digest_cpu'] == c['ledger_digest_gpu'], c['gpu_kernel_ms'], c['cpu_ms']))
    sv = os.path.join(d, 'SWARM-VALIDATION.txt')
    if os.path.exists(sv):
        shutil.copy(sv, os.path.join(R, 'SWARM-VALIDATION.txt'))
        L.append('')
        L.append('Swarm validation: SWARM-VALIDATION.txt (copied from %s)' % sv)
    open(os.path.join(R, 'RESULTS.txt'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main(sys.argv[1])
