#!/usr/bin/env python3
"""Summarize a qualification receipt into RESULTS.txt (and copy SWARM-VALIDATION.txt) in the component directory.

usage: summarize_receipt.py OUTPUT_LABEL_DIR [ANALYSIS_DIR]
v3.2 additions: swarm re-validation (SWARM-REVALIDATION.txt copied), cost per collision (cost stage and, if present,
the CPU benchmark bench-cost.txt), and the ion-beam / breakdown-mix / validation analyses of ANALYSIS_DIR.
"""
import json
import os
import shutil
import sys

R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main(d, analysis=None):
    rc = json.load(open(os.path.join(d, 'receipt.json')))
    L = ['Qualification results of water chemistry v3.2 (receipt %s)' % os.path.join(d, 'receipt.json'), '']
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
    rv = os.path.join(d, 'SWARM-REVALIDATION.txt')
    if os.path.exists(sv):
        note = ('Note (v3.2, decision of 2026-10-31): the agreement bands v3 stated (SCOPE.md section 8) were single-seed values, and at\n'
                'E/N between computed points they were interpolated (log-linear between 300 and 500 Td etc.); repeated seeds show a scatter\n'
                'of 0.5-3 % and the interpolation shifted the bulk drift at 396-497 Td by about +5 points. The v3 records are corrected\n'
                'to six-seed means run directly at the measured E/N (second table below, old records side by side); v3.2 is judged\n'
                'against them: each v3.2 mean within the v3 multi-seed mean +- max(2 SD, 1 percentage point). The first table is the\n'
                'single-seed, interpolated comparison of all sets as in v3 (kept for continuity).\n\n')
        body = open(sv).read()
        if os.path.exists(rv):
            body += '\n\n' + open(rv).read()
        open(os.path.join(R, 'SWARM-VALIDATION.txt'), 'w').write(note + body)
        L.append('')
        L.append('Swarm validation: SWARM-VALIDATION.txt (from %s, with the multi-seed re-validation table appended)' % sv)
    if os.path.exists(rv):
        L.append('Swarm re-validation of the changed electron set: ' + open(rv).read().splitlines()[2])
    cost = next((st for st in rc['stages'] if st['stage'] == 'cost'), None)
    if cost:
        L.append('')
        L.append('Sampling cost per collision (GPU kernel and single-thread CPU reference of tests/gpu_parity.cpp, 25000 carriers, 3 step cases):')
        for k, v in cost['cost_per_collision'].items():
            L.append('    %-24s %7d events  gpu %7.1f ns/event  cpu %6.2f us/event (single run)' % (k, v['events'], v['gpu_ns_per_event'], v['cpu_us_per_event']))
    bc = os.path.join(d, 'bench', 'bench-cost.txt')
    if os.path.exists(bc):
        best = {}
        for line in open(bc):
            cfg, case, js = line.split(' ', 2)
            j = json.loads(js)
            k = (cfg, case)
            if k not in best or j['min_ms'] < best[k]['min_ms']:
                best[k] = j
        rows = {}
        for (cfg, case), j in best.items():
            r = rows.setdefault(cfg, {'ms': 0.0, 'ev': 0, 'cand': 0})
            r['ms'] += j['min_ms']
            r['ev'] += j['events']
            r['cand'] += j['candidates']
        L.append('CPU benchmark on identical inputs (tests/bench_step.cpp; minimum over interleaved rounds of 5 runs each on one shared core; 3 cases, 25000 carriers each):')
        for cfg, r in rows.items():
            L.append('    %-20s %7d events %9d candidates  %8.1f ms  %6.2f us/event  %6.1f ns/candidate' % (cfg, r['ev'], r['cand'], r['ms'], 1e3 * r['ms'] / r['ev'],
                                                                                                      1e6 * r['ms'] / r['cand']))
    if analysis:
        for name in ('ion-beam-share.txt', 'breakdown-mix.txt', 'breakdown-mix-v3.1.txt'):
            p = os.path.join(analysis, name)
            if os.path.exists(p):
                L.append('')
                L.append('%s (%s):' % (name, p))
                L.extend('    ' + x[:200] for x in open(p).read().splitlines())
        p = os.path.join(analysis, 'validation-checks.json')
        if os.path.exists(p):
            v = json.load(open(p))
            r1 = [x['ratio'] for x in v['free_collision_vs_H_minus_sigma_m10']]
            r2 = [x['ratio'] for x in v['free_collision_vs_He_plus_sigma12'] if x['ratio']]
            o = v['ionization_6p6keV_O_plus']
            L.append('')
            L.append('Model checks (tools/validation_checks.py): free-collision estimate / H- sigma-10 (1-130 keV) %.2f-%.2f; / He+ sigma12 '
                     '(3000-10000 km/s) %.2f-%.2f; 6.6 keV O+: measured fragment-ion production %.2g m2, He+ rule %.2g, proton rule %.2g'
                     % (min(r1), max(r1), min(r2), max(r2), o['measured_total_positive_fragment_ion_production_m2'], o['He_plus_sigma11_equal_speed_m2'],
                        o['proton_sigma11_equal_speed_m2']))
    L.append('')
    L.append('Pending corrections of the corrective v3.1.1 (pic-water-chemistry-v3.1.1-20261031, in qualification, no receipt at the')
    L.append('v3.2 seal; v3.2 is sealed on v3.1 and is not the production binding -- v3.2.1 will fold them in and requalify):')
    for x in ('chem_event.h cover_support(): coverage exposure and bounds over the truncated thermal support (max of sigma_bound v over the',
              '  intersected bins, max per region), input validation of the step configuration and cells, exact energy of cold drifting',
              '  targets; raises every B_r (B/A tables and mix checks to be recomputed in v3.2.1).',
              'chem_host.hpp: neutral-parcel reservoir (pfc-neutral-parcels-v1, release API change), GAS-return non-debit, strict',
              '  unknown-channel handling (debit and release tests to be requalified; native v8 INTEGRATION note for the parcel API).',
              'make_coverage_policy.py strict production default: no unbounded region admitted below 10 MeV, placeholder regions unbounded',
              '  (ruling of 2026-10-31: adopted in v3.2.1, with the v3.2 placeholders below 1 keV, e + OH+ DR and H3O+ + OH- mutual',
              '  neutralization replaced by sampled estimates; a diagnostic opt-in may reproduce the v3.2 admissions for comparison).'):
        L.append('    ' + x)
    open(os.path.join(R, 'RESULTS.txt'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main(*sys.argv[1:3])
