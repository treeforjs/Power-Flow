#!/usr/bin/python -I
"""Compare two Hall v3.2 pre-run series on the same landing grid (test_v32 prerun series.jsonl) and their end-state
comparison (test_v32 compare compare.json): field energies, terminal voltage, foil temperatures, surface fields, Gauss
ledger ratios. Bounds of the v3.1 step-size evidence (100 ps window, 100 fs against 2.5 fs): E energy-norm error 1.2e-6,
terminal voltage 0.05 % of max |V|.

  compare_series.py COARSE_SERIES FINE_SERIES COMPARE_JSON OUT_JSON
"""
import json, sys


def load(p):
    return {r['m']: r for r in map(json.loads, open(p))}


def main(coarse, fine, cmp_path, out):
    a, b = load(coarse), load(fine)
    common = sorted(set(a) & set(b))
    if not common or any(a[m]['time_s'] != b[m]['time_s'] for m in common):
        raise SystemExit('series are not on the same landing grid')
    rel = lambda x, y: abs(x - y) / max(abs(y), 1e-300)
    vmax = max(abs(b[m]['terminal_V']) for m in common)
    dv = max(abs(a[m]['terminal_V'] - b[m]['terminal_V']) for m in common)
    r = {'schema': 'hall-v3.2-crosscheck-v1', 'coarse': coarse, 'fine': fine, 'samples': len(common),
         'time_range_s': [b[common[0]]['time_s'], b[common[-1]]['time_s']],
         'W_E_max_rel_diff': max(rel(a[m]['W_E_J'], b[m]['W_E_J']) for m in common),
         'W_B_max_rel_diff': max(rel(a[m]['W_B_J'], b[m]['W_B_J']) for m in common),
         'terminal_V_max_abs_diff_V': dv, 'terminal_V_max_abs_V': vmax, 'terminal_V_max_diff_over_max_V': dv / vmax}
    rise = {f: max(b[m][f]['max_T_K'] - 300 for m in common) for f in ('lower_foil', 'upper_foil')}
    r['foil_max_T_max_abs_diff_K'] = {f: max(abs(a[m][f]['max_T_K'] - b[m][f]['max_T_K']) for m in common) for f in rise}
    r['foil_max_T_rise_K'] = rise
    r['foil_max_T_diff_over_rise'] = {f: (r['foil_max_T_max_abs_diff_K'][f] / rise[f] if rise[f] > 0 else 0.) for f in rise}
    faces = list(b[common[-1]]['surface_field'])
    r['surface_field_max_rel_diff'] = {f: max(rel(a[m]['surface_field'][f]['max_abs_E_normal_V_m'], b[m]['surface_field'][f]['max_abs_E_normal_V_m'])
                                              for m in common if b[m]['surface_field'][f]['max_abs_E_normal_V_m'] > 0) for f in faces}
    g = lambda s, k: max(s[m]['gauss'][k] for m in common)
    r['gauss'] = {k: {'coarse_max': g(a, k), 'fine_max': g(b, k)} for k in
                  ('guard_ratio_max_abs_D_minus_L_over_bound', 'v31_form_ratio_max_abs_D_over_bound', 'implied_ratio_max_abs_D_over_bound_plus_abs_L',
                   'max_abs_L_C')}
    r['gauss_end'] = {'coarse': a[common[-1]]['gauss'], 'fine': b[common[-1]]['gauss']}
    r['divB_over_bound_max'] = {'coarse': max(a[m]['divB_over_bound_max_max'] for m in common), 'fine': max(b[m]['divB_over_bound_max_max'] for m in common)}
    c = json.load(open(cmp_path))
    r['end_state'] = c
    r['bounds'] = {'E_energy_norm': 1.2e-6, 'terminal_V_fraction': 5e-4}
    r['within_bounds'] = c['E_energy_norm_rel_diff'] <= 1.2e-6 and r['terminal_V_max_diff_over_max_V'] <= 5e-4
    open(out, 'w').write(json.dumps(r, indent=1) + '\n')
    print(json.dumps({k: r[k] for k in ('samples', 'W_E_max_rel_diff', 'W_B_max_rel_diff', 'terminal_V_max_diff_over_max_V', 'foil_max_T_diff_over_rise',
                                        'within_bounds')}))
    print(json.dumps({'E_energy_norm_rel_diff': c['E_energy_norm_rel_diff'], 'B_energy_norm_rel_diff': c['B_energy_norm_rel_diff']}))


if __name__ == '__main__':
    if len(sys.argv) != 5:
        raise SystemExit(__doc__)
    main(*sys.argv[1:])
