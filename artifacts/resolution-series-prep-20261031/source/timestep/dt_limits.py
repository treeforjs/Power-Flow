#!/usr/bin/python -I
"""Time-step limits per mesh of the resolution series (PLAN.md section 8, coordinator v7 controller constants).

Inputs: the mesh report and the measured field-only runs (B_max per ampere of drive, CG iterations versus dt).
Output: dt-limits.json and dt-limits.txt in the given output directory.

  dt_limits.py OUTDIR MESH_REPORT.json GPU_VACUUM_JSON... [--dtscan GPU_DTSCAN_JSON...]

Criteria (coordinator v7 step_controller_v7.hpp defaults; Hall guards in brackets):
  displacement  dt v <= c_disp h_min                  c_disp = 0.1   [Hall push guard 0.2 cell]
  gyro          omega_c,max dt <= c_gyro              c_gyro = 0.025 [Hall push guard 0.05]
  plasma        omega_p,max dt <= c_plasma            c_plasma = 0.03 [Hall motion guard 0.1]
  particle cap  dt <= 10 fs while carriers exist (Hall 8dddec0 neutral stream and the electrode-heat v3 source owner keep the
                1e-14 literal; coordinator v7 b5)
  field         particle-free steps: field accuracy (third-difference LTE, lte_tol) and the CG iteration budget (maxit 200
                per solve in the resident path)
"""
import json, math, sys
from pathlib import Path

E = 1.602176634e-19
ME = 9.1093837015e-31
EPS0 = 8.8541878128e-12
C = 299792458.0
C_DISP, C_GYRO, C_PLASMA = 0.1, 0.025, 0.03
G_DISP, G_GYRO, G_PLASMA = 0.2, 0.05, 0.1
PARTICLE_CAP = 1e-14
# C19 currents of the windows considered (shot 17133 trace): 0.75 ns 6.0 kA; ~11.3 ns 50 kA; ~18.5 ns 100 kA; ~24.6 ns 150 kA;
# peak 820.5 kA at 119.2 ns.
CURRENTS = {'0.75ns_6.0kA': 6.012e3, '11.3ns_50kA': 5.0e4, '18.5ns_100kA': 1.0e5, '24.6ns_150kA': 1.5e5, 'peak_119ns_820kA': 8.205e5}
ELECTRON_EV = [100.0, 1e3, 1e4, 5e4]
DENSITIES = [1e16, 1e18, 1e20, 1e22, 1e24]
TE_EV = [1.0, 3.0, 10.0]


def v_of_ev(ev):
    g = 1 + ev * E / (ME * C * C)
    return C * math.sqrt(1 - 1 / (g * g))


def omega_p(n):
    return math.sqrt(n * E * E / (EPS0 * ME))


def main(argv):
    out = Path(argv[0])
    mesh = json.load(open(argv[1]))
    files = argv[2:]
    dtscan = []
    if '--dtscan' in files:
        k = files.index('--dtscan')
        dtscan, files = files[k + 1:], files[:k]
    runs = {}
    for f in files:
        d = json.load(open(f))
        if len(d['spec']) > 3 and abs(d['spec'][3] - 0.009) > 1e-12:
            continue   # other box widths (separate study)
        runs[tuple(d['spec'][:3])] = d
    scans = {}
    for f in dtscan:
        d = json.load(open(f))
        scans[tuple(d['spec'][:3])] = d
    rows = []
    finest = None   # (h, B per A) of the finest measured transverse mesh, used where a mesh has no run
    for key, run in runs.items():
        h = 0.009 / key[0]
        bpa = run['B_bounds_final']['neighbourhood_bound_max_T'] / run['final_current_A']
        if finest is None or h < finest[0]:
            finest = (h, bpa, run['control']['B_max_T'] / run['final_current_A'], key)
    for m in mesh['meshes']:
        s = m['spec']
        key = (s['transverse'], s['axial'], s['air'])
        h = m['spacing']['dx_m']['value']
        dz = m['spacing']['dz_m']['value']
        vnode = h * h * dz   # interior node dual volume
        r = {'spec': key, 'h_m': h, 'dz_m': dz, 'node_volume_m3': vnode}
        r['courant_at_10fs'] = C * 1e-14 / h
        # displacement
        r['displacement'] = {'c_light': C_DISP * h / C, **{'%g_eV' % ev: C_DISP * h / v_of_ev(ev) for ev in ELECTRON_EV},
                             'guard_c_light': G_DISP * h / C}
        # gyro (measured |B| bound per ampere at the start of the window, scaled with the drive current)
        run = runs.get(key)
        if run:
            bnb = run['B_bounds_final']['neighbourhood_bound_max_T']
            bcc = run['control']['B_max_T']
            amp = run['final_current_A']
            r['B_per_A'] = {'neighbourhood_bound_T_per_A': bnb / amp, 'cell_centre_T_per_A': bcc / amp, 'at_current_A': amp, 'source': 'measured'}
        elif finest is not None:
            # no field run on this mesh (not admitted): |B| bound per ampere of the finest measured mesh (the measured
            # values converge with refinement: +9 % from 200 to 100 um, +2.5 % from 100 to 50 um)
            run = {'_extrapolated': True}
            bnb, bcc, amp = finest[1], finest[2], 1.0
            r['B_per_A'] = {'neighbourhood_bound_T_per_A': bnb, 'cell_centre_T_per_A': bcc, 'at_current_A': None,
                            'source': 'finest measured mesh %s' % (finest[3],)}
        if 'B_per_A' in r:
            r['gyro'] = {}
            for name, cur in CURRENTS.items():
                b = r['B_per_A']['neighbourhood_bound_T_per_A'] * cur
                wc = E * b / ME
                r['gyro'][name] = {'B_bound_T': b, 'dt_c_gyro_s': C_GYRO / wc, 'dt_guard_s': G_GYRO / wc}
        # plasma: physical density -> dt; one macro-carrier of weight w in one node
        r['plasma'] = {'%.0e_m-3' % n: C_PLASMA / omega_p(n) for n in DENSITIES}
        r['single_carrier_node_density_m-3'] = {'w=%g' % w: w / vnode for w in (1e4, 1e5, 1e6)}
        # Debye resolution: density up to which lambda_D >= h, and dt from c_plasma at that density
        r['debye'] = {}
        for te in TE_EV:
            n_res = EPS0 * te / (E * h * h)
            r['debye']['Te_%g_eV' % te] = {'n_max_resolved_m-3': n_res, 'dt_c_plasma_at_n_max_s': C_PLASMA / omega_p(n_res),
                                           'lambda_D_at_1e16_m': math.sqrt(EPS0 * te / (E * 1e16)),
                                           'lambda_D_at_1e18_m': math.sqrt(EPS0 * te / (E * 1e18))}
        # mesh-implied particle-step dt (relativistic electrons, breakdown window currents 50-150 kA)
        if run:
            cand = {'particle_cap': PARTICLE_CAP, 'displacement_c': r['displacement']['c_light'],
                    'displacement_50keV': r['displacement']['50000_eV']}
            imp = {}
            for name in CURRENTS:
                lim = dict(cand)
                lim['gyro'] = r['gyro'][name]['dt_c_gyro_s']
                k = min(lim, key=lim.get)
                imp[name] = {'dt_s': lim[k], 'limited_by': k, 'dt_with_50keV_bound_s': min(PARTICLE_CAP, r['displacement']['50000_eV'], lim['gyro'])}
            r['implied_particle_dt'] = imp
        # field-only: CG iterations versus dt (sealed resident path)
        sc = scans.get(key)
        if sc:
            pts = []
            for row in sc['rows']:
                st = row['stats']
                pts.append({'dt_s': row['dt_s'], 'courant': row['courant_transverse'], 'steps_done': row['steps_done'],
                            'refused': bool(row['refusal']), 'iterations_median': st.get('iterations', {}).get('median'),
                            'iterations_max': st.get('iterations', {}).get('max'), 'step_ms_median': st.get('step_ms', {}).get('median'),
                            'ms_per_simulated_ps': (st.get('step_ms', {}).get('median') or 0) / (row['dt_s'] / 1e-12) if st.get('step_ms') else None})
            r['field_only_dtscan'] = pts
            ok = [p for p in pts if not p['refused'] and p['steps_done'] == sc['steps_per_dt']]
            r['field_only_max_admitted_dt_s'] = max(p['dt_s'] for p in ok) if ok else None
            best = min((p for p in ok if p['ms_per_simulated_ps']), key=lambda p: p['ms_per_simulated_ps'], default=None)
            r['field_only_cheapest_dt_s'] = best['dt_s'] if best else None
        rows.append(r)
    res = {'schema': 'resolution-series-dt-limits-v1', 'constants': {'c_disp': C_DISP, 'c_gyro': C_GYRO, 'c_plasma': C_PLASMA,
           'guards': {'disp': G_DISP, 'gyro': G_GYRO, 'plasma': G_PLASMA}, 'particle_cap_s': PARTICLE_CAP}, 'currents_A': CURRENTS, 'meshes': rows}
    out.mkdir(parents=True, exist_ok=True)
    (out / 'dt-limits.json').write_text(json.dumps(res, indent=1) + '\n')
    lines = ['Time-step limits per mesh (s); controller factors c_disp 0.1, c_gyro 0.025, c_plasma 0.03; particle cap 10 fs', '']
    for r in rows:
        lines.append('mesh %s  h = %.0f um  dz = %.3g mm  Courant(10 fs) = %.3f' % (r['spec'], r['h_m'] * 1e6, r['dz_m'] * 1e3, r['courant_at_10fs']))
        d = r['displacement']
        lines.append('  displacement (c_disp h / v): v=c %.3g, 50 keV %.3g, 10 keV %.3g, 1 keV %.3g, 100 eV %.3g' %
                     (d['c_light'], d['50000_eV'], d['10000_eV'], d['1000_eV'], d['100_eV']))
        if 'gyro' in r:
            lines.append('  |B| bound per A %.3e T/A (cell centre %.3e); gyro dt (c_gyro): ' % (r['B_per_A']['neighbourhood_bound_T_per_A'], r['B_per_A']['cell_centre_T_per_A']) +
                         ', '.join('%s %.3g (B %.3g T)' % (k, v['dt_c_gyro_s'], v['B_bound_T']) for k, v in r['gyro'].items()))
        lines.append('  plasma dt (c_plasma): ' + ', '.join('%s %.3g' % (k, v) for k, v in r['plasma'].items()))
        lines.append('  one carrier in one node: ' + ', '.join('%s n=%.3g' % (k, v) for k, v in r['single_carrier_node_density_m-3'].items()))
        lines.append('  Debye resolved up to: ' + ', '.join('%s n=%.3g m^-3' % (k, v['n_max_resolved_m-3']) for k, v in r['debye'].items()))
        if 'implied_particle_dt' in r:
            lines.append('  implied particle-step dt: ' + ', '.join('%s %.3g (%s)' % (k, v['dt_s'], v['limited_by']) for k, v in r['implied_particle_dt'].items()))
        if 'field_only_dtscan' in r:
            lines.append('  field-only CG (sealed Jacobi PCG, maxit 200): ' + ', '.join(
                '%.0e: %s its%s' % (p['dt_s'], p['iterations_median'], ' REFUSED' if p['refused'] else '') for p in r['field_only_dtscan']))
            lines.append('  field-only max admitted dt %s, cheapest per simulated time %s' % (r['field_only_max_admitted_dt_s'], r['field_only_cheapest_dt_s']))
        lines.append('')
    (out / 'dt-limits.txt').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main(sys.argv[1:])
