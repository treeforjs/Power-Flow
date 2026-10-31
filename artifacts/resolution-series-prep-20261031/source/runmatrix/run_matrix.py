#!/usr/bin/python -I
"""Run matrix of the resolution series (mesh x macroparticle count x dt) with wall-time and memory estimates.

  run_matrix.py OUTDIR INPUTS.json

INPUTS.json names the measured files (mesh report, field-only vacuum runs, dt scans, coupled Hall runs, GMG runs, dt limits);
every number of the estimate is either read from those files or listed in ASSUMPTIONS below (written into the output).
Outputs: run-matrix.json, run-matrix.txt.
"""
import json, math, sys
from pathlib import Path

ASSUMPTIONS = {
    'window': {'t0_s': 7.49788055492e-10, 'T_on_s': 1.0e-8, 'coupled_length_s': 2.0e-9,
               'note': 'Hall-only particle-free pre-run 0.75 ns -> T_on, then a fresh native genesis and the coupled window'},
    'species_per_member': 3,   # electrons, H2O+, neutral H2O carriers (the pilot census has all three)
    'plasma_cells_note': 'particles per cell over the 5 mm x 5 mm x 25 mm gap volume (cells transverse x 4 axial)',
    'coordinator_current': {
        'non_hall_fixed_ms_200um': 22.3,   # coordinator v6.1 pilot GPU step 72.8 ms minus Hall proposal 41.9 and Hall preflight 8.8
        'mesh_dependent_ms_200um': 2.8,     # 3 extra total-field digests (0.58 ms each) + source field build (1.1 ms), scale with ne+nb
        'non_hall_per_carrier_us': 15.0,    # coordinator v7 README: controller 4.4 ms = ~1 % of the step at 1.9e4 carriers (~20 us/carrier
                                            # in total), minus the Hall per-carrier host work measured here
        'source': 'Power-Flow/artifacts/coordinator-v6.1-20261031/source/README.txt; coordinator-v7-20261031/source/README.txt line 179'},
    'projected': {
        'per_carrier_us': 1.0,   # target after native v8 binary exact sums / Merkle motion digests (native v6 binary admission
                                 # measured 0.17-0.31 us/carrier), Hall device deposits and digests, chemistry vectorized
        'field_floor': 'measured resident vacuum step minus the host independent residual (assumed moved to a device-side '
                       'compensated residual) and minus the host total-field digests (assumed device-side block digests)',
        'non_hall_fixed_ms': 15.0},
    'station': {'gpu': 'RTX PRO 4000 Blackwell (sm_120), 24 GB, 672 GB/s (public specification, not measured here)',
                'device_time_factor': 912.0 / 672.0, 'host_time_factor': 1.0, 'concurrent_members_max': 6,
                'note': 'device-bound parts scaled by the memory-bandwidth ratio to the RTX 3080 Ti (912 GB/s); host single-thread speed '
                        'assumed equal to the Ryzen 7 9850X3D (range 1.0-1.3x slower); members run concurrently up to the 24 GB device '
                        'memory and a device-busy sum <= 1'},
    'zenobia': {'gpu': 'RTX 3080 Ti, 12 GB (shared with other work; about 1.6 GB in use by others at measurement time)',
                'concurrent_members_max': 2, 'note': 'Power-Flow cores 6-7, one member per core'},
    'hall_v31_record': {'coupled2_over_vacuum': 10.2 / 2.80, 'per_carrier_us': 5.3,
                        'source': 'Hall v3.1 RESULTS sections 5 and 9: coupled 2 carriers 10.2 ms vs vacuum 2.80 ms per step at 200 um; '
                                  '1e4 -> 1e5 carriers slope 5.3 us per carrier (GPU push)'},
    'box_extra_ms_per_edge_face': 15e-6,   # measured: uniform 200 um boxes 9 -> 35 -> 45 mm, 2.85 -> 21.8 -> 33.1 ms per vacuum
                                           # step for 0.08 -> 1.20 -> 2.08 M edges+faces (16.9 and 12.8 ms per million); metal cells
                                           # (thermal owner) are unchanged by the box
    'box_extra_host_fraction': 0.4,        # share of the extra cost that is host O(edges) work (residual, digests) removed in the projection
    'device_bytes_per_carrier': 1400,   # Hall v3.1 GPU push engine buffers (hall_gpu_push_v3.cpp: 944 B fixed per carrier,
                                        # 8 CIC start/end records, sinks, flights) plus edge records (~20 per carrier x 24 B)
    'host_bytes_per_carrier': 4096 + 1024,   # Hall CarrierLimits::bytes_per_carrier (4096) + native/coordinator rows (estimate)
    'cuda_context_bytes': 400 * 2**20,
}
C = 299792458.0


def load(p):
    return json.load(open(p))


def main(argv):
    out = Path(argv[0]); inp = load(argv[1])
    mesh = {tuple(m['spec'][k] for k in ('transverse', 'axial', 'air')): m for m in load(inp['mesh_report'])['meshes']}
    vac = {tuple(d['spec'][:3]): d for d in map(load, inp['gpu_vacuum'])}
    scans = {tuple(d['spec'][:3]): d for d in map(load, inp['gpu_dtscan'])}
    coupled = {}
    for f in inp['gpu_coupled']:
        d = load(f)
        if d.get('steps_done', 0) > 0:
            coupled.setdefault(tuple(d['spec'][:3]), {})[d['carriers']] = d
    gmg = {}
    for f in inp.get('gmg', []):
        d = load(f)
        if d.get('options', {}).get('nu', 1) == 1 and 'rows' in d:
            gmg[(tuple(d['spec'][:3]), d['rhs'])] = d
    dtl = {tuple(r['spec']): r for r in load(inp['dt_limits'])['meshes']}
    proxy = {tuple(k): tuple(v) for k, v in inp.get('proxy', {}).items()} if isinstance(inp.get('proxy'), dict) else {}
    series = [(45, 4, 1), (90, 4, 1), (180, 4, 1), (360, 4, 1)]
    proxy_for = {(360, 4, 1): (180, 16, 4)}
    A = ASSUMPTIONS
    rows = []
    for key in series:
        m = mesh[key]; h = m['spacing']['dx_m']['value']
        ne, nb, nn = m['counts']['edges'], m['counts']['faces'], m['counts']['nodes']
        src = key if key in vac else proxy_for.get(key)
        if src not in vac:   # no proxy run: scale the finest measured series mesh by the edge+face count
            src = (180, 4, 1)
        scale = (ne + nb) / (mesh[src]['counts']['edges'] + mesh[src]['counts']['faces']) if src != key else 1.0
        v = vac[src]
        ph = v['phases']
        vac_ms = ph['step.total_ms']['median'] * scale
        host_res_ms = ph['prepare.independent_residual']['median'] * scale
        # coupled Hall step: fixed part (2 carriers) and slope per carrier from the largest run
        cp = coupled.get(src, {})
        if 2 in cp:
            fixed_ms = cp[2]['phases']['step.total_ms']['median'] * scale
        else:
            fixed_ms = None
        slope_us = None
        big = max((c for c in cp if c > 2), default=None)
        if big and fixed_ms is not None:
            slope_us = (cp[big]['phases']['step.total_ms']['median'] - cp[2]['phases']['step.total_ms']['median']) / big * 1e3
        iters = cp.get(2, {}).get('phases', {}).get('proposal.iterations', {}).get('median') if cp else None
        coupled_source = 'measured (gpu-coupled on %s)' % (src,)
        if fixed_ms is None:   # no coupled run on this mesh yet: Hall v3.1 RESULTS (200 um, build a5) ratios
            fixed_ms = vac_ms * A['hall_v31_record']['coupled2_over_vacuum']
            coupled_source = 'estimated: vacuum step x %.2f (Hall v3.1 RESULTS, 200 um)' % A['hall_v31_record']['coupled2_over_vacuum']
        if slope_us is None:
            slope_us = A['hall_v31_record']['per_carrier_us']
            coupled_source += '; per-carrier slope %.1f us (Hall v3.1 RESULTS section 5)' % slope_us
        # memory
        dev_field = m['device_memory_model']['resident_field_bytes']
        # dt
        d = dtl[key]
        dt_eff = min(1e-14, d['displacement']['c_light'], d['gyro']['11.3ns_50kA']['dt_c_gyro_s'])
        limiter = 'displacement (c)' if dt_eff == d['displacement']['c_light'] else ('gyro' if dt_eff == d['gyro']['11.3ns_50kA']['dt_c_gyro_s'] else 'particle cap')
        # pre-run
        W = A['window']
        pre_span = W['T_on_s'] - W['t0_s']
        pre = {}
        pre['dt_10fs'] = {'dt_s': 1e-14, 'ms_per_step': vac_ms, 'hours': pre_span / 1e-14 * vac_ms / 3.6e6}
        sc = scans.get(src)
        if sc:
            best = None
            for r in sc['rows']:
                if r['refusal'] or r['steps_done'] < sc['steps_per_dt']:
                    continue
                ms = r['stats']['step_ms']['median'] * scale
                if src != key:   # Jacobi iterations double per halving of h: scale the CG part (step time above the 10 fs step)
                    base = sc['rows'][0]['stats']['step_ms']['median'] * scale
                    its = r['stats']['iterations']['median'] * (mesh[src]['spacing']['dx_m']['value'] / h)
                    if its > 200:
                        continue
                    ms = base + (ms - base) * (mesh[src]['spacing']['dx_m']['value'] / h)
                hours = pre_span / r['dt_s'] * ms / 3.6e6
                if best is None or hours < best['hours']:
                    best = {'dt_s': r['dt_s'], 'ms_per_step': ms, 'hours': hours, 'jacobi_iterations': r['stats']['iterations']['median']}
            pre['jacobi_best_admitted'] = best
        g = gmg.get((key, 'random'))
        if g:
            r1 = [r for r in g['rows'] if abs(r['dt_s'] - 1e-12) < 1e-20]
            r0 = [r for r in g['rows'] if abs(r['dt_s'] - 1e-14) < 1e-20]
            if r1 and r0:
                ms = vac_ms + r1[0]['gmg']['solve_ms']   # conservative: the GMG solve is added to the whole 10 fs step
                pre['gmg_1ps'] = {'dt_s': 1e-12, 'ms_per_step': ms, 'hours': pre_span / 1e-12 * ms / 3.6e6, 'gmg_iterations': r1[0]['gmg']['iterations'],
                                  'jacobi_iterations': r1[0]['jacobi']['iterations'], 'note': 'resident step at 10 fs plus the whole GMG solve at 1 ps (random rhs, prototype kernels)'}
        cells = (0.005 / h) ** 2 * 4
        for N in (1e4, 1e5, 1e6):
            P = A['species_per_member'] * N
            for dtf, dtname in ((1.0, 'dt'), (0.5, 'dt/2')):
                dt = dt_eff * dtf
                steps = W['coupled_length_s'] / dt
                cc = A['coordinator_current']
                coord_ms = cc['non_hall_fixed_ms_200um'] + cc['mesh_dependent_ms_200um'] * (ne + nb) / (41676 + 39015)
                cur_ms = None
                if fixed_ms is not None and slope_us is not None:
                    cur_ms = fixed_ms + slope_us * P / 1e3 + coord_ms + cc['non_hall_per_carrier_us'] * P / 1e3
                floor_proj = (vac_ms - host_res_ms) * (fixed_ms / vac_ms if fixed_ms else 2.0)
                proj_ms = floor_proj + A['projected']['non_hall_fixed_ms'] + A['projected']['per_carrier_us'] * P / 1e3
                dev_bytes = dev_field + P * A['device_bytes_per_carrier'] + A['cuda_context_bytes']
                host_bytes = (v.get('max_rss_kB', 0) * 1024 * scale) + P * A['host_bytes_per_carrier']
                row = {'mesh': key, 'h_um': h * 1e6, 'edges': ne, 'carriers_per_species': N, 'carriers_total': P, 'ppc_gap': N / cells,
                       'dt_name': dtname, 'dt_s': dt, 'dt_limiter': limiter, 'coupled_steps': steps,
                       'device_bytes': dev_bytes, 'fits_12GB': dev_bytes < 12e9 - 1.7e9, 'fits_24GB': dev_bytes < 24e9 * 0.95,
                       'host_bytes': host_bytes, 'fits_host_zenobia_14GB': host_bytes < 14e9,
                       'per_carrier_us_needed_for_24h': max(0.0, (24 * 3.6e6 / steps - floor_proj - A['projected']['non_hall_fixed_ms']) / P * 1e3),
                       'ms_per_step_current': cur_ms, 'hours_coupled_current': cur_ms * steps / 3.6e6 if cur_ms else None,
                       'ms_per_step_projected': proj_ms, 'hours_coupled_projected': proj_ms * steps / 3.6e6,
                       'station_hours_coupled_projected': proj_ms * steps / 3.6e6 * (0.5 * A['station']['device_time_factor'] + 0.5),
                       'pre_run': pre, 'field_only_ms_per_step_10fs': vac_ms, 'hall_fixed_point_iterations': iters,
                       'hall_coupled_fixed_ms': fixed_ms, 'hall_per_carrier_us': slope_us, 'coupled_cost_source': coupled_source,
                       'proxy_mesh': src if src != key else None}
                rows.append(row)
    # minimal series: mesh series at 1e5/species (dt), macroparticle series at 50 um (1e4, 1e5, 1e6), dt series at 50 um 1e5 (dt/2)
    def pick(h, N, dtn):
        return next(r for r in rows if abs(r['h_um'] - h) < 1 and r['carriers_per_species'] == N and r['dt_name'] == dtn)
    members = [pick(h, 1e5, 'dt') for h in (200, 100, 50, 25)] + [pick(50, 1e4, 'dt'), pick(50, 1e6, 'dt'), pick(50, 1e5, 'dt/2')]
    def total(key, conc, cap_bytes):
        hrs = [m[key] for m in members if m[key] is not None]
        if len(hrs) < len(members):
            return None
        # members run `conc` at a time (longest first) when their device memory fits together
        hrs = sorted(hrs, reverse=True); lanes = [0.0] * conc
        for x in hrs:
            i = lanes.index(min(lanes)); lanes[i] += x
        return max(lanes)
    pre = {}
    for key in series:
        p = next(r for r in rows if r['mesh'] == key)['pre_run']
        opts = [p['dt_10fs']['hours']] + ([p['jacobi_best_admitted']['hours']] if p.get('jacobi_best_admitted') else []) + \
               ([p['gmg_1ps']['hours']] if p.get('gmg_1ps') else [])
        pre[str(key)] = {'hours_10fs': p['dt_10fs']['hours'], 'hours_best': min(opts)}
    pre_total = sum(v['hours_best'] for v in pre.values())
    series_tot = {'members': [{'h_um': m['h_um'], 'N': m['carriers_per_species'], 'dt': m['dt_name']} for m in members],
                  'pre_runs_hours_best_total': pre_total, 'pre_runs': pre,
                  'coupled_hours_current_zenobia': total('hours_coupled_current', A['zenobia']['concurrent_members_max'], 12e9),
                  'coupled_hours_projected_zenobia': total('hours_coupled_projected', A['zenobia']['concurrent_members_max'], 12e9),
                  'coupled_hours_projected_station': total('station_hours_coupled_projected', A['station']['concurrent_members_max'], 24e9),
                  'note': 'longest-first packing on 2 lanes (zenobia, cores 6-7) or 6 lanes (station); pre-runs are per mesh and can '
                          'share the lanes; the 50 um / 1e5 member is shared by the three series'}
    # box variants (graded transverse meshes to a 35 / 45 mm box, core 9 mm at h, q = 1.1, cells <= 2 mm): sizes from
    # graded-sizes, step cost scaled from the 9 mm mesh of the same h by (edges + faces); uniform big boxes measured where run
    box_rows = []
    gs = load(inp['graded_sizes'])['rows'] if inp.get('graded_sizes') else []
    for key in series:
        base = next(r for r in rows if r['mesh'] == key and r['carriers_per_species'] == 1e4 and r['dt_name'] == 'dt')
        m = mesh[key]; h = m['spacing']['dx_m']['value']
        ne9, nb9 = m['counts']['edges'], m['counts']['faces']
        for g in gs:
            if abs(g['h_m'] - h) > 1e-12 or abs(g['q'] - 1.1) > 1e-12 or g['box_half_m'] < .01:
                continue
            dE = (g['edges'] + g['faces']) - (ne9 + nb9)   # extra edges+faces of the graded outer region (no metal there)
            a = A['box_extra_ms_per_edge_face']
            for N in (1e4, 1e5):
                r9 = next(r for r in rows if r['mesh'] == key and r['carriers_per_species'] == N and r['dt_name'] == 'dt')
                P = r9['carriers_total']
                vac9 = r9['field_only_ms_per_step_10fs']; ratio = r9['hall_coupled_fixed_ms'] / vac9
                vac_box = vac9 + a * dE
                cc = A['coordinator_current']
                coord_ms = cc['non_hall_fixed_ms_200um'] + cc['mesh_dependent_ms_200um'] * (g['edges'] + g['faces']) / (41676 + 39015)
                cur = (r9['hall_coupled_fixed_ms'] + ratio * a * dE) + r9['hall_per_carrier_us'] * P / 1e3 + coord_ms + cc['non_hall_per_carrier_us'] * P / 1e3
                proj = r9['ms_per_step_projected'] + ratio * a * dE * (1 - A['box_extra_host_fraction'])
                dev = g['resident_field_bytes'] + P * A['device_bytes_per_carrier'] + A['cuda_context_bytes']
                box_rows.append({'mesh_h_um': h * 1e6, 'box_mm': 2e3 * g['box_half_m'], 'n_transverse': g['nx'], 'edges': g['edges'], 'faces': g['faces'],
                                 'nodes': g['nodes'], 'largest_cell_mm': g['largest_cell_m'] * 1e3, 'carriers_per_species': N,
                                 'resident_field_bytes': g['resident_field_bytes'], 'device_bytes': dev, 'fits_12GB': dev < 12e9 - 1.7e9,
                                 'field_only_ms_per_step_10fs': vac_box, 'ms_per_step_current': cur,
                                 'hours_coupled_current': cur * r9['coupled_steps'] / 3.6e6, 'ms_per_step_projected': proj,
                                 'hours_coupled_projected': proj * r9['coupled_steps'] / 3.6e6,
                                 'pre_run_hours_10fs': r9['pre_run']['dt_10fs']['hours'] * vac_box / vac9,
                                 'hall_v31_uniform_box_admissible': (round(2 * g['box_half_m'] / h) <= 500)})
    res = {'schema': 'resolution-series-run-matrix-v1', 'assumptions': A, 'inputs': inp, 'rows': rows, 'series_total': series_tot, 'box_variants': box_rows}
    out.mkdir(parents=True, exist_ok=True)
    (out / 'run-matrix.json').write_text(json.dumps(res, indent=1) + '\n')
    L = ['Run matrix (coupled window %.1f ns from T_on = %.1f ns; pre-run from 0.75 ns); hours per member' % (A['window']['coupled_length_s'] * 1e9, A['window']['T_on_s'] * 1e9),
         '%-6s %-7s %-5s %8s %7s %9s %8s %8s %10s %10s %10s %8s %8s %9s' % ('h(um)', 'N/sp', 'dt', 'dt(fs)', 'ppc', 'steps', 'dev(GB)', 'host(GB)', 'ms/st cur', 'h cur',
                                                                       'ms/st proj', 'h proj', 'h stat', 'us/c 24h')]
    for r in rows:
        L.append('%-6.0f %-7.0e %-5s %8.2f %7.1f %9.3g %8.2f %8.2f %10s %10s %10.1f %8.2f %8.2f %9.3f' % (
            r['h_um'], r['carriers_per_species'], r['dt_name'], r['dt_s'] * 1e15, r['ppc_gap'], r['coupled_steps'], r['device_bytes'] / 1e9,
            r['host_bytes'] / 1e9,
            '%.0f' % r['ms_per_step_current'] if r['ms_per_step_current'] else '-', '%.3g' % r['hours_coupled_current'] if r['hours_coupled_current'] else '-',
            r['ms_per_step_projected'], r['hours_coupled_projected'], r['station_hours_coupled_projected'], r['per_carrier_us_needed_for_24h']))
    L.append('')
    L.append('Pre-run (particle-free, 0.75 -> %.0f ns), hours:' % (A['window']['T_on_s'] * 1e9))
    for key in series:
        r = next(x for x in rows if x['mesh'] == key)
        p = r['pre_run']
        L.append('  %s: 10 fs %.2f h (%.1f ms/step)' % (key, p['dt_10fs']['hours'], p['dt_10fs']['ms_per_step']) +
                 ('; Jacobi best dt %.0e: %.2f h (%.1f ms/step)' % (p['jacobi_best_admitted']['dt_s'], p['jacobi_best_admitted']['hours'], p['jacobi_best_admitted']['ms_per_step']) if p.get('jacobi_best_admitted') else '') +
                 ('; GMG 1 ps: %.2f h (%.1f ms/step)' % (p['gmg_1ps']['hours'], p['gmg_1ps']['ms_per_step']) if p.get('gmg_1ps') else ''))
    L.append('')
    L.append('Minimal series (7 members: mesh 200/100/50/25 at 1e5/species; 50 um at 1e4 and 1e6; 50 um 1e5 at dt/2):')
    st = res['series_total']
    fmt = lambda x: ('%.1f h' % x) if x is not None else 'n/a'
    L.append('  pre-runs (best admitted dt per mesh): %.1f h in total' % st['pre_runs_hours_best_total'])
    L.append('  coupled windows: current code on zenobia %s; projected code on zenobia %s, on the station %s' % (
        fmt(st['coupled_hours_current_zenobia']), fmt(st['coupled_hours_projected_zenobia']), fmt(st['coupled_hours_projected_station'])))
    if box_rows:
        L.append('')
        L.append('Box variants (graded transverse mesh, uniform core 9 mm at h, growth 1.1 per cell up to ~2 mm cells; cost = 9 mm mesh + 15 ms per million extra edges+faces, measured on uniform 200 um boxes):')
        L.append('%-6s %-6s %-5s %10s %8s %-7s %9s %10s %9s %10s %9s %9s' % ('h(um)', 'box', 'n', 'edges', 'field GB', 'N/sp', 'dev(GB)', 'vac ms/st', 'ms cur', 'h cur', 'ms proj', 'h proj'))
        for b in box_rows:
            L.append('%-6.0f %-6.0f %-5d %10d %8.2f %-7.0e %9.2f %10.1f %9.0f %10.3g %9.1f %9.2f' % (
                b['mesh_h_um'], b['box_mm'], b['n_transverse'], b['edges'], b['resident_field_bytes'] / 1e9, b['carriers_per_species'], b['device_bytes'] / 1e9,
                b['field_only_ms_per_step_10fs'], b['ms_per_step_current'], b['hours_coupled_current'], b['ms_per_step_projected'], b['hours_coupled_projected']))
    (out / 'run-matrix.txt').write_text('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main(sys.argv[1:])
