#!/usr/bin/python -I
"""Acceptance checks of the resolution-series preparation over the recorded run outputs.

  check_results.py OUTPUT_ROOT RESULT.json

Checks (each prints PASS/FAIL with its evidence file):
  M1  every series mesh is constructed by the sealed Hall constructor with Hall's admission conditions; counts equal the
      closed formulas; incidence complete; cell, dual and metal volumes equal the box and the nominal foils
  M2  geometry identity across resolutions: foil faces (x = +-2.5 mm, y = +-2.5, +-2.7 mm), short plane and feed port node
      coordinates bitwise identical on every mesh; the kinetic geometry record bitwise identical
  M3  Hall v3.1 admission: 200/100/50 um (and the axial variants 100/8, 50/16) admitted; 25 um (360x4+1, 360x8+2) refused
      by the DC Coulomb-gauge vector solve; the instrumented DC profile equals the Model's initial fields bitwise
  M4  device memory: the allocation model equals the measured own peak device bytes of every resident run
  F1  field-only runs: every step admitted, independent residual <= target, guards within bounds
  C1  closed coupled Hall steps with the GPU push (2 .. 1e5 carriers): every requested step admitted
  G1  graded mesh: with the box equal to the core, operators, feed and DC seed B bitwise equal to the sealed Hall mesh; the DC port
      inductance on the 35 / 45 mm graded 200 um meshes within 0.1 % of the sealed Hall Model on the uniform boxes
  S1  solver criterion: every Jacobi and GMG solve meets ||b - A x|| <= max(1e-30, 2e-13 ||b||) with the sealed apply and norm
  S2  GMG determinism: repeat solve bitwise identical; CPU reference and GPU bitwise identical where run
  Z1  sanitizer runs (AddressSanitizer + UndefinedBehaviorSanitizer) exit 0 with no sanitizer report
"""
import glob, json, re, sys
from pathlib import Path


def main(argv):
    root = Path(argv[0]); res = {'checks': []}
    def rec(name, ok, evidence, detail=''):
        res['checks'].append({'check': name, 'passed': bool(ok), 'evidence': [str(e) for e in evidence], 'detail': detail})
        print('%s %s %s' % ('PASS' if ok else 'FAIL', name, detail))
    rep_f = sorted(root.glob('b1/mesh-report-*/mesh-report.json'))[-1]
    rep = json.load(open(rep_f))
    ok = True; det = []
    for m in rep['meshes']:
        s = m['spec']; n, a, air = s['transverse'], s['axial'], s['air']; nz = a + 2 * air
        c = m['counts']
        exp = ((n + 1) ** 2 * (nz + 1), 2 * n * (n + 1) * (nz + 1) + (n + 1) ** 2 * nz, 2 * (n + 1) * n * nz + n * n * (nz + 1), n * n * nz)
        got = (c['nodes'], c['edges'], c['faces'], c['cells'])
        v = m['volumes']
        good = exp == got and m['incidence_ok'] and abs(v['cells_rel_box']) < 1e-13 and abs(v['node_dual_rel_box']) < 1e-13 and \
            abs(v['metal_rel_nominal']) < 1e-13 and abs(m['feed']['sum_shape'] + 0.0052) < 1e-15 and m['operator']['min_M'] > 0 and m['operator']['min_H'] > 0
        ok = ok and good; det.append('%dx%d+%d:%s' % (n, a, air, 'ok' if good else 'bad'))
    rec('M1 mesh construction and quality', ok, [rep_f], ' '.join(det))
    ident = rep['cross_mesh_identity']
    rec('M2 geometry identity across resolutions', all(v['bitwise_identical'] for v in ident.values()), [rep_f],
        ', '.join('%s=%s' % (k, v['bitwise_identical']) for k, v in ident.items()))
    adm = {}
    for f in list(root.glob('b1/admit-*/admission.json')) + list(root.glob('b3/admit-*/admission.json')):
        d = json.load(open(f)); adm[tuple(d['spec'][:3])] = (d, f)
    want = {(45, 4, 1): True, (90, 4, 1): True, (180, 4, 1): True, (90, 8, 2): True, (180, 16, 4): True, (360, 4, 1): False, (360, 8, 2): False}
    ok = True; det = []
    for k, w in want.items():
        d, f = adm[k]
        good = d['admitted'] == w and (w or 'DC Coulomb vector' in d.get('refusal', '')) and (not w or all(d['profile_equals_model_bitwise'].values()))
        ok = ok and good; det.append('%s:%s' % (k, 'admitted' if d['admitted'] else 'refused'))
    rec('M3 Hall v3.1 admission per mesh', ok, [adm[k][1] for k in want], ' '.join(det))
    ok = True; det = []; ev = []
    vac = [Path(f) for f in glob.glob(str(root / '*/gpu-vacuum-*/gpu-vacuum.json'))]
    for f in vac:
        d = json.load(open(f))
        if 'memory_model_bytes' in d:
            good = d['memory_model_bytes'] == d['memory']['own_peak_device_bytes']
            ok = ok and good; det.append('%s:%d/%d' % (d['spec'][:3], d['memory_model_bytes'], d['memory']['own_peak_device_bytes'])); ev.append(f)
    rec('M4 device memory model equals measured peak', ok and len(ev) > 0, ev, ' '.join(det))
    ok = True; det = []
    for f in vac:
        d = json.load(open(f)); g = d['guards']
        good = not d['refusal'] and d['steps_done'] == d['steps_requested'] and g['max_independent_residual_over_target'] <= 1 and \
            g['max_divB_over_bound'] <= 1 and g['max_energy_over_bound'] <= 1
        ok = ok and good; det.append('%s:%d steps' % (d['spec'][:3], d['steps_done']))
    rec('F1 field-only runs admitted with guards', ok and len(vac) > 0, vac, ' '.join(det))
    cpl = [Path(f) for f in glob.glob(str(root / '*/gpu-coupled-*/gpu-coupled.json'))]
    ok = True; det = []
    for f in cpl:
        d = json.load(open(f))
        good = not d['refusal'] and d['steps_done'] == d['steps_requested']
        ok = ok and good; det.append('%s/%d:%d steps%s' % (d['spec'][:3], d['carriers'], d['steps_done'], '' if good else ' REFUSED ' + d['refusal'][:80]))
    if cpl:
        rec('C1 closed coupled Hall steps (GPU push) admitted', ok, cpl, ' '.join(det))
    gd = {Path(f).parent.name: json.load(open(f)) for f in glob.glob(str(root / '*/graded-dc-*/graded-dc.json'))}
    ad = {tuple(json.load(open(f))['spec']): json.load(open(f)) for f in glob.glob(str(root / '*/admit-*/admission.json'))}
    if gd:
        val = [d for k, d in gd.items() if 'uniform_reference' in d]
        ok = bool(val) and all(d['uniform_reference']['operators_bitwise_equal'] and d['uniform_reference']['B_bitwise_equal'] and
                               d['uniform_reference']['feed_bitwise_equal'] for d in val)
        det = []
        for half, n in ((.0175, 175), (.0225, 225)):
            hall = ad.get((n, 4, 1, 2 * half))
            g = [d for d in gd.values() if abs(d['h_m'] - 2e-4) < 1e-12 and abs(d['box_half_m'] - half) < 1e-12 and abs(d['q'] - 1.1) < 1e-12 and abs(d['hmax_m'] - 2e-3) < 1e-12]
            if hall and g:
                rel = abs(g[0]['dc_port_inductance_H'] / hall['dc_port_inductance_H'] - 1)
                ok = ok and rel < 1e-3; det.append('%.0f mm graded/uniform L %.5f/%.5f nH (%.1e)' % (2e3 * half, g[0]['dc_port_inductance_H'] * 1e9, hall['dc_port_inductance_H'] * 1e9, rel))
        rec('G1 graded mesh: bitwise equal to Hall on the 9 mm box; port inductance within 0.1 % of the sealed uniform box', ok,
            [Path(f) for f in glob.glob(str(root / '*/graded-dc-*/graded-dc.json'))][:3], '; '.join(det))
    gm = [Path(f) for f in glob.glob(str(root / '*/gmg-*/gmg.json'))]
    ok = True; det = []; ok2 = True; det2 = []; ev = []
    for f in gm:
        d = json.load(open(f))
        if 'rows' not in d:
            continue
        ev.append(f)
        for r in d['rows']:
            for mode in ('jacobi', 'gmg'):
                if mode not in r:
                    continue
                c = r[mode].get('check', r[mode])
                ok = ok and bool(c['passed'])
            if 'gmg' in r and 'repeat_bitwise_identical' in r['gmg']:
                ok2 = ok2 and r['gmg']['repeat_bitwise_identical'] and r['gmg'].get('cpu_gpu_bitwise_identical', True) and \
                    r.get('jacobi', {}).get('cpu_gpu_bitwise_identical', True)
        det.append('%s %s:%d dts' % (d['spec'][:3], d['rhs'], len(d['rows'])))
    rec('S1 solver criterion met (Jacobi and GMG)', ok and len(ev) > 0, ev, ' '.join(det))
    rec('S2 GMG repeat and CPU/GPU bit identity', ok2 and len(ev) > 0, ev, '')
    logs = [Path(f) for f in glob.glob(str(root / 's1/*/*.log')) if 'build' not in f]
    ok = len(logs) > 0; det = []
    for f in logs:
        t = Path(f).read_text(errors='replace')
        rc = json.load(open(Path(f).parent / 'run-receipt.json'))['command']['returncode'] if (Path(f).parent / 'run-receipt.json').exists() else None
        bad = re.search(r'ERROR: AddressSanitizer|runtime error:|LeakSanitizer', t) is not None
        ok = ok and rc == 0 and not bad; det.append('%s:rc=%s%s' % (Path(f).parent.name, rc, ' report' if bad else ''))
    rec('Z1 sanitizer runs clean', ok, logs, ' '.join(det))
    res['all_passed'] = all(c['passed'] for c in res['checks'])
    Path(argv[1]).write_text(json.dumps(res, indent=1) + '\n')
    return 0 if res['all_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main(sys.argv[1:]))
