#!/usr/bin/python -I
"""Physical state of the Hall v3.4 20 ns vacuum pre-run against the sealed v3.3 pre-run (build s1, gpu-prerun-20ns-130711,
same arguments): at every 1 ns checkpoint the compensated field words (E, E_low, B, B_low), revision and time, the charge
ledger arrays, the admitted Gauss ledger, the kinetic ledgers and the heat slices byte for byte. Identity-bearing values
(configuration, solver descriptor, state digests, checkpoint schema) differ by design (README section 4) and are listed,
not compared; v3.4 adds the divB gross ledger (divb_gross_T_m), reported with its maximum.
Usage: compare_prerun_v33.py V33_RUN_DIR V34_RUN_DIR OUT_JSON"""
import json, sys
from pathlib import Path

old, new, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
PHYSICAL = ('field', 'charge', 'admitted_gauss_C', 'bank_residue_C', 'residue_present', 'kinetic_ledgers', 'dt_max_s',
            'neutral_stream_dt_max_s', 'motion_forms_native_v7', 'motion_mesh_sha256', 'solver_name')
IDENTITY = ('configuration_sha256', 'solver_descriptor_sha256', 'field_sha256', 'heat_sha256', 'kinetic_ledgers_sha256', 'checkpoint', 'schema')


def ckpts(d):
    return {round(c['time_s'] * 1e12): Path(c['dir']) for c in json.loads((d / 'checkpoints.json').read_text())}


a, b = ckpts(old), ckpts(new)
rows, ok = [], True
for t in sorted(set(a) | set(b)):
    r = {'time_ps': t}
    if t not in a or t not in b:
        r['present'] = [t in a, t in b]; ok = False; rows.append(r); continue
    ja, jb = (json.loads((x[t] / 'checkpoint.json').read_text()) for x in (a, b))
    same = {k: ja.get(k) == jb.get(k) for k in PHYSICAL}
    ha, hb = a[t] / 'heat', b[t] / 'heat'
    names = sorted(p.name for p in ha.iterdir())
    same['heat_slices'] = names == sorted(p.name for p in hb.iterdir()) and all((ha / f).read_bytes() == (hb / f).read_bytes() for f in names)
    r['identical'] = same
    r['identity_values_differ'] = sorted(k for k in IDENTITY if ja.get(k) != jb.get(k))
    g = jb.get('charge', {}).get('divb_gross_T_m') or jb.get('divb_gross_T_m') or []
    r['divb_gross_max_T_m'] = max(g) if g else 0.0
    r['added_keys'] = sorted(set(jb) - set(ja)) + ['charge/' + k for k in sorted(set(jb.get('charge', {})) - set(ja.get('charge', {})))]
    rows.append(r)
    ok = ok and all(same.values())
res = {'schema': 'hall-v3.4-prerun-vs-v3.3-v1', 'v33_run': str(old), 'v34_run': str(new), 'checkpoints': len(rows),
       'physical_state_bit_identical': ok, 'rows': rows}
out.write_text(json.dumps(res, indent=1) + '\n')
print(json.dumps({'checkpoints': len(rows), 'physical_state_bit_identical': ok,
                  'identity_values_differ': rows[-1].get('identity_values_differ') if rows else None,
                  'added_keys': rows[-1].get('added_keys') if rows else None}))
