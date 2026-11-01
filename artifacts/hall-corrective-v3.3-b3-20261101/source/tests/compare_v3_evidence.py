#!/usr/bin/python -I
"""Compare the Hall v3 component outputs produced on the v3.1 sources with the Hall v3 evidence of record
(b7/cpu-211226): component outcomes, physical-state parity scenarios, motion records and sink records.

  compare_v3_evidence.py OUTFILE COMPONENTS_DIR PARITY_JSON MOTION_JSON SINK_JSON
"""
import json, sys
from pathlib import Path

REF = Path('/home/treefor/projects/.workspace/outputs/hall-corrective-v3-20261030/b7/cpu-211226')


def load(p):
    return json.loads(Path(p).read_text())


def main(out, comp, parity, motion, sink):
    r = {'schema': 'hall-v3.1-v3-evidence-comparison-v1', 'reference': str(REF)}
    c, cref = load(Path(comp) / 'components.json'), load(REF / 'components/components.json')
    r['components'] = {'passed': sum(1 for x in c['cases'] if x['passed']), 'cases': len(c['cases']), 'failed': c['failed'],
                       'same_case_names_and_outcomes': [(x['group'], x['name'], x['passed']) for x in c['cases']] ==
                       [(x['group'], x['name'], x['passed']) for x in cref['cases']]}
    keys = ('physical_sha256', 'motion_record_text_sha256', 'motion_seal_v2', 'motion_record_v3_sha256', 'motion_arrays_sha256', 'iterations',
            'combined_balance_J', 'field_heat_balance_J', 'field_receipt', 'refusal')
    p, pref = load(parity), load(REF / 'parity-v3.json')
    rows = [{'scenario': a['scenario'], 'step': a['step'], 'identical': all(a.get(k) == b.get(k) for k in keys),
             'differences': [k for k in keys if a.get(k) != b.get(k)]} for a, b in zip(pref, p)]
    r['parity_scenarios'] = {'identical': len(p) == len(pref) and all(x['identical'] for x in rows), 'steps': rows}
    m, mref = load(motion), load(REF / 'motion-dump-v3.json')
    strip = lambda rows: [{k: v for k, v in x.items() if k != 'wall_s'} for x in rows]   # wall-clock timing is not a record value
    r['motion_records'] = {'identical': strip(m) == strip(mref), 'entries': len(m), 'excluded_keys': ['wall_s']}
    s, sref = load(sink), load(REF / 'sink-dump-cpu.json')
    r['sink_records'] = {'identical': s == sref, 'fixtures': len(s)}
    r['all_identical'] = r['components']['same_case_names_and_outcomes'] and r['components']['failed'] == 0 and r['parity_scenarios']['identical'] \
        and r['motion_records']['identical'] and r['sink_records']['identical']
    Path(out).write_text(json.dumps(r, indent=1) + '\n')
    print(json.dumps({k: (v if k != 'parity_scenarios' else {'identical': v['identical']}) for k, v in r.items()}, indent=1))
    return 0 if r['all_identical'] else 1


if __name__ == '__main__':
    if len(sys.argv) != 6:
        raise SystemExit(__doc__)
    raise SystemExit(main(*sys.argv[1:]))
