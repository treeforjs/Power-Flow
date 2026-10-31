#!/usr/bin/python -I
"""Compare the coordinator v5 T1 step at the Stage 1 positions (z = 1 mm) with the Stage 1 (v4) receipts.

  compare_stage1.py STAGE1_T1_DIR V5_T1_DIR OUTPUT_JSON [--native-v5]

With --native-v5 (coordinator v5 linked against native v5, whose own identities and export format changed
by design) the native save export and native token audit are not compared structurally; the comparison
covers Hall field words and charge ledger, heat slices, the Hall motion record (endpoints, paths, budget),
the chemistry checkpoint values and the T1 numeric evidence.

Physical state must be bit for bit identical: Hall field words E, B, E_low, B_low, clock and revision,
the complete charge ledger (free net words and gross, conductor and exterior banks), every heat slice
file, every native particle, packet, count and ledger value, the chemistry checkpoint values, the motion
record numbers (endpoints, paths, budget) and the T1 receipt numbers. Every remaining difference must be a
64-hex digest string; such differences form the old -> new identity mapping, which must be a consistent
one-to-one map (the same old digest always maps to the same new digest). A short list of record schema
texts changed intentionally and is reported separately.
"""
import json, re, sys
from pathlib import Path

HEX = re.compile(r'^[0-9a-f]{64}$')
INTENTIONAL_TEXT = {('hall-compensated-field-owner-v1', 'hall-compensated-field-owner-v2'),
                    ('closed-parent-corrective-common-step-v4', 'closed-parent-corrective-common-step-v5'),
                    ('closed-parent-joint-genesis-v4', 'closed-parent-joint-genesis-v5')}
NEW_KEYS = {'field_solve_backend', 'field_solve_backend_descriptor_sha256', 'motion_provider', 'transport_record_v3_sha256'}
EMBEDDED_JSON = {'chemistry_checkpoint'}


class Diff:
    def __init__(self):
        self.mapping, self.roles, self.issues, self.texts, self.added = {}, {}, [], set(), set()

    def walk(self, a, b, path):
        if isinstance(a, dict) and isinstance(b, dict):
            for k in a:
                if k not in b:
                    self.issues.append(['missing key', path + '/' + k])
                elif k in EMBEDDED_JSON and isinstance(a[k], str) and isinstance(b[k], str):
                    self.walk(json.loads(a[k]), json.loads(b[k]), path + '/' + k)
                else:
                    self.walk(a[k], b[k], path + '/' + k)
            for k in b:
                if k not in a:
                    if k in NEW_KEYS:
                        self.added.add(k)
                    else:
                        self.issues.append(['added key', path + '/' + k])
        elif isinstance(a, list) and isinstance(b, list):
            if len(a) != len(b):
                self.issues.append(['length', path, len(a), len(b)])
                return
            for i, (x, y) in enumerate(zip(a, b)):
                self.walk(x, y, '%s[%d]' % (path, i))
        elif isinstance(a, str) and isinstance(b, str) and HEX.match(a) and HEX.match(b):
            if a != b:
                self.mapping.setdefault(a, set()).add(b)
                self.roles.setdefault(a, set()).add(re.sub(r'\[\d+\]', '[]', path))
        elif type(a) is not type(b) or a != b:
            if isinstance(a, str) and (a, b) in INTENTIONAL_TEXT:
                self.texts.add((a, b))
            else:
                self.issues.append(['value', path, repr(a)[:200], repr(b)[:200]])


def load(p):
    return json.loads(Path(p).read_text())


def main():
    old, new, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    native_v5 = '--native-v5' in sys.argv[4:]
    if native_v5:
        NEW_KEYS.add('row_sha256')   # native v5 typed journal row digest in the exposure audits
    d = Diff()
    report = {'schema': 'stage1-v5-identity-mapping-v1', 'stage1_dir': str(old), 'v5_dir': str(new), 'populations': {},
              'native': 'native v5 (native identities and export not compared structurally)' if native_v5 else 'candidate003'}
    ok = True
    for o_label, n_label in (('P1', 'P1-stage1-z1mm'), ('P2', 'P2-stage1-z1mm')):
        pop = {}
        ho, hn = old / (o_label + '-after-owner-bytes/hall-save'), new / (n_label + '-after-owner-bytes/hall-save')
        co, cn = load(ho / 'checkpoint.json'), load(hn / 'checkpoint.json')
        pop['hall_field_words_identical'] = co['field'] == cn['field']
        pop['hall_charge_ledger_identical'] = co['charge'] == cn['charge']
        heat = sorted(p.name for p in (ho / 'heat').iterdir())
        pop['heat_slices_byte_identical'] = heat == sorted(p.name for p in (hn / 'heat').iterdir()) and all(
            (ho / 'heat' / f).read_bytes() == (hn / 'heat' / f).read_bytes() for f in heat)
        pop['heat_slices'] = len(heat)
        for key in ('configuration_sha256', 'field_sha256', 'heat_sha256', 'motion_mesh_sha256', 'schema'):
            d.walk(co[key], cn[key], n_label + '/hall-checkpoint/' + key)
        if not native_v5:
            d.walk(load(old / (o_label + '-after-owner-bytes/native-save.json')), load(new / (n_label + '-after-owner-bytes/native-save.json')),
                   n_label + '/native-save')
            d.walk(load(old / (o_label + '-trial-audit.json')), load(new / (n_label + '-trial-audit.json')), n_label + '/trial-audit')
        else:
            ao, an = load(old / (o_label + '-trial-audit.json')), load(new / (n_label + '-trial-audit.json'))
            for key in ('transport', 'first_exposure', 'second_exposure', 'hall_iterations', 'field_heat_balance_J', 'motion_midpoint_sha256',
                        'midpoint_role_sha256', 'candidate_field_sha256', 'candidate_heat_sha256', 'begin_s', 'end_s', 'authoritative_dt_s'):
                d.walk(ao[key], an[key], n_label + '/trial-audit/' + key)
        ro, rn = load(old / 't1-receipt.json')[o_label], load(new / 't1-receipt.json')[n_label]
        keys = ('genesis', 'before', 'after', 'observed_before', 'observed_after', 'motion_budget', 'midpoint_compensation',
                  'combined_field_motion_heat_balance_J', 'hall_iterations', 'parent_events', 'publish', 'candidate_field_sha256',
                  'candidate_heat_sha256', 'motion_midpoint_sha256', 'midpoint_role_seal', 'trial_audit_sha256', 'position_m')
        if native_v5:
            keys = tuple(k for k in keys if k not in ('observed_before', 'observed_after', 'trial_audit_sha256'))
        for k in keys:
            if k in ro:
                if k == 'genesis':
                    ro_g = dict(ro[k])
                    rn_g = dict(rn[k])
                    if native_v5:
                        for x in (ro_g, rn_g):
                            x.pop('native_state_prefix', None)
                    d.walk(ro_g, rn_g, n_label + '/t1-receipt/' + k)
                elif k in ('before', 'after'):
                    a, b = dict(ro[k]), dict(rn[k])
                    for x in (a, b):
                        x.pop('hall_save_sha256', None)   # Hall save file digests: checkpoint text changed (schema, digests)
                        if native_v5:
                            for key in ('native_save_bytes', 'native_save_sha256', 'native_boundary', 'native_state_prefix', 'digest'):
                                x.pop(key, None)
                    d.walk(a, b, n_label + '/t1-receipt/' + k)
                elif k == 'publish':
                    a = {x: y for x, y in ro[k].items() if not x.startswith('allocations_') or x == 'allocations_between_first_and_last_swap'}
                    b = {x: y for x, y in rn[k].items() if not x.startswith('allocations_') or x == 'allocations_between_first_and_last_swap'}
                    d.walk(a, b, n_label + '/t1-receipt/' + k)
                else:
                    d.walk(ro[k], rn[k], n_label + '/t1-receipt/' + k)
        ok = ok and pop['hall_field_words_identical'] and pop['hall_charge_ledger_identical'] and pop['heat_slices_byte_identical']
        report['populations'][n_label] = pop
    inconsistent = {k: sorted(v) for k, v in d.mapping.items() if len(v) != 1}
    inverse = {}
    for k, v in d.mapping.items():
        for x in v:
            inverse.setdefault(x, set()).add(k)
    collisions = {k: sorted(v) for k, v in inverse.items() if len(v) != 1}
    report['identity_mapping'] = sorted(({'old': k, 'new': sorted(v)[0], 'roles': sorted(d.roles[k])[:12], 'role_count': len(d.roles[k])}
                                         for k, v in d.mapping.items()), key=lambda r: r['roles'][0])
    report['mapping_size'] = len(d.mapping)
    report['mapping_inconsistent'] = inconsistent
    report['mapping_collisions'] = collisions
    report['intentional_text_changes'] = sorted(d.texts)
    report['added_keys'] = sorted(d.added)
    report['non_digest_differences'] = d.issues[:200]
    report['non_digest_difference_count'] = len(d.issues)
    ok = ok and not d.issues and not inconsistent and not collisions
    report['physical_state_bit_identical_and_mapping_consistent'] = ok
    out.write_text(json.dumps(report, indent=1) + '\n')
    print(json.dumps({'ok': ok, 'mapping_size': len(d.mapping), 'non_digest_differences': len(d.issues),
                      'populations': report['populations'], 'first_issues': d.issues[:5]}, indent=1))
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main())
