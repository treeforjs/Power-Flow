"""Strict input completeness for the observed v3.2 swarm comparison.

This validates the existing scientific rule; it does not authorize that rule or
qualify the underlying cross-section models. No input/output run files mutate.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile


class Invalid(ValueError):
    pass


def need(condition, message):
    if not condition:
        raise Invalid(message)


def pairs(items):
    result = {}
    for key, value in items:
        need(key not in result, 'duplicate JSON key: ' + key)
        result[key] = value
    return result


def reject_constant(value):
    raise Invalid('nonfinite JSON constant: ' + value)


def finite_tree(value):
    if isinstance(value, float):
        need(math.isfinite(value), 'nonfinite number')
    elif isinstance(value, dict):
        for x in value.values():
            finite_tree(x)
    elif isinstance(value, list):
        for x in value:
            finite_tree(x)


def decode(raw):
    result = json.loads(raw, object_pairs_hook=pairs, parse_constant=reject_constant)
    finite_tree(result)
    return result


def stable_json(path):
    """Detect concurrent replacement/truncation during the read; never accept a prefix."""
    path = Path(path)
    with path.open('rb') as stream:
        before = os.fstat(stream.fileno())
        raw = stream.read()
        after = os.fstat(stream.fileno())
    final = path.stat()
    identity = lambda s: (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns)
    need(identity(before) == identity(after) == identity(final), str(path) + ': changed while reading')
    return decode(raw), hashlib.sha256(raw).hexdigest()


def number(x, name):
    need(type(x) in (int, float) and math.isfinite(x), name + ': expected finite number')
    return x


def integer(x, name):
    need(type(x) is int and x >= 0, name + ': expected nonnegative integer')
    return x


def statistics(values):
    try:
        for value in values:
            number(value, 'relative sample')
        mean = number(sum(values)/len(values), 'sample mean')
        sd = number(math.sqrt(sum((v-mean)**2 for v in values)/(len(values)-1)), 'sample SD')
        return mean, sd
    except (OverflowError, ZeroDivisionError) as exc:
        raise Invalid('unrepresentable sample statistics') from exc


class Contract:
    def __init__(self, path, root):
        self.data, self.sha256 = stable_json(path)
        self.root = Path(root)
        c = self.data
        need(c['schema'] == 'swarm-completeness-contract-v1', 'contract schema')
        need(c['seeds'] == list(range(20261101, 20261107)), 'contract seed plan')
        need(c['projectile'] == 'e' and c['target'] == 'H2O', 'contract species')
        need(c['criterion'] == {'formula': 'abs(new_mean-ref_mean) <= max(2*ref_sample_SD,0.01)',
                              'reference_sd_multiplier': 2.0, 'absolute_relative_floor': 0.01}, 'criterion changed')
        need(c['raw_energy_format'] == '.6g', 'raw energy serialization changed')
        self.measured = {}
        for rel, expected in c['measured_files'].items():
            data, actual = stable_json(self.root / rel)
            need(actual == expected, 'measured-source hash differs: ' + rel)
            need(len(data['x']) == len(data['y']) and len(data['x']) == len(set(data['x'])), 'measured rows')
            if 'mobility-N' in rel:
                for q in ('W_bulk_m_s', 'W_flux_m_s'):
                    self.measured[q] = {x: y*x*1e-21 for x, y in zip(data['x'], data['y'])}
            else:
                q = 'alpha_minus_eta_over_N_m2' if 'alpha-minus-eta' in rel else 'DT_over_mu_eV'
                self.measured[q] = dict(zip(data['x'], data['y']))
        need(set(self.measured) == set(METRICS), 'metric source set')
        self.energies = sorted({x for values in self.measured.values() for x in values})
        self.serialized = {float(format(x, c['raw_energy_format'])): x for x in self.energies}
        need(len(self.energies) == 36 and len(self.serialized) == 36, 'ambiguous/incomplete energy plan')
        self.rows = {(q, x) for q, values in self.measured.items() for x in values}
        self.judged = sum(judged(q, x) for q, x in self.rows)
        need(self.judged > 0, 'zero judged rows')
        self.processes = {}
        for label, spec in c['sets'].items():
            data, actual = stable_json(self.root / spec['spec_path'])
            need(actual == spec['spec_sha256'], 'spec-source hash differs: ' + label)
            groups = [g for g in data['groups'] if g['projectile'] == c['projectile'] and g['target'] == c['target']]
            need(len(groups) == 1, 'species group missing or duplicated')
            names = [p['id'] for p in groups[0]['processes']]
            need(len(names) == len(set(names)), 'duplicate process identity')
            self.processes[label] = set(names)
        need(set(self.processes) == {'v3', 'v32'}, 'set identities')
        for rel, expected in c['source_files'].items():
            need(hashlib.sha256((self.root / rel).read_bytes()).hexdigest() == expected,
                 'writer source differs: ' + rel)


METRICS = ('W_bulk_m_s', 'W_flux_m_s', 'alpha_minus_eta_over_N_m2', 'DT_over_mu_eV')
POINT_NUMBERS = (*METRICS, 'mean_energy_eV', 'nu_ion_s', 'nu_att_s', 'alpha_over_N_m2', 'eta_over_N_m2',
                 'k_ion_m3_s', 'k_att_m3_s', 'mobility_N_1_m_V_s', 'Ecap_eV', 'nu_max_s', 't_relax_s', 't_measure_s')
POINT_COUNTS = ('candidates', 'real_events', 'candidates_above_Ecap', 'candidates_below_grid')


def judged(quantity, energy):
    return quantity in ('W_bulk_m_s', 'W_flux_m_s') or (quantity == 'alpha_minus_eta_over_N_m2' and energy >= 198.0)


def validate_run(run, label, seed, contract):
    c = contract.data
    need(type(run) is dict, 'run must be object')
    need(run.get('blob_sha256') == c['sets'][label]['blob_sha256'], 'wrong blob identity')
    need(integer(run.get('seed'), 'seed') == seed, 'wrong seed identity')
    for key in ('N_m3', 'T_gas_K', 'electrons'):
        value = integer(run.get(key), key) if key == 'electrons' else number(run.get(key), key)
        need(value == c[key], 'wrong ' + key)
    need(type(run.get('points')) is list and len(run['points']) == len(contract.energies), 'incomplete/extra energy rows')
    points = {}
    for row in run['points']:
        energy = number(row.get('E_over_N_Td'), 'E/N')
        need(energy in contract.serialized, 'unexpected serialized E/N: ' + str(energy))
        need(energy not in points, 'duplicate E/N: ' + str(energy))
        for q in POINT_NUMBERS:
            number(row.get(q), q)
        for q in POINT_COUNTS:
            integer(row.get(q), q)
        counts = row.get('process_events')
        need(type(counts) is dict and set(counts) == contract.processes[label], 'species/process set mismatch')
        for name, count in counts.items():
            integer(count, name)
        need(sum(counts.values()) == row['real_events'], 'process total mismatch')
        need(row['real_events'] <= row['candidates'], 'event count exceeds candidates')
        points[energy] = row
    need(set(points) == set(contract.serialized), 'missing expected E/N')
    return {contract.serialized[e]: row for e, row in points.items()}


def summarize(directory, label, contract):
    need(label in contract.data['sets'], 'unknown set')
    directory = Path(directory)
    expected = {f'{label}-{s}.json': s for s in contract.data['seeds']}
    # Completion receipts are separate artifacts, not additional seed runs.
    actual = {p.name for p in directory.glob(label + '-*.json') if not p.name.endswith('.json.complete.json')}
    need(actual == set(expected), 'missing/extra seed files: ' + repr(sorted(actual ^ set(expected))))
    runs, evidence = [], []
    for name, seed in expected.items():
        run, digest = stable_json(directory / name)
        runs.append(validate_run(run, label, seed, contract))
        evidence.append({'file': name, 'sha256': digest})
    rows = []
    for q, values in contract.measured.items():
        for energy, measured in sorted(values.items()):
            need(number(measured, 'measured') != 0, 'zero measured denominator')
            vals = [r[energy][q]/measured-1 for r in runs]
            mean, sd = statistics(vals)
            rows.append({'quantity': q, 'E_over_N_Td': energy, 'measured': measured, 'n': len(vals),
                         'mean_rel': mean, 'sd_rel': sd, 'values_rel': vals})
    result = {'schema': 'swarm-complete-summary-v1', 'contract_sha256': contract.sha256, 'set': label,
              'projectile': contract.data['projectile'], 'target': contract.data['target'],
              'blobs': [contract.data['sets'][label]['blob_sha256']], 'seeds': contract.data['seeds'],
              'input_files': evidence, 'rows': rows, 'complete': True,
              'source_note': 'Raw files do not bind executable hash or successful process exit; separate run receipt required.'}
    validate_summary(result, label, contract)
    return result


def validate_summary(summary, label, contract):
    need(summary.get('schema') == 'swarm-complete-summary-v1' and summary.get('complete') is True, 'unqualified summary schema')
    need(summary.get('contract_sha256') == contract.sha256 and summary.get('set') == label, 'summary contract/set mismatch')
    need(summary.get('projectile') == 'e' and summary.get('target') == 'H2O', 'summary species mismatch')
    need(summary.get('seeds') == contract.data['seeds'], 'summary seed set mismatch')
    need(summary.get('blobs') == [contract.data['sets'][label]['blob_sha256']], 'summary blob mismatch')
    inputs = summary.get('input_files', [])
    need(len(inputs) == len(contract.data['seeds']), 'summary input-file count')
    need({v['file'] for v in inputs} == {f'{label}-{s}.json' for s in contract.data['seeds']}, 'summary input files')
    for entry in inputs:
        need(len(entry['sha256']) == 64 and all(c in '0123456789abcdef' for c in entry['sha256']), 'summary input hash')
    actual = set()
    for row in summary['rows']:
        key = row['quantity'], row['E_over_N_Td']
        need(key in contract.rows and key not in actual, 'unexpected or duplicate summary row')
        actual.add(key)
        need(row['measured'] == contract.measured[key[0]][key[1]], 'summary measured data mismatch')
        vals = row['values_rel']
        need(type(vals) is list and len(vals) == len(contract.data['seeds']) and row['n'] == len(vals), 'summary sample completeness')
        for value in vals:
            number(value, 'summary sample')
        mean, sd = statistics(vals)
        need(number(row['mean_rel'], 'summary mean') == mean and number(row['sd_rel'], 'summary SD') == sd, 'summary statistics mismatch')
    need(actual == contract.rows and sum(judged(*k) for k in actual) == contract.judged, 'incomplete/zero judged rows')


def validate_comparison(result, new, ref, contract):
    """Recompute the scientific-rule fields from current, complete raw summaries."""
    validate_summary(new, 'v32', contract)
    validate_summary(ref, 'v3', contract)
    need(result.get('contract_sha256') == contract.sha256 and result.get('complete') is True, 'comparison contract')
    need(result.get('raw_input_files') == {'v3':ref['input_files'],'v32':new['input_files']}, 'comparison raw-output identities differ')
    labels = {'W_bulk_m_s':'bulk drift','W_flux_m_s':'flux drift',
              'alpha_minus_eta_over_N_m2':'(alpha - eta)/N','DT_over_mu_eV':'D_T/mu'}
    refs = {(r['quantity'],r['E_over_N_Td']):r for r in ref['rows']}
    expected = {}
    failing = 0
    for row in new['rows']:
        key = row['quantity'],row['E_over_N_Td']; other = refs[key]
        tolerance = max(2*other['sd_rel'],0.01)
        difference = row['mean_rel']-other['mean_rel']; judge = judged(*key)
        ok = abs(difference) <= tolerance
        failing += int(judge and not ok)
        expected[(labels[key[0]],key[1])] = {'measured':row['measured'],'judged':judge,
            'v3_mean_rel':other['mean_rel'],'v3_sd_rel':other['sd_rel'],'v3_n':other['n'],
            'v32_mean_rel':row['mean_rel'],'v32_sd_rel':row['sd_rel'],'v32_n':row['n'],
            'tolerance':tolerance,'difference':difference,'pass':ok if judge else None}
    rows = result.get('rows',[])
    need(len(rows) == len(expected), 'comparison row count')
    seen = set()
    for row in rows:
        key = row['quantity'],row['E_over_N_Td']
        need(key in expected and key not in seen, 'comparison duplicate/unknown row')
        seen.add(key)
        need(all(row.get(k) == v for k,v in expected[key].items()), 'comparison values differ from current raw outputs')
    need(result.get('judged') == contract.judged and result.get('failing') == failing and
         result.get('passed') is (failing == 0), 'comparison disposition differs')


def atomic_json(path, value):
    """Publish a complete validation artifact. Input simulation files are untouched."""
    path = Path(path)
    fd, temporary = tempfile.mkstemp(prefix='.'+path.name+'.', dir=path.parent)
    try:
        with os.fdopen(fd, 'w') as stream:
            json.dump(value, stream, indent=1, allow_nan=False)
            stream.write('\n')
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
        directory_fd = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def inspect(directory, contract):
    result = {'schema': 'swarm-completeness-inspection-v1', 'contract_sha256': contract.sha256,
              'scientific_criterion_status': contract.data['scientific_criterion_status'], 'sets': {}}
    for label in contract.data['sets']:
        try:
            summary = summarize(directory, label, contract)
            result['sets'][label] = {'status': 'complete', 'runs': len(summary['seeds']), 'rows': len(summary['rows'])}
        except (Invalid, OSError, ValueError, KeyError, TypeError) as exc:
            result['sets'][label] = {'status': 'incomplete_or_invalid', 'reason': str(exc)}
    result['complete'] = all(v['status'] == 'complete' for v in result['sets'].values())
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('directory')
    parser.add_argument('--contract', required=True)
    parser.add_argument('--source-root', required=True)
    args = parser.parse_args()
    answer = inspect(args.directory, Contract(args.contract, args.source_root))
    print(json.dumps(answer, indent=2, allow_nan=False))
    raise SystemExit(0 if answer['complete'] else 2)
