"""Completion evidence for swarm jobs. No acceptance from filename existence."""
import hashlib
import importlib.util
import json
from pathlib import Path

_spec = importlib.util.spec_from_file_location('swarm_job_validation', Path(__file__).with_name('swarm_validation.py'))
V = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(V)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def identity(job, binding):
    a = job['argv']
    V.need(len(a) == 8 and a[2] == job['output'], 'unsupported swarm command/output shape')
    return {'argv': a, 'binary_sha256': sha(a[0]), 'blob_sha256': sha(a[1]),
            'qualification_binding_sha256': binding}


def validate_output(path, job):
    data, digest = V.stable_json(path)
    a = job['argv']
    V.need(data.get('blob_sha256') == sha(a[1]), 'job blob mismatch')
    for key, expected in [('seed', int(a[7])), ('N_m3', float(a[4])),
                          ('T_gas_K', float(a[5])), ('electrons', int(a[6]))]:
        value = V.integer(data.get(key), key) if key in ('seed', 'electrons') else V.number(data.get(key), key)
        V.need(value == expected, 'job configuration mismatch: '+key)
    expected = [float(format(float(e), '.6g')) for e in a[3].split(',')]
    actual = [r['E_over_N_Td'] for r in data.get('points', [])]
    V.need(len(expected) == len(set(expected)) and actual == expected, 'job E/N coverage/order mismatch')
    for row in data['points']:
        for key in V.POINT_NUMBERS:
            V.number(row.get(key), key)
        for key in V.POINT_COUNTS:
            V.integer(row.get(key), key)
        counts = row.get('process_events')
        V.need(type(counts) is dict and counts, 'missing process counts')
        for value in counts.values(): V.integer(value, 'process count')
        V.need(sum(counts.values()) == row['real_events'], 'process total mismatch')
        V.need(row['real_events'] <= row['candidates'], 'event count exceeds candidates')
    return digest


def completed(job, binding):
    """Missing/invalid proof is incomplete, never a successful resume."""
    try:
        proof, _ = V.stable_json(job['output'] + '.complete.json')
        V.need(proof.get('schema') == 'swarm-job-completion-v1' and type(proof.get('exit')) is int and proof['exit'] == 0,
               'missing successful completion')
        V.need(proof.get('identity') == identity(job, binding), 'job provenance mismatch')
        V.need(proof.get('output_sha256') == validate_output(job['output'], job), 'output changed')
        return True
    except (OSError, ValueError, KeyError, TypeError):
        return False


def publish(temp_path, job, before, binding):
    """Called only after child exit 0 and under the output's exclusive claim."""
    V.need(before == identity(job, binding), 'binary/blob/command changed while job ran')
    digest = validate_output(temp_path, job)
    # Atomic rename publishes only a complete, verified JSON document.
    with open(temp_path, 'rb') as stream:
        __import__('os').fsync(stream.fileno())
    Path(temp_path).replace(job['output'])
    V.atomic_json(job['output'] + '.complete.json', {'schema':'swarm-job-completion-v1',
                  'identity':before, 'exit':0, 'output_sha256':digest})
