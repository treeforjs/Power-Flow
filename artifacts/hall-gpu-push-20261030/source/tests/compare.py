"""Bitwise comparison of HGPENC1 encodings (tests/encoding.hpp).

compare(reference, candidate, subset=False) compares every array of the reference
with the candidate byte for byte; with subset=False the entry sets must also be
equal. The first differing element of each array is reported with both values
(decimal and IEEE bit pattern) and, for interleaved (x, y, z) arrays, the carrier
rank and component.
"""
import hashlib
import struct
import sys
import json

INTERLEAVED = {'endpoint.position_m', 'endpoint.proper_velocity_m_s', 'path.field_start_m', 'path.field_displacement_m'}
SIZE = {'d': 8, 'q': 8, 'i': 4, 's': 1}


def read(path):
    out = {}
    with open(path, 'rb') as f:
        data = f.read()
    if data[:8] != b'HGPENC1\0':
        raise ValueError('bad magic: ' + path)
    pos = 8
    (count,) = struct.unpack_from('<Q', data, pos)
    pos += 8
    for _ in range(count):
        (length,) = struct.unpack_from('<Q', data, pos)
        pos += 8
        name = data[pos:pos + length].decode()
        pos += length
        kind = chr(data[pos])
        pos += 1
        (elements,) = struct.unpack_from('<Q', data, pos)
        pos += 8
        size = SIZE[kind] * elements
        out[name] = (kind, elements, data[pos:pos + size])
        pos += size
    if pos != len(data):
        raise ValueError('trailing bytes in ' + path)
    return out


def _value(kind, raw, index):
    if kind == 'd':
        (bits,) = struct.unpack_from('<Q', raw, 8 * index)
        (v,) = struct.unpack_from('<d', raw, 8 * index)
        return {'value': repr(v), 'bits': '%016x' % bits}
    if kind == 'q':
        return {'value': struct.unpack_from('<Q', raw, 8 * index)[0]}
    if kind == 'i':
        return {'value': struct.unpack_from('<i', raw, 4 * index)[0]}
    return {'value': raw[index:index + 1].decode(errors='replace')}


def first_difference(name, a, b):
    kind_a, count_a, raw_a = a
    kind_b, count_b, raw_b = b
    if kind_a != kind_b or count_a != count_b:
        return {'kind': [kind_a, kind_b], 'count': [count_a, count_b]}
    if kind_a == 's':
        return {'reference': raw_a.decode(errors='replace'), 'candidate': raw_b.decode(errors='replace')}
    size = SIZE[kind_a]
    differing = 0
    first = None
    for i in range(count_a):
        if raw_a[size * i:size * (i + 1)] != raw_b[size * i:size * (i + 1)]:
            differing += 1
            if first is None:
                first = i
    report = {'index': first, 'differing_elements': differing, 'reference': _value(kind_a, raw_a, first),
              'candidate': _value(kind_b, raw_b, first)}
    if name in INTERLEAVED:
        report['carrier_rank'] = first // 3
        report['component'] = 'xyz'[first % 3]
    return report


def compare(reference_path, candidate_path, subset=False):
    ref = read(reference_path)
    cand = read(candidate_path)
    result = {'reference': reference_path, 'candidate': candidate_path, 'subset': subset, 'arrays': {}, 'missing': [], 'extra': []}
    for name, entry in ref.items():
        if name not in cand:
            result['missing'].append(name)
            continue
        same = entry == cand[name]
        item = {'identical': same, 'type': entry[0], 'elements': entry[1], 'sha256': hashlib.sha256(entry[2]).hexdigest()}
        if not same:
            item['candidate_sha256'] = hashlib.sha256(cand[name][2]).hexdigest()
            item['first_difference'] = first_difference(name, entry, cand[name])
        result['arrays'][name] = item
    if not subset:
        result['extra'] = sorted(set(cand) - set(ref))
    result['identical'] = not result['missing'] and not result['extra'] and all(v['identical'] for v in result['arrays'].values())
    result['differing_arrays'] = sorted(k for k, v in result['arrays'].items() if not v['identical'])
    return result


def manifest(path):
    return {name: {'type': e[0], 'elements': e[1], 'sha256': hashlib.sha256(e[2]).hexdigest()} for name, e in read(path).items()}


if __name__ == '__main__':
    subset = '--subset' in sys.argv
    args = [a for a in sys.argv[1:] if a != '--subset']
    print(json.dumps(compare(args[0], args[1], subset), indent=1))
