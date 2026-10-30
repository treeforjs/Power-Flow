"""CPU-only analysis after a parent-controlled native diagnostic replay."""
import array
import json
import os
from pathlib import Path
import sys

assert os.sched_getaffinity(0) == {6} and os.getpriority(os.PRIO_PROCESS, 0) >= 15
prefix = sys.argv[1]
def load(suffix, kind):
    values = array.array(kind)
    with open(prefix + suffix, 'rb') as stream:
        values.frombytes(stream.read())
    return values
rows = load('.rows.i32', 'i')
cols = load('.cols.i32', 'i')
permutation = load('.permutation.i32', 'i')
n = len(rows) - 1
assert len(permutation) == n and sorted(permutation) == list(range(n))
assert rows[0] == 0 and rows[-1] == len(cols)
inverse = [0] * n
for new, old in enumerate(permutation):
    inverse[old] = new
bandwidth = distance_sum = edges = 0
for old in range(n):
    for k in range(rows[old], rows[old + 1]):
        other = cols[k]
        if other == old:
            continue
        distance = abs(inverse[old] - inverse[other])
        bandwidth = max(bandwidth, distance)
        distance_sum += distance
        edges += 1
record = json.loads(Path(prefix + '.json').read_text())
result = {'rows': n, 'stored_upper_nnz': len(cols), 'permutation_valid': True,
          'permuted_index_bandwidth': bandwidth,
          'mean_permuted_index_distance_over_upper_offdiagonals': distance_sum / edges,
          'original_index_bandwidth': record['original_index_bandwidth'],
          'library_lu_nnz': record['lu_nnz'], 'library_lu_nnz_status': record['lu_nnz_status'],
          'limitations': ['Index bandwidth is not elimination-tree depth or measured solve parallelism.',
                          'Library LU_NNZ is retained with its native SPD semantics; no Cholesky-count conversion is assumed.']}
Path(prefix + '.permutation-analysis.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
