#!/usr/bin/env python3
"""Read-only CPU verification. Never grants GPU dispatch or touches production."""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
V2 = HERE.parent
FIXTURE = Path('/home/treefor/sim-workspace/powerflow/outputs/pic-profile/cuda-gather-investigation/actual-midpoint-fixture')
REVIEW = Path('/home/treefor/sim-workspace/powerflow/outputs/pic-profile/cuda-gather-independent-review/backend-v2/disposition.json')

def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()

def check(mapping, root=None):
    failures = []
    for name, expected in mapping.items():
        path = Path(name) if root is None else root / name
        try:
            actual = digest(path)
        except OSError as error:
            failures.append({'path': str(path), 'error': str(error)})
            continue
        if actual != expected:
            failures.append({'path': str(path), 'expected': expected, 'actual': actual})
    return failures

def main():
    review = json.loads(REVIEW.read_text())
    errors = {
        'harness_and_dependencies': check(json.loads((HERE / 'artifact-sha256.json').read_text())),
        'accepted_v2_source': check(review['accepted_source_hashes'], V2),
        'v1_preservation': check(json.loads((V2 / 'v1-baseline-sha256.json').read_text())),
        'actual_fixture': check(json.loads((FIXTURE / 'sha256.json').read_text()), FIXTURE),
    }
    if review['status'] != 'accepted_for_parent_preparation_of_bounded_GPU_lifecycle_qualification':
        errors['core_review'] = ['Unexpected review disposition']
    if digest(V2 / 'REVIEW-REQUEST.json') != review['review_request_sha256']:
        errors['core_review_request'] = ['Accepted core review request changed']
    passed = not any(errors.values())
    print(json.dumps({'passed': passed, 'hash_mismatches': errors,
                      'GPU_executed': False, 'GPU_dispatch_authorized_by_this_check': False,
                      'production_actions': False}, indent=2))
    return 0 if passed else 1

if __name__ == '__main__':
    raise SystemExit(main())
