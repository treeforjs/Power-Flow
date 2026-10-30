#!/usr/bin/env python3
"""Generate the coverage-v1 port of the sealed v3 parent regression test.

The sealed file is read, never written. Every substitution is line-local, has
an exact expected count, and keeps the line count, so every check name and
check order of the sealed test are preserved. With the default coverage policy
the port must reproduce the sealed regression (run_parent_regression.py
semantics): 442 identical check results, identical recorded event physics and
the same 73/128 first-event count.
"""
import hashlib
import sys
from pathlib import Path

SEALED = Path('/home/treefor/projects/runtimes/pic-water-common-commit-v3-20261030/tests/test_native_parent_subflow.cpp')
PORT = Path(__file__).resolve().parent / 'tests/test_v3_parent_regression_port.cpp'

# (old, new, expected occurrences)
SUBSTITUTIONS = [
    ('#include "native_parent_subflow.hpp"', '#include "parent_coverage_subflow.hpp"', 1),
    ('namespace w=raikou::pic_water;', 'namespace w=raikou::pic_water;namespace cv=raikou::native_water::coverage_v1;', 1),
    ('nw::ParentSubflow', 'cv::ParentSubflow', 12),
    ('w::ParentBridge bridge(data,w::parent_model,b.dataset_sha256,b.kernel_sha256,mn);',
     'w::ParentBridge bridge(data,w::parent_model,b.dataset_sha256,b.kernel_sha256,mn);'
     'const cv::CoverageScope scope(data,b.dataset_sha256,bridge);', 1),
    ('"native_water_parent_subflow/v1"', 'cv::ParentSubflow::integration_owner()', 1),
    ('{"physical_trajectory_advanced",false},{"full_coupled_run",false}};',
     '{"physical_trajectory_advanced",false},{"full_coupled_run",false},'
     '{"coverage_scope_sha256",scope.identity()},{"coverage_scope",scope.description()}};', 1),
]
# Constructor calls gain the scope argument right after the bridge.
CONSTRUCTOR_LINES = 11


def generate(text: str) -> str:
    for old, new, count in SUBSTITUTIONS:
        found = text.count(old)
        if found != count:
            raise SystemExit(f'expected {count} occurrences of {old!r}, found {found}')
        text = text.replace(old, new)
    lines = text.split('\n')
    changed = 0
    for i, line in enumerate(lines):
        if 'cv::ParentSubflow ' in line and ',bridge,' in line:
            if line.count(',bridge,') != 1:
                raise SystemExit(f'ambiguous constructor line {i + 1}')
            lines[i] = line.replace(',bridge,', ',bridge,scope,')
            changed += 1
    if changed != CONSTRUCTOR_LINES:
        raise SystemExit(f'expected {CONSTRUCTOR_LINES} constructor lines, found {changed}')
    return '\n'.join(lines)


def main() -> None:
    sealed = SEALED.read_bytes()
    digest = hashlib.sha256(sealed).hexdigest()
    pinned = (Path(__file__).resolve().parent / 'V3-TEST-SHA256').read_text().strip()
    if digest != pinned:
        raise SystemExit(f'sealed v3 test changed: {digest}')
    port = generate(sealed.decode())
    if port.count('\n') != sealed.decode().count('\n'):
        raise SystemExit('port changed the line count')
    if len(sys.argv) > 1 and sys.argv[1] == '--check':
        if PORT.read_text() != port:
            raise SystemExit('committed port differs from the generated port')
        print('port matches generator')
        return
    PORT.write_text(port)
    print(f'wrote {PORT} from sealed {digest}')


if __name__ == '__main__':
    main()
