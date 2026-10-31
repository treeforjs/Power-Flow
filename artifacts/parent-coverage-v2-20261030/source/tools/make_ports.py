#!/usr/bin/env python3
"""Generate the coverage-v2 ports of the sealed v3 parent regression and of the
coverage-v1 component test.

Sealed files are read, never written. Every substitution is line-local, has an
exact expected count, and keeps the line count, so check names and order of
the sealed tests are preserved.

  test_v3_parent_regression_port_v2a.cpp  v3 regression through the v2 subflow,
      scope without high-energy bands (v1 coverage map). Must reproduce the
      sealed baseline: 442 identical check results, identical recorded event
      physics, 73/128 first events.
  test_v3_parent_regression_port_v2b.cpp  v3 regression through the v2 subflow
      with all v2 bands and the composite channel registration. The 8000 eV
      entry of the unsupported-energy loop is removed (that energy is covered
      in v2); the other 440 check results must equal the baseline's, with
      identical recorded event physics and first events.
  test_v1_coverage_port_v2.cpp  the coverage-v1 component test (171 checks)
      against the v2 subflow and scope without high-energy bands.
"""
import hashlib
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
V3_TEST = Path('/home/treefor/projects/runtimes/pic-water-common-commit-v3-20261030/tests/test_native_parent_subflow.cpp')
V1_TEST = Path('/home/treefor/projects/Power-Flow/artifacts/parent-coverage-v1-20261030/source/tests/test_parent_coverage.cpp')
OUTER3 = '/home/treefor/projects/runtimes/h2o-collision-data-20261029/data/collisions/h2o-gas-20261029-outer3.json'
OUTER3_SHA = 'c436020d3f2d07902fb059cfc0a269b2e914a3cbb2789300093d3631df896551'
KERNEL_BINDING = '/home/treefor/projects/runtimes/pic-water-parent-coverage-v2-20261030/kernel-binding-v2.json'

V3_COMMON = [
    ('#include "native_parent_subflow.hpp"', '#include "parent_coverage_v2_subflow.hpp"', 1),
    ('namespace w=raikou::pic_water;', 'namespace w=raikou::pic_water;namespace cv=raikou::native_water::coverage_v2;', 1),
    ('nw::ParentSubflow', 'cv::ParentSubflow', 12),
    ('"native_water_parent_subflow/v1"', 'cv::ParentSubflow::integration_owner()', 1),
    ('{"physical_trajectory_advanced",false},{"full_coupled_run",false}};',
     '{"physical_trajectory_advanced",false},{"full_coupled_run",false},'
     '{"coverage_scope_sha256",scope.identity()},{"coverage_scope",scope.description()}};', 1),
]
BRIDGE = 'w::ParentBridge bridge(data,w::parent_model,b.dataset_sha256,b.kernel_sha256,mn);'
V3_A = V3_COMMON + [
    (BRIDGE, BRIDGE + 'const cv::CoverageScope scope(data,b.dataset_sha256,bridge);', 1),
]
CHANNELS = ('std::vector<m::ArenaChannel> channels(const nw::Binding& b){return {{w::parent_channel,"electron","H2O",'
            '{"electron","electron","H2O+"},b.dataset_sha256,b.kernel_sha256,b.kernel_binding_file_sha256,w::parent_model}};}')
V3_B = V3_COMMON + [
    (BRIDGE, BRIDGE + 'auto bands=std::make_shared<const cv::HighEnergyBands>(cv::HighEnergyBands::read("' + OUTER3 + '","' +
     OUTER3_SHA + '",*data,b.dataset_sha256,"' + KERNEL_BINDING + '",b.kernel_sha256,b.kernel_binding_file_sha256));'
     'const cv::CoverageScope scope(data,b.dataset_sha256,bridge,bands);'
     'registered=scope.arena_binding(b.kernel_binding_file_sha256);', 1),
    (CHANNELS, 'cv::ArenaBinding registered;std::vector<m::ArenaChannel> channels(const nw::Binding&){return {{w::parent_channel,'
     '"electron","H2O",{"electron","electron","H2O+"},registered.dataset_sha256,registered.kernel_manifest_sha256,'
     'registered.kernel_binding_file_sha256,registered.model}};}', 1),
    ('for(double E:{13.,8000.}){', 'for(double E:{13.}){', 1),
]
V1_PORT = [
    ('#include "parent_coverage_subflow.hpp"', '#include "parent_coverage_v2_subflow.hpp"', 1),
    ('namespace cv = raikou::native_water::coverage_v1;', 'namespace cv = raikou::native_water::coverage_v2;', 1),
]
CONSTRUCTOR_LINES = 11


def substitute(text, table):
    for old, new, count in table:
        found = text.count(old)
        if found != count:
            raise SystemExit(f'expected {count} occurrences of {old[:70]!r}, found {found}')
        text = text.replace(old, new)
    return text


def v3_port(text, table):
    text = substitute(text, table)
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


def pinned(path, pin_file):
    data = path.read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    if digest != (HERE / pin_file).read_text().strip():
        raise SystemExit(f'sealed test changed: {path}: {digest}')
    return data.decode()


def outputs():
    v3 = pinned(V3_TEST, 'V3-TEST-SHA256')
    v1 = pinned(V1_TEST, 'V1-TEST-SHA256')
    out = {
        HERE / 'tests/test_v3_parent_regression_port_v2a.cpp': v3_port(v3, V3_A),
        HERE / 'tests/test_v3_parent_regression_port_v2b.cpp': v3_port(v3, V3_B),
        HERE / 'tests/test_v1_coverage_port_v2.cpp': substitute(v1, V1_PORT),
    }
    for path, text in out.items():
        source = v1 if 'v1_coverage' in path.name else v3
        if text.count('\n') != source.count('\n'):
            raise SystemExit(f'port changed the line count: {path.name}')
    return out


def main():
    out = outputs()
    if len(sys.argv) > 1 and sys.argv[1] == '--check':
        for path, text in out.items():
            if path.read_text() != text:
                raise SystemExit(f'committed port differs from the generated port: {path.name}')
        print('ports match generator')
        return
    for path, text in out.items():
        path.write_text(text)
        print(f'wrote {path}')


if __name__ == '__main__':
    main()
