"""Apply the reserved legacy-solver retirement and freeze its exact scope."""
from pathlib import Path
import ast
import hashlib
import json
import re

ROOT = Path('/workspace/Power-Flow')
WORK = Path('/workspace/research_work/mitl-github-push-20261020T0601/mhd-retirement')
REMOVED = (
    'mhdlab/mhd.py', 'mhdlab/materials.py', 'mhdlab/conductivity_tables.py',
    'MHD_heating.py', 'tools/build_conductivity_table.py',
    'configs/mykonos.yaml', 'configs/mykonos_preview.yaml',
    'tests/test_mhd.py', 'tests/test_smoke_run.py', 'tests/test_conductivity_tables.py',
    'data/conductivity/README.md', 'data/conductivity/ss304_knoepfel_seed.csv',
)
MODIFIED = ('mhdlab/runner.py', 'mhdlab/__init__.py', 'tests/test_cr_data.py')
CREATED = ('configs/cr_reference_library.yaml', 'tests/test_legacy_mhd_retirement.py')
HELPERS = ('_grid_shape_from_config', '_thermionic_electron_energy_ev')


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def inventory(paths):
    return {str(path.relative_to(ROOT)): digest(path.read_bytes()) for path in sorted(set(paths))
            if path.is_file() and '__pycache__' not in path.parts}


assert not (WORK/'source-freeze.json').exists()
assert all((ROOT/path).is_file() for path in REMOVED + MODIFIED)
assert all(not (ROOT/path).exists() for path in CREATED)
before = {path: (ROOT/path).read_bytes() for path in REMOVED + MODIFIED}
for path, raw in before.items():
    destination = WORK/'before-files'/path
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(raw)

owned = set(REMOVED + MODIFIED + CREATED)
watched = []
for directory in ('mhdlab', 'data', 'src/mhd_core', 'tools', 'tests', 'configs'):
    watched.extend(path for path in (ROOT/directory).rglob('*')
                   if path.is_file() and str(path.relative_to(ROOT)) not in owned)
watched.append(ROOT/'CMakeLists.txt')
preserved_before = inventory(watched)
(WORK/'preserved-before.json').write_text(json.dumps(preserved_before,indent=2,sort_keys=True)+'\n')

old_runner = before['mhdlab/runner.py'].decode()
old_tree = ast.parse(old_runner)
old_helpers = {node.name: node for node in old_tree.body
               if isinstance(node, ast.FunctionDef) and node.name in HELPERS}
assert set(old_helpers) == set(HELPERS)
new_runner = '''"""Compatibility helpers for geometry and thermionic boundary calculations."""

from __future__ import annotations

import numpy as np

from .constants import KB, QE
from .geometry import Geometry


'''+ '\n\n\n'.join(ast.get_source_segment(old_runner, old_helpers[name]) for name in HELPERS)+'\n'
(ROOT/'mhdlab/runner.py').write_text(new_runner)
(ROOT/'mhdlab/__init__.py').write_text('''"""Surface-source research interfaces and external FLASH interchange."""

from .config import load_config

__all__ = ["load_config"]
''')

old_config = before['configs/mykonos.yaml'].decode()
match = re.search(r'^cr_model:\n.*?(?=^\S|\Z)', old_config, flags=re.M|re.S)
assert match is not None
cr_subsection = match.group()
preview_match = re.search(r'^cr_model:\n.*?(?=^\S|\Z)', before['configs/mykonos_preview.yaml'].decode(), flags=re.M|re.S)
assert preview_match is not None and preview_match.group() == cr_subsection
(ROOT/'configs/cr_reference_library.yaml').write_text(cr_subsection)
cr_test = before['tests/test_cr_data.py'].decode()
assert cr_test.count('configs/mykonos.yaml') == 1
(ROOT/'tests/test_cr_data.py').write_text(cr_test.replace('configs/mykonos.yaml', 'configs/cr_reference_library.yaml'))

for path in REMOVED:
    (ROOT/path).unlink()

(ROOT/'tests/test_legacy_mhd_retirement.py').write_text('''"""Import-boundary checks for retirement of the legacy production solver."""

from pathlib import Path
import subprocess
import sys
import textwrap


ROOT = Path(__file__).resolve().parents[1]


def _fresh_python(source):
    completed = subprocess.run(
        [sys.executable, "-I", "-c", textwrap.dedent(source), str(ROOT)],
        cwd=ROOT, capture_output=True, text=True, check=False, timeout=60,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr


def test_package_import_exposes_config_without_eager_scientific_imports():
    _fresh_python("""
        import importlib.abc
        import pathlib
        import sys
        sys.path.insert(0, sys.argv[1])
        blocked = ("mhdlab.runner", "mhdlab.mhd", "mhdlab.materials",
                   "mhdlab.conductivity_tables", "mhdlab.neutrals", "mhdlab.bolsig",
                   "mhdlab.chemistry", "mhdlab.cr_data", "numpy", "scipy")
        class RejectEagerImports(importlib.abc.MetaPathFinder):
            def find_spec(self, fullname, path=None, target=None):
                if any(fullname == name or fullname.startswith(name + ".") for name in blocked):
                    raise AssertionError("unexpected eager import: " + fullname)
        sys.meta_path.insert(0, RejectEagerImports())
        import mhdlab
        assert pathlib.Path(mhdlab.__file__).resolve() == pathlib.Path(sys.argv[1]) / "mhdlab/__init__.py"
        assert mhdlab.__all__ == ["load_config"]
        assert callable(mhdlab.load_config)
        assert not hasattr(mhdlab, "run_from_config")
        assert not any(name in sys.modules for name in blocked)
    """)


def test_legacy_mhd_modules_entrypoints_and_configs_are_absent():
    _fresh_python("""
        import importlib.util
        import pathlib
        import sys
        sys.path.insert(0, sys.argv[1])
        for name in ("mhdlab.mhd", "mhdlab.materials", "mhdlab.conductivity_tables", "MHD_heating"):
            assert importlib.util.find_spec(name) is None, name
        root = pathlib.Path(sys.argv[1])
        removed = (
            "mhdlab/mhd.py", "mhdlab/materials.py", "mhdlab/conductivity_tables.py",
            "MHD_heating.py", "tools/build_conductivity_table.py",
            "configs/mykonos.yaml", "configs/mykonos_preview.yaml",
            "tests/test_mhd.py", "tests/test_smoke_run.py", "tests/test_conductivity_tables.py",
            "data/conductivity/README.md", "data/conductivity/ss304_knoepfel_seed.csv",
        )
        assert all(not (root / name).exists() for name in removed)
        # The native extension contains arithmetic helpers and remains available to build.
        assert (root / "src/mhd_core/reduced.cpp").is_file()
        assert (root / "src/mhd_core/reduced.hpp").is_file()
        assert (root / "src/mhd_core/bindings.cpp").is_file()
        assert (root / "CMakeLists.txt").is_file()
    """)


def test_runner_preserves_only_the_two_compatibility_helpers():
    _fresh_python("""
        import importlib.abc
        import sys
        sys.path.insert(0, sys.argv[1])
        class RejectLegacySolverImports(importlib.abc.MetaPathFinder):
            def find_spec(self, fullname, path=None, target=None):
                if fullname in ("mhdlab.mhd", "mhdlab.materials", "mhdlab.conductivity_tables",
                                "mhdlab.bolsig", "mhdlab.neutrals", "mhdlab.chemistry", "mhdlab.cr_data"):
                    raise AssertionError("unexpected solver import: " + fullname)
        sys.meta_path.insert(0, RejectLegacySolverImports())
        from mhdlab import runner
        functions = {name for name, value in vars(runner).items()
                     if callable(value) and getattr(value, "__module__", None) == runner.__name__}
        assert functions == {"_grid_shape_from_config", "_thermionic_electron_energy_ev"}
        assert not hasattr(runner, "run_from_config")
    """)


def test_flash_heating_hall_and_source_reference_interfaces_remain_importable():
    _fresh_python("""
        import importlib
        import importlib.abc
        import sys
        sys.path.insert(0, sys.argv[1])
        class RejectLegacySolverImports(importlib.abc.MetaPathFinder):
            def find_spec(self, fullname, path=None, target=None):
                if fullname in ("mhdlab.mhd", "mhdlab.materials", "mhdlab.conductivity_tables"):
                    raise AssertionError("reference depends on retired solver: " + fullname)
        sys.meta_path.insert(0, RejectLegacySolverImports())
        interfaces = {
            "mhdlab.flash_surface": ("FlashSurfaceHistory", "load_surface_history"),
            "mhdlab.heating_benchmarks": ("HalfSpaceHeating",),
            "mhdlab.hall_induction_benchmark": ("PeriodicHallInduction",),
            "mhdlab.thermal_release": ("BindingReservoir", "release_step"),
        }
        for name, symbols in interfaces.items():
            module = importlib.import_module(name)
            assert all(callable(getattr(module, symbol)) for symbol in symbols)
    """)
''')

new_tree = ast.parse((ROOT/'mhdlab/runner.py').read_text())
new_helpers = {node.name: node for node in new_tree.body if isinstance(node, ast.FunctionDef)}
assert set(new_helpers) == set(HELPERS)
helper_evidence = {name: {
    'before_AST_sha256': digest(ast.dump(old_helpers[name], include_attributes=False).encode()),
    'after_AST_sha256': digest(ast.dump(new_helpers[name], include_attributes=False).encode()),
    'AST_unchanged': ast.dump(old_helpers[name], include_attributes=False) == ast.dump(new_helpers[name], include_attributes=False),
    'source_unchanged': ast.get_source_segment(old_runner, old_helpers[name]) == ast.get_source_segment(new_runner, new_helpers[name]),
} for name in HELPERS}
assert all(row['AST_unchanged'] and row['source_unchanged'] for row in helper_evidence.values())
assert (ROOT/'configs/cr_reference_library.yaml').read_bytes() == cr_subsection.encode()
assert (ROOT/'tests/test_cr_data.py').read_text().replace('configs/cr_reference_library.yaml','configs/mykonos.yaml') == cr_test
for path in MODIFIED + ('tests/test_legacy_mhd_retirement.py',):
    compile((ROOT/path).read_bytes(), str(ROOT/path), 'exec')
preserved_after = inventory(watched)
assert preserved_after == preserved_before
owned_evidence = {path: {'before_sha256': digest(before[path]) if path in before else None,
                         'after_sha256': digest((ROOT/path).read_bytes()) if (ROOT/path).is_file() else None}
                  for path in REMOVED + MODIFIED + CREATED}
record = {
    'schema': 'mhdlab-legacy-mhd-retirement-source-freeze-v1',
    'status': 'source_and_focused_tests_frozen_before_target_execution',
    'removed_paths': list(REMOVED), 'modified_paths': list(MODIFIED), 'created_paths': list(CREATED),
    'owned_file_hashes': owned_evidence, 'preserved_file_count': len(preserved_before),
    'preserved_inventory_sha256': digest((WORK/'preserved-before.json').read_bytes()),
    'all_preserved_files_byte_identical': True,
    'helper_evidence': helper_evidence,
    'cr_model_subsection_sha256': digest(cr_subsection.encode()), 'cr_model_subsection_byte_identical': True,
    'cr_metadata_test_only_filename_changed': True, 'syntax_only_author_checks_passed': True,
    'target_test_runs_before_freeze': 0, 'solver_or_science_suite_runs': 0,
    'recommended_minimal_tests': [
        'tests/test_legacy_mhd_retirement.py', 'tests/test_cr_data.py', 'tests/test_geometry.py',
        'tests/test_neutrals.py::test_thermionic_electron_energy_floor_uses_surface_temperature',
    ],
}
(WORK/'source-freeze.json').write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
print(json.dumps({'source_freeze':str(WORK/'source-freeze.json'), 'source_freeze_sha256':digest((WORK/'source-freeze.json').read_bytes()),
                  'removed_files':len(REMOVED),'preserved_files':len(preserved_before),'helper_evidence':helper_evidence},sort_keys=True))
