"""Import-boundary checks for retirement of the legacy production solver."""

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
