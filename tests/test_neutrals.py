import numpy as np
import pytest

from mhdlab.constants import SPECIES_MASS_KG
from mhdlab.constants import QE
from mhdlab.cross_sections import CrossSectionLibrary, CrossSectionTable
from mhdlab.geometry import Geometry, boundary_mask
from mhdlab.neutrals import VelocityGrid, half_range_flux_moment, shift_fractional_no_wrap, thermal_weights
from mhdlab.neutrals import KineticNeutralSolver
from mhdlab.runner import _thermionic_electron_energy_ev


def test_half_range_flux_moment_positive_and_scales_with_temperature():
    mass = SPECIES_MASS_KG["H2O"]
    cold = half_range_flux_moment(300.0, mass)
    hot = half_range_flux_moment(1200.0, mass)
    assert cold > 0.0
    assert np.isclose(hot / cold, 2.0, rtol=1e-12)


def test_thermal_velocity_weights_are_normalized():
    grid = VelocityGrid.polar(max_speed_m_s=5000.0, n_speed=4, n_angle=8)
    weights = thermal_weights(grid, SPECIES_MASS_KG["H2O"], 600.0)
    assert np.isclose(weights.sum(), 1.0)
    assert np.all(weights >= 0.0)


def test_surface_emission_deposits_neutrals_into_neighboring_vacuum():
    geom = Geometry.from_json("examples/geometry/mykonos_foil_gap_5mm.json")
    raster = geom.rasterize(nx=40, ny=48)
    grid = VelocityGrid.polar(max_speed_m_s=5000.0, n_speed=3, n_angle=6)
    solver = KineticNeutralSolver(raster, ["H2O"], grid, {"wall_sticking": 0.0})
    state = solver.initial_state()
    material = raster.mask_by_kind("cathode") | raster.mask_by_kind("anode")
    surface = boundary_mask(material, raster.mask_by_kind("vacuum") | ~material)
    source = np.zeros(raster.shape)
    source[surface] = 1.0e20
    temperature = np.full(raster.shape, 600.0)

    next_state = solver.step(
        state,
        dt_s=1.0e-9,
        surface_source_m2_s=source,
        surface_temperature_k=temperature,
        source_species="H2O",
        surface_mask=surface,
    )

    density = next_state.density("H2O")
    assert density.sum() > 0.0
    assert np.all(density[~solver.vacuum_mask] == 0.0)


def test_fractional_shift_moves_density_without_rounding_to_zero():
    arr = np.zeros((5, 5))
    arr[2, 2] = 1.0

    shifted = shift_fractional_no_wrap(arr, sx=0.5, sy=0.0)

    assert np.isclose(shifted.sum(), 1.0)
    assert np.isclose(shifted[2, 2], 0.5)
    assert np.isclose(shifted[2, 3], 0.5)


def test_neutral_state_memory_guard_fails_before_allocation():
    geom = Geometry.from_json("examples/geometry/mykonos_foil_gap_5mm.json")
    raster = geom.rasterize(nx=40, ny=48)
    grid = VelocityGrid.polar(max_speed_m_s=5000.0, n_speed=3, n_angle=6)
    solver = KineticNeutralSolver(
        raster,
        ["H2O", "H2"],
        grid,
        {"state_dtype": "float32", "max_working_set_gb": 1.0e-9},
    )

    try:
        solver.initial_state()
    except MemoryError as exc:
        assert "neutral velocity-space state is too large" in str(exc)
        assert "max_working_set_gb" in str(exc)
    else:
        raise AssertionError("expected MemoryError")


def test_neutral_solver_can_step_on_gpu_when_available():
    cp = pytest.importorskip("cupy")
    try:
        cp.cuda.runtime.getDeviceCount()
    except Exception as exc:
        pytest.skip(f"CUDA device unavailable: {exc}")

    geom = Geometry.from_json("examples/geometry/mykonos_foil_gap_5mm.json")
    raster = geom.rasterize(nx=16, ny=16)
    grid = VelocityGrid.polar(max_speed_m_s=5000.0, n_speed=2, n_angle=4)
    solver = KineticNeutralSolver(
        raster,
        ["H2O"],
        grid,
        {"backend": "cuda", "state_dtype": "float32", "max_working_set_gb": 1.0},
    )
    state = solver.initial_state()
    assert solver.backend.name == "cupy"
    assert isinstance(state.f["H2O"], cp.ndarray)

    material = raster.mask_by_kind("cathode") | raster.mask_by_kind("anode")
    surface = boundary_mask(material, raster.mask_by_kind("vacuum") | ~material)
    source = np.zeros(raster.shape)
    source[surface] = 1.0e19
    temperature = np.full(raster.shape, 600.0)
    next_state = solver.step(
        state,
        dt_s=1.0e-9,
        surface_source_m2_s=source,
        surface_temperature_k=temperature,
        source_species="H2O",
        surface_mask=surface,
    )
    density = solver.asnumpy(next_state.density("H2O"))
    assert np.isfinite(density).all()


def test_projectile_energy_comes_from_velocity_distribution():
    geom = Geometry.from_json("examples/geometry/mykonos_foil_gap_5mm.json")
    raster = geom.rasterize(nx=8, ny=8)
    grid = VelocityGrid.polar(max_speed_m_s=1000.0, n_speed=1, n_angle=1)
    solver = KineticNeutralSolver(raster, ["H"], grid, {"state_dtype": "float32"})
    state = solver.initial_state()
    state.f["H"][:] = 1.0
    energy = solver._mean_projectile_energy_ev(state.f["H"], "H")
    expected = 0.5 * SPECIES_MASS_KG["H"] * grid.speed[0] ** 2 / QE
    assert np.allclose(solver.asnumpy(energy), expected)


def test_electron_cross_sections_use_dynamic_bolsig_energy():
    geom = Geometry.from_json("examples/geometry/mykonos_foil_gap_5mm.json")
    raster = geom.rasterize(nx=8, ny=8)
    grid = VelocityGrid.polar(max_speed_m_s=1000.0, n_speed=1, n_angle=4)
    solver = KineticNeutralSolver(raster, ["H2O", "H2O+"], grid, {"state_dtype": "float32"})
    state = solver.initial_state()
    state.f["H2O"][:, solver.vacuum_mask] = 1.0e12
    library = CrossSectionLibrary(
        [
            CrossSectionTable(
                name="fake_e_H2O",
                incident="e",
                target="H2O",
                products={"e": 2.0, "H2O+": 1.0},
                process="test",
                reference="test",
                energy_ev=np.asarray([1.0, 10.0]),
                cross_section_m2=np.asarray([0.0, 1.0e-14]),
                source_file="test",
            )
        ]
    )
    electron_energy = np.full(raster.shape, 10.0)
    next_state = solver.step(
        state,
        dt_s=1.0e-9,
        surface_source_m2_s=np.zeros(raster.shape),
        surface_temperature_k=np.full(raster.shape, 300.0),
        source_species="H2O",
        surface_mask=np.zeros(raster.shape, dtype=bool),
        cross_sections=library,
        incident_energies_ev={"e": 0.0},
        electron_energy_ev=electron_energy,
    )
    assert next_state.density("H2O+").sum() > 0.0


def test_thermionic_electron_energy_floor_uses_surface_temperature():
    temperature = np.full((3, 3), 300.0)
    temperature[1, 1] = 3000.0
    surface = np.zeros((3, 3), dtype=bool)
    surface[1, 1] = True
    vacuum = ~surface
    energy = _thermionic_electron_energy_ev(
        temperature,
        surface,
        vacuum,
        flux_mean_factor=2.0,
        min_temperature_k=300.0,
    )
    expected = 2.0 * 3000.0 * 8.617333262e-5
    assert np.isclose(energy[1, 0], expected, rtol=1.0e-6)
    assert energy[0, 0] == 0.0
