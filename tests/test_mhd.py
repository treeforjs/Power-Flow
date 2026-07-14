import numpy as np

from mhdlab.geometry import Geometry
from mhdlab.materials import material_from_config
from mhdlab.mhd import ReducedMHDSolver


def test_out_of_plane_current_density_integrates_to_drive_current():
    geom = Geometry.from_json("examples/geometry/mykonos_foil_gap_5mm.json")
    raster = geom.rasterize(nx=80, ny=96)
    solver = ReducedMHDSolver(raster, material_from_config({"preset": "SS304"}), {"backend": "numpy"})

    current_a = 850.0e3
    state = solver.initial_state()
    jz, _az = solver._solve_inductive_current(state, current_a, dt_s=5.0e-10)
    cell_area = raster.dx * raster.dy

    assert np.isclose(jz[solver.cathode_mask].sum() * cell_area, current_a)
    assert np.isclose(jz[solver.anode_mask].sum() * cell_area, -current_a)


def test_inductive_current_solve_can_develop_nonuniform_jz():
    geom = Geometry.from_json("examples/geometry/mykonos_foil_gap_5mm.json")
    raster = geom.rasterize(nx=120, ny=144)
    solver = ReducedMHDSolver(
        raster,
        material_from_config({"preset": "SS304"}),
        {"backend": "numpy", "poisson_iterations": 50, "induction_iterations": 4},
    )

    current_a = 850.0e3
    state = solver.initial_state()
    jz0, az0 = solver._solve_inductive_current(state, current_a, dt_s=5.0e-10)
    state.az_wb_m = az0
    state.jz_a_m2 = jz0
    jz1, _az1 = solver._solve_inductive_current(state, current_a, dt_s=5.0e-10)
    cathode_j = np.abs(jz1[solver.cathode_mask])

    assert cathode_j.std() / cathode_j.mean() > 1.0e-5


def test_open_boundary_impedance_current_crowds_lateral_edges():
    geom = Geometry.from_json("examples/geometry/mykonos_foil_gap_5mm.json")
    raster = geom.rasterize(nx=160, ny=192)
    solver = ReducedMHDSolver(
        raster,
        material_from_config({"preset": "SS304"}),
        {
            "backend": "numpy",
            "induction_solver": "open_boundary_impedance",
            "induction_max_unknowns": 256,
            "enforce_unidirectional_region_current": True,
        },
    )

    current_a = 850.0e3
    state = solver.initial_state()
    jz, _az = solver._solve_inductive_current(state, current_a, dt_s=5.0e-10)
    active_columns = solver.cathode_mask.any(axis=0)
    profile = np.where(solver.cathode_mask[:, active_columns], np.abs(jz[:, active_columns]), np.nan)
    profile = np.nanmean(profile, axis=0)
    edge_mean = 0.5 * (profile[0] + profile[-1])
    center = profile[profile.size // 2]

    assert edge_mean / center > 2.0
    assert np.isclose(profile[0], profile[-1], rtol=1.0e-2)
