import numpy as np

from mhdlab.conductivity_tables import ConductivityTable
from mhdlab.geometry import Geometry
from mhdlab.materials import material_from_config
from mhdlab.mhd import ReducedMHDSolver


def test_conductivity_table_interpolates_log_space(tmp_path):
    path = tmp_path / "sigma.csv"
    path.write_text(
        "\n".join(
            [
                "ensemble_index,density_kg_m3,temperature_k,conductivity_s_m",
                "0,1000,100,100",
                "0,1000,1000,10",
                "0,10000,100,1000",
                "0,10000,1000,100",
            ]
        ),
        encoding="utf-8",
    )
    table = ConductivityTable.from_file(path)
    value = table.interpolate(np.asarray([[3162.27766017]]), np.asarray([[316.227766017]]))
    assert np.isclose(value[0, 0], 100.0, rtol=1.0e-6)


def test_conductivity_table_ensemble_selection(tmp_path):
    path = tmp_path / "sigma.csv"
    path.write_text(
        "\n".join(
            [
                "ensemble_index,density_kg_m3,temperature_k,conductivity_s_m",
                "0,1000,100,100",
                "0,1000,1000,100",
                "0,10000,100,100",
                "0,10000,1000,100",
                "1,1000,100,200",
                "1,1000,1000,200",
                "1,10000,100,200",
                "1,10000,1000,200",
            ]
        ),
        encoding="utf-8",
    )
    table = ConductivityTable.from_file(path)
    rho = np.asarray([[5000.0]])
    temp = np.asarray([[500.0]])
    assert np.isclose(table.interpolate(rho, temp, ensemble_index=1)[0, 0], 200.0)
    assert np.isclose(table.interpolate(rho, temp, statistic="mean")[0, 0], 150.0)


def test_mhd_solver_can_use_conductivity_table(tmp_path):
    path = tmp_path / "sigma.csv"
    path.write_text(
        "\n".join(
            [
                "density_kg_m3,temperature_k,conductivity_s_m",
                "1000,100,500000",
                "1000,1000,500000",
                "10000,100,500000",
                "10000,1000,500000",
            ]
        ),
        encoding="utf-8",
    )
    geom = Geometry.from_json("examples/geometry/mykonos_foil_gap_5mm.json")
    raster = geom.rasterize(nx=24, ny=24)
    solver = ReducedMHDSolver(
        raster,
        material_from_config({"preset": "SS304"}),
        {"backend": "numpy", "conductivity": {"model": "table", "file": str(path), "minimum_s_m": 1.0}},
    )
    state = solver.initial_state()
    sigma = solver._electrical_conductivity(state.temperature_k, state.density_kg_m3, state.jz_a_m2)
    assert np.isclose(sigma[solver.cathode_mask].mean(), 500000.0)
    assert np.all(sigma[solver.vacuum_mask] == 0.0)
