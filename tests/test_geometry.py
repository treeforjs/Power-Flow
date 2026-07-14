import numpy as np

from mhdlab.geometry import Geometry, boundary_mask
from mhdlab.runner import _grid_shape_from_config


def test_rasterizes_parallel_plate_geometry():
    geom = Geometry.from_json("examples/geometry/parallel_plate_gap.json")
    raster = geom.rasterize(nx=20, ny=12)
    assert raster.mask_by_kind("cathode").sum() > 0
    assert raster.mask_by_kind("anode").sum() > 0
    assert raster.mask_by_kind("vacuum").sum() > 0
    surface = boundary_mask(
        raster.mask_by_kind("cathode") | raster.mask_by_kind("anode"),
        raster.mask_by_kind("vacuum"),
    )
    assert surface.sum() > 0


def test_grid_shape_can_refine_by_material_thickness():
    geom = Geometry.from_json("examples/geometry/mykonos_foil_gap_5mm.json")
    nx, ny, note = _grid_shape_from_config(
        geom,
        {
            "min_cells_per_material_thickness": 16,
            "max_cells": 300000,
        },
    )
    assert nx >= 400
    assert ny >= 400
    assert nx * ny <= 300000
    assert "target cell" in note
