import numpy as np

from mhdlab.geometry import Geometry, boundary_mask


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
