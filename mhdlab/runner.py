"""Compatibility helpers for geometry and thermionic boundary calculations."""

from __future__ import annotations

import numpy as np

from .constants import KB, QE
from .geometry import Geometry


def _grid_shape_from_config(geometry: Geometry, grid_cfg: dict) -> tuple[int, int, str]:
    xmin, xmax, ymin, ymax = geometry.bounds
    width = xmax - xmin
    height = ymax - ymin
    nx = grid_cfg.get("nx")
    ny = grid_cfg.get("ny")
    note = "explicit grid"
    target_cell_m = None

    if "target_cell_size_m" in grid_cfg:
        target_cell_m = float(grid_cfg["target_cell_size_m"])
    elif "target_cell_size_um" in grid_cfg:
        target_cell_m = float(grid_cfg["target_cell_size_um"]) * 1.0e-6

    min_cells = grid_cfg.get("min_cells_per_material_thickness", grid_cfg.get("min_cells_per_foil_thickness"))
    material_scales = geometry.material_length_scales()
    if min_cells is not None and material_scales:
        material_cell = min(material_scales) / max(float(min_cells), 1.0)
        target_cell_m = material_cell if target_cell_m is None else min(target_cell_m, material_cell)

    if target_cell_m is not None:
        if target_cell_m <= 0.0:
            raise ValueError("grid target cell size must be positive")
        nx = int(np.ceil(width / target_cell_m))
        ny = int(np.ceil(height / target_cell_m))
        note = f"target cell <= {target_cell_m * 1.0e6:.2f} um"

    nx = int(nx if nx is not None else 64)
    ny = int(ny if ny is not None else 48)
    max_cells = grid_cfg.get("max_cells")
    if max_cells is not None and nx * ny > int(max_cells):
        scale = np.sqrt((nx * ny) / max(int(max_cells), 1))
        nx = max(4, int(np.floor(nx / scale)))
        ny = max(4, int(np.floor(ny / scale)))
        note += f", capped at {int(max_cells)} cells"
    return nx, ny, note


def _thermionic_electron_energy_ev(
    temperature_k: np.ndarray,
    surface_mask: np.ndarray,
    vacuum_mask: np.ndarray,
    flux_mean_factor: float = 2.0,
    min_temperature_k: float = 0.0,
) -> np.ndarray:
    surface_energy = flux_mean_factor * KB * np.maximum(temperature_k, min_temperature_k) / QE
    total = np.zeros_like(temperature_k, dtype=float)
    count = np.zeros_like(temperature_k, dtype=float)
    shifts = (
        (slice(0, -1), slice(None), slice(1, None), slice(None)),
        (slice(1, None), slice(None), slice(0, -1), slice(None)),
        (slice(None), slice(0, -1), slice(None), slice(1, None)),
        (slice(None), slice(1, None), slice(None), slice(0, -1)),
    )
    for dst_y, dst_x, src_y, src_x in shifts:
        adjacent = surface_mask[src_y, src_x] & vacuum_mask[dst_y, dst_x]
        total[dst_y, dst_x] += np.where(adjacent, surface_energy[src_y, src_x], 0.0)
        count[dst_y, dst_x] += adjacent.astype(float)
    return np.divide(total, count, out=np.zeros_like(total), where=count > 0.0)
