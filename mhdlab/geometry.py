"""2D polygon geometry and rasterization."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


@dataclass(frozen=True)
class Region:
    name: str
    kind: str
    points: np.ndarray


@dataclass
class Raster:
    x: np.ndarray
    y: np.ndarray
    dx: float
    dy: float
    masks: dict[str, np.ndarray]
    kinds: dict[str, str]

    @property
    def nx(self) -> int:
        return int(self.x.size)

    @property
    def ny(self) -> int:
        return int(self.y.size)

    @property
    def shape(self) -> tuple[int, int]:
        return self.y.size, self.x.size

    @property
    def x_min(self) -> float:
        return float(self.x[0] - 0.5 * self.dx)

    @property
    def x_max(self) -> float:
        return float(self.x[-1] + 0.5 * self.dx)

    @property
    def y_min(self) -> float:
        return float(self.y[0] - 0.5 * self.dy)

    @property
    def y_max(self) -> float:
        return float(self.y[-1] + 0.5 * self.dy)

    @property
    def xx(self) -> np.ndarray:
        return np.broadcast_to(self.x[None, :], self.shape)

    @property
    def yy(self) -> np.ndarray:
        return np.broadcast_to(self.y[:, None], self.shape)

    def mask_by_kind(self, kind: str) -> np.ndarray:
        out = np.zeros(self.shape, dtype=bool)
        for name, mask in self.masks.items():
            if self.kinds.get(name) == kind:
                out |= mask
        return out


@dataclass
class Geometry:
    bounds: tuple[float, float, float, float]
    regions: list[Region]
    los: list[dict]

    @classmethod
    def from_json(cls, path: str | Path) -> "Geometry":
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        regions = [
            Region(
                name=item["name"],
                kind=item["kind"],
                points=np.asarray(item["points"], dtype=float),
            )
            for item in data.get("regions", [])
        ]
        return cls(bounds=tuple(data["bounds"]), regions=regions, los=data.get("los", []))

    def rasterize(self, nx: int, ny: int) -> Raster:
        xmin, xmax, ymin, ymax = self.bounds
        dx = (xmax - xmin) / nx
        dy = (ymax - ymin) / ny
        x = xmin + (np.arange(nx) + 0.5) * dx
        y = ymin + (np.arange(ny) + 0.5) * dy
        xx, yy = np.meshgrid(x, y)
        masks = {r.name: points_in_poly(xx, yy, r.points) for r in self.regions}
        kinds = {r.name: r.kind for r in self.regions}
        return Raster(x=x, y=y, dx=dx, dy=dy, masks=masks, kinds=kinds)

    def material_length_scales(self) -> list[float]:
        """Return bounding-box material scales useful for mesh resolution checks."""
        scales: list[float] = []
        for region in self.regions:
            if region.kind not in {"cathode", "anode", "material"}:
                continue
            span = np.ptp(region.points, axis=0)
            positive = span[span > 0.0]
            if positive.size:
                scales.append(float(np.min(positive)))
        return scales


def points_in_poly(xx: np.ndarray, yy: np.ndarray, poly: np.ndarray) -> np.ndarray:
    """Vectorized ray-casting point-in-polygon test."""
    x = xx.ravel()
    y = yy.ravel()
    inside = np.zeros(x.size, dtype=bool)
    px = poly[:, 0]
    py = poly[:, 1]
    j = len(poly) - 1
    for i in range(len(poly)):
        yi = py[i]
        yj = py[j]
        xi = px[i]
        xj = px[j]
        crosses = ((yi > y) != (yj > y)) & (
            x < (xj - xi) * (y - yi) / (yj - yi + 1e-300) + xi
        )
        inside ^= crosses
        j = i
    return inside.reshape(xx.shape)


def boundary_mask(material_mask: np.ndarray, vacuum_mask: np.ndarray) -> np.ndarray:
    """Cells in material adjacent to vacuum."""
    adjacent_vacuum = np.zeros_like(material_mask, dtype=bool)
    for shifted in neighbors4(vacuum_mask):
        adjacent_vacuum |= shifted
    return material_mask & adjacent_vacuum


def vacuum_neighbors(surface_mask: np.ndarray, vacuum_mask: np.ndarray) -> np.ndarray:
    out = np.zeros_like(surface_mask, dtype=bool)
    for shifted in neighbors4(surface_mask):
        out |= shifted
    return out & vacuum_mask


def neighbors4(mask: np.ndarray) -> Iterable[np.ndarray]:
    up = np.zeros_like(mask)
    up[:-1, :] = mask[1:, :]
    down = np.zeros_like(mask)
    down[1:, :] = mask[:-1, :]
    left = np.zeros_like(mask)
    left[:, :-1] = mask[:, 1:]
    right = np.zeros_like(mask)
    right[:, 1:] = mask[:, :-1]
    return up, down, left, right
