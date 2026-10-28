"""Explicitly configured scalar extraction from structured FLASH HDF5 blocks.

Only Cartesian, cell-centered blocks without guard cells are supported. The
caller supplies dataset names, storage axes, active physical axes, node codes,
units, time and a stationary material-side plane. No field name, geometry,
temperature floor, material surface or simulation time is inferred.

The result uses the surface-history interchange as a transport container, but
its manifest identifies a depth-sampled plane diagnostic. It is not a boundary
temperature extrapolation. Two-dimensional data require an explicit extrusion
interval in the current direction; three-dimensional data are area averaged
over that interval. Scalar cells are piecewise constant.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
import hashlib
import math
from numbers import Integral, Real
from pathlib import Path
from types import MappingProxyType
from typing import Any

import h5py
import numpy as np

from .flash_surface import FlashSurfaceHistory


EXTRACTION_SCHEMA = "mhdlab-flash-plane-extraction-v1"
_RAW_AXES = {"x", "y", "z"}
_PHYSICAL_AXES = {"width", "normal", "current"}


def _text(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be an explicit nonempty string")
    return value


def _number(value, name, *, positive=False):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number")
    result = float(value)
    if not math.isfinite(result) or (positive and result <= 0):
        raise ValueError(f"{name} must be finite" + (" and positive" if positive else ""))
    return result


def _integer(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral):
        raise ValueError(f"{name} must be an integer")
    return int(value)


def _numeric_array(value, name):
    array = np.asarray(value)
    if array.dtype.kind not in "iufb":
        raise ValueError(f"{name} must contain real numeric values")
    return array.astype(float)


@dataclass(frozen=True, kw_only=True)
class FlashLayout:
    bounding_box_dataset: str
    node_type_dataset: str
    temperature_dataset: str
    material_fraction_dataset: str
    bounding_box_axes: tuple[str, ...]
    scalar_axes: tuple[str, ...]
    physical_axes: Mapping[str, str]
    leaf_node_value: int
    nonleaf_node_values: tuple[int, ...]
    coordinate_unit: str
    coordinate_scale_to_m: float
    temperature_unit: str
    temperature_scale_to_k: float
    temperature_offset_k: float
    min_material_fraction: float
    geometry: str = "cartesian"
    axis_signs: Mapping[str, int] | None = None
    axis_offsets_m: Mapping[str, float] | None = None
    data_location: str = "cell_centered_no_guard_cells"

    def __post_init__(self):
        names = [self.bounding_box_dataset, self.node_type_dataset,
                 self.temperature_dataset, self.material_fraction_dataset]
        for name in names:
            _text(name, "dataset name")
        if len(set(names)) != 4:
            raise ValueError("the four configured datasets must be distinct")
        bbox_axes, scalar_axes = tuple(self.bounding_box_axes), tuple(self.scalar_axes)
        if not 2 <= len(bbox_axes) <= 3 or len(set(bbox_axes)) != len(bbox_axes) or not set(bbox_axes) <= _RAW_AXES:
            raise ValueError("bounding_box_axes must explicitly list two or three distinct raw x/y/z axes")
        if (not 3 <= len(scalar_axes) <= 4 or len(set(scalar_axes)) != len(scalar_axes)
                or scalar_axes.count("block") != 1 or not set(scalar_axes) <= _RAW_AXES | {"block"}):
            raise ValueError("scalar_axes must explicitly identify block and two or three raw storage axes")
        if not isinstance(self.physical_axes, Mapping):
            raise ValueError("physical_axes must map active raw axes to physical axes")
        axes = dict(self.physical_axes)
        if any(not isinstance(raw, str) or not isinstance(physical, str) for raw, physical in axes.items()):
            raise ValueError("physical_axes keys and values must be explicit axis names")
        if (len(axes) not in (2, 3) or not set(axes) <= set(bbox_axes)
                or not set(axes) <= set(scalar_axes) or len(set(axes.values())) != len(axes)
                or set(axes.values()) not in ({"width", "normal"}, _PHYSICAL_AXES)):
            raise ValueError("active physical_axes must bijectively identify width/normal, and current only in 3D")
        if not set(scalar_axes) - {"block"} <= set(bbox_axes):
            raise ValueError("storage axes must have explicitly named bounding-box columns")
        signs = {axis: 1 for axis in axes} if self.axis_signs is None else dict(self.axis_signs)
        offsets = {axis: 0.0 for axis in axes} if self.axis_offsets_m is None else dict(self.axis_offsets_m)
        if set(signs) != set(axes) or set(offsets) != set(axes):
            raise ValueError("axis_signs and axis_offsets_m must cover exactly the active raw axes")
        signs = {axis: _integer(value, "axis sign") for axis, value in signs.items()}
        if any(value not in (-1, 1) for value in signs.values()):
            raise ValueError("axis signs must be -1 or +1")
        offsets = {axis: _number(value, "axis offset") for axis, value in offsets.items()}
        leaf = _integer(self.leaf_node_value, "leaf_node_value")
        nonleaf = tuple(_integer(value, "nonleaf node code") for value in self.nonleaf_node_values)
        if not nonleaf or len(set(nonleaf)) != len(nonleaf) or leaf in nonleaf:
            raise ValueError("explicit distinct leaf/nonleaf node codes are required")
        for name in ("coordinate_unit", "temperature_unit"):
            _text(getattr(self, name), name)
        for name in ("coordinate_scale_to_m", "temperature_scale_to_k"):
            object.__setattr__(self, name, _number(getattr(self, name), name, positive=True))
        offset = _number(self.temperature_offset_k, "temperature_offset_k")
        purity = _number(self.min_material_fraction, "min_material_fraction")
        if not 0.99 <= purity <= 1.0:
            raise ValueError("min_material_fraction must be an explicit near-pure threshold in [0.99,1]")
        if self.geometry != "cartesian" or self.data_location != "cell_centered_no_guard_cells":
            raise ValueError("only declared Cartesian cell-centered scalars without guard cells are supported")
        object.__setattr__(self, "bounding_box_axes", bbox_axes)
        object.__setattr__(self, "scalar_axes", scalar_axes)
        object.__setattr__(self, "physical_axes", MappingProxyType(axes))
        object.__setattr__(self, "axis_signs", MappingProxyType(signs))
        object.__setattr__(self, "axis_offsets_m", MappingProxyType(offsets))
        object.__setattr__(self, "leaf_node_value", leaf)
        object.__setattr__(self, "nonleaf_node_values", nonleaf)
        object.__setattr__(self, "temperature_offset_k", offset)
        object.__setattr__(self, "min_material_fraction", purity)


@dataclass(frozen=True)
class PlanarFace:
    x_edges_m: tuple[float, ...]
    face_normal_coordinate_m: float
    outward_sign: int
    material_sample_depth_m: float
    current_interval_m: tuple[float, float]
    geometry_atol_m: float = 0.0

    def __post_init__(self):
        edges = tuple(_number(item, "x edge") for item in self.x_edges_m)
        current = tuple(_number(item, "current interval") for item in self.current_interval_m)
        if len(edges) < 2 or any(right <= left for left, right in zip(edges, edges[1:])):
            raise ValueError("x_edges_m must contain strictly increasing full-width edges")
        if len(current) != 2 or current[1] <= current[0]:
            raise ValueError("current_interval_m must explicitly specify a positive extent")
        face = _number(self.face_normal_coordinate_m, "face_normal_coordinate_m")
        depth = _number(self.material_sample_depth_m, "material_sample_depth_m", positive=True)
        sign = _integer(self.outward_sign, "outward_sign")
        tolerance = _number(self.geometry_atol_m, "geometry_atol_m")
        if sign not in (-1, 1) or tolerance < 0:
            raise ValueError("outward_sign must be +/-1 and geometry tolerance nonnegative")
        plane = face - sign * depth
        if not math.isfinite(plane) or plane == face:
            raise ValueError("sample depth must resolve a finite distinct material-side plane")
        lengths = [right - left for left, right in zip(edges, edges[1:])] + [current[1] - current[0]]
        if any(not math.isfinite(length) or length <= 2 * tolerance for length in lengths):
            raise ValueError("projection intervals must be finite and resolvable at the geometry tolerance")
        object.__setattr__(self, "x_edges_m", edges)
        object.__setattr__(self, "current_interval_m", current)
        object.__setattr__(self, "face_normal_coordinate_m", face)
        object.__setattr__(self, "material_sample_depth_m", depth)
        object.__setattr__(self, "outward_sign", sign)
        object.__setattr__(self, "geometry_atol_m", tolerance)

    @property
    def sample_normal_coordinate_m(self):
        return self.face_normal_coordinate_m - self.outward_sign * self.material_sample_depth_m


@dataclass(frozen=True)
class FlashSnapshot:
    path: str | Path
    time_value: float
    time_unit: str
    time_scale_to_s: float
    time_offset_s: float = 0.0

    def __post_init__(self):
        object.__setattr__(self, "path", Path(self.path))
        object.__setattr__(self, "time_value", _number(self.time_value, "time_value"))
        object.__setattr__(self, "time_scale_to_s", _number(self.time_scale_to_s, "time_scale_to_s", positive=True))
        object.__setattr__(self, "time_offset_s", _number(self.time_offset_s, "time_offset_s"))
        _text(self.time_unit, "time_unit")
        if not math.isfinite(self.time_s) or self.time_s < 0:
            raise ValueError("explicitly converted snapshot time must be finite and nonnegative")

    @property
    def time_s(self):
        return self.time_value * self.time_scale_to_s + self.time_offset_s


def _dataset(handle, path, expected_unit=None):
    # External/virtual storage would need additional source hashes and is not
    # covered by this single-file snapshot contract.
    parent = handle
    parts = path.strip("/").split("/")
    for index, part in enumerate(parts):
        if not isinstance(parent.get(part, getlink=True), h5py.HardLink):
            raise ValueError(f"required dataset path cannot use soft/external links: {path}")
        item = parent[part]
        if index < len(parts) - 1:
            if not isinstance(item, h5py.Group):
                raise ValueError(f"invalid dataset path: {path}")
            parent = item
    if not isinstance(item, h5py.Dataset) or item.is_virtual or item.external:
        raise ValueError(f"required field must be a nonvirtual internal HDF5 dataset: {path}")
    if expected_unit is not None and "units" in item.attrs:
        unit = item.attrs["units"]
        if isinstance(unit, bytes):
            unit = unit.decode("utf-8")
        if unit != expected_unit:
            raise ValueError(f"{path} unit attribute conflicts with the explicit conversion")
    if item.dtype.kind not in "iufb":
        raise ValueError(f"{path} must contain real numeric scalar data")
    return item


def _sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _signature(path):
    stat = path.stat()
    return (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)


def _coverage(rectangles, face):
    """Require exact pointwise topology, even below the area-sum tolerance.

    An integral-area tolerance cannot license unresolved gaps or overlaps:
    alternating missing and duplicate slabs can have the correct total area.
    Adjacent raw block edges must therefore agree after the declared conversion.
    """
    events = defaultdict(list)
    for index, (x0, x1, _z0, _z1, _temperature) in enumerate(rectangles):
        events[x0].append((1, index))
        events[x1].append((-1, index))
    events[face.x_edges_m[0]]
    events[face.x_edges_m[-1]]
    positions = sorted(events)
    active = set()
    for index, left in enumerate(positions[:-1]):
        for change, rectangle in events[left]:
            if change == 1:
                active.add(rectangle)
            else:
                active.remove(rectangle)
        right = positions[index + 1]
        if right <= left:
            raise ValueError("projected width sweep is numerically unresolved")
        cursor = face.current_interval_m[0]
        for lower, upper in sorted((rectangles[item][2], rectangles[item][3]) for item in active):
            if lower > cursor:
                raise ValueError("leaf-cell projection has a hole in width/current coverage")
            if lower < cursor:
                raise ValueError("leaf-cell projection overlaps; ambiguous leaf AMR coverage")
            cursor = upper
        if cursor != face.current_interval_m[1]:
            raise ValueError("leaf-cell projection has incomplete current-direction coverage")


def _extract_snapshot(snapshot, layout, face):
    path = snapshot.path.resolve()
    before = _signature(path)
    rectangles = []
    materials = []
    selected_blocks = set()
    with h5py.File(path, "r") as handle:
        bbox_ds = _dataset(handle, layout.bounding_box_dataset, layout.coordinate_unit)
        node_ds = _dataset(handle, layout.node_type_dataset)
        temperature_ds = _dataset(handle, layout.temperature_dataset, layout.temperature_unit)
        material_ds = _dataset(handle, layout.material_fraction_dataset, "1")
        if bbox_ds.dtype.kind == "b" or temperature_ds.dtype.kind == "b":
            raise ValueError("bounding boxes and temperature must be numeric scalars, not Boolean masks")
        boxes = _numeric_array(bbox_ds[...], "bounding boxes")
        nodes = node_ds[...]
        if nodes.dtype.kind not in "iu" or nodes.ndim != 1 or not nodes.size:
            raise ValueError("node types must be a nonempty one-dimensional integer array")
        block_count = nodes.size
        if boxes.shape != (block_count, len(layout.bounding_box_axes), 2) or not np.all(np.isfinite(boxes)):
            raise ValueError("bounding-box shape must be (blocks, configured raw axes, 2) with finite bounds")
        if not np.all(np.isin(nodes, (layout.leaf_node_value, *layout.nonleaf_node_values))):
            raise ValueError("node type contains an unconfigured leaf/nonleaf code")
        if len(temperature_ds.shape) != len(layout.scalar_axes) or temperature_ds.shape != material_ds.shape:
            raise ValueError("temperature/material scalar shapes must match the configured storage axes")
        sizes = dict(zip(layout.scalar_axes, temperature_ds.shape))
        if sizes["block"] != block_count or any(size <= 0 for size in sizes.values()):
            raise ValueError("scalar dimensions do not match the block count")
        for raw in set(layout.scalar_axes) - {"block"} - set(layout.physical_axes):
            if sizes[raw] != 1:
                raise ValueError("inactive storage axes must be explicitly singleton, not a hidden third dimension")
        raw_for = {physical: raw for raw, physical in layout.physical_axes.items()}
        for block in range(block_count):
            coordinates = {}
            for raw in layout.physical_axes:
                lo, hi = boxes[block, layout.bounding_box_axes.index(raw)]
                if hi <= lo:
                    raise ValueError("active bounding-box extents must be positive")
                values = (np.linspace(lo, hi, sizes[raw] + 1) * layout.coordinate_scale_to_m
                          * layout.axis_signs[raw] + layout.axis_offsets_m[raw])
                if layout.axis_signs[raw] == -1:
                    values = values[::-1]
                if not np.all(np.isfinite(values)) or np.any(np.diff(values) <= 0):
                    raise ValueError("converted cell edges are nonfinite or numerically unresolved")
                coordinates[layout.physical_axes[raw]] = values
            if nodes[block] != layout.leaf_node_value:
                continue
            normal = coordinates["normal"]
            plane = face.sample_normal_coordinate_m
            # Physical normal cells are [lower, upper). Exact internal edges
            # select the higher-normal cell, including across block boundaries.
            if not normal[0] <= plane < normal[-1]:
                continue
            normal_index = int(np.searchsorted(normal, plane, side="right") - 1)
            width = coordinates["width"]
            current = coordinates.get("current", np.asarray(face.current_interval_m))
            for i in range(len(width) - 1):
                x0, x1 = max(width[i], face.x_edges_m[0]), min(width[i + 1], face.x_edges_m[-1])
                if x1 <= x0:
                    continue
                for k in range(len(current) - 1):
                    z0, z1 = max(current[k], face.current_interval_m[0]), min(current[k + 1], face.current_interval_m[1])
                    if z1 <= z0:
                        continue
                    physical_indices = {"width": i, "normal": normal_index, "current": k}
                    selection = []
                    for axis in layout.scalar_axes:
                        if axis == "block":
                            selection.append(block)
                        elif axis not in layout.physical_axes:
                            selection.append(0)
                        else:
                            at = physical_indices[layout.physical_axes[axis]]
                            selection.append(at if layout.axis_signs[axis] == 1 else sizes[axis] - 1 - at)
                    sample = tuple(selection)
                    material = float(material_ds[sample])
                    if not math.isfinite(material) or not 0 <= material <= 1 or material < layout.min_material_fraction:
                        raise ValueError("sampled cell is vacuum/mixed/invalid for the configured material-fraction threshold")
                    temperature = float(temperature_ds[sample]) * layout.temperature_scale_to_k + layout.temperature_offset_k
                    if not math.isfinite(temperature) or temperature <= 0:
                        raise ValueError("sampled material temperature must convert to positive finite kelvin")
                    rectangles.append((float(x0), float(x1), float(z0), float(z1), temperature))
                    materials.append(material)
                    selected_blocks.add(block)
        shapes = {"bounding_box": list(bbox_ds.shape), "node_type": list(node_ds.shape),
                  "temperature": list(temperature_ds.shape), "material_fraction": list(material_ds.shape)}
    if not rectangles:
        raise ValueError("no leaf material cells cover the explicitly sampled plane")
    _coverage(rectangles, face)
    values, covered_areas = [], []
    current_length = face.current_interval_m[1] - face.current_interval_m[0]
    for lo, hi in zip(face.x_edges_m, face.x_edges_m[1:]):
        expected = (hi - lo) * current_length
        contributions = [(max(0.0, min(hi, r[1]) - max(lo, r[0])) * (r[3] - r[2]), r[4])
                         for r in rectangles]
        covered = math.fsum(area for area, _ in contributions)
        allowed = 1e-12 * expected + face.geometry_atol_m * ((hi - lo) + current_length)
        if not math.isfinite(expected) or expected <= 0 or abs(covered - expected) > allowed:
            raise ValueError("projected bin area is nonfinite, unresolved or incompletely covered")
        mean = math.fsum((area / covered) * temperature for area, temperature in contributions)
        if not math.isfinite(mean) or mean <= 0:
            raise ValueError("projected material temperature is invalid")
        values.append(mean)
        covered_areas.append(covered)
    checksum = _sha256(path)
    if _signature(path) != before:
        raise ValueError("source snapshot changed during extraction/checksum; use completed immutable output")
    return values, {
        "path": str(path), "sha256": checksum, "file_size_bytes": before[2],
        "time_value": snapshot.time_value, "time_unit": snapshot.time_unit,
        "time_scale_to_s": snapshot.time_scale_to_s, "time_offset_s": snapshot.time_offset_s,
        "time_s": snapshot.time_s, "dataset_shapes": shapes,
        "selected_leaf_blocks": sorted(selected_blocks), "selected_cell_count": len(rectangles),
        "sampled_material_fraction_min": min(materials),
        "covered_bin_area_m2": covered_areas,
    }


def extract_surface_history(
    snapshots: Iterable[FlashSnapshot], layout: FlashLayout, face: PlanarFace, *,
    source_kind: str, run_id: str, source_reference: str, surface_region: str,
    flash_version: str = "",
) -> tuple[FlashSurfaceHistory, dict[str, Any]]:
    """Extract one depth-sampled stationary face and a recoverable manifest.

    Times and source classification are caller assertions, not authenticated
    run metadata. Raw fields and scalar storage are never auto-detected. Supply
    source_kind='synthetic' for manufactured files. Parent AMR blocks are never
    a fallback for incomplete leaf coverage. No exporter writes or solver runs
    occur here.
    """
    if not isinstance(layout, FlashLayout) or not isinstance(face, PlanarFace):
        raise ValueError("layout and face must be validated explicit configuration records")
    inputs = tuple(snapshots)
    if len(inputs) < 2 or any(not isinstance(item, FlashSnapshot) for item in inputs):
        raise ValueError("at least two explicit FlashSnapshot records are required")
    times = [item.time_s for item in inputs]
    # Validate interchange metadata and time ordering before reading raw files.
    FlashSurfaceHistory(times, face.x_edges_m, np.ones((len(times), len(face.x_edges_m) - 1)),
                        source_kind, run_id, source_reference, surface_region, flash_version,
                        temperature_location="material_sample_plane",
                        material_sample_depth_m=face.material_sample_depth_m,
                        boundary_temperature_reference="",
                        current_averaging_interval_m=face.current_interval_m if len(layout.physical_axes) == 3 else None,
                        current_uniformity_reference="")
    paths = [item.path.resolve() for item in inputs]
    if len(set(paths)) != len(paths):
        raise ValueError("each snapshot must refer to a distinct completed source file")
    temperatures, sources = [], []
    try:
        for snapshot in inputs:
            temperature, manifest = _extract_snapshot(snapshot, layout, face)
            temperatures.append(temperature)
            sources.append(manifest)
    except (OSError, KeyError, TypeError, UnicodeDecodeError) as error:
        raise ValueError(f"could not extract the configured FLASH scalar layout: {error}") from error
    history = FlashSurfaceHistory(times, face.x_edges_m, temperatures,
                                  source_kind, run_id, source_reference, surface_region, flash_version,
                                  temperature_location="material_sample_plane",
                                  material_sample_depth_m=face.material_sample_depth_m,
                                  boundary_temperature_reference="",
                                  current_averaging_interval_m=face.current_interval_m if len(layout.physical_axes) == 3 else None,
                                  current_uniformity_reference="")
    manifest = {
        **history.attributes, "schema": EXTRACTION_SCHEMA,
        "surface_history_schema": history.attributes["schema"],
        "diagnostic_kind": "material_side_depth_sampled_plane_not_boundary_extrapolation",
        "sampling": {"face_normal_coordinate_m": face.face_normal_coordinate_m,
                     "outward_sign": face.outward_sign, "material_sample_depth_m": face.material_sample_depth_m,
                     "sample_normal_coordinate_m": face.sample_normal_coordinate_m,
                     "normal_cell_convention": "physical-normal half-open [lower,upper); exact edges select higher-normal cell",
                     "current_interval_m": list(face.current_interval_m), "x_edges_m": list(face.x_edges_m),
                     "geometry_atol_m": face.geometry_atol_m,
                     "coverage_topology": "strict pointwise coverage; geometry_atol_m only bounds integral-area roundoff",
                     "projection": "exact rectangle overlap area mean of piecewise-constant scalar cells",
                     "current_model": "explicit_uniform_extrusion_2d" if len(layout.physical_axes) == 2 else "resolved_current_interval_average_3d"},
        "layout": {"geometry": layout.geometry, "data_location": layout.data_location,
                   "bounding_box_dataset": layout.bounding_box_dataset, "node_type_dataset": layout.node_type_dataset,
                   "temperature_dataset": layout.temperature_dataset, "material_fraction_dataset": layout.material_fraction_dataset,
                   "bounding_box_axes": list(layout.bounding_box_axes), "scalar_axes": list(layout.scalar_axes),
                   "physical_axes": dict(layout.physical_axes), "axis_signs": dict(layout.axis_signs),
                   "axis_offsets_m": dict(layout.axis_offsets_m), "leaf_node_value": layout.leaf_node_value,
                   "nonleaf_node_values": list(layout.nonleaf_node_values),
                   "coordinate_unit": layout.coordinate_unit, "coordinate_scale_to_m": layout.coordinate_scale_to_m,
                   "temperature_unit": layout.temperature_unit, "temperature_scale_to_k": layout.temperature_scale_to_k,
                   "temperature_offset_k": layout.temperature_offset_k, "material_fraction_unit": "1",
                   "min_material_fraction": layout.min_material_fraction},
        "snapshots": sources,
        "limitations": ["stationary declared planar face; no surface reconstruction or expansion tracking",
                        "depth-sampled material scalar; no extrapolation to physical boundary temperature",
                        "scalar location/units/geometry and snapshot times require independent run verification",
                        "averaging must be converged before nonlinear surface-release models are applied"],
    }
    return history, manifest
