"""Condition original planar source components on disjoint emission windows.

Selections use offsets inside their original frozen-temperature interval.
Canonical source validation replays its declared law/energy checks but never
evolves a new history or mutates an inventory. Positive original inventories
remain components even when a partial represented quantity is unqualified.
"""
from __future__ import annotations

from collections.abc import Sequence
from dataclasses import InitVar, dataclass, field, fields
import math
from numbers import Integral, Real
from types import MappingProxyType

import numpy as np

from .planar_release import _edges, _plain
from .temkin_planar import TemkinPlanarReleaseHistory, _planar_attributes
from .temkin_release import TemkinReservoir
from .temkin_source_window import FrozenTemkinSourceWindow, STATUSES


SCHEMA = "mhdlab-temkin-planar-source-windows-v1"
SKIP_REASON = "exact_zero_canonical_initial_inventory"
MAX_SELECTED_COMPONENTS = 200_000
MAX_SELECTIONS = 10_000
_JUST_VALIDATED = object()


def _index(value, name):
    if np.ma.isMaskedArray(value) or isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer, not boolean or masked")
    return int(value)


def _real(value, name):
    if np.ma.isMaskedArray(value) or isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real scalar, not boolean or masked")
    try:
        result = float(value)
    except (ValueError, TypeError, OverflowError) as error:
        raise ValueError(f"{name} must be finite and representable") from error
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _text(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a nonempty explicit string")


def _bounds(value, name):
    if np.ma.isMaskedArray(value):
        raise ValueError(f"{name} cannot contain masked bounds")
    try:
        values = tuple(value)
    except TypeError as error:
        raise ValueError(f"{name} requires two finite increasing coordinates") from error
    if len(values) != 2:
        raise ValueError(f"{name} requires two finite increasing coordinates")
    result = tuple(_real(item, name) for item in values)
    if result[0] >= result[1]:
        raise ValueError(f"{name} requires increasing coordinates")
    return result


@dataclass(frozen=True, slots=True)
class SourceWindowSelection:
    interval_index: int
    start_s: float
    stop_s: float

    def __post_init__(self):
        object.__setattr__(self, "interval_index", _index(self.interval_index, "interval_index"))
        object.__setattr__(self, "start_s", _real(self.start_s, "start_s"))
        object.__setattr__(self, "stop_s", _real(self.stop_s, "stop_s"))
        if not 0 <= self.start_s < self.stop_s:
            raise ValueError("source selection requires 0<=start_s<stop_s in its original interval")

    def manifest(self):
        return {"interval_index": self.interval_index, "start_s": self.start_s, "stop_s": self.stop_s}


def _selections(value, time_edges):
    if np.ma.isMaskedArray(value) or isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ValueError("selections must be a finite Sequence of SourceWindowSelection records")
    count = len(value)
    if not 0 < count <= MAX_SELECTIONS:
        raise ValueError("selections must be nonempty and within the selection-count work limit")
    selections = tuple(value[index] for index in range(count))
    if not selections or any(not isinstance(item, SourceWindowSelection) for item in selections):
        raise ValueError("selections must be a nonempty sequence of SourceWindowSelection records")
    # Reconstruct the compact records to reject forged critical selection fields.
    selections = tuple(SourceWindowSelection(item.interval_index, item.start_s, item.stop_s) for item in selections)
    grouped = {}
    for item in selections:
        if item.interval_index >= len(time_edges)-1:
            raise ValueError("selection interval_index is outside the original source history")
        duration = float(time_edges[item.interval_index+1]-time_edges[item.interval_index])
        if item.stop_s > duration:
            raise ValueError("selection offsets must remain inside their original source interval")
        grouped.setdefault(item.interval_index, []).append(item)
    for windows in grouped.values():
        ordered = sorted(windows, key=lambda item: (item.start_s, item.stop_s))
        if any(right.start_s < left.stop_s for left, right in zip(ordered[:-1], ordered[1:])):
            raise ValueError("selections overlap or duplicate within an original source interval")
    return selections


@dataclass(frozen=True, slots=True, eq=False)
class TemkinPlanarWindowComponent:
    selection_index: int
    interval_index: int
    reservoir_index: int
    current_index: int
    width_index: int
    interval_origin_s: float
    x_bounds_m: tuple[float, float]
    z_bounds_m: tuple[float, float]
    source_window: FrozenTemkinSourceWindow

    def __post_init__(self):
        for name in ("selection_index", "interval_index", "reservoir_index", "current_index", "width_index"):
            object.__setattr__(self, name, _index(getattr(self, name), name))
        origin = _real(self.interval_origin_s, "interval_origin_s")
        if origin < 0:
            raise ValueError("interval_origin_s must be nonnegative")
        object.__setattr__(self, "interval_origin_s", origin)
        object.__setattr__(self, "x_bounds_m", _bounds(self.x_bounds_m, "x_bounds_m"))
        object.__setattr__(self, "z_bounds_m", _bounds(self.z_bounds_m, "z_bounds_m"))
        if not isinstance(self.source_window, FrozenTemkinSourceWindow):
            raise ValueError("source_window must be a qualified FrozenTemkinSourceWindow association")
        area = (self.x_bounds_m[1]-self.x_bounds_m[0])*(self.z_bounds_m[1]-self.z_bounds_m[0])
        if not math.isfinite(area) or area <= 0 or area != self.source_window.patch_area_m2:
            raise ValueError("source window area must equal its original x/z patch area once")

    def manifest(self):
        return {**{name: getattr(self, name) for name in ("selection_index", "interval_index", "reservoir_index", "current_index", "width_index", "interval_origin_s")},
                "x_bounds_m": list(self.x_bounds_m), "z_bounds_m": list(self.z_bounds_m),
                "source_window": self.source_window.manifest()}


@dataclass(frozen=True, slots=True, eq=False)
class TemkinPlanarSourceWindows:
    components: tuple[TemkinPlanarWindowComponent, ...]
    selections: tuple[SourceWindowSelection, ...]
    x_edges_m: np.ndarray
    z_edges_m: np.ndarray
    time_edges_s: np.ndarray
    surface_attributes: object
    parameter_provenance: str
    status: str
    provenance: str
    log_mass_tolerance: float
    source_reservoirs: tuple[TemkinReservoir, ...]
    skipped_zero_initials: tuple[object, ...]
    _source_history: TemkinPlanarReleaseHistory = field(repr=False)
    _factory_validation_token: InitVar[object] = None

    def __post_init__(self, _factory_validation_token):
        if not isinstance(self._source_history, TemkinPlanarReleaseHistory):
            raise ValueError("conditioned sources require their authoritative canonical source history")
        original = self._source_history
        if _factory_validation_token is not _JUST_VALIDATED:
            original = TemkinPlanarReleaseHistory(**{item.name: getattr(original, item.name) for item in fields(TemkinPlanarReleaseHistory)})
        object.__setattr__(self, "_source_history", original)
        for name in ("x_edges_m", "z_edges_m", "time_edges_s"):
            object.__setattr__(self, name, _edges(getattr(self, name), name, time=name == "time_edges_s"))
            if not np.array_equal(getattr(self, name), getattr(original, name)):
                raise ValueError("conditioned source geometry or time edges change the canonical source")
        object.__setattr__(self, "surface_attributes", _planar_attributes(self.surface_attributes, self.x_edges_m, self.z_edges_m))
        if _plain(self.surface_attributes) != _plain(original.surface_attributes):
            raise ValueError("conditioned source metadata changes the canonical source association")
        object.__setattr__(self, "selections", _selections(self.selections, self.time_edges_s))
        _text(self.parameter_provenance, "original parameter_provenance")
        if self.parameter_provenance != original.parameter_provenance:
            raise ValueError("parameter_provenance must retain the canonical source declaration")
        _text(self.provenance, "window provenance")
        if not isinstance(self.status, str) or self.status not in STATUSES:
            raise ValueError("status must be an explicit supported source-window label")
        tolerance = _real(self.log_mass_tolerance, "log_mass_tolerance")
        if tolerance <= 0:
            raise ValueError("log_mass_tolerance must be positive")
        object.__setattr__(self, "log_mass_tolerance", tolerance)
        if not isinstance(self.source_reservoirs, Sequence) or isinstance(self.source_reservoirs, (str, bytes)):
            raise ValueError("source_reservoirs must be a finite sequence of original pools")
        if len(self.source_reservoirs) != len(original.reservoirs):
            raise ValueError("source_reservoirs must retain the original pool count")
        pools = tuple(self.source_reservoirs[index] for index in range(len(self.source_reservoirs)))
        if not pools or any(not isinstance(pool, TemkinReservoir) for pool in pools) or len({pool.name for pool in pools}) != len(pools):
            raise ValueError("source_reservoirs must retain the unique original Temkin pools")
        if pools != original.reservoirs:
            raise ValueError("conditioned source pools change the canonical source declaration")
        object.__setattr__(self, "source_reservoirs", pools)
        if len(self.selections)*len(pools)*(len(self.z_edges_m)-1)*(len(self.x_edges_m)-1) > MAX_SELECTED_COMPONENTS:
            raise ValueError("conditioned source exceeds the selected-component memory/work limit")
        if not isinstance(self.components, Sequence) or isinstance(self.components, (str, bytes)) or len(self.components) > MAX_SELECTED_COMPONENTS:
            raise ValueError("components must be a bounded finite sequence")
        supplied_components = tuple(self.components[index] for index in range(len(self.components)))
        components = []
        order = []
        for component in supplied_components:
            if not isinstance(component, TemkinPlanarWindowComponent):
                raise ValueError("components must contain TemkinPlanarWindowComponent records")
            component = TemkinPlanarWindowComponent(**{item.name: getattr(component, item.name) for item in fields(TemkinPlanarWindowComponent)})
            key = (component.selection_index, component.reservoir_index, component.current_index, component.width_index)
            if (component.selection_index >= len(self.selections) or component.reservoir_index >= len(pools)
                    or component.current_index >= len(self.z_edges_m)-1 or component.width_index >= len(self.x_edges_m)-1):
                raise ValueError("conditioned component indices disagree with original selections or geometry")
            selection = self.selections[component.selection_index]
            window = component.source_window
            t, z, x = selection.interval_index, component.current_index, component.width_index
            if (component.interval_index != t or component.interval_origin_s != self.time_edges_s[t]
                    or component.x_bounds_m != tuple(self.x_edges_m[x:x+2])
                    or component.z_bounds_m != tuple(self.z_edges_m[z:z+2])
                    or window.parent_duration_s != float(self.time_edges_s[t+1]-self.time_edges_s[t])
                    or window.start_s != selection.start_s or window.stop_s != selection.stop_s
                    or window.reservoir != pools[component.reservoir_index]
                    or window.initial_coverage_m2 <= 0
                    or window.initial_coverage_m2 != original.coverage_m2[t, component.reservoir_index, z, x]
                    or window.temperature_k != original.temperature_k[t, z, x]
                    or window.status != self.status or window.provenance != self.provenance
                    or window.log_mass_tolerance != tolerance):
                raise ValueError("conditioned component changes original law, time/geometry or window association")
            order.append(key)
            components.append(component)
        if order != sorted(set(order)):
            raise ValueError("components must be unique and ordered by selection/reservoir/current/width")
        object.__setattr__(self, "components", tuple(components))
        skipped, skip_keys = [], []
        if not isinstance(self.skipped_zero_initials, Sequence) or isinstance(self.skipped_zero_initials, (str, bytes)) or len(self.skipped_zero_initials) > MAX_SELECTED_COMPONENTS:
            raise ValueError("zero-inventory skips must be a bounded finite sequence")
        expected = {"selection_index", "interval_index", "reservoir_index", "current_index", "width_index", "reason"}
        for skip_index in range(len(self.skipped_zero_initials)):
            record = self.skipped_zero_initials[skip_index]
            if not hasattr(record, "keys") or set(record) != expected or record["reason"] != SKIP_REASON:
                raise ValueError("zero-inventory skips require explicit original component indices and reason")
            copied = {name: _index(record[name], name) for name in expected-{ "reason" }}
            s, r, z, x = (copied[name] for name in ("selection_index", "reservoir_index", "current_index", "width_index"))
            if (s >= len(self.selections) or r >= len(pools) or z >= len(self.z_edges_m)-1 or x >= len(self.x_edges_m)-1
                    or copied["interval_index"] != self.selections[s].interval_index):
                raise ValueError("zero-inventory skip indices disagree with original source selections")
            if original.coverage_m2[copied["interval_index"], r, z, x] != 0:
                raise ValueError("a positive canonical initial inventory cannot be declared a zero skip")
            skip_keys.append((s, r, z, x))
            skipped.append(MappingProxyType({**copied, "reason": SKIP_REASON}))
        if skip_keys != sorted(set(skip_keys)) or set(skip_keys) & set(order):
            raise ValueError("zero-inventory skips must be unique, ordered and disjoint from positive components")
        expected_keys = {(s, r, z, x) for s in range(len(self.selections)) for r in range(len(pools))
                         for z in range(len(self.z_edges_m)-1) for x in range(len(self.x_edges_m)-1)}
        if set(order) | set(skip_keys) != expected_keys:
            raise ValueError("every selected original source component must be retained or explicitly skipped as zero")
        object.__setattr__(self, "skipped_zero_initials", tuple(skipped))

    def validated_copy(self):
        """Recheck authoritative history and component associations; no evolution."""
        return TemkinPlanarSourceWindows(**{item.name: getattr(self, item.name) for item in fields(type(self)) if item.init})

    def manifest(self):
        from dataclasses import asdict
        return {"schema": SCHEMA, "selections": [item.manifest() for item in self.selections],
                "components": [component.manifest() for component in self.components],
                "skipped_zero_initials": [dict(record) for record in self.skipped_zero_initials],
                "source_reservoirs": [asdict(pool) for pool in self.source_reservoirs],
                "x_edges_m": self.x_edges_m.tolist(), "z_edges_m": self.z_edges_m.tolist(), "time_edges_s": self.time_edges_s.tolist(),
                "surface_attributes": _plain(self.surface_attributes), "parameter_provenance": self.parameter_provenance,
                "status": self.status, "provenance": self.provenance, "log_mass_tolerance": self.log_mass_tolerance,
                "host_budget_accepted": False,
                "time_coordinate": "interval_origin_s remains separate; selection start/stop are original interval offsets; emission quantiles are local delays from selection start",
                "limits": ["Positive original initial inventories are retained regardless of individual partial-quantity representation flags.",
                           "Exact-zero canonical initial densities are skipped explicitly; no cross-interval mathematical survival is reconstructed from a rounded zero.",
                           "No source history is evolved, no inventory is changed and no particle or energy transaction is accepted.",
                           "Status and provenance are caller declarations; original surface attributes, pool records and parameter provenance remain separate.",
                           "Source-emission restrictions are not photon observation gates; no transport, collisions, chemistry or spectra are evaluated."]}


def condition_temkin_planar_source(history, selections, *, log_mass_tolerance, status, provenance):
    """Validate a canonical history and condition its original local components."""
    if not isinstance(history, TemkinPlanarReleaseHistory):
        raise ValueError("history must be a canonical TemkinPlanarReleaseHistory")
    # Replay only the established immutable record validation, not the evolution
    # driver. In particular, forged source moments/diagnostics are not accepted.
    history = TemkinPlanarReleaseHistory(**{item.name: getattr(history, item.name) for item in fields(TemkinPlanarReleaseHistory)})
    selections = _selections(selections, history.time_edges_s)
    if len(selections)*len(history.reservoirs)*(len(history.z_edges_m)-1)*(len(history.x_edges_m)-1) > MAX_SELECTED_COMPONENTS:
        raise ValueError("conditioned source exceeds the selected-component memory/work limit")
    tolerance = _real(log_mass_tolerance, "log_mass_tolerance")
    if tolerance <= 0:
        raise ValueError("log_mass_tolerance must be positive")
    if not isinstance(status, str) or status not in STATUSES:
        raise ValueError("status must be an explicit supported source-window label")
    _text(provenance, "window provenance")
    components, skipped = [], []
    areas = history.patch_areas_m2
    for s, selection in enumerate(selections):
        t = selection.interval_index
        duration = float(history.time_edges_s[t+1]-history.time_edges_s[t])
        for r, pool in enumerate(history.reservoirs):
            for z in range(len(history.z_edges_m)-1):
                for x in range(len(history.x_edges_m)-1):
                    initial = float(history.coverage_m2[t, r, z, x])
                    indices = {"selection_index": s, "interval_index": t, "reservoir_index": r, "current_index": z, "width_index": x}
                    if initial == 0:
                        skipped.append({**indices, "reason": SKIP_REASON})
                        continue
                    window = FrozenTemkinSourceWindow(pool, initial, float(history.temperature_k[t, z, x]), duration,
                                                      selection.start_s, selection.stop_s, float(areas[z, x]),
                                                      log_mass_tolerance=tolerance, status=status, provenance=provenance)
                    components.append(TemkinPlanarWindowComponent(**indices, interval_origin_s=float(history.time_edges_s[t]),
                                      x_bounds_m=tuple(history.x_edges_m[x:x+2]), z_bounds_m=tuple(history.z_edges_m[z:z+2]), source_window=window))
    return TemkinPlanarSourceWindows(tuple(components), selections, history.x_edges_m, history.z_edges_m, history.time_edges_s,
                                    history.surface_attributes, history.parameter_provenance, status, provenance, tolerance,
                                    history.reservoirs, tuple(skipped), history, _factory_validation_token=_JUST_VALIDATED)
