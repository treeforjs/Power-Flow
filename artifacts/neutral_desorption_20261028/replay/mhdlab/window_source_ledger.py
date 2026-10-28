"""Conservative diagnostics of individually qualified source-window amounts.

Aggregation is exact on represented component values, using Fractions without
cutoffs. Mathematical source logs and independently qualified number/energy
fields remain separate. A missing positive contributor makes that quantity's
complete aggregate unavailable; its qualified subtotal is still reported.

The matched uniform-width control changes width only. It redistributes retained
current/state/time/temperature/cost mixtures by exact geometric width fractions,
without materializing a component-by-target-width expansion. These are source
diagnostics, not host acceptance, wall feedback, transport or photon gates.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict, fields
from fractions import Fraction
import json
import math

import numpy as np

from .constants import KB
from .planar_release import _plain
from .temkin_emission_window import ReducedTemkinEmissionWindow
from .temkin_release import TemkinReservoir
from .temkin_source_window import FrozenTemkinSourceWindow
from .temkin_windowed_source import MAX_SELECTED_COMPONENTS, TemkinPlanarSourceWindows


SCHEMA = "mhdlab-window-source-ledger-v1"
QUANTITIES = ("released_number", "translational_energy_j", "internal_energy_j",
              "binding_energy_j", "total_wall_energy_debit_j")
MAX_LEDGER_GROUPS = 250_000
_ZERO = Fraction(0)
_TINY = float(np.finfo(float).tiny)


class WindowSourceLedgerRangeError(ValueError):
    """The diagnostic exceeds its explicitly bounded work or input domain."""


def _q(value):
    return Fraction.from_float(float(value))


def _rational(value):
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


def _number(value):
    try:
        represented = float(value)
    except OverflowError:
        represented = None
    if represented is not None and not math.isfinite(represented):
        represented = None
    residual = None if represented is None else _q(represented)-value
    return {"exact_sum": _rational(value), "value": represented,
            "projection": {
                "exact_zero": value == 0, "positive": value > 0, "negative": value < 0,
                "finite": represented is not None,
                "underflowed": value != 0 and represented == 0.,
                "overflowed": represented is None,
                "rounded": residual is None or residual != 0,
                "subnormal": represented is not None and 0 < abs(represented) < _TINY,
                "exact_rounding_error": None if residual is None else _rational(residual),
            }}


def _amount(subtotal, missing, area, *, scale=Fraction(1), source_ref=None):
    subtotal *= scale
    complete = not missing
    aggregate = _number(subtotal) if complete else {"exact_sum": None, "value": None, "projection": None}
    density = _number(subtotal/area) if complete else {"exact_sum": None, "value": None, "projection": None}
    evidence = {"count": len(missing)}
    if source_ref is None:
        evidence["catalog_indices"] = list(missing)
    else:
        evidence["source_group_ref"] = source_ref
    return {"available": complete, **aggregate,
            "qualified_subtotal": _number(subtotal),
            "per_area": {"available": complete, **density, "qualified_subtotal": _number(subtotal/area)},
            "unqualified": evidence}


class _Accumulator:
    """Trusted internal input: already replayed per-component representations."""
    def __init__(self):
        self.subtotals = {name: _ZERO for name in QUANTITIES}
        self.missing = {name: [] for name in QUANTITIES}
        self.count_times_2kbt = _ZERO
        self.components = 0

    def add(self, values, missing_ids, temperature_k):
        for name in QUANTITIES:
            value = values[name]
            if value is None:
                self.missing[name].append(missing_ids[name])
            else:
                self.subtotals[name] += value
        if values["released_number"] is not None:
            self.count_times_2kbt += values["released_number"]*2*_q(KB)*_q(temperature_k)
        self.components += 1


def _area(exact, canonical, definition):
    projected = _number(exact)
    return {"exact_m2": projected["exact_sum"], "value_m2": projected["value"],
            "projection": projected["projection"],
            "canonical_patch_area_sum": _number(canonical),
            "canonical_minus_exact_residual_m2": _number(canonical-exact),
            "definition": definition}


def _residual(available, value, unavailable):
    return {"available": available,
            **(_number(value) if available else {"exact_sum": None, "value": None, "projection": None}),
            "unavailable_quantities": unavailable}


def _group(accumulator, key, exact_area, canonical_area, definition, *, scale=Fraction(1), source_ref=None):
    quantities = {name: _amount(accumulator.subtotals[name], accumulator.missing[name], exact_area,
                                scale=scale, source_ref=None if source_ref is None else source_ref+"/quantities/"+name)
                  for name in QUANTITIES}
    energy_names = QUANTITIES[1:]
    absent_energy = [name for name in energy_names if accumulator.missing[name]]
    wall_residual = accumulator.subtotals["total_wall_energy_debit_j"]-sum(
        (accumulator.subtotals[name] for name in energy_names[:-1]), _ZERO)
    thermal_absent = [name for name in ("translational_energy_j", "released_number") if accumulator.missing[name]]
    thermal_residual = accumulator.subtotals["translational_energy_j"]-accumulator.count_times_2kbt
    count_energy = _amount(accumulator.count_times_2kbt, accumulator.missing["released_number"], exact_area,
                           scale=scale, source_ref=None if source_ref is None else source_ref+"/residuals/count_times_2kbt_j")
    result = {"key": dict(key), "area": _area(exact_area, canonical_area, definition),
              "component_count": accumulator.components,
              "quantities": quantities,
              "residuals": {
                  "wall_minus_parts_j": _residual(not absent_energy, scale*wall_residual, absent_energy),
                  "translation_minus_count_2kbt_j": _residual(not thermal_absent, scale*thermal_residual, thermal_absent),
                  "count_times_2kbt_j": count_energy,
              }, "host_budget_accepted": False}
    if source_ref is not None:
        result["redistribution"] = {"source_group_ref": source_ref, "fraction": _rational(scale),
                                    "interpretation": "change width only; preserve original current/state/time/temperature/cost mixture"}
    return result


def _replayed_window(window):
    if not isinstance(window, FrozenTemkinSourceWindow):
        raise ValueError("components require FrozenTemkinSourceWindow records")
    pool = TemkinReservoir(**{item.name: getattr(window.reservoir, item.name) for item in fields(TemkinReservoir)})
    arguments = {item.name: getattr(window, item.name) for item in fields(FrozenTemkinSourceWindow) if item.init}
    arguments["reservoir"] = pool
    replay = FrozenTemkinSourceWindow(**arguments)
    # JSON distinguishes boolean flags from equal-valued integers. Compare the
    # private law caches too; the manifest deliberately does not expose all of
    # the query implementation fields.
    original_json = json.dumps(window.manifest(), sort_keys=True, allow_nan=False)
    replay_manifest = replay.manifest()
    if original_json != json.dumps(replay_manifest, sort_keys=True, allow_nan=False):
        raise ValueError("source-window quantity flags or physical associations disagree with constructor replay")
    for item in fields(ReducedTemkinEmissionWindow):
        original, expected = getattr(window.time_law, item.name), getattr(replay.time_law, item.name)
        if type(original) is not type(expected) or original != expected:
            raise ValueError("source-window timing cache disagrees with constructor replay")
    return replay_manifest


def _geometry(source):
    x, z = tuple(map(float, source.x_edges_m)), tuple(map(float, source.z_edges_m))
    dx = tuple(_q(right)-_q(left) for left, right in zip(x[:-1], x[1:]))
    dz = tuple(_q(right)-_q(left) for left, right in zip(z[:-1], z[1:]))
    width, length = _q(x[-1])-_q(x[0]), _q(z[-1])-_q(z[0])
    canonical = tuple(tuple(_q(float((z[j+1]-z[j])*(x[i+1]-x[i]))) for i in range(len(dx)))
                      for j in range(len(dz)))
    whole = sum((sum(row, _ZERO) for row in canonical), _ZERO)
    widths = tuple(sum((row[i] for row in canonical), _ZERO) for i in range(len(dx)))
    currents = tuple(sum(row, _ZERO) for row in canonical)
    return x, z, dx, dz, width, length, canonical, whole, widths, currents


def window_source_ledger(source: TemkinPlanarSourceWindows, *, include_components=True):
    """Return detached exact-represented ledgers and a matched width control.

    ``available`` means all contributing source fields were individually
    qualified. Aggregate/density float projections have their own range and
    rounding flags; a complete exact sum can overflow its SI display. Qualified
    zeros are distinct from positive unrepresented amounts. No threshold,
    population repair, energy adjustment or host acceptance is performed.
    """
    if type(include_components) is not bool:
        raise ValueError("include_components must be an explicit boolean")
    if not isinstance(source, TemkinPlanarSourceWindows):
        raise ValueError("source must be a TemkinPlanarSourceWindows record")
    try:
        source = source.validated_copy()
        pools = source.source_reservoirs
        states = tuple(dict.fromkeys(pool.state_id for pool in pools))
        nx, nz, nr, ns = len(source.x_edges_m)-1, len(source.z_edges_m)-1, len(pools), len(source.selections)
        rows = 2*(1+nx+nz+len(states)+nr+ns+nz*len(states)*nx)+nz*len(states)
        if rows > MAX_LEDGER_GROUPS or len(source.components) > MAX_SELECTED_COMPONENTS:
            raise WindowSourceLedgerRangeError("source ledger exceeds its group/component work limit")
        x, z, dx, dz, width, length, patch_canonical, whole_canonical, width_canonical, current_canonical = _geometry(source)
        total = _Accumulator()
        by_width = [_Accumulator() for _ in range(nx)]
        by_current = [_Accumulator() for _ in range(nz)]
        by_state = {state: _Accumulator() for state in states}
        by_reservoir = [_Accumulator() for _ in pools]
        by_selection = [_Accumulator() for _ in source.selections]
        patch_state = {(j, state, i): _Accumulator() for j in range(nz) for state in states for i in range(nx)}
        mixture = {(j, state): _Accumulator() for j in range(nz) for state in states}
        components, unavailable = [], []
        for index, component in enumerate(source.components):
            manifest = _replayed_window(component.source_window)
            j, i, r, s = component.current_index, component.width_index, component.reservoir_index, component.selection_index
            state = pools[r].state_id
            association = {"component_index": index, "selection_index": s, "interval_index": component.interval_index,
                           "reservoir_index": r, "reservoir_name": pools[r].name, "state_id": state,
                           "current_index": j, "width_index": i, "interval_origin_s": component.interval_origin_s,
                           "x_bounds_m": list(component.x_bounds_m), "z_bounds_m": list(component.z_bounds_m)}
            values, missing = {}, {}
            for name in QUANTITIES:
                record = manifest["quantities"][name]
                if record["qualified"]:
                    values[name] = _q(record["value"])
                else:
                    values[name] = None
                    missing[name] = len(unavailable)
                    unavailable.append({"catalog_index": len(unavailable), **association, "quantity_name": name,
                                        "positive": not record["exact_zero"], "representation": deepcopy(record)})
            for accumulator in (total, by_width[i], by_current[j], by_state[state], by_reservoir[r],
                                by_selection[s], patch_state[j, state, i], mixture[j, state]):
                accumulator.add(values, missing, manifest["temperature_k"])
            if include_components:
                components.append({**association, "source_window": manifest,
                                   "geometric_patch_area_m2": _number(dx[i]*dz[j]),
                                   "canonical_minus_geometric_patch_area_m2": _number(_q(manifest["patch_area_m2"])-dx[i]*dz[j])})
        full_area = width*length
        resolved = {
            "total": _group(total, {}, full_area, whole_canonical, "full physical face once"),
            "by_width": [_group(value, {"width_index": i}, dx[i]*length, width_canonical[i], "width bin times full current length once")
                         for i, value in enumerate(by_width)],
            "by_current": [_group(value, {"current_index": j}, width*dz[j], current_canonical[j], "full width times current bin once")
                           for j, value in enumerate(by_current)],
            "by_state": [_group(by_state[state], {"state_id": state}, full_area, whole_canonical, "full physical face once per declared exact state ID")
                         for state in states],
            "by_reservoir": [_group(value, {"reservoir_index": r, "reservoir_name": pools[r].name, "state_id": pools[r].state_id},
                                     full_area, whole_canonical, "full physical face once per declared reservoir")
                             for r, value in enumerate(by_reservoir)],
            "by_selection": [_group(value, {"selection_index": s}, full_area, whole_canonical, "full physical face once per original selection")
                             for s, value in enumerate(by_selection)],
            "patch_state": [_group(patch_state[j, state, i], {"current_index": j, "state_id": state, "width_index": i},
                                    dx[i]*dz[j], patch_canonical[j][i], "one physical patch once per declared exact state ID")
                            for j in range(nz) for state in states for i in range(nx)],
        }
        bases = [_group(mixture[j, state], {"current_index": j, "state_id": state}, width*dz[j], current_canonical[j],
                        "full width times one current bin, retaining original component mixture")
                 for j in range(nz) for state in states]
        base_indices = {(j, state): index for index, (j, state) in enumerate((j, state) for j in range(nz) for state in states)}
        uniform = {name: deepcopy(resolved[name]) for name in ("total", "by_current", "by_state", "by_reservoir", "by_selection")}
        uniform["by_width"] = [_group(total, {"width_index": i}, dx[i]*length, width_canonical[i],
                                      "width bin times full current length once", scale=dx[i]/width, source_ref="resolved/total")
                               for i in range(nx)]
        uniform["patch_state"] = [_group(mixture[j, state], {"current_index": j, "state_id": state, "width_index": i},
                                          dx[i]*dz[j], patch_canonical[j][i], "one physical patch once per declared exact state ID",
                                          scale=dx[i]/width, source_ref=f"uniform_width_bases/{base_indices[j, state]}")
                                    for j in range(nz) for state in states for i in range(nx)]
        result = {
            "schema": SCHEMA, "resolved": resolved, "matched_uniform_width": uniform,
            "uniform_width_bases": bases, "unqualified_contributors": unavailable,
            "skipped_zero_initials": [dict(value) for value in source.skipped_zero_initials],
            "source_association": {"status": source.status, "provenance": source.provenance,
                                   "parameter_provenance": source.parameter_provenance,
                                   "log_mass_tolerance": source.log_mass_tolerance,
                                   "surface_attributes": _plain(source.surface_attributes),
                                   "source_reservoirs": [asdict(pool) for pool in pools],
                                   "selections": [selection.manifest() for selection in source.selections]},
            "geometry": {"x_edges_m": list(x), "z_edges_m": list(z), "time_edges_s": source.time_edges_s.tolist(),
                         "width_m": _number(width), "current_length_m": _number(length)},
            "quantity_units": {name: "parent" if name == "released_number" else "J" for name in QUANTITIES},
            "density_units": {name: "parent m^-2" if name == "released_number" else "J m^-2" for name in QUANTITIES},
            "component_count": len(source.components), "components_included": include_components,
            "host_budget_accepted": False,
            "interpretation": "exact aggregation of qualified represented component amounts; SI projection flags are independent of source qualification",
            "limits": [
                "Unqualified positive contributor amounts remain unavailable; qualified subtotals are not complete totals.",
                "Exact-zero source costs and empty declared bins are distinct from unrepresented positive amounts.",
                "SI projections can underflow, overflow or round coarsely; exact rational amounts and signed projection errors remain authoritative diagnostics, not accepted budgets.",
                "Density divides the unique exact geometric area once; canonical float patch areas used by source quantities are retained through signed area residuals.",
                "Matched uniform width changes only width position; current/state/time/temperature and supplied independent energy-cost mixtures remain those of the original components.",
                "Translation is expected outgoing flux energy 2*kB*T, not sampled particle energy; binding/removal cost is independent of the kinetic barrier.",
                "This ledger accepts no host transaction, changes no source history, and evolves no gap, molecular chemistry, radiation, photon gate or thermal feedback.",
            ],
            "work_limits": {"max_ledger_groups": MAX_LEDGER_GROUPS, "max_components": MAX_SELECTED_COMPONENTS},
        }
        if include_components:
            result["components"] = components
        return result
    except (AttributeError, TypeError) as error:
        raise ValueError("source-window ledger associations are incomplete or modified") from error
