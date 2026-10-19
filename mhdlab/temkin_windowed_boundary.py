"""Paired thermal packets from qualified partial planar source components.

The two hypotheses change width position only. Original interval origin,
selection start and conditional local delay remain separate, including when
their sum loses timing precision. No random numbers, host transactions or gap
evolution are generated here.
"""

from __future__ import annotations

from bisect import bisect_right
from collections.abc import Mapping
from dataclasses import dataclass, fields
from fractions import Fraction
import json
import math
from numbers import Integral, Real
from typing import ClassVar

import numpy as np

from .constants import KB, QE
from .planar_boundary import _validate_extracted_face
from .temkin_release import LAW_IDENTITY, TemkinReservoir
from .temkin_source_window import FrozenTemkinSourceWindow
from .temkin_windowed_source import TemkinPlanarSourceWindows, TemkinPlanarWindowComponent


SCHEMA = "mhdlab-temkin-windowed-boundary-packets-v1"
ALLOCATION_MODEL = "one_per_conditioned_positive_component_then_exact_count_cdf_refinement"
VELOCITY_MODEL = "nonrelativistic_outward_flux_weighted_maxwellian"
_EPS = float(np.finfo(float).eps)
_MAX_SPEED = .05 * 299_792_458.
_ZERO = Fraction(0)


class TemkinWindowedBoundaryRangeError(ValueError):
    """A packet amount, coordinate or local time lacks representation support."""


def _real(value, name, *, positive=False):
    if np.ma.isMaskedArray(value) or isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be an explicit finite real scalar")
    try:
        value = float(value)
    except (ValueError, OverflowError) as error:
        raise ValueError(f"{name} must be finite") from error
    if not math.isfinite(value) or (positive and value <= 0):
        raise ValueError(f"{name} must be finite{' and positive' if positive else ''}")
    return value


def _immutable(value):
    array = np.asarray(value)
    return np.frombuffer(array.tobytes(order="C"), dtype=array.dtype).reshape(array.shape)


def _fraction_log(value):
    if value <= 0:
        raise TemkinWindowedBoundaryRangeError("a positive mathematical quantity is required")
    return math.log(value.numerator) - math.log(value.denominator)


def _finite_fraction(value, name):
    try:
        result = float(value)
    except OverflowError as error:
        raise TemkinWindowedBoundaryRangeError(f"{name} exceeds finite float64 range") from error
    if not math.isfinite(result):
        raise TemkinWindowedBoundaryRangeError(f"{name} exceeds finite float64 range")
    if value and result == 0:
        raise TemkinWindowedBoundaryRangeError(f"{name} is nonzero but below representable float64 range")
    return result


def _positive_fraction(value, name, tolerance):
    result = _finite_fraction(value, name)
    if result <= 0:
        raise TemkinWindowedBoundaryRangeError(f"{name} must remain positive")
    ratio = float(Fraction(result) / value)
    rounding = abs(math.log(ratio))
    allowance = 64 * _EPS * (1 + abs(_fraction_log(value)))
    if rounding + allowance > tolerance:
        raise TemkinWindowedBoundaryRangeError(f"{name} float64 quantization exceeds the packet log allowance")
    return result


def _draws(value):
    if np.ma.isMaskedArray(value):
        raise ValueError("uniforms cannot contain masked variates")
    try:
        objects = np.asarray(value, dtype=object)
    except (ValueError, TypeError) as error:
        raise ValueError("uniforms require an (n,7) real array") from error
    if objects.ndim != 2 or objects.shape[1] != 7 or objects.shape[0] == 0:
        raise ValueError("uniforms require a nonempty (n,7) array")
    for value in objects.flat:
        if not 0 < _real(value, "uniform variate") < 1:
            raise ValueError("uniforms must lie strictly inside (0,1)")
    return np.array(objects, dtype=float)


def _revalidate_source(source):
    if not isinstance(source, TemkinPlanarSourceWindows):
        raise ValueError("source must be a TemkinPlanarSourceWindows record")
    try:
        source = source.validated_copy()
        pools = tuple(TemkinReservoir(**{item.name: getattr(pool, item.name) for item in fields(TemkinReservoir)})
                      for pool in source.source_reservoirs)
        components = []
        for component in source.components:
            if not isinstance(component, TemkinPlanarWindowComponent) or not isinstance(component.source_window, FrozenTemkinSourceWindow):
                raise ValueError("conditioned components require original physical source-window records")
            original = component.source_window
            pool = TemkinReservoir(**{item.name: getattr(original.reservoir, item.name) for item in fields(TemkinReservoir)})
            arguments = {item.name: getattr(original, item.name) for item in fields(FrozenTemkinSourceWindow) if item.init}
            arguments["reservoir"] = pool
            replay = FrozenTemkinSourceWindow(**arguments)
            if replay.manifest() != original.manifest():
                raise ValueError("physical source-window quantities or timing association were modified")
            arguments = {item.name: getattr(component, item.name) for item in fields(TemkinPlanarWindowComponent) if item.init}
            arguments["source_window"] = replay
            components.append(TemkinPlanarWindowComponent(**arguments))
        # validated_copy checks the retained canonical baseline, including the
        # source of every zero skip. The scalar replay above additionally guards
        # the cached physical/time records without integrating a new history.
        if pools != source.source_reservoirs:
            raise ValueError("original reservoir records were modified")
        return source
    except (AttributeError, TypeError) as error:
        raise ValueError("source-window record association is incomplete or modified") from error


def _position(left, right, uniform, tolerance, name):
    left, right, uniform = map(Fraction, (left, right, uniform))
    width = right - left
    local_exact = width * uniform
    exact = left + local_exact
    value = float(exact)
    offset = float(local_exact)
    if not math.isfinite(value):
        raise TemkinWindowedBoundaryRangeError(f"{name} exceeds finite float64 range")
    error = max(abs(Fraction(value) - exact), abs(Fraction(offset) - local_exact)) / width
    if error > Fraction(tolerance):
        raise TemkinWindowedBoundaryRangeError(f"{name} rounding is unresolved relative to its interval width")
    if not left <= Fraction(value) <= right:
        raise TemkinWindowedBoundaryRangeError(f"{name} left its original geometry")
    return value, offset


def _global_time(origin, start, delay, tolerance):
    exact = Fraction(origin) + Fraction(start) + Fraction(delay)
    value = _finite_fraction(exact, "best-effort global emission time")
    error = Fraction(value) - exact
    recovered_delay = Fraction(value) - Fraction(origin) - Fraction(start)
    qualified = recovered_delay > 0
    if qualified:
        try:
            relative_error = float(error / Fraction(delay))
            log_error = abs(math.log1p(relative_error)) if relative_error > -1 else math.inf
        except (OverflowError, ValueError):
            log_error = math.inf
        qualified = log_error + 64 * _EPS <= tolerance
    return value, bool(qualified), _finite_fraction(abs(error), "global-time absolute rounding error")


def _positive_exp(logarithm, name, tolerance):
    try:
        value = math.exp(logarithm)
    except OverflowError as error:
        raise TemkinWindowedBoundaryRangeError(f"{name} exceeds finite float64 range") from error
    if not math.isfinite(value) or value <= 0:
        raise TemkinWindowedBoundaryRangeError(f"positive {name} is outside representable float64 range")
    if abs(math.log(value) - logarithm) + 64 * _EPS * (1 + abs(logarithm)) > tolerance:
        raise TemkinWindowedBoundaryRangeError(f"{name} float64 projection exceeds the packet log allowance")
    return value


def _weighted_energy(indices, weights, masses, velocities):
    total = _ZERO
    for packet in indices:
        squares = sum((Fraction(float(value))**2 for value in velocities[packet]), _ZERO)
        total += Fraction(float(weights[packet])) * Fraction(float(masses[packet])) * squares / 2
    return total


_ARRAY_UNITS = {
    **{name: "1" for name in ("component_index", "selection_index", "interval_index", "reservoir_index", "source_current_index",
        "source_width_index", "global_time_precision_qualified", "component_selection_index", "component_interval_index",
        "component_reservoir_index", "component_current_index", "component_width_index", "component_probabilities",
        "component_probability_underflowed", "component_probability_representation_qualified")},
    **{name: "s" for name in ("interval_origin_s", "window_start_offset_s", "window_stop_offset_s", "local_time_delay_s",
                              "best_effort_time_s", "global_time_absolute_rounding_error_s", "time_edges_s")},
    "log_local_time_delay_s": "ln(delay / 1 s)", "component_log_probabilities": "ln(probability)",
    "component_log_released_particles": "ln(number / 1 parent)",
    **{name: "m" for name in ("resolved_x_m", "uniform_x_m", "resolved_position_m", "uniform_width_position_m", "x_edges_m", "z_edges_m",
                              "resolved_width_offset_m", "uniform_width_offset_m", "source_current_offset_m")},
    "temperature_k": "K", "mass_kg": "kg", "velocity_m_s": "m s^-1", "momentum_kg_m_s": "kg m s^-1",
    "particle_weight": "parent per packet", "internal_energy_ev": "eV per parent", "desorption_energy_ev": "eV per parent",
    "packet_counts": "packets", "component_packet_weight": "parent per packet",
    "component_released_particles": "parent", "component_particle_balance_residual": "parent",
    **{name: "J" for name in ("component_expected_translational_energy_j", "component_sampled_translational_energy_j",
        "component_translational_energy_sampling_residual_j", "component_expected_internal_energy_j", "component_weighted_internal_energy_j",
        "component_expected_binding_energy_j", "component_weighted_binding_energy_j", "component_expected_total_wall_energy_debit_j",
        "component_sampled_total_wall_energy_debit_j", "component_wall_energy_sampling_residual_j")},
}


@dataclass(frozen=True, slots=True, eq=False, init=False)
class TemkinWindowedEmissionPackets:
    component_index: np.ndarray
    selection_index: np.ndarray
    interval_index: np.ndarray
    reservoir_index: np.ndarray
    source_current_index: np.ndarray
    source_width_index: np.ndarray
    interval_origin_s: np.ndarray
    window_start_offset_s: np.ndarray
    window_stop_offset_s: np.ndarray
    local_time_delay_s: np.ndarray
    log_local_time_delay_s: np.ndarray
    best_effort_time_s: np.ndarray
    global_time_precision_qualified: np.ndarray
    global_time_absolute_rounding_error_s: np.ndarray
    resolved_x_m: np.ndarray
    uniform_x_m: np.ndarray
    resolved_position_m: np.ndarray
    uniform_width_position_m: np.ndarray
    resolved_width_offset_m: np.ndarray
    uniform_width_offset_m: np.ndarray
    source_current_offset_m: np.ndarray
    temperature_k: np.ndarray
    mass_kg: np.ndarray
    velocity_m_s: np.ndarray
    momentum_kg_m_s: np.ndarray
    particle_weight: np.ndarray
    internal_energy_ev: np.ndarray
    desorption_energy_ev: np.ndarray
    component_selection_index: np.ndarray
    component_interval_index: np.ndarray
    component_reservoir_index: np.ndarray
    component_current_index: np.ndarray
    component_width_index: np.ndarray
    component_released_particles: np.ndarray
    component_log_released_particles: np.ndarray
    component_probabilities: np.ndarray
    component_log_probabilities: np.ndarray
    component_probability_underflowed: np.ndarray
    component_probability_representation_qualified: np.ndarray
    packet_counts: np.ndarray
    component_packet_weight: np.ndarray
    component_particle_balance_residual: np.ndarray
    component_expected_translational_energy_j: np.ndarray
    component_sampled_translational_energy_j: np.ndarray
    component_translational_energy_sampling_residual_j: np.ndarray
    component_expected_internal_energy_j: np.ndarray
    component_weighted_internal_energy_j: np.ndarray
    component_expected_binding_energy_j: np.ndarray
    component_weighted_binding_energy_j: np.ndarray
    component_expected_total_wall_energy_debit_j: np.ndarray
    component_sampled_total_wall_energy_debit_j: np.ndarray
    component_wall_energy_sampling_residual_j: np.ndarray
    x_edges_m: np.ndarray
    z_edges_m: np.ndarray
    time_edges_s: np.ndarray
    state_ids: tuple[str, ...]
    reservoir_names: tuple[str, ...]
    required_minimum_n: int
    total_released_particles: float
    total_packet_particles: float
    total_particle_balance_residual: float
    expected_translational_energy_j: float
    sampled_translational_energy_j: float
    translational_energy_sampling_residual_j: float
    expected_internal_energy_j: float
    weighted_internal_energy_j: float
    expected_binding_energy_j: float
    weighted_binding_energy_j: float
    expected_total_wall_energy_debit_j: float
    sampled_total_wall_energy_debit_j: float
    wall_energy_sampling_residual_j: float
    face_normal_coordinate_m: float
    outward_normal_sign: int
    unequal_macro_weights: bool
    status: str
    provenance: str
    parameter_provenance: str
    _manifest_json: str
    schema: ClassVar[str] = SCHEMA
    host_budget_accepted: ClassVar[bool] = False

    def __init__(self, *args, **kwargs):
        raise TypeError("use sample_temkin_windowed_emission_packets to construct packet records")

    def manifest(self):
        return json.loads(self._manifest_json)


def _packet_record(payload, metadata):
    expected = {item.name for item in fields(TemkinWindowedEmissionPackets)} - {"_manifest_json"}
    if set(payload) != expected:
        raise RuntimeError("internal windowed packet record fields are incomplete")
    record = object.__new__(TemkinWindowedEmissionPackets)
    array_names = {name for name, value in payload.items() if isinstance(value, np.ndarray)}
    if array_names != set(_ARRAY_UNITS):
        raise RuntimeError("windowed packet array units are incomplete")
    for name, value in payload.items():
        if name in array_names:
            if not np.all(np.isfinite(value)):
                raise TemkinWindowedBoundaryRangeError(f"packet array {name} is not finite")
            value = _immutable(value)
        object.__setattr__(record, name, value)
    metadata["array_units"] = dict(_ARRAY_UNITS)
    object.__setattr__(record, "_manifest_json", json.dumps(metadata, allow_nan=False))
    return record


def sample_temkin_windowed_emission_packets(source, state_masses_kg, uniforms, *, face_normal_coordinate_m, outward_normal_sign):
    """Sample every conditioned positive component with exact-count refinement.

    Columns are component refinement, conditional delay, width position, normal
    Rayleigh variate, tangent radius, tangent angle and current position. Each
    component receives one mandatory packet. No internal generator is used.
    Qualified extensive released number and source energies are required;
    underflowed areal quantities do not erase an otherwise qualified component.
    """
    source = _revalidate_source(source)
    draws = _draws(uniforms)
    coordinate = _real(face_normal_coordinate_m, "face_normal_coordinate_m")
    if isinstance(outward_normal_sign, (bool, np.bool_)) or not isinstance(outward_normal_sign, Integral) or outward_normal_sign not in (-1, 1):
        raise ValueError("outward_normal_sign must be explicitly -1 or +1")
    sign = int(outward_normal_sign)
    _validate_extracted_face(source, coordinate, sign)
    if not isinstance(state_masses_kg, Mapping):
        raise ValueError("state_masses_kg must map exact source states to positive masses")
    components = source.components
    required, number = len(components), len(draws)
    if required == 0:
        raise ValueError("conditioned source has no positive initial components to sample")
    if number < required:
        raise ValueError(f"at least {required} packets are required for every positive conditioned component")
    quantities, population, expected = [], [], {name: [] for name in ("translational", "internal", "binding", "total_wall")}
    mass_by_state = {}
    for component in components:
        window = component.source_window
        records = window.manifest()["quantities"]
        names = ("released_number", "translational_energy_j", "internal_energy_j", "binding_energy_j", "total_wall_energy_debit_j")
        if any(not records[name]["qualified"] for name in names):
            raise TemkinWindowedBoundaryRangeError("every positive conditioned component requires qualified extensive release and energy; no component is omitted")
        count = records["released_number"]["value"]
        if count is None or count <= 0:
            raise TemkinWindowedBoundaryRangeError("conditioned emitted number must be positive and represented")
        population.append(float(count))
        quantities.append(records)
        for key, name in zip(expected, names[1:]):
            expected[key].append(float(records[name]["value"]))
        state = window.reservoir.state_id
        if state not in state_masses_kg:
            raise ValueError(f"missing exact-state mass for {state!r}; species fallback is forbidden")
        mass_by_state[state] = _real(state_masses_kg[state], "exact-state mass_kg", positive=True)
    fractions = [Fraction(value) for value in population]
    total_exact = sum(fractions, _ZERO)
    total = _positive_fraction(total_exact, "total released particles", source.log_mass_tolerance)
    cumulative, subtotal = [], _ZERO
    for count in fractions:
        subtotal += count
        cumulative.append(subtotal)
    assigned = np.empty(number, dtype=np.int64)
    assigned[:required] = np.arange(required)
    for packet in range(required, number):
        assigned[packet] = bisect_right(cumulative, Fraction(float(draws[packet, 0])) * total_exact)
    counts = np.bincount(assigned, minlength=required)
    weights = np.array([_positive_fraction(count / int(n), "component packet weight", source.log_mass_tolerance)
                        for count, n in zip(fractions, counts)])
    probabilities, log_probabilities, probability_underflowed, probability_qualified = [], [], [], []
    for count in fractions:
        fraction = count / total_exact
        probability = float(fraction)
        probabilities.append(probability)
        log_probabilities.append(_fraction_log(fraction))
        probability_underflowed.append(probability == 0)
        error = abs(math.log(float(Fraction(probability) / fraction))) if probability > 0 else math.inf
        probability_qualified.append(error + 64 * _EPS * (1 + abs(_fraction_log(fraction))) <= source.log_mass_tolerance)
    names = ("selection_index", "interval_index", "reservoir_index", "current_index", "width_index")
    component_labels = {name: np.array([getattr(component, name) for component in components], dtype=np.int64) for name in names}
    temperatures = np.array([component.source_window.temperature_k for component in components])[assigned]
    state_ids = tuple(components[index].source_window.reservoir.state_id for index in assigned)
    reservoir_names = tuple(components[index].source_window.reservoir.name for index in assigned)
    masses = np.array([mass_by_state[state] for state in state_ids])
    packet_weights = weights[assigned]
    origins, starts, stops, delays, log_delays, times, global_flags, global_errors = (np.empty(number) for _ in range(8))
    velocities, momenta = np.empty((number, 3)), np.empty((number, 3))
    resolved, uniform = np.empty((number, 3)), np.empty((number, 3))
    resolved_width_offsets, uniform_width_offsets, current_offsets = (np.empty(number) for _ in range(3))
    for packet, selected in enumerate(assigned):
        component = components[int(selected)]
        window = component.source_window
        origin, start, stop = component.interval_origin_s, window.start_s, window.stop_s
        log_delay = float(window.time_law.log_quantile(float(draws[packet, 1])))
        # A required local delay must satisfy the timing law's own projection
        # qualification; a finite logarithm alone is insufficient for packets.
        delay = float(window.time_law.quantile(float(draws[packet, 1])))
        time, flag, error = _global_time(origin, start, delay, source.log_mass_tolerance)
        origins[packet], starts[packet], stops[packet] = origin, start, stop
        delays[packet], log_delays[packet], times[packet] = delay, log_delay, time
        global_flags[packet], global_errors[packet] = flag, error
        sigma_log = .5 * math.fsum((math.log(KB), math.log(float(temperatures[packet])), -math.log(float(masses[packet]))))
        normal = _positive_exp(sigma_log + .5 * math.log(-2 * math.log1p(-float(draws[packet, 3]))), "outgoing normal speed", source.log_mass_tolerance)
        radius = _positive_exp(sigma_log + .5 * math.log(-2 * math.log(float(draws[packet, 4]))), "tangential radius", source.log_mass_tolerance)
        angle = 2 * math.pi * float(draws[packet, 5])
        for axis, trig in ((0, math.cos(angle)), (2, math.sin(angle))):
            product = Fraction(radius) * Fraction(trig)
            velocities[packet, axis] = (math.copysign(_positive_fraction(abs(product), "tangential velocity component", source.log_mass_tolerance), product)
                                        if product else 0.)
        velocities[packet, 1] = sign * normal
        if not np.all(np.isfinite(velocities[packet])) or math.hypot(*velocities[packet]) > _MAX_SPEED:
            raise TemkinWindowedBoundaryRangeError("thermal sample exceeds the nonrelativistic speed limit 0.05c")
        for axis in range(3):
            product = Fraction(float(masses[packet])) * Fraction(float(velocities[packet, axis]))
            momenta[packet, axis] = (math.copysign(_positive_fraction(abs(product), "sampled momentum", source.log_mass_tolerance), product)
                                     if product else 0.)
        x, resolved_width_offsets[packet] = _position(*component.x_bounds_m, float(draws[packet, 2]), source.log_mass_tolerance, "resolved width position")
        comparison_x, uniform_width_offsets[packet] = _position(float(source.x_edges_m[0]), float(source.x_edges_m[-1]), float(draws[packet, 2]), source.log_mass_tolerance, "uniform-width position")
        z, current_offsets[packet] = _position(*component.z_bounds_m, float(draws[packet, 6]), source.log_mass_tolerance, "current position")
        resolved[packet], uniform[packet] = (x, coordinate, z), (comparison_x, coordinate, z)
    internal_ev = np.array([components[index].source_window.reservoir.internal_energy_ev for index in assigned])
    binding_ev = np.array([components[index].source_window.reservoir.desorption_energy_ev for index in assigned])
    component_number_exact = [Fraction(float(weight)) * int(count) for weight, count in zip(weights, counts)]
    component_residual = [_finite_fraction(actual - target, "component number closure residual") for actual, target in zip(component_number_exact, fractions)]
    sampled, weighted_internal, weighted_binding, sampled_total = [], [], [], []
    sampled_exact, internal_exact, binding_exact = [], [], []
    for index, component in enumerate(components):
        indices = np.flatnonzero(assigned == index)
        thermal = _weighted_energy(indices, packet_weights, masses, velocities)
        pool = component.source_window.reservoir
        internal = component_number_exact[index] * Fraction(QE) * Fraction(pool.internal_energy_ev)
        binding = component_number_exact[index] * Fraction(QE) * Fraction(pool.desorption_energy_ev)
        combined = thermal + internal + binding
        sampled_exact.append(thermal)
        internal_exact.append(internal)
        binding_exact.append(binding)
        sampled.append(_positive_fraction(thermal, "sampled component translation", source.log_mass_tolerance))
        weighted_internal.append(_positive_fraction(internal, "weighted component internal energy", source.log_mass_tolerance) if internal else 0.)
        weighted_binding.append(_positive_fraction(binding, "weighted component binding energy", source.log_mass_tolerance) if binding else 0.)
        sampled_total.append(_positive_fraction(combined, "sampled component total wall debit", source.log_mass_tolerance))
    expected = {key: np.array(value) for key, value in expected.items()}
    sampled, weighted_internal, weighted_binding, sampled_total = map(np.array, (sampled, weighted_internal, weighted_binding, sampled_total))
    total_packet_exact = sum(component_number_exact, _ZERO)
    totals = {
        "total_released_particles": total, "total_packet_particles": _positive_fraction(total_packet_exact, "total packet particles", source.log_mass_tolerance),
        "total_particle_balance_residual": _finite_fraction(total_packet_exact - total_exact, "total number closure residual"),
    }
    for key, array in (("expected_translational_energy_j", expected["translational"]), ("expected_internal_energy_j", expected["internal"]),
                       ("expected_binding_energy_j", expected["binding"]), ("expected_total_wall_energy_debit_j", expected["total_wall"])):
        exact = sum((Fraction(float(value)) for value in array), _ZERO)
        totals[key] = _positive_fraction(exact, key, source.log_mass_tolerance) if exact else 0.
    for key, exact in (("sampled_translational_energy_j", sum(sampled_exact, _ZERO)), ("weighted_internal_energy_j", sum(internal_exact, _ZERO)),
                       ("weighted_binding_energy_j", sum(binding_exact, _ZERO)),
                       ("sampled_total_wall_energy_debit_j", sum(sampled_exact + internal_exact + binding_exact, _ZERO))):
        totals[key] = _positive_fraction(exact, key, source.log_mass_tolerance) if exact else 0.
    totals["translational_energy_sampling_residual_j"] = math.fsum((totals["sampled_translational_energy_j"], -totals["expected_translational_energy_j"]))
    totals["wall_energy_sampling_residual_j"] = math.fsum((totals["sampled_total_wall_energy_debit_j"], -totals["expected_total_wall_energy_debit_j"]))
    payload = dict(
        component_index=assigned, selection_index=component_labels["selection_index"][assigned], interval_index=component_labels["interval_index"][assigned],
        reservoir_index=component_labels["reservoir_index"][assigned], source_current_index=component_labels["current_index"][assigned], source_width_index=component_labels["width_index"][assigned],
        interval_origin_s=origins, window_start_offset_s=starts, window_stop_offset_s=stops, local_time_delay_s=delays, log_local_time_delay_s=log_delays,
        best_effort_time_s=times, global_time_precision_qualified=global_flags.astype(bool), global_time_absolute_rounding_error_s=global_errors,
        resolved_x_m=resolved[:, 0], uniform_x_m=uniform[:, 0], resolved_position_m=resolved, uniform_width_position_m=uniform,
        resolved_width_offset_m=resolved_width_offsets, uniform_width_offset_m=uniform_width_offsets, source_current_offset_m=current_offsets,
        temperature_k=temperatures, mass_kg=masses, velocity_m_s=velocities, momentum_kg_m_s=momenta, particle_weight=packet_weights,
        internal_energy_ev=internal_ev, desorption_energy_ev=binding_ev,
        component_selection_index=component_labels["selection_index"], component_interval_index=component_labels["interval_index"],
        component_reservoir_index=component_labels["reservoir_index"], component_current_index=component_labels["current_index"], component_width_index=component_labels["width_index"],
        component_released_particles=np.array(population), component_log_released_particles=np.array([records["released_number"]["log_value"] for records in quantities]),
        component_probabilities=np.array(probabilities), component_log_probabilities=np.array(log_probabilities),
        component_probability_underflowed=np.array(probability_underflowed), component_probability_representation_qualified=np.array(probability_qualified),
        packet_counts=counts, component_packet_weight=weights, component_particle_balance_residual=np.array(component_residual),
        component_expected_translational_energy_j=expected["translational"], component_sampled_translational_energy_j=sampled,
        component_translational_energy_sampling_residual_j=np.array([math.fsum((float(a), -float(b))) for a, b in zip(sampled, expected["translational"])]),
        component_expected_internal_energy_j=expected["internal"], component_weighted_internal_energy_j=weighted_internal,
        component_expected_binding_energy_j=expected["binding"], component_weighted_binding_energy_j=weighted_binding,
        component_expected_total_wall_energy_debit_j=expected["total_wall"], component_sampled_total_wall_energy_debit_j=sampled_total,
        component_wall_energy_sampling_residual_j=np.array([math.fsum((float(a), -float(b))) for a, b in zip(sampled_total, expected["total_wall"])]),
        x_edges_m=np.asarray(source.x_edges_m), z_edges_m=np.asarray(source.z_edges_m), time_edges_s=np.asarray(source.time_edges_s),
        state_ids=state_ids, reservoir_names=reservoir_names, required_minimum_n=required, **totals,
        face_normal_coordinate_m=coordinate, outward_normal_sign=sign, unequal_macro_weights=bool(np.any(packet_weights != packet_weights[0])),
        status=source.status, provenance=source.provenance, parameter_provenance=source.parameter_provenance,
    )
    metadata = {
        "schema": SCHEMA, "law_identity": LAW_IDENTITY, "allocation_model": ALLOCATION_MODEL, "velocity_model": VELOCITY_MODEL,
        "source": source.manifest(), "status": source.status, "provenance": source.provenance, "parameter_provenance": source.parameter_provenance,
        "packets": number, "required_minimum_n": required, "host_budget_accepted": False,
        "face_normal_coordinate_m": coordinate, "outward_normal_sign": sign, "state_masses_kg": mass_by_state,
        "totals": totals, "global_time_precision_unqualified_packets": int(np.count_nonzero(~global_flags.astype(bool))),
        "components": [{"component_index": index, **{name: int(component_labels[name][index]) for name in names},
                        "released_particles": population[index], "packet_count": int(counts[index]), "packet_weight": float(weights[index]),
                        "particle_balance_residual": component_residual[index], "probability": probabilities[index],
                        "log_probability": log_probabilities[index], "probability_underflowed": probability_underflowed[index],
                        "probability_representation_qualified": probability_qualified[index],
                        "expected_translational_energy_j": float(expected["translational"][index]), "sampled_translational_energy_j": float(sampled[index]),
                        "expected_internal_energy_j": float(expected["internal"][index]), "weighted_internal_energy_j": float(weighted_internal[index]),
                        "expected_binding_energy_j": float(expected["binding"][index]), "weighted_binding_energy_j": float(weighted_binding[index]),
                        "expected_total_wall_energy_debit_j": float(expected["total_wall"][index]), "sampled_total_wall_energy_debit_j": float(sampled_total[index])}
                       for index in range(required)],
        "time_coordinate": "original interval_origin_s + window_start_offset_s + local_time_delay_s; retain all three and log_local_time_delay_s separately",
        "limits": ["Global times are best-effort rounded sums; flags qualify recovery of the retained local delay, not a shared experimental timing origin.",
                   "All positive conditioned components are represented or the call fails; zero probabilities are explicitly flagged and exact-count refinement retains their mandatory packet.",
                   "Packet rounding allowance is additional to the source-window numerical qualification; neither is a certified bound or host acceptance.",
                   "Position offsets and global coordinates are rounded independently with error qualified relative to each original interval width; unsupported geometry projection is rejected.",
                   "Resolved and uniform hypotheses share the complete time/state/temperature/current/velocity/momentum/weight mixture and change only width position.",
                   "Sampled translation and weighted state/removal costs remain separate from the source expectation; signed residuals are retained without velocity or inventory repair.",
                   "Stationary prescribed-temperature emission only; no field advance, collisions, radiation, gap transport, photon gating, wall feedback or host transaction."],
    }
    return _packet_record(payload, metadata)
