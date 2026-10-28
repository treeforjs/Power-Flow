"""Manufactured zero-barrier, finite prescribed-power single-patch reference.

The parameter records imported below are reused as immutable data only. This
module advances its own joint inventory/temperature RK4 states and gross energy
channels; it never calls the zero-power evolution. Each constant-power segment
owns all four stages of its final step. No stage straddles a power jump.
"""
from __future__ import annotations

from collections.abc import Sequence
from dataclasses import asdict, dataclass, fields
from fractions import Fraction
import math
from numbers import Real
from types import MappingProxyType
from typing import Any, Mapping

from .thermal_feedback_reference import FeedbackProblem, FeedbackReservoir

KB = 1.380649e-23
QE = 1.602176634e-19
_KB = Fraction("1.380649e-23")
_QE = Fraction("1.602176634e-19")
MAX_TOTAL_STEPS = 1_000_000
MAX_SEGMENTS = 1024
MODEL_IDENTITY = "manufactured-zero-barrier-piecewise-power-joint-rk4-v1"


class UnsupportedPowerRangeError(ValueError):
    """Unsupported barrier, domain, time-grid or positive-transfer representation."""


def _label(value: Any, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a nonempty string")
    return value


def _real(value: Any, name: str, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real scalar, not a boolean")
    if value < 0:
        raise ValueError(f"{name} must be nonnegative before floating-point conversion")
    try:
        result = float(value)
    except (OverflowError, ValueError) as exc:
        raise ValueError(f"{name} is outside finite floating-point range") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    if value > 0 and result == 0:
        raise UnsupportedPowerRangeError(f"positive {name} underflows its representation")
    if positive and result <= 0:
        raise ValueError(f"{name} must be positive")
    return result


def _time(value: Any, name: str, *, positive: bool = False) -> float:
    result = _real(value, name, positive=positive)
    # The library's coordinates are already represented floats. A caller with
    # exact Real coordinates must not round away an interval or boundary gap.
    # A JSON reader may instead validate raw decimal support first, then supply
    # explicitly projected canonical floats. Decimal power coefficients need
    # not be exactly binary-representable.
    if value != result:
        raise UnsupportedPowerRangeError(f"{name} is not an exactly represented time coordinate")
    return result


def _f(value: float) -> Fraction:
    return Fraction.from_float(value)


def _project(value: Fraction, name: str) -> float:
    if not value:
        return 0.0
    try:
        result = float(value)
    except OverflowError as exc:
        raise UnsupportedPowerRangeError(f"{name} overflows its representation") from exc
    if not math.isfinite(result) or result == 0:
        raise UnsupportedPowerRangeError(f"nonzero {name} is not representable")
    return result


@dataclass(frozen=True, slots=True)
class PowerSegment:
    start_s: float
    stop_s: float
    power_w_m2: float
    steps: int
    power_provenance: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "start_s", _time(self.start_s, "start_s"))
        object.__setattr__(self, "stop_s", _time(self.stop_s, "stop_s", positive=True))
        object.__setattr__(self, "power_w_m2", _real(self.power_w_m2, "power_w_m2"))
        if not self.start_s < self.stop_s:
            raise ValueError("each power segment requires start_s < stop_s")
        if type(self.steps) is not int or not 1 <= self.steps <= MAX_TOTAL_STEPS:
            raise ValueError(f"segment steps must be an integer in [1,{MAX_TOTAL_STEPS}]")
        _label(self.power_provenance, "power_provenance")


@dataclass(frozen=True, slots=True)
class PowerSchedule:
    duration_s: float
    segments: Sequence[PowerSegment]
    parameter_status: str
    parameter_provenance: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "duration_s", _time(self.duration_s, "duration_s", positive=True))
        _label(self.parameter_status, "parameter_status")
        _label(self.parameter_provenance, "parameter_provenance")
        if self.parameter_status != "manufactured":
            raise ValueError("power schedule requires explicit parameter_status='manufactured'")
        if isinstance(self.segments, (str, bytes)) or not isinstance(self.segments, Sequence):
            raise ValueError("segments must be a finite nonempty sequence of PowerSegment records")
        if not 1 <= len(self.segments) <= MAX_SEGMENTS:
            raise ValueError(f"schedule must contain between 1 and {MAX_SEGMENTS} segments")
        segments = tuple(self.segments)
        if any(not isinstance(segment, PowerSegment) for segment in segments):
            raise ValueError("every segment must be a PowerSegment")
        if segments[0].start_s != 0 or segments[-1].stop_s != self.duration_s:
            raise ValueError("power schedule must cover exactly [0,duration_s]")
        if any(left.stop_s != right.start_s for left, right in zip(segments, segments[1:])):
            raise ValueError("power schedule must have exact adjacency without a gap or overlap")
        if sum(segment.steps for segment in segments) > MAX_TOTAL_STEPS:
            raise ValueError(f"total segment steps exceed {MAX_TOTAL_STEPS}")
        object.__setattr__(self, "segments", segments)


@dataclass(frozen=True, slots=True)
class PowerSample:
    segment_index: int
    segment_step: int
    time_s: float
    coverage_m2: float
    temperature_k: float
    released_m2: float
    input_energy_j_m2: float
    prescribed_input_energy_j_m2: float
    input_energy_residual_j_m2: float
    binding_energy_j_m2: float
    internal_energy_j_m2: float
    translational_energy_j_m2: float
    total_wall_energy_debit_j_m2: float
    gas_energy_j_m2: float
    wall_thermal_change_j_m2: float
    adsorbate_reference_energy_j_m2: float
    adsorbate_reference_change_j_m2: float
    adsorbate_reference_inventory_change_j_m2: float
    inventory_residual_m2: float
    adsorbate_reference_residual_j_m2: float
    wall_debit_residual_j_m2: float
    wall_channel_residual_j_m2: float
    binding_allocation_residual_j_m2: float
    internal_allocation_residual_j_m2: float
    total_energy_residual_j_m2: float
    extensive: Mapping[str, float]

    def __post_init__(self) -> None:
        object.__setattr__(self, "extensive", MappingProxyType(dict(self.extensive)))

    def to_dict(self) -> dict[str, Any]:
        return {field.name: (dict(self.extensive) if field.name == "extensive" else getattr(self, field.name))
                for field in fields(self)}


@dataclass(frozen=True, slots=True)
class PowerHistory:
    problem: FeedbackProblem
    schedule: PowerSchedule
    samples: tuple[PowerSample, ...]

    def to_dict(self) -> dict[str, Any]:
        return {"schema": "mhdlab-thermal-power-reference-v1", "model_identity": MODEL_IDENTITY,
                "method": "segment_local_fixed_step_RK4_joint_coverage_temperature_and_gross_channels",
                "problem": asdict(self.problem), "schedule": asdict(self.schedule),
                "samples": [sample.to_dict() for sample in self.samples],
                "time_representation": "Segment coordinates are canonical floats; logical step times are exact rational start+(stop-start)*step/steps. time_s is their nearest float; segment_index/segment_step reconstruct each logical time.",
                "boundary_state_rule": "One sample at each shared boundary belongs to the segment just completed. The initial sample has segment_index=0 and segment_step=0.",
                "stage_power_rule": "All stages of a segment use its power, including its final stage. The next segment uses its new power from its first stage.",
                "zero_activation_barrier_only": True, "physical_parameters_validated": False,
                "host_budget_accepted": False, "owner_implementation_verified": False,
                "limits": ["Manufactured stationary one-patch constant areal heat capacity and finite nonnegative prescribed power only.",
                           "No temperature-dependent activation kinetics, conduction, radiation, phase change, motion, gap/chemistry evolution or host coupling.",
                           "Fixed-step RK4 errors need independent trajectory/channel and refinement checks; no adaptive or rigorous global error certificate.",
                           "Unrepresentable positive transfers, lost nonzero state/cumulative changes and unsupported domain/time-grid stages are rejected.",
                           "Exact signed net thermal balance may be zero while gross input and wall channels remain positive.",
                           "Input power is explicitly manufactured; it is not inferred from a published current trace or physical foil closure."]}


def _problem(problem: Any) -> FeedbackProblem:
    if not isinstance(problem, FeedbackProblem):
        raise ValueError("problem must be an explicit validated FeedbackProblem")
    if problem.parameter_status != "manufactured":
        raise ValueError("this reference requires manufactured problem parameters")
    if problem.reservoir.activation_energy_ev != 0:
        raise UnsupportedPowerRangeError("only zero activation barrier is qualified for prescribed-power evolution")
    return problem


def _state(problem: FeedbackProblem, coverage: float, temperature: float) -> None:
    if not 0 <= coverage <= problem.reservoir.max_coverage_m2:
        raise UnsupportedPowerRangeError("coverage is outside the declared supported domain")
    if not problem.temperature_min_k <= temperature <= problem.temperature_max_k:
        raise UnsupportedPowerRangeError("temperature is outside the declared supported domain")


def _increment(problem: FeedbackProblem, coverage: float, temperature: float, dt: Fraction, power: float) -> tuple[Fraction, ...]:
    _state(problem, coverage, temperature)
    pool = problem.reservoir
    # Zero-barrier kinetics need no exponential/logarithmic rate approximation.
    # Keep this exact stage amount through every SI product and net balance.
    amount = dt*_f(pool.attempt_frequency_s)*_f(coverage)
    input_energy = dt*_f(power)
    binding = amount*_f(pool.desorption_energy_ev)*_QE
    internal = amount*_f(pool.internal_energy_ev)*_QE
    translation = amount*2*_KB*_f(temperature)
    wall = binding+internal+translation
    thermal = (input_energy-wall)/_f(problem.areal_heat_capacity_j_m2_k)
    return amount, thermal, input_energy, binding, internal, translation, wall


def power_rhs(problem: FeedbackProblem, *, coverage_m2: float, temperature_k: float, power_w_m2: float) -> dict[str, float]:
    """Instantaneous signed thermal derivative and separate gross positive rates."""
    problem = _problem(problem)
    n = _real(coverage_m2, "coverage_m2")
    t = _real(temperature_k, "temperature_k", positive=True)
    p = _real(power_w_m2, "power_w_m2")
    dn, thermal, input_energy, binding, internal, translation, wall = _increment(problem, n, t, Fraction(1), p)
    return {"release_flux_m2_s": _project(dn, "release flux"),
            "coverage_rate_m2_s": _project(-dn, "coverage rate"),
            "temperature_rate_k_s": _project(thermal, "signed thermal rate"),
            "input_power_j_m2_s": _project(input_energy, "input power"),
            "binding_power_j_m2_s": _project(binding, "binding power"),
            "internal_power_j_m2_s": _project(internal, "internal power"),
            "translational_power_j_m2_s": _project(translation, "translation power"),
            "total_wall_power_j_m2_s": _project(wall, "wall power")}


def _advance(problem: FeedbackProblem, coverage: float, temperature: float, increment: tuple[Fraction, ...], scale: Fraction) -> tuple[float, float]:
    dn, thermal = increment[0]*scale, increment[1]*scale
    exact_n = _f(coverage)-dn
    if exact_n < 0 or (coverage > 0 and exact_n == 0):
        raise UnsupportedPowerRangeError("RK4 coverage stage exhausts finite inventory; reduce the step or duration")
    n = _project(exact_n, "coverage state")
    if dn > 0 and n >= coverage:
        raise UnsupportedPowerRangeError("positive coverage debit cannot change its floating-point state")
    exact_t = _f(temperature)+thermal
    if not _f(problem.temperature_min_k) <= exact_t <= _f(problem.temperature_max_k):
        raise UnsupportedPowerRangeError("RK4 temperature stage leaves its declared supported domain")
    t = _project(exact_t, "temperature state")
    if (thermal > 0 and t <= temperature) or (thermal < 0 and t >= temperature):
        raise UnsupportedPowerRangeError("nonzero net thermal increment cannot change its floating-point state")
    _state(problem, n, t)
    return n, t


def _weighted(k1: tuple[Fraction, ...], k2: tuple[Fraction, ...], k3: tuple[Fraction, ...], k4: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    return tuple((a+2*b+2*c+d)/6 for a, b, c, d in zip(k1, k2, k3, k4))


def _increase(current: float, amount: Fraction, name: str) -> float:
    result = _project(_f(current)+amount, name)
    if amount > 0 and result <= current:
        raise UnsupportedPowerRangeError(f"positive {name} increment cannot change its cumulative representation")
    return result


def _prescribed_input(schedule: PowerSchedule, time: Fraction) -> Fraction:
    if not 0 <= time <= _f(schedule.duration_s):
        raise UnsupportedPowerRangeError("prescribed power may not be extrapolated beyond its finite support")
    result = Fraction(0)
    for segment in schedule.segments:
        length = min(time, _f(segment.stop_s))-_f(segment.start_s)
        if length > 0:
            result += _f(segment.power_w_m2)*length
    return result


def _sample(problem: FeedbackProblem, schedule: PowerSchedule, segment_index: int, segment_step: int, time: Fraction, coverage: float, temperature: float, cumulative: tuple[float, ...]) -> PowerSample:
    released, input_energy, binding, internal, translation, wall = cumulative
    d = _f(problem.reservoir.desorption_energy_ev)*_QE
    i = _f(problem.reservoir.internal_energy_ev)*_QE
    n, m, ein = _f(coverage), _f(released), _f(input_energy)
    thermal = _f(problem.areal_heat_capacity_j_m2_k)*(_f(temperature)-_f(problem.initial_temperature_k))
    ads_change = d*m
    ads_inventory_change = d*(_f(problem.initial_coverage_m2)-n)
    gas = _f(internal)+_f(translation)
    prescribed = _prescribed_input(schedule, time)
    quantities = {
        "coverage_m2": coverage, "temperature_k": temperature, "released_m2": released,
        "input_energy_j_m2": input_energy,
        "prescribed_input_energy_j_m2": _project(prescribed, "prescribed input energy"),
        "input_energy_residual_j_m2": _project(ein-prescribed, "input energy residual"),
        "binding_energy_j_m2": binding, "internal_energy_j_m2": internal,
        "translational_energy_j_m2": translation, "total_wall_energy_debit_j_m2": wall,
        "gas_energy_j_m2": _project(gas, "gas energy"),
        "wall_thermal_change_j_m2": _project(thermal, "wall thermal change"),
        "adsorbate_reference_energy_j_m2": _project(-d*n, "adsorbate reference energy"),
        "adsorbate_reference_change_j_m2": _project(ads_change, "adsorbate reference change"),
        "adsorbate_reference_inventory_change_j_m2": _project(ads_inventory_change, "inventory adsorbate reference change"),
        "inventory_residual_m2": _project(n+m-_f(problem.initial_coverage_m2), "inventory residual"),
        "adsorbate_reference_residual_j_m2": _project(ads_inventory_change-ads_change, "adsorbate reference residual"),
        "wall_debit_residual_j_m2": _project(thermal-ein+_f(wall), "wall debit residual"),
        "wall_channel_residual_j_m2": _project(_f(wall)-_f(binding)-gas, "wall channel residual"),
        "binding_allocation_residual_j_m2": _project(_f(binding)-d*m, "binding allocation residual"),
        "internal_allocation_residual_j_m2": _project(_f(internal)-i*m, "internal allocation residual"),
        "total_energy_residual_j_m2": _project(thermal+ads_inventory_change+gas-ein, "total energy residual"),
    }
    extensive = {}
    for key, value in quantities.items():
        if key.endswith("_j_m2"):
            name = key[:-3]
        elif key in ("coverage_m2", "released_m2", "inventory_residual_m2"):
            name = key[:-3]+"_number"
        else:
            continue
        extensive[name] = _project(_f(value)*_f(problem.patch_area_m2), f"extensive {name}")
    return PowerSample(segment_index=segment_index, segment_step=segment_step,
                       time_s=_project(time, "sample time"), **quantities, extensive=extensive)


def evolve_power(problem: FeedbackProblem, schedule: PowerSchedule) -> PowerHistory:
    """Joint N/T RK4 with segment-local power and six gross cumulative channels.

    Exact rational stage amounts and energy products are retained until state
    and cumulative projections. Nonzero net temperature increments must change
    the represented state. Exactly zero thermal balance is allowed and never
    removes the separately accumulated positive input or desorption channels.
    """
    problem = _problem(problem)
    if not isinstance(schedule, PowerSchedule):
        raise ValueError("schedule must be an explicit validated PowerSchedule")
    coverage, temperature = problem.initial_coverage_m2, problem.initial_temperature_k
    cumulative = (0.0,)*6
    samples = [_sample(problem, schedule, 0, 0, Fraction(0), coverage, temperature, cumulative)]
    for segment_index, segment in enumerate(schedule.segments):
        start = _f(segment.start_s)
        dt = (_f(segment.stop_s)-start)/segment.steps
        _project(dt, "positive segment time step")
        for step in range(1, segment.steps+1):
            time = start+step*dt
            if _project(time, "sample time") <= samples[-1].time_s:
                raise UnsupportedPowerRangeError("sample time cannot advance its floating-point representation")
            # Segment power is deliberately a constant argument to all stages.
            k1 = _increment(problem, coverage, temperature, dt, segment.power_w_m2)
            n2, t2 = _advance(problem, coverage, temperature, k1, Fraction(1, 2))
            k2 = _increment(problem, n2, t2, dt, segment.power_w_m2)
            n3, t3 = _advance(problem, coverage, temperature, k2, Fraction(1, 2))
            k3 = _increment(problem, n3, t3, dt, segment.power_w_m2)
            n4, t4 = _advance(problem, coverage, temperature, k3, Fraction(1))
            k4 = _increment(problem, n4, t4, dt, segment.power_w_m2)
            increment = _weighted(k1, k2, k3, k4)
            coverage, temperature = _advance(problem, coverage, temperature, increment, Fraction(1))
            additions = (increment[0], increment[2], increment[3], increment[4], increment[5], increment[6])
            cumulative = tuple(_increase(old, amount, "gross cumulative channel") for old, amount in zip(cumulative, additions))
            samples.append(_sample(problem, schedule, segment_index, step, time, coverage, temperature, cumulative))
    return PowerHistory(problem, schedule, tuple(samples))
