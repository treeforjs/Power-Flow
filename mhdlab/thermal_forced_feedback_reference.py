"""Manufactured positive-barrier, prescribed-power joint thermal reference.

Only immutable parameter/schedule records are reused. This module advances its
own inventory, temperature and gross channels. A bounded 50-digit Decimal
exponential is a numerical rate approximation, not an exact exponential law or
certified enclosure. Its tiny positive factor is retained before population and
energy scaling, so it is never silently replaced by a binary64 zero hazard.
"""
from __future__ import annotations

from collections.abc import Sequence
from dataclasses import asdict, dataclass, fields
from decimal import Context, Decimal, DecimalException, DivisionByZero, InvalidOperation, Overflow, ROUND_HALF_EVEN, localcontext
from fractions import Fraction
import math
from numbers import Real
from types import MappingProxyType
from typing import Any, Mapping

from .thermal_feedback_reference import FeedbackProblem, FeedbackReservoir
from .thermal_power_reference import PowerSchedule, PowerSegment

KB = 1.380649e-23
QE = 1.602176634e-19
_KB = Fraction("1.380649e-23")
_QE = Fraction("1.602176634e-19")
ARRHENIUS_DECIMAL_PRECISION = 50
ARRHENIUS_DECIMAL_EMIN = -10_000
ARRHENIUS_DECIMAL_EMAX = 10_000
MAX_DIMENSIONLESS_BARRIER = 4096
MAX_TOTAL_STEPS = 1_000_000
MAX_SEGMENTS = 1024
MODEL_IDENTITY = "manufactured-positive-barrier-piecewise-power-joint-rk4-v1"


class UnsupportedForcedFeedbackRangeError(ValueError):
    """Unsupported kinetic domain, state, time grid or transfer representation."""


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
        raise UnsupportedForcedFeedbackRangeError(f"positive {name} underflows its representation")
    if positive and result <= 0:
        raise ValueError(f"{name} must be positive")
    return result


def _f(value: float) -> Fraction:
    return Fraction.from_float(value)


def _project(value: Fraction, name: str) -> float:
    if not value:
        return 0.0
    try:
        result = float(value)
    except OverflowError as exc:
        raise UnsupportedForcedFeedbackRangeError(f"{name} overflows its representation") from exc
    if not math.isfinite(result) or result == 0:
        raise UnsupportedForcedFeedbackRangeError(f"nonzero {name} is not representable")
    return result


def _problem_snapshot(problem: Any) -> FeedbackProblem:
    if not isinstance(problem, FeedbackProblem):
        raise ValueError("problem must be an explicit FeedbackProblem data record")
    # Revalidate every parameter as data; no prior numerical method is called.
    result = FeedbackProblem(**{field.name: getattr(problem, field.name) for field in fields(FeedbackProblem)})
    if result.reservoir.activation_energy_ev <= 0:
        raise UnsupportedForcedFeedbackRangeError("strictly positive base activation energy is required")
    return result


def _schedule_snapshot(schedule: Any) -> PowerSchedule:
    if not isinstance(schedule, PowerSchedule):
        raise ValueError("schedule must be an explicit PowerSchedule data record")
    raw = schedule.segments
    if isinstance(raw, (str, bytes)) or not isinstance(raw, Sequence) or not 1 <= len(raw) <= MAX_SEGMENTS:
        raise ValueError("schedule requires a finite nonempty bounded segment sequence")
    if any(not isinstance(segment, PowerSegment) for segment in raw):
        raise ValueError("every schedule segment must be a PowerSegment data record")
    segments = tuple(PowerSegment(**{field.name: getattr(segment, field.name) for field in fields(PowerSegment)})
                     for segment in raw)
    return PowerSchedule(duration_s=schedule.duration_s, segments=segments,
                         parameter_status=schedule.parameter_status,
                         parameter_provenance=schedule.parameter_provenance)


@dataclass(frozen=True, slots=True)
class ForcedFeedbackSample:
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
class ForcedFeedbackHistory:
    problem: FeedbackProblem
    schedule: PowerSchedule
    samples: tuple[ForcedFeedbackSample, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": "mhdlab-thermal-forced-feedback-reference-v1", "model_identity": MODEL_IDENTITY,
            "method": "segment_local_fixed_step_RK4_joint_coverage_temperature_and_gross_channels",
            "problem": asdict(self.problem), "schedule": asdict(self.schedule),
            "samples": [sample.to_dict() for sample in self.samples],
            "time_representation": "Canonical float segment endpoints; logical steps use exact rational start+(stop-start)*step/steps. time_s is the nearest float and segment_index/segment_step reconstruct the logical time.",
            "boundary_state_rule": "One shared sample belongs to the segment just completed; the initial sample has segment_index=0 and segment_step=0.",
            "stage_power_rule": "All stages of the old segment, including its last stage, use old power; the next segment uses new power from its first stage.",
            "arrhenius_evaluation": {
                "effective_barrier": "Ea0*(1-alpha*N/Nref), with fixed supplied Nref and exact rational complement",
                "base_activation_energy_strictly_positive": True,
                "supported_dimensionless_barrier": [0, MAX_DIMENSIONLESS_BARRIER],
                "range_applies_to": "evaluated nonempty states; empty inventory has exact zero release",
                "precision_decimal_digits": ARRHENIUS_DECIMAL_PRECISION,
                "rounding": "ROUND_HALF_EVEN", "context_Emin": ARRHENIUS_DECIMAL_EMIN,
                "context_Emax": ARRHENIUS_DECIMAL_EMAX,
                "context_traps": ["InvalidOperation", "DivisionByZero", "Overflow"],
                "effective_zero_barrier_factor": "exact one",
                "positive_suppression_must_resolve": True,
                "factor_and_products": "Decimal exponential approximation retained as Fraction before nu*N*dt, all SI costs, net thermal changes and RK weights",
                "qualification": "A 50-digit numerical exponential approximation, not an exact rational/transcendental law, certified enclosure or error certificate."
            },
            "physical_parameters_validated": False, "host_budget_accepted": False,
            "owner_implementation_verified": False,
            "limits": ["Manufactured one stationary patch, one parent reservoir, fixed positive areal heat capacity and finite prescribed nonnegative power only.",
                       "Strictly positive base activation barrier; no physical activation-energy cap is inferred from the supported numerical beta range.",
                       "N and T advance independently; no frozen-T, maintained-T, zero-barrier or unforced analytic thermal trajectory is substituted.",
                       "RHS derivatives and finite-step evolution have different representation limits: a finite tiny rate need not yield a resolvable N/T step.",
                       "Unsupported nonzero projections, lost state/cumulative changes and domain/time-grid violations reject without clipping, floors, retries or repairs.",
                       "Exact net thermal balance may vanish while positive gross channels remain distinct.",
                       "Fixed-step/refinement and independent reference evidence are numerical checks, not physical validation or rigorous error bounds.",
                       "No actual-water, material, phase, source-history, spectral, host, gap, owner or MHD coupling is supplied."]
        }


def _state(problem: FeedbackProblem, coverage: float, temperature: float) -> None:
    if not 0 <= coverage <= problem.reservoir.max_coverage_m2:
        raise UnsupportedForcedFeedbackRangeError("coverage is outside the declared supported domain")
    if not problem.temperature_min_k <= temperature <= problem.temperature_max_k:
        raise UnsupportedForcedFeedbackRangeError("temperature is outside the declared supported domain")


def _arrhenius_factor(problem: FeedbackProblem, coverage: float, temperature: float) -> Fraction:
    pool = problem.reservoir
    complement = 1-_f(pool.coverage_coefficient)*_f(coverage)/_f(pool.reference_coverage_m2)
    if complement < 0:
        raise UnsupportedForcedFeedbackRangeError("effective kinetic barrier is negative")
    beta = _f(pool.activation_energy_ev)*_QE*complement/(_KB*_f(temperature))
    if beta < 0 or beta > MAX_DIMENSIONLESS_BARRIER:
        raise UnsupportedForcedFeedbackRangeError("dimensionless kinetic barrier is outside [0,4096]")
    if beta == 0:
        return Fraction(1)
    context = Context(prec=ARRHENIUS_DECIMAL_PRECISION, rounding=ROUND_HALF_EVEN,
                      Emin=ARRHENIUS_DECIMAL_EMIN, Emax=ARRHENIUS_DECIMAL_EMAX,
                      clamp=0, traps=[InvalidOperation, DivisionByZero, Overflow])
    try:
        with localcontext(context):
            exponent = Decimal(beta.numerator)/Decimal(beta.denominator)
            factor = (-exponent).exp()
    except (DecimalException, OverflowError, ValueError) as exc:
        raise UnsupportedForcedFeedbackRangeError("Arrhenius approximation could not be represented") from exc
    if not factor.is_finite() or not 0 < factor < 1:
        raise UnsupportedForcedFeedbackRangeError("positive kinetic suppression is unresolvable at the declared Decimal precision")
    return Fraction(factor)


def _increment(problem: FeedbackProblem, coverage: float, temperature: float, dt: Fraction, power: float) -> tuple[Fraction, ...]:
    _state(problem, coverage, temperature)
    pool = problem.reservoir
    factor = _arrhenius_factor(problem, coverage, temperature) if coverage > 0 else Fraction(0)
    amount = dt*_f(pool.attempt_frequency_s)*_f(coverage)*factor
    input_energy = dt*_f(power)
    binding = amount*_f(pool.desorption_energy_ev)*_QE
    internal = amount*_f(pool.internal_energy_ev)*_QE
    translation = amount*2*_KB*_f(temperature)
    wall = binding+internal+translation
    thermal = (input_energy-wall)/_f(problem.areal_heat_capacity_j_m2_k)
    return amount, thermal, input_energy, binding, internal, translation, wall


def forced_feedback_rhs(problem: FeedbackProblem, *, coverage_m2: float, temperature_k: float, power_w_m2: float) -> dict[str, float]:
    """Temperature-dependent release and signed thermal/gross energy rates."""
    problem = _problem_snapshot(problem)
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
        raise UnsupportedForcedFeedbackRangeError("RK4 coverage stage exhausts finite inventory; reduce the step or duration")
    n = _project(exact_n, "coverage state")
    if dn > 0 and n >= coverage:
        raise UnsupportedForcedFeedbackRangeError("positive coverage debit cannot change its floating-point state")
    exact_t = _f(temperature)+thermal
    if not _f(problem.temperature_min_k) <= exact_t <= _f(problem.temperature_max_k):
        raise UnsupportedForcedFeedbackRangeError("RK4 temperature stage leaves its declared supported domain")
    t = _project(exact_t, "temperature state")
    if (thermal > 0 and t <= temperature) or (thermal < 0 and t >= temperature):
        raise UnsupportedForcedFeedbackRangeError("nonzero net thermal increment cannot change its floating-point state")
    _state(problem, n, t)
    return n, t


def _weighted(k1: tuple[Fraction, ...], k2: tuple[Fraction, ...], k3: tuple[Fraction, ...], k4: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    return tuple((a+2*b+2*c+d)/6 for a, b, c, d in zip(k1, k2, k3, k4))


def _increase(current: float, amount: Fraction, name: str) -> float:
    result = _project(_f(current)+amount, name)
    if amount > 0 and result <= current:
        raise UnsupportedForcedFeedbackRangeError(f"positive {name} increment cannot change its cumulative representation")
    return result


def _prescribed_input(schedule: PowerSchedule, time: Fraction) -> Fraction:
    if not 0 <= time <= _f(schedule.duration_s):
        raise UnsupportedForcedFeedbackRangeError("prescribed power may not be extrapolated beyond finite support")
    result = Fraction(0)
    for segment in schedule.segments:
        length = min(time, _f(segment.stop_s))-_f(segment.start_s)
        if length > 0:
            result += _f(segment.power_w_m2)*length
    return result


def _sample(problem: FeedbackProblem, schedule: PowerSchedule, segment_index: int, segment_step: int, time: Fraction, coverage: float, temperature: float, cumulative: tuple[float, ...]) -> ForcedFeedbackSample:
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
    return ForcedFeedbackSample(segment_index=segment_index, segment_step=segment_step,
                                time_s=_project(time, "sample time"), **quantities, extensive=extensive)


def evolve_forced_feedback(problem: FeedbackProblem, schedule: PowerSchedule) -> ForcedFeedbackHistory:
    """Direct coupled N/T RK4 with positive-barrier rates and gross channels."""
    problem = _problem_snapshot(problem)
    schedule = _schedule_snapshot(schedule)
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
                raise UnsupportedForcedFeedbackRangeError("sample time cannot advance its floating-point representation")
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
    return ForcedFeedbackHistory(problem, schedule, tuple(samples))
