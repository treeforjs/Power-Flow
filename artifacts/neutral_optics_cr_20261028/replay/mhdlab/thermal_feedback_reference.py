"""Manufactured zero-power single-patch desorption/cooling reference.

Coverage and temperature are separate RK4 state variables.  The kinetic barrier
is independent of the removal cost.  No analytic temperature/coverage relation,
existing source solver, material law, host transaction or physical default is
used.  This deliberately bounded floating-point reference rejects unsupported
positive increments and stages rather than clipping or silently dropping them.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, fields
from fractions import Fraction
import math
from numbers import Real
from types import MappingProxyType
from typing import Any, Mapping

KB = 1.380649e-23
QE = 1.602176634e-19
_KB = Fraction("1.380649e-23")
_QE = Fraction("1.602176634e-19")
MAX_STEPS = 1_000_000
MODEL_IDENTITY = "manufactured-zero-power-single-patch-joint-rk4-v1"


class UnsupportedFeedbackRangeError(ValueError):
    """The declared state, time step or representation is unsupported."""


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
    if not math.isfinite(result) or result < 0 or (positive and result <= 0):
        raise ValueError(f"{name} must be finite and {'positive' if positive else 'nonnegative'}")
    if value > 0 and result == 0:
        raise UnsupportedFeedbackRangeError(f"positive {name} underflows its representation")
    return result


def _f(value: float) -> Fraction:
    return Fraction.from_float(value)


def _project(value: Fraction, name: str) -> float:
    if not value:
        return 0.0
    try:
        result = float(value)
    except OverflowError as exc:
        raise UnsupportedFeedbackRangeError(f"{name} overflows its representation") from exc
    if not math.isfinite(result) or result == 0:
        raise UnsupportedFeedbackRangeError(f"nonzero {name} is not representable")
    return result


@dataclass(frozen=True, slots=True)
class FeedbackReservoir:
    """Explicit manufactured Temkin parameters; all fields are required.

    Ea(N)=Ea0*(1-alpha*N/Nref). Nref is fixed.  Ea0 controls kinetics;
    desorption_energy_ev is the separate constant removal cost D.
    """

    name: str
    state_id: str
    attempt_frequency_s: float
    activation_energy_ev: float
    desorption_energy_ev: float
    internal_energy_ev: float
    coverage_coefficient: float
    reference_coverage_m2: float
    max_coverage_m2: float
    kinetic_provenance: str
    energy_provenance: str

    def __post_init__(self) -> None:
        for name in ("name", "state_id", "kinetic_provenance", "energy_provenance"):
            _label(getattr(self, name), name)
        for name in ("attempt_frequency_s", "reference_coverage_m2", "max_coverage_m2"):
            object.__setattr__(self, name, _real(getattr(self, name), name, positive=True))
        for name in ("activation_energy_ev", "desorption_energy_ev", "internal_energy_ev", "coverage_coefficient"):
            object.__setattr__(self, name, _real(getattr(self, name), name))
        if self.activation_energy_ev > 0:
            effect = _f(self.coverage_coefficient)*_f(self.max_coverage_m2)/_f(self.reference_coverage_m2)
            if effect > 1:
                raise ValueError("alpha*max_coverage_m2/reference_coverage_m2 exceeds one")


def _copy_reservoir(value: Any) -> FeedbackReservoir:
    # A structural copy also accepts an existing TemkinReservoir as data only.
    # No type import, source evolution, inversion or method call is performed.
    try:
        data = {field.name: getattr(value, field.name) for field in fields(FeedbackReservoir)}
    except AttributeError as exc:
        raise ValueError("reservoir must supply every explicit FeedbackReservoir parameter") from exc
    return FeedbackReservoir(**data)


@dataclass(frozen=True, slots=True)
class FeedbackProblem:
    reservoir: Any
    initial_coverage_m2: float
    initial_temperature_k: float
    areal_heat_capacity_j_m2_k: float
    patch_area_m2: float
    temperature_min_k: float
    temperature_max_k: float
    parameter_status: str
    parameter_provenance: str
    energy_reference: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "reservoir", _copy_reservoir(self.reservoir))
        object.__setattr__(self, "initial_coverage_m2", _real(self.initial_coverage_m2, "initial_coverage_m2"))
        for name in ("initial_temperature_k", "areal_heat_capacity_j_m2_k", "patch_area_m2", "temperature_min_k", "temperature_max_k"):
            object.__setattr__(self, name, _real(getattr(self, name), name, positive=True))
        for name in ("parameter_status", "parameter_provenance", "energy_reference"):
            _label(getattr(self, name), name)
        if self.parameter_status != "manufactured":
            raise ValueError("this reference requires explicit parameter_status='manufactured'")
        if not self.temperature_min_k < self.temperature_max_k:
            raise ValueError("temperature_min_k must be below temperature_max_k")
        if not self.temperature_min_k <= self.initial_temperature_k <= self.temperature_max_k:
            raise ValueError("initial temperature is outside the declared supported range")
        if self.initial_coverage_m2 > self.reservoir.max_coverage_m2:
            raise ValueError("initial coverage exceeds the declared reservoir domain")


@dataclass(frozen=True, slots=True)
class FeedbackSample:
    time_s: float
    coverage_m2: float
    temperature_k: float
    released_m2: float
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
class FeedbackHistory:
    problem: FeedbackProblem
    duration_s: float
    steps: int
    samples: tuple[FeedbackSample, ...]

    def to_dict(self) -> dict[str, Any]:
        return {"schema": "mhdlab-thermal-feedback-reference-v1", "model_identity": MODEL_IDENTITY,
                "method": "fixed_step_RK4_joint_coverage_temperature_and_cumulative_channels",
                "problem": asdict(self.problem), "duration_s": self.duration_s, "steps": self.steps,
                "samples": [sample.to_dict() for sample in self.samples],
                "zero_input_power": True, "physical_parameters_validated": False,
                "host_budget_accepted": False, "owner_implementation_verified": False,
                "limits": ["Manufactured stationary one-patch constant areal heat capacity only.",
                           "No conduction, radiation, phase change, wall motion, external power or gap/chemistry evolution.",
                           "Finite-step RK4 errors require independent oracle/refinement checks; no adaptive error certificate.",
                           "Unsupported positive-flow underflow, state stagnation, exhaustion and temperature-domain stages are rejected.",
                           "Translation is the continuum flux-Maxwellian expectation 2*k_B*T per parent; no sampled packet or host acceptance."]}


def _state(problem: FeedbackProblem, coverage: float, temperature: float) -> None:
    if not 0 <= coverage <= problem.reservoir.max_coverage_m2:
        raise UnsupportedFeedbackRangeError("coverage is outside the declared supported domain")
    if not problem.temperature_min_k <= temperature <= problem.temperature_max_k:
        raise UnsupportedFeedbackRangeError("temperature is outside the declared supported domain")


def _stage_increment(problem: FeedbackProblem, coverage: float, temperature: float, dt: float) -> tuple[float, ...]:
    _state(problem, coverage, temperature)
    if coverage == 0:
        return (0.0,)*6
    pool = problem.reservoir
    if pool.activation_energy_ev == 0:
        beta = 0.0
    else:
        complement = 1-_f(pool.coverage_coefficient)*_f(coverage)/_f(pool.reference_coverage_m2)
        beta = _project(_f(pool.activation_energy_ev)*_QE*complement/(_KB*_f(temperature)), "dimensionless kinetic barrier")
    log_amount = math.fsum((math.log(dt), math.log(pool.attempt_frequency_s), math.log(coverage), -beta))
    try:
        released = math.exp(log_amount)
    except OverflowError as exc:
        raise UnsupportedFeedbackRangeError("positive stage release overflows its representation") from exc
    if not math.isfinite(released) or released <= 0:
        raise UnsupportedFeedbackRangeError("mathematically positive stage release is not representable")
    amount = _f(released)
    binding = amount*_f(pool.desorption_energy_ev)*_QE
    internal = amount*_f(pool.internal_energy_ev)*_QE
    translation = amount*2*_KB*_f(temperature)
    wall = binding+internal+translation
    cooling = wall/_f(problem.areal_heat_capacity_j_m2_k)
    return (released, _project(cooling, "positive stage cooling"),
            _project(binding, "stage binding energy"), _project(internal, "stage internal energy"),
            _project(translation, "positive stage translation"), _project(wall, "positive stage wall debit"))


def feedback_rhs(problem: FeedbackProblem, *, coverage_m2: float, temperature_k: float) -> dict[str, float]:
    """Instantaneous derivatives, with the same explicit representation guards."""
    if not isinstance(problem, FeedbackProblem):
        raise ValueError("problem must be a FeedbackProblem")
    coverage = _real(coverage_m2, "coverage_m2")
    temperature = _real(temperature_k, "temperature_k", positive=True)
    dn, cooling, binding, internal, translation, wall = _stage_increment(problem, coverage, temperature, 1.0)
    return {"release_flux_m2_s": dn, "coverage_rate_m2_s": -dn, "temperature_rate_k_s": -cooling,
            "binding_power_j_m2_s": binding, "internal_power_j_m2_s": internal,
            "translational_power_j_m2_s": translation, "total_wall_power_j_m2_s": wall}


def _decrease(current: float, debit: Fraction, name: str, lower: Fraction, *, strict_lower: bool = False) -> float:
    value = _f(current)-debit
    if value < lower or (strict_lower and value == lower):
        raise UnsupportedFeedbackRangeError(f"RK4 {name} stage leaves its supported domain; reduce the step or duration")
    result = _project(value, name)
    if debit > 0 and result >= current:
        raise UnsupportedFeedbackRangeError(f"positive {name} debit cannot change its floating-point state")
    return result


def _stage_state(problem: FeedbackProblem, coverage: float, temperature: float, increment: tuple[float, ...], scale: Fraction) -> tuple[float, float]:
    n = _decrease(coverage, _f(increment[0])*scale, "coverage", Fraction(0), strict_lower=coverage > 0)
    t = _decrease(temperature, _f(increment[1])*scale, "temperature", _f(problem.temperature_min_k))
    _state(problem, n, t)
    return n, t


def _weighted(k1: tuple[float, ...], k2: tuple[float, ...], k3: tuple[float, ...], k4: tuple[float, ...]) -> tuple[float, ...]:
    return tuple(_project((_f(a)+2*_f(b)+2*_f(c)+_f(d))/6, "RK4 weighted increment")
                 for a, b, c, d in zip(k1, k2, k3, k4))


def _increase(current: float, amount: float, name: str) -> float:
    result = _project(_f(current)+_f(amount), name)
    if amount > 0 and result <= current:
        raise UnsupportedFeedbackRangeError(f"positive {name} increment cannot change its cumulative representation")
    return result


def _sample(problem: FeedbackProblem, time: float, coverage: float, temperature: float, cumulative: tuple[float, ...]) -> FeedbackSample:
    released, binding, internal, translation, wall = cumulative
    pool = problem.reservoir
    d = _f(pool.desorption_energy_ev)*_QE
    i = _f(pool.internal_energy_ev)*_QE
    n, m = _f(coverage), _f(released)
    thermal = _f(problem.areal_heat_capacity_j_m2_k)*(_f(temperature)-_f(problem.initial_temperature_k))
    ads_change = d*m
    ads_inventory_change = d*(_f(problem.initial_coverage_m2)-n)
    gas = _f(internal)+_f(translation)
    quantities = {
        "coverage_m2": coverage, "temperature_k": temperature, "released_m2": released,
        "binding_energy_j_m2": binding, "internal_energy_j_m2": internal,
        "translational_energy_j_m2": translation, "total_wall_energy_debit_j_m2": wall,
        "gas_energy_j_m2": _project(gas, "gas energy"),
        "wall_thermal_change_j_m2": _project(thermal, "wall thermal change"),
        "adsorbate_reference_energy_j_m2": _project(-d*n, "adsorbate reference energy"),
        "adsorbate_reference_change_j_m2": _project(ads_change, "adsorbate reference change"),
        "adsorbate_reference_inventory_change_j_m2": _project(ads_inventory_change, "inventory adsorbate reference change"),
        "inventory_residual_m2": _project(n+m-_f(problem.initial_coverage_m2), "inventory residual"),
        "adsorbate_reference_residual_j_m2": _project(ads_inventory_change-ads_change, "adsorbate reference residual"),
        "wall_debit_residual_j_m2": _project(thermal+_f(wall), "wall debit residual"),
        "wall_channel_residual_j_m2": _project(_f(wall)-_f(binding)-gas, "wall channel residual"),
        "binding_allocation_residual_j_m2": _project(_f(binding)-d*m, "binding allocation residual"),
        "internal_allocation_residual_j_m2": _project(_f(internal)-i*m, "internal allocation residual"),
        "total_energy_residual_j_m2": _project(thermal+ads_inventory_change+gas, "total energy residual"),
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
    return FeedbackSample(time_s=time, **quantities, extensive=extensive)


def evolve_feedback(problem: FeedbackProblem, *, duration_s: float, steps: int) -> FeedbackHistory:
    """Evolve N and T jointly with fixed-step RK4 and separate positive channels.

    Every intermediate stage must remain inside the declared thermal/coverage
    domain. Positive stage amounts, final debits and cumulative additions must
    be representable and change the relevant floating-point state. Failure is
    explicit; this routine never reduces a time step, pads an endpoint or clips.
    """
    if not isinstance(problem, FeedbackProblem):
        raise ValueError("problem must be a FeedbackProblem")
    duration = _real(duration_s, "duration_s", positive=True)
    if type(steps) is not int or not 1 <= steps <= MAX_STEPS:
        raise ValueError(f"steps must be an integer in [1,{MAX_STEPS}]")
    dt = _project(_f(duration)/steps, "positive fixed time step")
    coverage, temperature = problem.initial_coverage_m2, problem.initial_temperature_k
    cumulative = (0.0,)*5
    samples = [_sample(problem, 0.0, coverage, temperature, cumulative)]
    for index in range(1, steps+1):
        k1 = _stage_increment(problem, coverage, temperature, dt)
        n2, t2 = _stage_state(problem, coverage, temperature, k1, Fraction(1, 2))
        k2 = _stage_increment(problem, n2, t2, dt)
        n3, t3 = _stage_state(problem, coverage, temperature, k2, Fraction(1, 2))
        k3 = _stage_increment(problem, n3, t3, dt)
        n4, t4 = _stage_state(problem, coverage, temperature, k3, Fraction(1))
        k4 = _stage_increment(problem, n4, t4, dt)
        increment = _weighted(k1, k2, k3, k4)
        coverage, temperature = _stage_state(problem, coverage, temperature, increment, Fraction(1))
        additions = (increment[0], increment[2], increment[3], increment[4], increment[5])
        cumulative = tuple(_increase(old, amount, "cumulative channel") for old, amount in zip(cumulative, additions))
        time = _project(_f(duration)*index/steps, "sample time")
        if time <= samples[-1].time_s:
            raise UnsupportedFeedbackRangeError("sample time cannot advance its floating-point representation")
        samples.append(_sample(problem, time, coverage, temperature, cumulative))
    return FeedbackHistory(problem, duration, steps, tuple(samples))
