"""Reservoir, area and diagnostic-energy association for one source window.

Mathematical log amounts remain distinct from representable float64 budgets.
No inventory is committed and no gap, chemistry or photon gate is advanced.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
import math
from numbers import Real
from types import MappingProxyType

import numpy as np

from .constants import KB, QE
from .temkin_emission_window import ReducedTemkinEmissionWindow
from .temkin_release import LAW_IDENTITY, TemkinReservoir, _thermal_parameters
from .thermal_release import _log_hazard

SCHEMA = "mhdlab-temkin-source-window-v1"
STATUSES = frozenset({"manufactured", "reconstructed", "sourced", "unvalidated", "owner_supplied"})
LD = np.longdouble
EPS = float(np.finfo(float).eps)


def _positive(value, name):
    if np.ma.isMaskedArray(value) or isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be an explicit finite positive real scalar")
    try:
        value = float(value)
    except (ValueError, OverflowError) as error:
        raise ValueError(f"{name} must be representable") from error
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be an explicit finite positive real scalar")
    return value


@dataclass(frozen=True, slots=True)
class _Representation:
    value: float | None
    log_value: float | None
    exact_zero: bool
    qualified: bool
    reason: str
    log_rounding_error: float | None
    log_rounding_allowance: float


def _representation(log_value, operation_scale, tolerance):
    if log_value is None:
        return _Representation(0., None, True, True, "exact_zero_supplied_energy_cost", 0., 0.)
    logarithm = LD(log_value)
    allowance = float(32*LD(EPS)*(1+LD(operation_scale)))
    if not np.isfinite(logarithm) or abs(logarithm) > np.finfo(float).max:
        return _Representation(None, None, False, False, "log_amount_outside_float64_range", None, allowance)
    reported_log = float(logarithm)
    if logarithm > np.log(LD(np.finfo(float).max)):
        return _Representation(None, reported_log, False, False, "positive_amount_overflow", None, allowance)
    with np.errstate(under="ignore", over="ignore", invalid="ignore"):
        value = float(np.exp(logarithm))
    if value == 0:
        return _Representation(0., reported_log, False, False, "positive_amount_underflow", None, allowance)
    if not math.isfinite(value):
        return _Representation(None, reported_log, False, False, "positive_amount_overflow", None, allowance)
    rounding = float(abs(np.log(LD(value))-logarithm))
    qualified = rounding+allowance <= tolerance
    return _Representation(value, reported_log, False, qualified,
                           "represented_within_declared_log_allowance" if qualified else "log_rounding_allowance_exceeded",
                           rounding, allowance)


@dataclass(frozen=True, slots=True, eq=False)
class FrozenTemkinSourceWindow:
    """Nonmutating source-window amounts for one supplied stationary patch.

Initial density, temperature and actual patch area are strictly positive.
Empty inventories/windows must be handled explicitly by the source caller.
Half the supplied log allowance is allocated to time-law/mass evaluation and
half to number/energy representation. These are numerical qualifications with
stated floating-point assumptions, not rigorous interval certificates.
"""

    reservoir: TemkinReservoir
    initial_coverage_m2: float
    temperature_k: float
    parent_duration_s: float
    start_s: float
    stop_s: float
    patch_area_m2: float
    log_mass_tolerance: float = field(kw_only=True)
    status: str = field(kw_only=True)
    provenance: str = field(kw_only=True)
    time_law: ReducedTemkinEmissionWindow = field(init=False, repr=False)
    _quantities: object = field(init=False, repr=False)

    def __post_init__(self):
        if not isinstance(self.reservoir, TemkinReservoir):
            raise ValueError("reservoir must be a supplied TemkinReservoir")
        for name in ("initial_coverage_m2", "temperature_k", "parent_duration_s", "patch_area_m2", "log_mass_tolerance"):
            object.__setattr__(self, name, _positive(getattr(self, name), name))
        if not isinstance(self.status, str) or self.status not in STATUSES:
            raise ValueError("status must be an explicit supported input label")
        if not isinstance(self.provenance, str) or not self.provenance.strip():
            raise ValueError("provenance must be a nonempty source association")
        if self.initial_coverage_m2 > self.reservoir.max_coverage_m2:
            raise ValueError("initial density exceeds the reservoir's supported coverage domain")
        if self.reservoir.coverage_coefficient == 0 or self.reservoir.activation_energy_ev == 0:
            coupling = 0.
            log_rate = _log_hazard(self.reservoir._energy_record(), self.temperature_k, 1.)
        else:
            coupling, barrier = _thermal_parameters(self.reservoir, self.initial_coverage_m2, self.temperature_k)
            log_rate = math.fsum((math.log(self.reservoir.attempt_frequency_s), -barrier))
        law = ReducedTemkinEmissionWindow(coupling, log_rate, self.parent_duration_s, self.start_s, self.stop_s,
                                          log_mass_tolerance=self.log_mass_tolerance/2)
        object.__setattr__(self, "start_s", law.start_s)
        object.__setattr__(self, "stop_s", law.stop_s)
        object.__setattr__(self, "time_law", law)
        initial_log = np.log(LD(self.initial_coverage_m2))
        area_log = np.log(LD(self.patch_area_m2))
        start = LD(law.start_depletion)
        delta = LD(law.local_depletion)
        source_log = LD(law.log_source_fraction)
        density_logs = {
            "initial_coverage_m2": (initial_log, abs(initial_log)),
            "remaining_start_m2": (initial_log-start, abs(initial_log)+start),
            "remaining_stop_m2": (initial_log-start-delta, abs(initial_log)+start+delta),
            "released_m2": (initial_log+source_log, abs(initial_log)+abs(source_log)),
        }
        number_keys = {"initial_coverage_m2":"initial_number", "remaining_start_m2":"remaining_start_number",
                       "remaining_stop_m2":"remaining_stop_number", "released_m2":"released_number"}
        records = {}
        half_allowance = self.log_mass_tolerance/2
        for key, (log_amount, scale) in density_logs.items():
            records[key] = _representation(log_amount, scale, half_allowance)
            records[number_keys[key]] = _representation(log_amount+area_log, scale+abs(area_log), half_allowance)
        # Declared gas internal energy and wall removal cost are independent of
        # the kinetic activation barrier. Translation is the flux mean 2*kB*T.
        # Form positive costs in logs before scaling by the released population.
        # A platform's long double may share float64's exponent range; a tiny
        # per-parent cost can still give a representable macroscopic energy.
        translation_log = np.log(LD(2))+np.log(LD(KB))+np.log(LD(self.temperature_k))
        internal_log = (np.log(LD(QE))+np.log(LD(self.reservoir.internal_energy_ev))
                        if self.reservoir.internal_energy_ev > 0 else None)
        binding_log = (np.log(LD(QE))+np.log(LD(self.reservoir.desorption_energy_ev))
                       if self.reservoir.desorption_energy_ev > 0 else None)
        positive_logs = [value for value in (translation_log,internal_log,binding_log) if value is not None]
        pivot = max(positive_logs)
        with np.errstate(under="ignore"):
            total_log = pivot+np.log(sum(np.exp(value-pivot) for value in positive_logs))
        cost_logs = {"translational_energy":translation_log, "internal_energy":internal_log,
                     "binding_energy":binding_log, "total_wall_energy_debit":total_log}
        release_log, release_scale = density_logs["released_m2"]
        for name, cost_log in cost_logs.items():
            energy_log = release_log+cost_log if cost_log is not None else None
            scale = release_scale+(abs(cost_log) if cost_log is not None else 0)
            records[name+"_j_m2"] = _representation(energy_log, scale, half_allowance)
            records[name+"_j"] = _representation(energy_log+area_log if energy_log is not None else None,
                                                 scale+abs(area_log), half_allowance)
        # A nonrepresentable initial parent budget cannot be an accepted source
        # association, even if a later mathematical partial amount is finite.
        if not records["initial_coverage_m2"].qualified or not records["initial_number"].qualified:
            raise ValueError("initial parent density/area number is outside the requested representation qualification")
        object.__setattr__(self, "_quantities", MappingProxyType(records))

    @property
    def all_budget_values_qualified(self):
        return all(record.qualified for record in self._quantities.values())

    def _balance(self, suffix):
        names = ("remaining_stop_"+suffix, "remaining_start_"+suffix, "released_"+suffix)
        values = [self._quantities[name].value for name in names]
        if any(value is None for value in values):
            return None
        return math.fsum((values[0], -values[1], values[2]))

    def _energy_balance(self, suffix):
        names = ("total_wall_energy_debit_"+suffix,"translational_energy_"+suffix,
                 "internal_energy_"+suffix,"binding_energy_"+suffix)
        values = [self._quantities[name].value for name in names]
        if any(value is None for value in values):
            return None
        return math.fsum((values[0],-values[1],-values[2],-values[3]))

    def _units(self):
        return {name:("J m^-2" if name.endswith("_j_m2") else "J" if name.endswith("_j")
                      else "parent" if name.endswith("_number") else "parent m^-2")
                for name in self._quantities}

    def manifest(self):
        """Return detached associations, logs, represented values and limits."""
        return {
            "schema":SCHEMA, "status":self.status, "provenance":self.provenance,
            "reservoir":asdict(self.reservoir), "law_identity":LAW_IDENTITY,
            **{name:getattr(self,name) for name in ("initial_coverage_m2", "temperature_k", "parent_duration_s",
                                                   "start_s", "stop_s", "patch_area_m2", "log_mass_tolerance")},
            "time_window":self.time_law.manifest(),
            "quantities":{name:asdict(record) for name,record in self._quantities.items()},
            "quantity_units":self._units(),
            "all_budget_values_qualified":self.all_budget_values_qualified,
            "inventory_balance_residual_m2":self._balance("m2"),
            "inventory_balance_residual_number":self._balance("number"),
            "energy_balance_residual_j_m2":self._energy_balance("j_m2"),
            "energy_balance_residual_j":self._energy_balance("j"),
            "host_budget_accepted":False,
            "numerical_allowance_allocation":{"time_law_log_mass":self.log_mass_tolerance/2,
                                               "representation_log_amount":self.log_mass_tolerance/2,
                                               "assumption":"Extended intermediates plus conservative float64 elementary-operation allowance; not a certified interval enclosure."},
            "limits":[
                "Partial source emission in one stationary frozen-temperature patch; no photons, transport or chemistry are evolved.",
                "Density, actual patch number and diagnostic energy each retain their own representation flags; a finite log is not an accepted budget.",
                "Initial independent pool/state inventory remains attached to the supplied reservoir; no new parent, state branch or daughter yield is created.",
                "Signed inventory rounding residuals are retained and require the owner ledger's explicit treatment; no automatic residual repair or transaction is performed.",
                "Translation is expected outgoing flux energy 2*kB*T. The supplied internal/removal costs are separate from the kinetic barrier; no wall cooling feedback or individual velocity sample is calculated.",
                "Source start and local delay remain separate from any global timestamp or recorded photon gate.",
            ],
        }

    def represented_budget(self):
        """Return qualified represented fields; this performs no host acceptance."""
        if not self.all_budget_values_qualified:
            rejected = [name for name,record in self._quantities.items() if not record.qualified]
            raise ValueError("source-window quantities are not represented within the declared allowance: "+", ".join(rejected))
        return {"schema":"mhdlab-temkin-window-represented-budget-v1",
                "values":{name:record.value for name,record in self._quantities.items()},
                "quantity_units":self._units(),
                "inventory_balance_residual_m2":self._balance("m2"),
                "inventory_balance_residual_number":self._balance("number"),
                "energy_balance_residual_j_m2":self._energy_balance("j_m2"),
                "energy_balance_residual_j":self._energy_balance("j"),
                "source_association":{"reservoir":asdict(self.reservoir),"temperature_k":self.temperature_k,
                                      "patch_area_m2":self.patch_area_m2,"start_s":self.start_s,"stop_s":self.stop_s,
                                      "status":self.status,"provenance":self.provenance},
                "host_budget_accepted":False}
