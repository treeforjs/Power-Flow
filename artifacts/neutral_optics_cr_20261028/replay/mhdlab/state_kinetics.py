"""Validated, state-resolved inputs for kinetic reference calculations.

This module contains no default water mechanism or population closure. Every
reaction names explicit state IDs; energies include formation and excitation
relative to ONE common reference chosen by the supplying dataset. An energy
measured from each species' separate ground state is insufficient for chemistry.

Quantum assignments are source-defined, including polyatomic mode tuples,
symmetry and parity. State IDs, rather than approximate assignments, identify
levels. Cross sections are integral channel cross sections versus relative
kinetic energy. A kinematics label and its provenance do not by themselves
supply angular distributions or product energy-sharing data.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
import hashlib
import json
import math
from numbers import Integral, Real
import re
from types import MappingProxyType
from typing import Any

import numpy as np


ELEMENTARY_CHARGE_C = 1.602176634e-19
RELATIVE_ENERGY_AXIS = "relative_kinetic_energy_ev"
MASS_RELATIVE_TOLERANCE = 1e-12
_UNREADY = re.compile(
    r"(?<![a-z0-9])(?:placeholder|provisional|todo|tbd|pending|fabricated)(?![a-z0-9])", re.I)


class CrossSectionDomainError(ValueError):
    """A physically possible query lies outside the supplied cross-section data."""


def _name(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value.strip() or value != value.strip():
        raise ValueError(f"{label} must be a nonempty string without surrounding whitespace")
    return value


def _finite(value: Any, label: str, *, positive: bool = False) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{label} must be a finite real number")
    number = float(value)
    if not math.isfinite(number) or (positive and number <= 0):
        qualifier = "positive finite" if positive else "finite"
        raise ValueError(f"{label} must be {qualifier}")
    return number


def _freeze(value: Any, label: str) -> Any:
    """Copy JSON-like assignments/metadata into recursively immutable values."""
    if isinstance(value, Mapping):
        return MappingProxyType({_name(key, f"{label} key"): _freeze(item, label)
                                 for key, item in value.items()})
    if isinstance(value, (list, tuple)):
        return tuple(_freeze(item, label) for item in value)
    if value is None or isinstance(value, (str, bool)):
        return value
    if isinstance(value, Integral):
        return int(value)
    if isinstance(value, Real):
        return _finite(value, label)
    raise ValueError(f"{label} values must be finite, JSON-like data")


def _plain(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {key: _plain(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_plain(item) for item in value]
    return value


def _unready(value: Any) -> bool:
    if isinstance(value, Mapping):
        return any(_unready(item) for item in value.values())
    if isinstance(value, tuple):
        return any(_unready(item) for item in value)
    return isinstance(value, str) and bool(_UNREADY.search(value))


@dataclass(frozen=True, slots=True)
class Species:
    """A chemical/isotopic species; mass convention must match all reactions.

    ``nuclei`` maps source-defined nuclear identities, e.g. ``H``/``O`` or
    ``1H``/``16O``, to positive counts. Electrons have an empty mapping.
    """

    name: str
    mass_kg: float
    charge_number: int
    nuclei: Mapping[str, int]

    def __post_init__(self) -> None:
        _name(self.name, "species name")
        object.__setattr__(self, "mass_kg", _finite(self.mass_kg, "mass_kg", positive=True))
        if isinstance(self.charge_number, (bool, np.bool_)) or not isinstance(self.charge_number, Integral):
            raise ValueError("charge_number must be an integer")
        object.__setattr__(self, "charge_number", int(self.charge_number))
        if not isinstance(self.nuclei, Mapping):
            raise ValueError("nuclei must be a mapping")
        nuclei = {}
        for key, count in self.nuclei.items():
            _name(key, "nuclear identity")
            if isinstance(count, (bool, np.bool_)) or not isinstance(count, Integral) or count <= 0:
                raise ValueError("nuclear counts must be positive integers")
            nuclei[key] = int(count)
        object.__setattr__(self, "nuclei", MappingProxyType(nuclei))

    @property
    def charge_c(self) -> float:
        return self.charge_number * ELEMENTARY_CHARGE_C


@dataclass(frozen=True, slots=True)
class InternalState:
    """One identified state; energy_ev includes formation AND excitation."""

    state_id: str
    species: str
    energy_ev: float
    quantum_numbers: Mapping[str, Any]

    def __post_init__(self) -> None:
        _name(self.state_id, "state_id")
        _name(self.species, "state species")
        object.__setattr__(self, "energy_ev", _finite(self.energy_ev, "state energy_ev"))
        if not isinstance(self.quantum_numbers, Mapping):
            raise ValueError("quantum_numbers must be a mapping")
        object.__setattr__(self, "quantum_numbers", _freeze(self.quantum_numbers, "quantum_numbers"))


def _indexed(values, kind, key):
    if isinstance(values, Mapping):
        result = dict(values)
        for name, item in result.items():
            if not isinstance(item, kind) or name != getattr(item, key):
                raise ValueError(f"catalog mapping keys must match {kind.__name__}.{key}")
        return result
    result = {}
    for item in values:
        if not isinstance(item, kind):
            raise ValueError(f"catalog requires {kind.__name__} records")
        name = getattr(item, key)
        if name in result:
            raise ValueError(f"duplicate {kind.__name__} identifier: {name}")
        result[name] = item
    return result


@dataclass(frozen=True, slots=True, init=False)
class StateCatalog:
    """Immutable catalog with exact state lookup and a content fingerprint.

    ``species`` and ``states`` accept iterables or mappings keyed by identifier.
    Each supplied species needs at least one supplied state. Molecular states
    require explicit assignments, without prescribing a diatomic label scheme.
    Construction validates structure, not the accuracy of the source physics.
    """

    species: Mapping[str, Species]
    states: Mapping[str, InternalState]
    fingerprint: str

    def __init__(self, species: Iterable[Species] | Mapping[str, Species],
                 states: Iterable[InternalState] | Mapping[str, InternalState]):
        species_map = _indexed(species, Species, "name")
        state_map = _indexed(states, InternalState, "state_id")
        if not species_map or not state_map:
            raise ValueError("catalog requires nonempty species and states")
        used = set()
        for state in state_map.values():
            if state.species not in species_map:
                raise ValueError(f"unknown species {state.species!r} for state {state.state_id!r}")
            if sum(species_map[state.species].nuclei.values()) > 1 and not state.quantum_numbers:
                raise ValueError(f"molecular state {state.state_id!r} requires quantum assignments")
            used.add(state.species)
        if missing := species_map.keys() - used:
            raise ValueError(f"species without explicit states: {sorted(missing)}")
        content = {
            "species": [{"name": item.name, "mass_kg": item.mass_kg,
                         "charge_number": item.charge_number, "nuclei": dict(item.nuclei)}
                        for _, item in sorted(species_map.items())],
            "states": [{"state_id": item.state_id, "species": item.species,
                        "energy_ev": item.energy_ev, "quantum_numbers": _plain(item.quantum_numbers)}
                       for _, item in sorted(state_map.items())],
        }
        fingerprint = hashlib.sha256(json.dumps(content, sort_keys=True, separators=(",", ":"),
                                                allow_nan=False).encode()).hexdigest()
        object.__setattr__(self, "species", MappingProxyType(species_map))
        object.__setattr__(self, "states", MappingProxyType(state_map))
        object.__setattr__(self, "fingerprint", fingerprint)

    def state(self, state_id: str) -> InternalState:
        try:
            return self.states[state_id]
        except KeyError:
            raise KeyError(f"unknown explicit state ID {state_id!r}; species fallback is forbidden") from None

    def __getitem__(self, state_id: str) -> InternalState:
        return self.state(state_id)

    def species_for_state(self, state_id: str) -> Species:
        return self.species[self.state(state_id).species]

    def energy_ev(self, state_id: str) -> float:
        return self.state(state_id).energy_ev

    def mass_kg(self, state_id: str) -> float:
        return self.species_for_state(state_id).mass_kg

    def charge_number(self, state_id: str) -> int:
        return self.species_for_state(state_id).charge_number

    def require_species(self, names: Iterable[str]) -> None:
        if missing := set(names) - self.species.keys():
            raise ValueError(f"required species absent from catalog: {sorted(missing)}")


@dataclass(frozen=True, slots=True)
class ReactionChannel:
    """Explicit chemical branch with finite cross-section support.

    ``energy_ev`` is relative kinetic energy, NOT incident laboratory energy.
    ``metadata['kinematics_provenance']`` is required for production use. Store
    angular/energy-sharing data or their explicit physical approximation there;
    no product velocities are inferred by this input contract.

    ``validate(production=False)`` allows clearly labeled reference fixtures.
    Disabled channels are still structurally/conservatively validated, but are
    never evaluated by ``sigma_at``. Partial ionic yields and total ionization
    observables are not automatically exclusive chemical branches.
    """

    name: str
    reactants: tuple[str, ...]
    products: tuple[str, ...]
    energy_ev: tuple[float, ...]
    sigma_m2: tuple[float, ...]
    kinematics_model: str
    provenance: str
    metadata: Mapping[str, Any] = field(default_factory=dict)
    energy_axis: str = RELATIVE_ENERGY_AXIS
    enabled: bool = True
    provisional: bool = False
    catalog_fingerprint: str | None = None

    def __post_init__(self) -> None:
        for key in ("name", "kinematics_model", "provenance"):
            _name(getattr(self, key), key)
        if self.energy_axis != RELATIVE_ENERGY_AXIS:
            raise ValueError(f"energy_axis must be {RELATIVE_ENERGY_AXIS!r}")
        for key in ("reactants", "products"):
            raw = getattr(self, key)
            if isinstance(raw, str):
                raise ValueError(f"{key} must be a sequence of explicit state IDs")
            entries = tuple(raw)
            if not entries:
                raise ValueError(f"{key} cannot be empty")
            for state_id in entries:
                _name(state_id, f"{key} state ID")
            object.__setattr__(self, key, entries)
        for key in ("enabled", "provisional"):
            if not isinstance(getattr(self, key), bool):
                raise ValueError(f"{key} must be boolean")
        if self.catalog_fingerprint is not None:
            if not isinstance(self.catalog_fingerprint, str) or not re.fullmatch(r"[0-9a-f]{64}", self.catalog_fingerprint):
                raise ValueError("catalog_fingerprint must be a lowercase SHA-256 digest")
        for key in ("energy_ev", "sigma_m2"):
            raw = getattr(self, key)
            if isinstance(raw, str):
                raise ValueError(f"{key} must be a numeric sequence")
            entries = tuple(_finite(value, key) for value in raw)
            if len(entries) < 2 or min(entries) < 0:
                raise ValueError(f"{key} needs at least two nonnegative finite samples")
            object.__setattr__(self, key, entries)
        if len(self.energy_ev) != len(self.sigma_m2):
            raise ValueError("energy_ev and sigma_m2 lengths differ")
        if any(right <= left for left, right in zip(self.energy_ev, self.energy_ev[1:])):
            raise ValueError("energy_ev must be strictly increasing")
        if not isinstance(self.metadata, Mapping):
            raise ValueError("metadata must be a mapping")
        object.__setattr__(self, "metadata", _freeze(self.metadata, "reaction metadata"))

    def energy_defect_ev(self, catalog: StateCatalog) -> float:
        """Products minus reactants on the supplied common energy reference."""
        try:
            defect = math.fsum(catalog.energy_ev(item) for item in self.products) - math.fsum(
                catalog.energy_ev(item) for item in self.reactants)
        except OverflowError as error:
            raise ValueError("reaction state-energy sum is not finite") from error
        return _finite(defect, "reaction energy defect")

    def threshold_ev(self, catalog: StateCatalog) -> float:
        return max(0.0, self.energy_defect_ev(catalog))

    def validate(self, catalog: StateCatalog, *, production: bool = True) -> None:
        if self.catalog_fingerprint is not None and self.catalog_fingerprint != catalog.fingerprint:
            raise ValueError("reaction catalog fingerprint does not match state catalog")
        # Lookup every ID, even when the same chemical species occurs on both sides.
        reactants = [catalog.species_for_state(item) for item in self.reactants]
        products = [catalog.species_for_state(item) for item in self.products]
        reactant_nuclei, product_nuclei = Counter(), Counter()
        for item in reactants:
            reactant_nuclei.update(item.nuclei)
        for item in products:
            product_nuclei.update(item.nuclei)
        if reactant_nuclei != product_nuclei:
            raise ValueError("reaction does not conserve nuclei")
        if sum(item.charge_number for item in reactants) != sum(item.charge_number for item in products):
            raise ValueError("reaction does not conserve charge")
        mass_in = math.fsum(item.mass_kg for item in reactants)
        mass_out = math.fsum(item.mass_kg for item in products)
        if not math.isclose(mass_in, mass_out, rel_tol=MASS_RELATIVE_TOLERANCE, abs_tol=0):
            raise ValueError("reaction does not conserve mass under the catalog mass convention")
        threshold = self.threshold_ev(catalog)
        if any(energy < threshold and sigma != 0 for energy, sigma in zip(self.energy_ev, self.sigma_m2)):
            raise ValueError("nonzero cross section below the state-energy reaction threshold")
        if self.enabled and production:
            if self.provisional or _unready((self.name, self.provenance, self.kinematics_model, self.metadata)):
                raise ValueError("enabled provisional/placeholder reaction is forbidden in production")
            if self.metadata.get("exclusive_reaction") is False or self.metadata.get("is_exclusive_reaction") is False:
                raise ValueError("nonexclusive yield data cannot be used as a chemical branch")
            kind = str(self.metadata.get("data_kind", self.metadata.get("process", ""))).lower()
            if "partial_ion" in kind or "total_ionization" in kind or "ionic_yield" in kind:
                raise ValueError("ionization yield/total data require validated exclusive branching")
            _name(self.metadata.get("kinematics_provenance"), "metadata kinematics_provenance")

    def sigma_at(self, relative_energy_ev: Any, catalog: StateCatalog):
        """Piecewise-linear sigma; zero only below the physical state threshold.

        Domain checking is independent of production readiness. Call
        ``validate(catalog, production=True)`` before enabling a production run.
        """
        self.validate(catalog, production=False)
        if not self.enabled:
            raise ValueError("cannot evaluate a disabled reaction channel")
        try:
            raw_query = np.asarray(relative_energy_ev)
            if raw_query.dtype.kind not in "iuf":
                raise ValueError("relative energy must contain real numbers")
            query = raw_query.astype(float)
        except (ValueError, TypeError) as error:
            raise ValueError("relative energy must be a finite nonnegative scalar or vector") from error
        if query.ndim > 1 or query.size == 0 or not np.all(np.isfinite(query)) or np.any(query < 0):
            raise ValueError("relative energy must be a finite nonnegative scalar or vector")
        threshold = self.threshold_ev(catalog)
        below = query < threshold
        outside = (~below) & ((query < self.energy_ev[0]) | (query > self.energy_ev[-1]))
        if np.any(outside):
            raise CrossSectionDomainError(
                f"{self.name}: cross section unsupported above threshold outside "
                f"[{self.energy_ev[0]}, {self.energy_ev[-1]}] eV")
        result = np.where(below, 0.0, np.interp(query, self.energy_ev, self.sigma_m2))
        return float(result) if result.ndim == 0 else result


@dataclass(frozen=True, slots=True)
class RadiativeTransition:
    """One bound-bound spontaneous transition; no inferred shell populations."""

    upper_state: str
    lower_state: str
    a_per_s: float
    provenance: str

    def __post_init__(self) -> None:
        _name(self.upper_state, "upper_state")
        _name(self.lower_state, "lower_state")
        _name(self.provenance, "radiative provenance")
        object.__setattr__(self, "a_per_s", _finite(self.a_per_s, "a_per_s", positive=True))

    def validate(self, catalog: StateCatalog) -> None:
        upper, lower = catalog.state(self.upper_state), catalog.state(self.lower_state)
        if upper.species != lower.species:
            raise ValueError("bound-bound radiative states must belong to the same species")
        _finite(upper.energy_ev - lower.energy_ev, "spontaneous photon energy", positive=True)

    def photon_energy_ev(self, catalog: StateCatalog) -> float:
        self.validate(catalog)
        return catalog.energy_ev(self.upper_state) - catalog.energy_ev(self.lower_state)
