"""Deterministic, nonrelativistic state-resolved reaction kinematics.

``react`` maps explicitly supplied Jacobi-mode energy fractions and directions
onto momenta for two incoming massive particles and N >= 2 outgoing particles.
It does not sample a differential cross section, select a chemical branch, or
invent fragment/secondary-electron distributions. The supplying collision kernel
must establish the validity, energy frame and provenance of those distributions.

For ordered products (0, ..., N-1), mode i is the relative motion of product i
against the center of mass of products (i+1, ..., N-1). In particular a mode's
energy fraction is NOT the kinetic-energy fraction of product i. Product order
is part of the kernel's definition. Ionization's indistinguishable-electron
labeling, exchange convention and secondary-energy law remain kernel inputs.

Masses use the catalog's mass-conserving nonrelativistic convention; state
energies include chemical formation and internal excitation on a common
reference. Photons, surface recoil, relativistic reactions and nuclear mass
defects are unsupported. A caller-controlled speed cap enforces the declared
nonrelativistic regime. A relativistic orbit pusher does not relax this cap.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import math
from types import MappingProxyType
from typing import Any

import numpy as np

from .state_kinetics import (
    ELEMENTARY_CHARGE_C,
    MASS_RELATIVE_TOLERANCE,
    ReactionChannel,
    StateCatalog,
)


SPEED_OF_LIGHT_M_S = 299_792_458.0
SUPPORTED_KINEMATICS_MODEL = "jacobi_energy_partition"
_NORMALIZATION_TOLERANCE = 1e-12
_ROUNDOFF_FACTOR = 64.0 * np.finfo(float).eps


class InaccessibleReactionError(ValueError):
    """Incoming relative energy cannot supply the specified product states."""


class NonrelativisticLimitError(ValueError):
    """An incoming or outgoing speed exceeds this constructor's validity cap."""


@dataclass(frozen=True, slots=True)
class ReactionOutcome:
    """One ordered event and its laboratory/COM energy and momentum ledger.

    Residuals mean after minus before. ``energy_residual_j`` is evaluated as a
    compensated sum of the kinetic-energy change and chemical energy defect,
    avoiding subtraction of potentially large common internal-energy offsets.
    ``threshold_roundoff_correction_j`` records any tiny negative available
    energy treated as zero; it is not an adjustable physical energy floor.
    Arrays are copied and made read-only by the constructor.
    """

    reactant_state_ids: tuple[str, ...]
    product_state_ids: tuple[str, ...]
    outgoing_momenta_kg_m_s: np.ndarray
    available_energy_j: float
    incoming_com_kinetic_energy_j: float
    outgoing_com_kinetic_energy_j: float
    kinetic_energy_before_j: float
    kinetic_energy_after_j: float
    internal_energy_before_j: float
    internal_energy_after_j: float
    total_energy_before_j: float
    total_energy_after_j: float
    energy_defect_ev: float
    energy_residual_j: float
    momentum_before_kg_m_s: np.ndarray
    momentum_after_kg_m_s: np.ndarray
    momentum_residual_kg_m_s: np.ndarray
    mass_residual_kg: float
    center_of_mass_velocity_m_s: np.ndarray
    jacobi_energy_fractions: np.ndarray
    jacobi_directions: np.ndarray
    maximum_incoming_beta: float
    maximum_outgoing_beta: float
    threshold_roundoff_correction_j: float
    metadata: Mapping[str, Any]


def _real_array(value: Any, shape: tuple[int, ...], label: str) -> np.ndarray:
    if np.iscomplexobj(value):
        raise ValueError(f"{label} must contain real numbers")
    try:
        result = np.array(value, dtype=float, copy=True)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{label} must be a finite real array of shape {shape}") from exc
    if result.shape != shape or not np.all(np.isfinite(result)):
        raise ValueError(f"{label} must be a finite real array of shape {shape}")
    return result


def _readonly(value: np.ndarray) -> np.ndarray:
    result = np.array(value, dtype=float, copy=True)
    result.flags.writeable = False
    return result


def _vector_sum(vectors: np.ndarray) -> np.ndarray:
    return np.array([math.fsum(float(item) for item in vectors[:, axis]) for axis in range(3)])


def _kinetic_energy(momenta: np.ndarray, masses: np.ndarray) -> float:
    return math.fsum(float(np.dot(momentum, momentum) / (2.0 * mass))
                     for momentum, mass in zip(momenta, masses))


def _maximum_beta(velocities: np.ndarray, maximum: float, label: str) -> float:
    # hypot is more resistant to intermediate overflow than sum(v*v).
    speed = max(math.hypot(*(float(component) for component in velocity))
                for velocity in velocities)
    beta = speed / SPEED_OF_LIGHT_M_S
    if not math.isfinite(beta) or beta > maximum:
        raise NonrelativisticLimitError(
            f"{label} speed beta={beta:.8g} exceeds max_nonrel_beta={maximum:.8g}; "
            "relativistic reaction kinematics are not implemented")
    return beta


def react(
    channel: ReactionChannel,
    catalog: StateCatalog,
    incoming_momenta_kg_m_s: Sequence[Sequence[float]] | np.ndarray,
    *,
    jacobi_energy_fractions: Sequence[float] | np.ndarray,
    jacobi_directions: Sequence[Sequence[float]] | np.ndarray,
    max_nonrel_beta: float = 0.05,
) -> ReactionOutcome:
    """Construct a conservative 2->N event from externally supplied kinematics.

    Fractions must be nonnegative and sum to one; each direction must be a unit
    vector. Normalization errors no larger than 1e-12 are accepted as roundoff
    and normalized before construction. There is no random or equipartition
    default. The outgoing tuple retains ``channel.products`` order.

    ``channel.validate(production=True)`` checks nuclei, charge, masses and
    provenance. The caller must separately select a physically possible event
    from ``channel.sigma_at`` and a supported differential kernel. This function
    does not test or extend the external cross-section/kernel energy domain.
    """
    channel.validate(catalog, production=True)
    if not channel.enabled:
        raise ValueError("disabled reaction channels cannot construct events")
    if channel.kinematics_model != SUPPORTED_KINEMATICS_MODEL:
        raise ValueError(
            f"unsupported kinematics_model {channel.kinematics_model!r}; "
            f"expected {SUPPORTED_KINEMATICS_MODEL!r}")
    reactants, products = tuple(channel.reactants), tuple(channel.products)
    if len(reactants) != 2 or len(products) < 2:
        raise ValueError("this constructor requires exactly two reactants and at least two massive products")
    if (isinstance(max_nonrel_beta, (bool, np.bool_))
            or not np.isscalar(max_nonrel_beta)):
        raise ValueError("max_nonrel_beta must be finite and strictly between zero and one")
    try:
        maximum_beta = float(max_nonrel_beta)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("max_nonrel_beta must be finite and strictly between zero and one") from exc
    if not math.isfinite(maximum_beta) or not 0.0 < maximum_beta < 1.0:
        raise ValueError("max_nonrel_beta must be finite and strictly between zero and one")

    incoming = _real_array(incoming_momenta_kg_m_s, (2, 3), "incoming_momenta_kg_m_s")
    mode_count = len(products) - 1
    fractions = _real_array(jacobi_energy_fractions, (mode_count,), "jacobi_energy_fractions")
    if np.any(fractions < 0):
        raise ValueError("jacobi_energy_fractions must be nonnegative")
    fraction_sum = math.fsum(float(item) for item in fractions)
    if abs(fraction_sum - 1.0) > _NORMALIZATION_TOLERANCE:
        raise ValueError("jacobi_energy_fractions must sum to one within 1e-12")
    fractions /= fraction_sum
    directions = _real_array(jacobi_directions, (mode_count, 3), "jacobi_directions")
    norms = np.array([math.hypot(*(float(component) for component in direction))
                      for direction in directions])
    if not np.all(np.abs(norms - 1.0) <= _NORMALIZATION_TOLERANCE):
        raise ValueError("jacobi_directions must be unit vectors within 1e-12")
    directions /= norms[:, None]

    masses_in = np.array([catalog.mass_kg(state) for state in reactants], dtype=float)
    masses_out = np.array([catalog.mass_kg(state) for state in products], dtype=float)
    if (not np.all(np.isfinite(masses_in)) or not np.all(np.isfinite(masses_out))
            or np.any(masses_in <= 0) or np.any(masses_out <= 0)):
        raise ValueError("all reactant/product masses must be positive finite massive-particle masses")
    mass_in, mass_out = math.fsum(masses_in), math.fsum(masses_out)
    mass_residual = mass_out - mass_in
    if abs(mass_residual) > MASS_RELATIVE_TOLERANCE * max(mass_in, mass_out):
        raise ValueError("reaction does not conserve the catalog's nonrelativistic additive mass")
    velocities_in = incoming / masses_in[:, None]
    incoming_beta = _maximum_beta(velocities_in, maximum_beta, "incoming")
    momentum_before = _vector_sum(incoming)
    velocity_com_in = momentum_before / mass_in
    velocity_com_out = momentum_before / mass_out
    relative_velocity = velocities_in[0] - velocities_in[1]
    reduced_mass_in = masses_in[0] * (masses_in[1] / mass_in)
    com_energy_before = 0.5 * reduced_mass_in * float(np.dot(relative_velocity, relative_velocity))
    defect_ev = float(channel.energy_defect_ev(catalog))
    if not math.isfinite(defect_ev):
        raise ValueError("reaction energy defect must be finite")
    defect_j = defect_ev * ELEMENTARY_CHARGE_C
    # Accepted floating-point mass differences are small, but correcting their
    # bulk-energy effect also conserves the laboratory ledger for moving COMs.
    # This expression avoids subtracting two large bulk kinetic energies.
    bulk_mass_roundoff_correction = (
        0.5 * mass_in * float(np.dot(velocity_com_in, velocity_com_in))
        * (mass_residual / mass_out)
    )
    available = math.fsum((com_energy_before, -defect_j, bulk_mass_roundoff_correction))
    threshold_tolerance = _ROUNDOFF_FACTOR * max(
        com_energy_before, abs(defect_j), abs(bulk_mass_roundoff_correction))
    if available < -threshold_tolerance:
        raise InaccessibleReactionError(
            "incoming relative energy cannot access product states: "
            f"available_energy_ev={available / ELEMENTARY_CHARGE_C:.12g}")
    threshold_roundoff = max(0.0, -available)
    available = max(0.0, available)

    # The orthogonal Jacobi modes sum to the available COM kinetic energy.
    # At each recursion, separate particle i from the remaining group, while
    # preserving that group's previously assigned center-of-mass velocity.
    velocities_com = np.zeros((len(products), 3), dtype=float)
    group_velocity = np.zeros(3, dtype=float)
    for index in range(mode_count):
        remaining_mass = math.fsum(masses_out[index + 1:])
        group_mass = masses_out[index] + remaining_mass
        reduced_mass = masses_out[index] * (remaining_mass / group_mass)
        mode_speed = math.sqrt((2.0 * available * fractions[index]) / reduced_mass)
        relative = mode_speed * directions[index]
        velocities_com[index] = group_velocity + (remaining_mass / group_mass) * relative
        group_velocity = group_velocity - (masses_out[index] / group_mass) * relative
    velocities_com[-1] = group_velocity
    velocities_out = velocities_com + velocity_com_out
    outgoing_beta = _maximum_beta(velocities_out, maximum_beta, "outgoing")
    outgoing = masses_out[:, None] * velocities_out
    com_momenta_out = masses_out[:, None] * velocities_com
    momentum_after = _vector_sum(outgoing)
    kinetic_before = _kinetic_energy(incoming, masses_in)
    kinetic_after = _kinetic_energy(outgoing, masses_out)
    internal_before = math.fsum(catalog.energy_ev(state) for state in reactants) * ELEMENTARY_CHARGE_C
    internal_after = math.fsum(catalog.energy_ev(state) for state in products) * ELEMENTARY_CHARGE_C
    total_before = math.fsum((kinetic_before, internal_before))
    total_after = math.fsum((kinetic_after, internal_after))
    energy_residual = math.fsum((kinetic_after, -kinetic_before, defect_j))
    if not all(math.isfinite(value) for value in (
        available, kinetic_before, kinetic_after, internal_before, internal_after,
        total_before, total_after, energy_residual,
    )):
        raise ValueError("nonfinite reaction energy ledger")
    return ReactionOutcome(
        reactant_state_ids=reactants,
        product_state_ids=products,
        outgoing_momenta_kg_m_s=_readonly(outgoing),
        available_energy_j=available,
        incoming_com_kinetic_energy_j=com_energy_before,
        outgoing_com_kinetic_energy_j=_kinetic_energy(com_momenta_out, masses_out),
        kinetic_energy_before_j=kinetic_before,
        kinetic_energy_after_j=kinetic_after,
        internal_energy_before_j=internal_before,
        internal_energy_after_j=internal_after,
        total_energy_before_j=total_before,
        total_energy_after_j=total_after,
        energy_defect_ev=defect_ev,
        energy_residual_j=energy_residual,
        momentum_before_kg_m_s=_readonly(momentum_before),
        momentum_after_kg_m_s=_readonly(momentum_after),
        momentum_residual_kg_m_s=_readonly(momentum_after - momentum_before),
        mass_residual_kg=mass_residual,
        center_of_mass_velocity_m_s=_readonly(velocity_com_out),
        jacobi_energy_fractions=_readonly(fractions),
        jacobi_directions=_readonly(directions),
        maximum_incoming_beta=incoming_beta,
        maximum_outgoing_beta=outgoing_beta,
        threshold_roundoff_correction_j=threshold_roundoff,
        metadata=MappingProxyType({
            "channel_name": channel.name,
            "kinematics_model": SUPPORTED_KINEMATICS_MODEL,
            "channel_metadata": channel.metadata,
            "catalog_fingerprint": catalog.fingerprint,
            "kinematic_distribution": "externally supplied ordered Jacobi fractions and unit directions",
            "cross_section_and_differential_kernel_validity": "caller responsibility; no extrapolation or physical sampling performed here",
            "max_nonrel_beta": maximum_beta,
            "bulk_mass_roundoff_correction_j": bulk_mass_roundoff_correction,
        }),
    )
