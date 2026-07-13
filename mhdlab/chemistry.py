"""YAML-driven reaction and collisional-radiative scaffolding."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np


@dataclass
class Reaction:
    name: str
    reactants: dict[str, float]
    products: dict[str, float]
    rate_s: float | None = None
    bolsig_rate: str | None = None
    reference: str | None = None
    provisional: bool = False


@dataclass
class SpectralLine:
    name: str
    species: str
    wavelength_nm: float
    upper: str
    lower: str
    a_s: float
    energy_ev: float


@dataclass
class Mechanism:
    species: list[str]
    reactions: list[Reaction] = field(default_factory=list)
    spectral_lines: list[SpectralLine] = field(default_factory=list)

    @classmethod
    def from_yaml(cls, path: str | Path) -> "Mechanism":
        import yaml

        data = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
        reactions = [
            Reaction(
                name=item["name"],
                reactants=item.get("reactants", {}),
                products=item.get("products", {}),
                rate_s=item.get("rate_s"),
                bolsig_rate=item.get("bolsig_rate"),
                reference=item.get("reference"),
                provisional=bool(item.get("provisional", False)),
            )
            for item in data.get("reactions", [])
        ]
        lines = [
            SpectralLine(
                name=item["name"],
                species=item["species"],
                wavelength_nm=float(item["wavelength_nm"]),
                upper=item.get("upper", ""),
                lower=item.get("lower", ""),
                a_s=float(item.get("a_s", 0.0)),
                energy_ev=float(item.get("energy_ev", 0.0)),
            )
            for item in data.get("spectral_lines", [])
        ]
        return cls(species=list(data.get("species", [])), reactions=reactions, spectral_lines=lines)


def local_reaction_rates(
    mechanism: Mechanism,
    bolsig_rates: dict[str, np.ndarray | float],
    allow_provisional: bool,
) -> dict[str, float]:
    rates: dict[str, float] = {}
    missing = []
    for reaction in mechanism.reactions:
        if reaction.provisional and not allow_provisional:
            missing.append(reaction.name)
            continue
        value: float | None = None
        if reaction.rate_s is not None:
            value = float(reaction.rate_s)
        elif reaction.bolsig_rate and reaction.bolsig_rate in bolsig_rates:
            raw = bolsig_rates[reaction.bolsig_rate]
            value = float(np.nanmean(raw))
        if value is None:
            missing.append(reaction.name)
        else:
            rates[reaction.name] = value
            if reaction.bolsig_rate:
                rates[reaction.bolsig_rate] = value
    if missing:
        rates["_missing_count"] = float(len(missing))
    return rates


def line_emissivity(
    lines: list[SpectralLine],
    electron_density_m3: np.ndarray,
    species_densities: dict[str, np.ndarray],
    excitation_scale: float = 1.0e-21,
) -> dict[str, np.ndarray]:
    """Produce simple optically thin emissivity maps for line-synthesis plumbing."""
    emiss = {}
    for line in lines:
        donor = species_densities.get(line.species)
        if donor is None:
            donor = np.zeros_like(electron_density_m3)
        emiss[line.name] = excitation_scale * line.a_s * electron_density_m3 * donor
    return emiss
