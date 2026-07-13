"""Energy-dependent collision cross-section tables and PIC-MC helpers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from .constants import QE


@dataclass(frozen=True)
class CrossSectionTable:
    name: str
    incident: str
    target: str
    products: dict[str, float]
    process: str
    reference: str
    energy_ev: np.ndarray
    cross_section_m2: np.ndarray
    source_file: str

    def sigma(self, energy_ev: np.ndarray | float) -> np.ndarray:
        energy = np.asarray(energy_ev, dtype=float)
        return np.interp(
            energy,
            self.energy_ev,
            self.cross_section_m2,
            left=0.0,
            right=0.0,
        )


class CrossSectionLibrary:
    def __init__(self, tables: list[CrossSectionTable]):
        self.tables = {table.name: table for table in tables}

    @classmethod
    def from_manifest(cls, manifest_path: str | Path) -> "CrossSectionLibrary":
        import yaml

        manifest_path = Path(manifest_path).resolve()
        data = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
        tables = []
        for item in data.get("cross_sections", []):
            file_name = item.get("file")
            if not file_name:
                continue
            table_path = (manifest_path.parent / file_name).resolve()
            if not table_path.exists():
                continue
            energy, sigma = load_cross_section_csv(table_path)
            tables.append(
                CrossSectionTable(
                    name=item["name"],
                    incident=item["incident"],
                    target=item["target"],
                    products=item.get("products", {}),
                    process=item.get("process", "unknown"),
                    reference=item.get("reference", ""),
                    energy_ev=energy,
                    cross_section_m2=sigma,
                    source_file=str(table_path),
                )
            )
        return cls(tables)

    def reaction_rates(self, incident_energies_ev: dict[str, float], target_densities_m3: dict[str, float]) -> dict[str, float]:
        """Return collision frequencies in s^-1 for available tables.

        This is the null-collision/PIC-MC scalar rate, nu = n_target sigma(E) v.
        Full particle sampling can use ``collision_probability`` directly.
        """
        rates = {}
        for table in self.tables.values():
            energy = incident_energies_ev.get(table.incident)
            density = target_densities_m3.get(table.target, 0.0)
            if energy is None or density <= 0.0:
                continue
            sigma = float(table.sigma(energy))
            speed = incident_speed_from_energy_ev(energy, table.incident)
            rates[f"{table.name}_s"] = density * sigma * speed
        return rates


def load_cross_section_csv(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    table = np.genfromtxt(path, delimiter=",", names=True)
    names = table.dtype.names or ()
    if "energy_ev" not in names or "cross_section_m2" not in names:
        raise ValueError(f"{path} must contain energy_ev and cross_section_m2 columns")
    energy = np.asarray(table["energy_ev"], dtype=float)
    sigma = np.asarray(table["cross_section_m2"], dtype=float)
    valid = np.isfinite(energy) & np.isfinite(sigma)
    order = np.argsort(energy[valid])
    return energy[valid][order], sigma[valid][order]


def collision_probability(target_density_m3: np.ndarray, sigma_m2: np.ndarray, relative_speed_m_s: np.ndarray, dt_s: float) -> np.ndarray:
    nu_dt = np.maximum(target_density_m3, 0.0) * np.maximum(sigma_m2, 0.0) * np.maximum(relative_speed_m_s, 0.0) * dt_s
    return 1.0 - np.exp(-nu_dt)


def sample_collisions(probability: np.ndarray, rng: np.random.Generator | None = None) -> np.ndarray:
    rng = rng or np.random.default_rng()
    return rng.random(probability.shape) < probability


def expected_collision_loss(population: np.ndarray, probability: np.ndarray) -> np.ndarray:
    return population * np.clip(probability, 0.0, 1.0)


def incident_speed_from_energy_ev(energy_ev: float, incident: str) -> float:
    # Electrons dominate velocity for electron-impact channels. Ion/neutral
    # masses should be supplied later through species metadata; this fallback
    # keeps table-based rates finite and explicit for first-use scaffolding.
    if incident == "e":
        mass_kg = 9.1093837139e-31
    else:
        mass_kg = 1.67262192595e-27
    return float(np.sqrt(2.0 * max(energy_ev, 0.0) * QE / mass_kg))
