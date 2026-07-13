"""Temkin-style thermal desorption surface source."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .constants import KB, QE


@dataclass
class TemkinDesorption:
    species: str = "H2O"
    attempt_frequency_s: float = 1.0e13
    activation_energy_ev: float = 0.95
    coverage_interaction_ev: float = 0.15
    initial_coverage_m2: float = 3.0e18
    min_temperature_k: float = 250.0

    def __post_init__(self) -> None:
        self.attempt_frequency_s = float(self.attempt_frequency_s)
        self.activation_energy_ev = float(self.activation_energy_ev)
        self.coverage_interaction_ev = float(self.coverage_interaction_ev)
        self.initial_coverage_m2 = float(self.initial_coverage_m2)
        self.min_temperature_k = float(self.min_temperature_k)

    def initial_coverage(self, shape: tuple[int, int]) -> np.ndarray:
        return np.full(shape, self.initial_coverage_m2, dtype=float)

    def rate(self, temperature_k: np.ndarray, coverage_m2: np.ndarray) -> np.ndarray:
        temp = np.maximum(temperature_k, self.min_temperature_k)
        theta = np.clip(coverage_m2 / max(self.initial_coverage_m2, 1.0), 0.0, 1.0)
        barrier_j = (self.activation_energy_ev - self.coverage_interaction_ev * theta) * QE
        return self.attempt_frequency_s * coverage_m2 * np.exp(-barrier_j / (KB * temp))

    def step(
        self,
        temperature_k: np.ndarray,
        coverage_m2: np.ndarray,
        surface_mask: np.ndarray,
        dt_s: float,
    ) -> tuple[np.ndarray, np.ndarray]:
        rate = np.zeros_like(temperature_k, dtype=float)
        local_rate = self.rate(temperature_k, coverage_m2)
        rate[surface_mask] = local_rate[surface_mask]
        emitted = rate * dt_s
        next_coverage = np.maximum(coverage_m2 - emitted, 0.0)
        return rate, next_coverage
