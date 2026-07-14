"""Material property presets."""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np


@dataclass(frozen=True)
class Material:
    name: str
    density_kg_m3: float
    heat_capacity_j_kg_k: float
    thermal_conductivity_w_m_k: float
    electrical_resistivity_ohm_m: float
    thermal_expansion_1_k: float
    bulk_modulus_pa: float
    initial_temperature_k: float = 300.0
    melting_temperature_k: float = 1670.0
    boiling_temperature_k: float = 3135.0
    liquid_heat_capacity_j_kg_k: float = 820.0
    vapor_heat_capacity_j_kg_k: float = 1000.0
    latent_heat_fusion_j_kg: float = 2.6e5
    latent_heat_vaporization_j_kg: float = 6.3e6
    max_temperature_k: float = 12000.0

    @property
    def thermal_diffusivity_m2_s(self) -> float:
        return self.thermal_conductivity_w_m_k / (
            self.density_kg_m3 * self.heat_capacity_j_kg_k
        )

    def specific_enthalpy_from_temperature(self, temperature_k) -> np.ndarray:
        """Piecewise sensible/latent enthalpy per mass in J/kg.

        This is a compact Knoepfel-style conductor heating surrogate: Joule
        energy advances enthalpy, while melting and vaporization are represented
        as constant-temperature latent heat intervals.
        """
        temp = np.asarray(temperature_k, dtype=float)
        tm = self.melting_temperature_k
        tb = self.boiling_temperature_k
        cp_s = self.heat_capacity_j_kg_k
        cp_l = self.liquid_heat_capacity_j_kg_k
        cp_v = self.vapor_heat_capacity_j_kg_k
        h_melt_start = cp_s * tm
        h_melt_end = h_melt_start + self.latent_heat_fusion_j_kg
        h_boil_start = h_melt_end + cp_l * max(tb - tm, 0.0)
        h_boil_end = h_boil_start + self.latent_heat_vaporization_j_kg

        h = cp_s * np.minimum(temp, tm)
        h = np.where(temp > tm, h_melt_end + cp_l * np.minimum(temp - tm, tb - tm), h)
        h = np.where(temp > tb, h_boil_end + cp_v * (temp - tb), h)
        return h

    def temperature_from_specific_enthalpy(self, enthalpy_j_kg) -> np.ndarray:
        h = np.asarray(enthalpy_j_kg, dtype=float)
        tm = self.melting_temperature_k
        tb = self.boiling_temperature_k
        cp_s = max(self.heat_capacity_j_kg_k, 1.0e-30)
        cp_l = max(self.liquid_heat_capacity_j_kg_k, 1.0e-30)
        cp_v = max(self.vapor_heat_capacity_j_kg_k, 1.0e-30)
        h_melt_start = cp_s * tm
        h_melt_end = h_melt_start + self.latent_heat_fusion_j_kg
        h_boil_start = h_melt_end + cp_l * max(tb - tm, 0.0)
        h_boil_end = h_boil_start + self.latent_heat_vaporization_j_kg

        temp = h / cp_s
        temp = np.where((h > h_melt_start) & (h <= h_melt_end), tm, temp)
        temp = np.where((h > h_melt_end) & (h <= h_boil_start), tm + (h - h_melt_end) / cp_l, temp)
        temp = np.where((h > h_boil_start) & (h <= h_boil_end), tb, temp)
        temp = np.where(h > h_boil_end, tb + (h - h_boil_end) / cp_v, temp)
        return np.minimum(temp, self.max_temperature_k)


SS304 = Material(
    name="SS304",
    density_kg_m3=8000.0,
    heat_capacity_j_kg_k=500.0,
    thermal_conductivity_w_m_k=16.2,
    electrical_resistivity_ohm_m=7.2e-7,
    thermal_expansion_1_k=17.3e-6,
    bulk_modulus_pa=160e9,
)


PRESETS = {"SS304": SS304, "stainless_304": SS304, "Stainless 304": SS304}


def material_from_config(data: dict | None) -> Material:
    data = data or {}
    preset = data.get("preset", "SS304")
    material = PRESETS.get(preset, SS304)
    overrides = {k: v for k, v in data.items() if hasattr(material, k)}
    return replace(material, **overrides)
