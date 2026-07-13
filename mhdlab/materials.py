"""Material property presets."""

from __future__ import annotations

from dataclasses import dataclass, replace


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

    @property
    def thermal_diffusivity_m2_s(self) -> float:
        return self.thermal_conductivity_w_m_k / (
            self.density_kg_m3 * self.heat_capacity_j_kg_k
        )


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
