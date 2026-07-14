"""Discrete-ordinate kinetic neutral transport."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .backend import ArrayBackend
from .constants import KB, PI, QE, SPECIES_MASS_KG
from .cross_sections import CrossSectionLibrary, collision_probability, incident_speed_from_energy_ev
from .geometry import neighbors4


@dataclass
class VelocityGrid:
    vx: np.ndarray
    vy: np.ndarray
    speed: np.ndarray
    weight: np.ndarray

    @classmethod
    def polar(cls, max_speed_m_s: float, n_speed: int, n_angle: int) -> "VelocityGrid":
        speeds = (np.arange(n_speed) + 0.5) * max_speed_m_s / n_speed
        angles = (np.arange(n_angle) + 0.5) * 2.0 * PI / n_angle
        speed, angle = np.meshgrid(speeds, angles, indexing="ij")
        speed = speed.ravel()
        angle = angle.ravel()
        vx = speed * np.cos(angle)
        vy = speed * np.sin(angle)
        # Area element in 2D velocity polar coordinates, normalized for isotropic reuse.
        weight = speed.copy()
        weight /= weight.sum()
        return cls(vx=vx, vy=vy, speed=speed, weight=weight)


@dataclass
class NeutralState:
    f: dict[str, np.ndarray]
    electron_density_m3: np.ndarray

    def density(self, species: str) -> np.ndarray:
        return self.f[species].sum(axis=0)

    def total_neutral_density(self) -> np.ndarray:
        total = None
        for arr in self.f.values():
            dens = arr.sum(axis=0)
            total = dens if total is None else total + dens
        return total if total is not None else np.array(0.0)


class KineticNeutralSolver:
    def __init__(self, raster, species: list[str], velocity_grid: VelocityGrid, config: dict):
        self.raster = raster
        self.species = species
        self.velocity_grid = velocity_grid
        self.backend = ArrayBackend.from_preference(config.get("backend", "numpy"))
        self.xp = self.backend.xp
        self.bgk_frequency_s = float(config.get("bgk_frequency_s", 0.0))
        self.wall_sticking = float(config.get("wall_sticking", 0.5))
        self.initial_electron_density_m3 = float(config.get("initial_electron_density_m3", 1.0e10))
        self.state_dtype = np.dtype(config.get("state_dtype", "float32"))
        self.max_working_set_gb = float(config.get("max_working_set_gb", 8.0))
        self.max_gpu_memory_fraction = float(config.get("max_gpu_memory_fraction", 0.75))
        self.material_mask = (
            raster.mask_by_kind("cathode")
            | raster.mask_by_kind("anode")
            | raster.mask_by_kind("material")
        )
        self.vacuum_mask = raster.mask_by_kind("vacuum") | ~self.material_mask
        self._material_mask_xp = None
        self._vacuum_mask_xp = None
        self._wall_adjacent = None
        self._wall_adjacent_xp = None
        self._velocity_weight_xp = None
        self._velocity_vx_xp = None
        self._velocity_vy_xp = None

    def initial_state(self) -> NeutralState:
        self._check_memory_budget()
        shape = self.raster.shape
        xp = self.xp
        f = {
            sp: xp.zeros((self.velocity_grid.vx.size, *shape), dtype=self.state_dtype)
            for sp in self.species
        }
        ne = xp.full(shape, self.initial_electron_density_m3, dtype=self.state_dtype)
        ne[~self._vacuum_mask()] = 0.0
        return NeutralState(f=f, electron_density_m3=ne)

    def estimated_working_set_gb(self) -> float:
        ny, nx = self.raster.shape
        bytes_per_value = self.state_dtype.itemsize
        cells = nx * ny
        f_bytes = len(self.species) * self.velocity_grid.vx.size * cells * bytes_per_value
        electron_bytes = cells * bytes_per_value
        # The transport step holds old and new distributions at once, plus
        # moments/reaction scratch arrays. Keep the estimate intentionally
        # conservative so a run fails before the machine starts paging.
        return (2.25 * f_bytes + 4.0 * electron_bytes) / (1024.0**3)

    def asnumpy(self, value):
        return self.backend.asnumpy(value)

    def _check_memory_budget(self) -> None:
        estimate_gb = self.estimated_working_set_gb()
        budget_gb = self.max_working_set_gb
        if self.backend.is_gpu:
            free_bytes, total_bytes = self.xp.cuda.runtime.memGetInfo()
            gpu_budget_gb = (free_bytes / (1024.0**3)) * min(max(self.max_gpu_memory_fraction, 0.05), 1.0)
            budget_gb = min(budget_gb, gpu_budget_gb)
        if estimate_gb <= budget_gb:
            return
        ny, nx = self.raster.shape
        backend_note = ""
        if self.backend.is_gpu:
            free_bytes, total_bytes = self.xp.cuda.runtime.memGetInfo()
            backend_note = (
                f" GPU memory free/total is {free_bytes / (1024.0**3):.2f}/"
                f"{total_bytes / (1024.0**3):.2f} GiB with "
                f"max_gpu_memory_fraction={self.max_gpu_memory_fraction:.2f}."
            )
        raise MemoryError(
            "neutral velocity-space state is too large: "
            f"estimated working set {estimate_gb:.2f} GiB exceeds "
            f"the active budget {budget_gb:.2f} GiB "
            f"(neutrals.max_working_set_gb={self.max_working_set_gb:.2f} GiB)."
            f"{backend_note} "
            f"for mesh {nx} x {ny}, {len(self.species)} species, "
            f"{self.velocity_grid.vx.size} velocity ordinates, dtype {self.state_dtype}. "
            "Coarsen grid.nx/grid.ny or grid.target_cell_size_um, reduce velocity_grid bins, "
            "reduce the species list, or raise neutrals.max_working_set_gb if the machine can handle it."
        )

    def step(
        self,
        state: NeutralState,
        dt_s: float,
        surface_source_m2_s: np.ndarray,
        surface_temperature_k: np.ndarray,
        source_species: str,
        surface_mask: np.ndarray,
        reaction_rates: dict[str, float] | None = None,
        cross_sections: CrossSectionLibrary | None = None,
        incident_energies_ev: dict[str, float] | None = None,
        electron_energy_ev: np.ndarray | float | None = None,
    ) -> NeutralState:
        next_f = {sp: self._advect_species(arr, dt_s) for sp, arr in state.f.items()}
        self._apply_walls(next_f)
        if source_species in next_f:
            self._emit_surface(
                next_f[source_species],
                surface_source_m2_s,
                surface_temperature_k,
                source_species,
                surface_mask,
                dt_s,
            )
        if self.bgk_frequency_s > 0.0:
            self._bgk_relax(next_f, dt_s)
        ne = state.electron_density_m3.copy()
        if reaction_rates:
            ne = self._apply_reaction_sources(next_f, ne, reaction_rates, dt_s)
        if cross_sections:
            ne = self._apply_cross_section_collisions(
                next_f,
                ne,
                cross_sections,
                incident_energies_ev or {},
                dt_s,
                electron_energy_ev=electron_energy_ev,
            )
        return NeutralState(f=next_f, electron_density_m3=ne)

    def moments(self, state: NeutralState) -> dict[str, dict[str, np.ndarray]]:
        out = {}
        vx = self._velocity_vx()[:, None, None]
        vy = self._velocity_vy()[:, None, None]
        for sp, arr in state.f.items():
            dens = arr.sum(axis=0)
            fx = (arr * vx).sum(axis=0)
            fy = (arr * vy).sum(axis=0)
            out[sp] = {"density_m3": dens, "flux_x_m2_s": fx, "flux_y_m2_s": fy}
        return out

    def _advect_species(self, arr: np.ndarray, dt_s: float) -> np.ndarray:
        xp = self.xp
        out = xp.zeros_like(arr)
        for k, (vx, vy) in enumerate(zip(self.velocity_grid.vx, self.velocity_grid.vy)):
            sx = vx * dt_s / self.raster.dx
            sy = vy * dt_s / self.raster.dy
            out[k] = self._shift_fractional_no_wrap(arr[k], sx=sx, sy=sy)
        out[:, ~self._vacuum_mask()] = 0.0
        return out

    def _apply_walls(self, f: dict[str, np.ndarray]) -> None:
        if self.wall_sticking <= 0.0:
            return
        wall_adjacent = self._wall_adjacent_mask()
        for arr in f.values():
            arr[:, wall_adjacent] *= max(0.0, 1.0 - self.wall_sticking)

    def _emit_surface(
        self,
        arr: np.ndarray,
        source_m2_s: np.ndarray,
        surface_temperature_k: np.ndarray,
        species: str,
        surface_mask: np.ndarray,
        dt_s: float,
    ) -> None:
        mass = SPECIES_MASS_KG[species]
        targets = emission_targets(surface_mask, self.vacuum_mask)
        if not targets.any():
            return
        local_source = np.zeros_like(source_m2_s)
        for shifted_mask, shifted_source in zip(neighbors4(surface_mask), neighbor_values(source_m2_s)):
            local_source += np.where(targets & shifted_mask, shifted_source, 0.0)
        local_temp = np.where(targets, neighbor_average(surface_temperature_k, surface_mask), 300.0)
        emission = np.zeros((self.velocity_grid.vx.size, *self.raster.shape), dtype=self.state_dtype)
        for j, i in zip(*np.nonzero(targets)):
            temp = max(float(local_temp[j, i]), 1.0)
            weights = thermal_weights(self.velocity_grid, mass, temp)
            # Convert molecules / m2 / s into molecules / m3 by depositing across one cell width.
            density_increment = local_source[j, i] * dt_s / max(min(self.raster.dx, self.raster.dy), 1e-30)
            emission[:, j, i] += density_increment * weights
        arr += self.xp.asarray(emission)

    def _bgk_relax(self, f: dict[str, np.ndarray], dt_s: float) -> None:
        alpha = min(max(self.bgk_frequency_s * dt_s, 0.0), 1.0)
        if alpha == 0.0:
            return
        weights = self._velocity_weight()[:, None, None]
        for arr in f.values():
            density = arr.sum(axis=0, keepdims=True)
            arr *= 1.0 - alpha
            arr += alpha * density * weights
            arr[:, ~self._vacuum_mask()] = 0.0

    def _apply_reaction_sources(
        self,
        f: dict[str, np.ndarray],
        electron_density_m3: np.ndarray,
        rates: dict[str, float],
        dt_s: float,
    ) -> np.ndarray:
        # Minimal physically transparent plumbing: H2O electron-impact dissociation
        # and ionization terms can be supplied by the CR/BOLSIG mechanism.
        xp = self.xp
        ne = electron_density_m3.copy()
        if "H2O_dissociation_s" in rates and "H2O" in f:
            h2o_density = f["H2O"].sum(axis=0)
            loss = xp.minimum(h2o_density, h2o_density * rates["H2O_dissociation_s"] * dt_s)
            remove_fraction = self._safe_divide(loss, h2o_density)
            f["H2O"] *= 1.0 - remove_fraction[None, :, :]
            for sp, stoich in (("H", 2.0), ("O", 1.0)):
                if sp in f:
                    f[sp] += stoich * loss[None, :, :] * self._velocity_weight()[:, None, None]
        if "impact_ionization_s" in rates:
            source = rates["impact_ionization_s"] * self.total_density_from_f(f) * dt_s
            ne += source
        ne[~self._vacuum_mask()] = 0.0
        return ne

    def _apply_cross_section_collisions(
        self,
        f: dict[str, np.ndarray],
        electron_density_m3: np.ndarray,
        library: CrossSectionLibrary,
        incident_energies_ev: dict[str, float],
        dt_s: float,
        electron_energy_ev: np.ndarray | float | None = None,
    ) -> np.ndarray:
        xp = self.xp
        ne = electron_density_m3.copy()
        velocity_weights = self._velocity_weight()[:, None, None]
        for table in library.tables.values():
            if not table.products or table.target not in f:
                continue
            if table.incident == "e":
                incident_density = ne
                consumed_electrons = 1.0
                energy = self._electron_energy_array(electron_energy_ev, incident_energies_ev)
            elif table.incident in f:
                incident_density = f[table.incident].sum(axis=0)
                consumed_electrons = 0.0
                energy = self._mean_projectile_energy_ev(f[table.incident], table.incident)
            else:
                continue
            if energy is None:
                continue

            sigma = self._sigma_for_energy(table, energy)
            speed = self._speed_from_energy_ev(energy, table.incident)
            if self.backend.is_gpu:
                probability = 1.0 - xp.exp(
                    -xp.maximum(incident_density, 0.0)
                    * xp.maximum(sigma, 0.0)
                    * xp.maximum(speed, 0.0)
                    * dt_s
                )
            else:
                probability = collision_probability(incident_density, sigma, speed, dt_s)
            target_density = f[table.target].sum(axis=0)
            loss = xp.minimum(target_density, target_density * probability)
            if not bool(self.backend.asnumpy(xp.any(loss > 0.0))):
                continue

            remove_fraction = self._safe_divide(loss, target_density)
            f[table.target] *= 1.0 - remove_fraction[None, :, :]
            if table.incident != "e" and table.incident in f and table.incident != table.target:
                incident_target = f[table.incident].sum(axis=0)
                incident_fraction = self._safe_divide(loss, incident_target)
                f[table.incident] *= 1.0 - xp.clip(incident_fraction[None, :, :], 0.0, 1.0)

            electron_delta = (float(table.products.get("e", 0.0)) - consumed_electrons) * loss
            ne += electron_delta
            for species, coeff in table.products.items():
                if species == "e":
                    continue
                if species in f:
                    f[species] += float(coeff) * loss[None, :, :] * velocity_weights
        ne[~self._vacuum_mask()] = 0.0
        return xp.maximum(ne, 0.0)

    def total_density_from_f(self, f: dict[str, np.ndarray]) -> np.ndarray:
        total = None
        for arr in f.values():
            dens = arr.sum(axis=0)
            total = dens if total is None else total + dens
        return total if total is not None else 0.0

    def _shift_fractional_no_wrap(self, arr, sx: float, sy: float):
        if self.backend.is_gpu:
            return shift_fractional_no_wrap_xp(arr, sx=sx, sy=sy, xp=self.xp)
        return shift_fractional_no_wrap(arr, sx=sx, sy=sy)

    def _safe_divide(self, numerator, denominator):
        xp = self.xp
        return xp.where(denominator > 0.0, numerator / xp.maximum(denominator, 1.0e-30), xp.zeros_like(numerator))

    def _electron_energy_array(
        self,
        electron_energy_ev: np.ndarray | float | None,
        fallback_energies_ev: dict[str, float],
    ):
        if electron_energy_ev is not None:
            return self.xp.asarray(electron_energy_ev, dtype=self.state_dtype)
        if "e" not in fallback_energies_ev:
            return None
        return float(fallback_energies_ev["e"])

    def _mean_projectile_energy_ev(self, arr, species: str):
        xp = self.xp
        mass = SPECIES_MASS_KG.get(species)
        if mass is None:
            return None
        density = arr.sum(axis=0)
        speed2 = self.xp.asarray(self.velocity_grid.speed**2, dtype=self.state_dtype)[:, None, None]
        kinetic_j_m3 = 0.5 * mass * (arr * speed2).sum(axis=0)
        return self._safe_divide(kinetic_j_m3, density) / QE

    def _sigma_for_energy(self, table, energy_ev):
        if np.isscalar(energy_ev):
            return float(table.sigma(float(energy_ev)))
        if self.backend.is_gpu:
            energy_axis = self.xp.asarray(table.energy_ev, dtype=self.state_dtype)
            sigma_axis = self.xp.asarray(table.cross_section_m2, dtype=self.state_dtype)
            return self.xp.interp(energy_ev, energy_axis, sigma_axis, left=0.0, right=0.0)
        return table.sigma(energy_ev)

    def _speed_from_energy_ev(self, energy_ev, incident: str):
        if np.isscalar(energy_ev):
            return incident_speed_from_energy_ev(float(energy_ev), incident)
        if incident == "e":
            mass_kg = 9.1093837139e-31
        else:
            mass_kg = SPECIES_MASS_KG.get(incident, 1.67262192595e-27)
        return self.xp.sqrt(2.0 * self.xp.maximum(energy_ev, 0.0) * QE / mass_kg)

    def _vacuum_mask(self):
        if not self.backend.is_gpu:
            return self.vacuum_mask
        if self._vacuum_mask_xp is None:
            self._vacuum_mask_xp = self.xp.asarray(self.vacuum_mask)
        return self._vacuum_mask_xp

    def _wall_adjacent_mask(self):
        if self._wall_adjacent is None:
            wall_adjacent = np.zeros_like(self.material_mask, dtype=bool)
            for shifted in neighbors4(self.material_mask):
                wall_adjacent |= shifted
            wall_adjacent &= self.vacuum_mask
            self._wall_adjacent = wall_adjacent
        if not self.backend.is_gpu:
            return self._wall_adjacent
        if self._wall_adjacent_xp is None:
            self._wall_adjacent_xp = self.xp.asarray(self._wall_adjacent)
        return self._wall_adjacent_xp

    def _velocity_weight(self):
        if not self.backend.is_gpu:
            return self.velocity_grid.weight.astype(self.state_dtype, copy=False)
        if self._velocity_weight_xp is None:
            self._velocity_weight_xp = self.xp.asarray(self.velocity_grid.weight, dtype=self.state_dtype)
        return self._velocity_weight_xp

    def _velocity_vx(self):
        if not self.backend.is_gpu:
            return self.velocity_grid.vx
        if self._velocity_vx_xp is None:
            self._velocity_vx_xp = self.xp.asarray(self.velocity_grid.vx, dtype=self.state_dtype)
        return self._velocity_vx_xp

    def _velocity_vy(self):
        if not self.backend.is_gpu:
            return self.velocity_grid.vy
        if self._velocity_vy_xp is None:
            self._velocity_vy_xp = self.xp.asarray(self.velocity_grid.vy, dtype=self.state_dtype)
        return self._velocity_vy_xp


def thermal_weights(grid: VelocityGrid, mass_kg: float, temperature_k: float) -> np.ndarray:
    weights = grid.weight * np.exp(-mass_kg * grid.speed**2 / (2.0 * KB * temperature_k))
    total = weights.sum()
    if total <= 0.0:
        return np.full_like(weights, 1.0 / weights.size)
    return weights / total


def half_range_flux_moment(temperature_k: float, mass_kg: float, density_m3: float = 1.0) -> float:
    return density_m3 * np.sqrt(KB * temperature_k / (2.0 * PI * mass_kg))


def emission_targets(surface_mask: np.ndarray, vacuum_mask: np.ndarray) -> np.ndarray:
    targets = np.zeros_like(surface_mask, dtype=bool)
    for shifted in neighbors4(surface_mask):
        targets |= shifted
    return targets & vacuum_mask


def neighbor_average(values: np.ndarray, mask: np.ndarray) -> np.ndarray:
    total = np.zeros_like(values, dtype=float)
    count = np.zeros_like(values, dtype=float)
    for shifted, shifted_values in zip(neighbors4(mask), neighbor_values(values)):
        total += np.where(shifted, shifted_values, 0.0)
        count += shifted.astype(float)
    return np.divide(total, count, out=np.zeros_like(total), where=count > 0)


def neighbor_values(values: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    up = np.zeros_like(values)
    up[:-1, :] = values[1:, :]
    down = np.zeros_like(values)
    down[1:, :] = values[:-1, :]
    left = np.zeros_like(values)
    left[:, :-1] = values[:, 1:]
    right = np.zeros_like(values)
    right[:, 1:] = values[:, :-1]
    return up, down, left, right


def shift_fractional_no_wrap(arr: np.ndarray, sx: float, sy: float) -> np.ndarray:
    ny, nx = arr.shape
    yy, xx = np.indices(arr.shape, dtype=float)
    src_x = xx - sx
    src_y = yy - sy
    x0 = np.floor(src_x).astype(int)
    y0 = np.floor(src_y).astype(int)
    x1 = x0 + 1
    y1 = y0 + 1
    wx = src_x - x0
    wy = src_y - y0

    out = np.zeros_like(arr, dtype=float)
    for ix, iy, weight in (
        (x0, y0, (1.0 - wx) * (1.0 - wy)),
        (x1, y0, wx * (1.0 - wy)),
        (x0, y1, (1.0 - wx) * wy),
        (x1, y1, wx * wy),
    ):
        valid = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
        out[valid] += arr[iy[valid], ix[valid]] * weight[valid]
    return out


def shift_fractional_no_wrap_xp(arr, sx: float, sy: float, xp):
    ny, nx = arr.shape
    yy, xx = xp.indices(arr.shape, dtype=float)
    src_x = xx - sx
    src_y = yy - sy
    x0 = xp.floor(src_x).astype(xp.int64)
    y0 = xp.floor(src_y).astype(xp.int64)
    x1 = x0 + 1
    y1 = y0 + 1
    wx = src_x - x0
    wy = src_y - y0

    out = xp.zeros_like(arr)
    for ix, iy, weight in (
        (x0, y0, (1.0 - wx) * (1.0 - wy)),
        (x1, y0, wx * (1.0 - wy)),
        (x0, y1, (1.0 - wx) * wy),
        (x1, y1, wx * wy),
    ):
        valid = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
        out[valid] += arr[iy[valid], ix[valid]] * weight[valid]
    return out


def shift_no_wrap(arr: np.ndarray, sx: int, sy: int) -> np.ndarray:
    out = np.zeros_like(arr)
    ny, nx = arr.shape
    src_x0 = max(0, -sx)
    src_x1 = min(nx, nx - sx)
    dst_x0 = max(0, sx)
    dst_x1 = min(nx, nx + sx)
    src_y0 = max(0, -sy)
    src_y1 = min(ny, ny - sy)
    dst_y0 = max(0, sy)
    dst_y1 = min(ny, ny + sy)
    if src_x0 < src_x1 and src_y0 < src_y1:
        out[dst_y0:dst_y1, dst_x0:dst_x1] = arr[src_y0:src_y1, src_x0:src_x1]
    return out
