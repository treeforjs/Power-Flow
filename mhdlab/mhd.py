"""Reduced 2D MHD and electrostatic field solver."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .backend import ArrayBackend
from .constants import MU0
from .geometry import Raster, boundary_mask
from .materials import Material


@dataclass
class MHDState:
    temperature_k: np.ndarray
    density_kg_m3: np.ndarray
    pressure_pa: np.ndarray
    magnetic_pressure_pa: np.ndarray
    az_wb_m: np.ndarray
    bx_t: np.ndarray
    by_t: np.ndarray
    jz_a_m2: np.ndarray
    phi_v: np.ndarray
    ex_v_m: np.ndarray
    ey_v_m: np.ndarray
    surface_displacement_m: np.ndarray


class ReducedMHDSolver:
    def __init__(self, raster: Raster, material: Material, config: dict):
        self.raster = raster
        self.material = material
        self.depth_m = float(config.get("out_of_plane_depth_m", 0.01))
        self.backend = ArrayBackend.from_preference(config.get("backend", "numpy"))
        self.poisson_iterations = int(config.get("poisson_iterations", 250))
        self.electrostatic_iterations = int(config.get("electrostatic_iterations", 350))
        self.induction_iterations = int(config.get("induction_iterations", 4))
        self.induction_relaxation = float(config.get("induction_relaxation", 0.7))
        self.temperature_floor_k = float(config.get("temperature_floor_k", 1.0))
        self.conductivity_config = dict(config.get("conductivity", {"model": "temperature"}))
        self.material_mask = (
            raster.mask_by_kind("cathode")
            | raster.mask_by_kind("anode")
            | raster.mask_by_kind("material")
        )
        self.vacuum_mask = raster.mask_by_kind("vacuum") | ~self.material_mask
        self.cathode_mask = raster.mask_by_kind("cathode")
        self.anode_mask = raster.mask_by_kind("anode")
        self.surface_mask = boundary_mask(self.material_mask, self.vacuum_mask)

    def initial_state(self) -> MHDState:
        shape = self.raster.shape
        temp = np.full(shape, self.material.initial_temperature_k, dtype=float)
        density = np.where(self.material_mask, self.material.density_kg_m3, 0.0)
        zero = np.zeros(shape, dtype=float)
        return MHDState(
            temperature_k=temp,
            density_kg_m3=density,
            pressure_pa=zero.copy(),
            magnetic_pressure_pa=zero.copy(),
            az_wb_m=zero.copy(),
            bx_t=zero.copy(),
            by_t=zero.copy(),
            jz_a_m2=zero.copy(),
            phi_v=zero.copy(),
            ex_v_m=zero.copy(),
            ey_v_m=zero.copy(),
            surface_displacement_m=zero.copy(),
        )

    def step(self, state: MHDState, current_a: float, voltage_v: float, dt_s: float) -> MHDState:
        if self.backend.is_gpu:
            return self._step_gpu(state, current_a=current_a, voltage_v=voltage_v, dt_s=dt_s)
        jz, az = self._solve_inductive_current(state, current_a, dt_s)
        bx = np.gradient(az, self.raster.dy, axis=0)
        by = -np.gradient(az, self.raster.dx, axis=1)
        b2 = bx * bx + by * by
        mag_pressure = b2 / (2.0 * MU0)

        conductivity = self._electrical_conductivity(state.temperature_k, state.density_kg_m3, jz)
        heat = jz * jz / np.maximum(conductivity, 1.0e-30)
        lap_t = laplacian(state.temperature_k, self.raster.dx, self.raster.dy)
        temp = state.temperature_k + dt_s * (
            heat / (self.material.density_kg_m3 * self.material.heat_capacity_j_kg_k)
            + self.material.thermal_diffusivity_m2_s * lap_t
        )
        temp = np.where(self.material_mask, np.maximum(temp, self.temperature_floor_k), state.temperature_k)

        delta_t = temp - self.material.initial_temperature_k
        volume_factor = np.maximum(1.0 + 3.0 * self.material.thermal_expansion_1_k * delta_t, 0.05)
        density = np.where(self.material_mask, self.material.density_kg_m3 / volume_factor, 0.0)
        thermal_pressure = np.where(
            self.material_mask,
            self.material.bulk_modulus_pa * self.material.thermal_expansion_1_k * np.maximum(delta_t, 0.0),
            0.0,
        )
        pressure = thermal_pressure + mag_pressure

        phi = solve_electrostatic(
            shape=self.raster.shape,
            cathode=self.cathode_mask,
            anode=self.anode_mask,
            voltage_v=voltage_v,
            dx=self.raster.dx,
            dy=self.raster.dy,
            iterations=self.electrostatic_iterations,
        )
        ey, ex = np.gradient(-phi, self.raster.dy, self.raster.dx)

        surface_speed = np.zeros_like(temp)
        surface_speed[self.surface_mask] = np.sqrt(
            np.maximum(pressure[self.surface_mask], 0.0)
            / np.maximum(density[self.surface_mask], 1.0)
        )
        surface_displacement = state.surface_displacement_m + surface_speed * dt_s

        return MHDState(
            temperature_k=temp,
            density_kg_m3=density,
            pressure_pa=pressure,
            magnetic_pressure_pa=mag_pressure,
            az_wb_m=az,
            bx_t=bx,
            by_t=by,
            jz_a_m2=jz,
            phi_v=phi,
            ex_v_m=ex,
            ey_v_m=ey,
            surface_displacement_m=surface_displacement,
        )

    def _step_gpu(self, state: MHDState, current_a: float, voltage_v: float, dt_s: float) -> MHDState:
        xp = self.backend.xp
        material_mask = xp.asarray(self.material_mask)
        surface_mask = xp.asarray(self.surface_mask)
        jz_np, az_np = self._solve_inductive_current(state, current_a, dt_s)
        conductivity_np = self._electrical_conductivity(state.temperature_k, state.density_kg_m3, jz_np)
        jz = xp.asarray(jz_np)
        conductivity = xp.asarray(conductivity_np)
        temp0 = xp.asarray(state.temperature_k)
        density0 = xp.asarray(state.density_kg_m3)
        displacement0 = xp.asarray(state.surface_displacement_m)

        az = xp.asarray(az_np)
        bx = xp.gradient(az, self.raster.dy, axis=0)
        by = -xp.gradient(az, self.raster.dx, axis=1)
        b2 = bx * bx + by * by
        mag_pressure = b2 / (2.0 * MU0)

        heat = jz * jz / xp.maximum(conductivity, 1.0e-30)
        lap_t = laplacian_xp(temp0, self.raster.dx, self.raster.dy, xp)
        temp = temp0 + dt_s * (
            heat / (self.material.density_kg_m3 * self.material.heat_capacity_j_kg_k)
            + self.material.thermal_diffusivity_m2_s * lap_t
        )
        temp = xp.where(material_mask, xp.maximum(temp, self.temperature_floor_k), temp0)

        delta_t = temp - self.material.initial_temperature_k
        volume_factor = xp.maximum(1.0 + 3.0 * self.material.thermal_expansion_1_k * delta_t, 0.05)
        density = xp.where(material_mask, self.material.density_kg_m3 / volume_factor, 0.0)
        thermal_pressure = xp.where(
            material_mask,
            self.material.bulk_modulus_pa * self.material.thermal_expansion_1_k * xp.maximum(delta_t, 0.0),
            0.0,
        )
        pressure = thermal_pressure + mag_pressure
        phi = solve_electrostatic_xp(
            shape=self.raster.shape,
            cathode=xp.asarray(self.cathode_mask),
            anode=xp.asarray(self.anode_mask),
            voltage_v=voltage_v,
            dx=self.raster.dx,
            dy=self.raster.dy,
            iterations=self.electrostatic_iterations,
            xp=xp,
        )
        ey, ex = xp.gradient(-phi, self.raster.dy, self.raster.dx)
        surface_speed = xp.zeros_like(temp)
        surface_speed = xp.where(
            surface_mask,
            xp.sqrt(xp.maximum(pressure, 0.0) / xp.maximum(density, 1.0)),
            surface_speed,
        )
        surface_displacement = displacement0 + surface_speed * dt_s

        return MHDState(
            temperature_k=self.backend.asnumpy(temp),
            density_kg_m3=self.backend.asnumpy(density),
            pressure_pa=self.backend.asnumpy(pressure),
            magnetic_pressure_pa=self.backend.asnumpy(mag_pressure),
            az_wb_m=self.backend.asnumpy(az),
            bx_t=self.backend.asnumpy(bx),
            by_t=self.backend.asnumpy(by),
            jz_a_m2=self.backend.asnumpy(jz),
            phi_v=self.backend.asnumpy(phi),
            ex_v_m=self.backend.asnumpy(ex),
            ey_v_m=self.backend.asnumpy(ey),
            surface_displacement_m=self.backend.asnumpy(surface_displacement),
        )

    def _solve_inductive_current(
        self,
        state: MHDState,
        current_a: float,
        dt_s: float,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Solve a reduced resistive-induction update for out-of-plane current.

        In each conductor, Ohm's law is approximated as
        Jz = sigma * (Ez_drive - dAz/dt).  The scalar Ez_drive is solved
        separately for the cathode and anode so each region carries the
        requested total current.  Az is then updated from Ampere's law.  This
        is still reduced/quasi-static, but current crowding is now produced by
        the coupled magnetic diffusion problem rather than by a prescribed
        spatial profile.
        """
        previous_az = np.asarray(state.az_wb_m, dtype=float)
        dt_s = max(float(dt_s), 1.0e-30)
        cell_area = self.raster.dx * self.raster.dy
        az = previous_az.copy()
        jz = np.asarray(state.jz_a_m2, dtype=float).copy()
        relaxation = min(max(self.induction_relaxation, 0.0), 1.0)

        for _ in range(max(self.induction_iterations, 1)):
            conductivity = self._electrical_conductivity(state.temperature_k, state.density_kg_m3, jz)
            jz.fill(0.0)
            self._fill_region_inductive_current(
                jz,
                az,
                previous_az,
                self.cathode_mask,
                target_current_a=float(current_a),
                conductivity=conductivity,
                cell_area=cell_area,
                dt_s=dt_s,
            )
            self._fill_region_inductive_current(
                jz,
                az,
                previous_az,
                self.anode_mask,
                target_current_a=-float(current_a),
                conductivity=conductivity,
                cell_area=cell_area,
                dt_s=dt_s,
            )
            next_az = solve_poisson(
                source=-MU0 * jz,
                dx=self.raster.dx,
                dy=self.raster.dy,
                iterations=self.poisson_iterations,
            )
            az = relaxation * next_az + (1.0 - relaxation) * az
        return jz, az

    def _electrical_conductivity(
        self,
        temperature_k: np.ndarray,
        density_kg_m3: np.ndarray,
        jz_a_m2: np.ndarray,
    ) -> np.ndarray:
        cfg = self.conductivity_config
        model = str(cfg.get("model", "constant")).lower()
        sigma0 = 1.0 / max(self.material.electrical_resistivity_ohm_m, 1.0e-30)
        rho_ratio = np.divide(
            density_kg_m3,
            max(self.material.density_kg_m3, 1.0e-30),
            out=np.zeros_like(density_kg_m3, dtype=float),
            where=density_kg_m3 > 0.0,
        )

        if model in {"constant", "uniform"}:
            sigma = np.full_like(temperature_k, sigma0, dtype=float)
        elif model in {"temperature", "knoepfel_like"}:
            alpha = float(cfg.get("temperature_coefficient_1_k", 9.4e-4))
            t0 = float(cfg.get("reference_temperature_k", self.material.initial_temperature_k))
            density_exponent = float(cfg.get("density_exponent", 0.0))
            denom = 1.0 + alpha * np.maximum(temperature_k - t0, 0.0)
            sigma = sigma0 * np.power(np.maximum(rho_ratio, 0.0), density_exponent) / np.maximum(denom, 1.0e-30)
        elif model in {"current_dependent_resistivity", "anomalous_resistivity"}:
            eta0 = float(cfg.get("base_resistivity_ohm_m", self.material.electrical_resistivity_ohm_m))
            j0 = max(float(cfg.get("current_density_scale_a_m2", 1.0e12)), 1.0e-30)
            exponent = float(cfg.get("exponent", 2.0))
            eta = eta0 * (1.0 + np.power(np.abs(jz_a_m2) / j0, exponent))
            sigma = 1.0 / np.maximum(eta, 1.0e-30)
        else:
            raise ValueError(f"unsupported conductivity model: {model}")

        minimum = float(cfg.get("minimum_s_m", 1.0e3))
        sigma = np.where(self.material_mask, np.maximum(sigma, minimum), 0.0)
        return sigma

    @staticmethod
    def _fill_region_inductive_current(
        jz: np.ndarray,
        az: np.ndarray,
        previous_az: np.ndarray,
        mask: np.ndarray,
        target_current_a: float,
        conductivity: np.ndarray,
        cell_area: float,
        dt_s: float,
    ) -> None:
        if not mask.any():
            return
        sigma_region = conductivity[mask]
        dadt = (az[mask] - previous_az[mask]) / dt_s
        sigma_area = float(sigma_region.sum()) * cell_area
        induction_term = float((sigma_region * dadt).sum()) * cell_area
        ez_drive = (target_current_a + induction_term) / max(sigma_area, 1.0e-30)
        jz[mask] = sigma_region * (ez_drive - dadt)


def laplacian(values: np.ndarray, dx: float, dy: float) -> np.ndarray:
    out = np.zeros_like(values, dtype=float)
    out[1:-1, 1:-1] = (
        (values[1:-1, :-2] - 2.0 * values[1:-1, 1:-1] + values[1:-1, 2:]) / dx**2
        + (values[:-2, 1:-1] - 2.0 * values[1:-1, 1:-1] + values[2:, 1:-1]) / dy**2
    )
    return out


def solve_poisson(source: np.ndarray, dx: float, dy: float, iterations: int) -> np.ndarray:
    u = np.zeros_like(source, dtype=float)
    dx2 = dx * dx
    dy2 = dy * dy
    denom = 2.0 * (dx2 + dy2)
    for _ in range(iterations):
        u_new = u.copy()
        u_new[1:-1, 1:-1] = (
            dy2 * (u[1:-1, :-2] + u[1:-1, 2:])
            + dx2 * (u[:-2, 1:-1] + u[2:, 1:-1])
            - source[1:-1, 1:-1] * dx2 * dy2
        ) / denom
        u = u_new
    return u


def solve_electrostatic(
    shape: tuple[int, int],
    cathode: np.ndarray,
    anode: np.ndarray,
    voltage_v: float,
    dx: float,
    dy: float,
    iterations: int,
) -> np.ndarray:
    ny, nx = shape
    phi = np.linspace(voltage_v, 0.0, ny)[:, None] * np.ones((1, nx))
    fixed = cathode | anode
    dx2 = dx * dx
    dy2 = dy * dy
    denom = 2.0 * (dx2 + dy2)
    for _ in range(iterations):
        old = phi
        phi = old.copy()
        phi[1:-1, 1:-1] = (
            dy2 * (old[1:-1, :-2] + old[1:-1, 2:])
            + dx2 * (old[:-2, 1:-1] + old[2:, 1:-1])
        ) / denom
        phi[cathode] = voltage_v
        phi[anode] = 0.0
        phi[0, :] = phi[1, :]
        phi[-1, :] = phi[-2, :]
        phi[:, 0] = phi[:, 1]
        phi[:, -1] = phi[:, -2]
        phi[fixed] = np.where(cathode[fixed], voltage_v, 0.0)
    return phi


def laplacian_xp(values, dx: float, dy: float, xp):
    out = xp.zeros_like(values)
    out[1:-1, 1:-1] = (
        (values[1:-1, :-2] - 2.0 * values[1:-1, 1:-1] + values[1:-1, 2:]) / dx**2
        + (values[:-2, 1:-1] - 2.0 * values[1:-1, 1:-1] + values[2:, 1:-1]) / dy**2
    )
    return out


def solve_poisson_xp(source, dx: float, dy: float, iterations: int, xp):
    u = xp.zeros_like(source)
    dx2 = dx * dx
    dy2 = dy * dy
    denom = 2.0 * (dx2 + dy2)
    for _ in range(iterations):
        u_new = u.copy()
        u_new[1:-1, 1:-1] = (
            dy2 * (u[1:-1, :-2] + u[1:-1, 2:])
            + dx2 * (u[:-2, 1:-1] + u[2:, 1:-1])
            - source[1:-1, 1:-1] * dx2 * dy2
        ) / denom
        u = u_new
    return u


def solve_electrostatic_xp(shape, cathode, anode, voltage_v: float, dx: float, dy: float, iterations: int, xp):
    ny, nx = shape
    phi = xp.linspace(voltage_v, 0.0, ny)[:, None] * xp.ones((1, nx))
    fixed = cathode | anode
    dx2 = dx * dx
    dy2 = dy * dy
    denom = 2.0 * (dx2 + dy2)
    for _ in range(iterations):
        old = phi
        phi = old.copy()
        phi[1:-1, 1:-1] = (
            dy2 * (old[1:-1, :-2] + old[1:-1, 2:])
            + dx2 * (old[:-2, 1:-1] + old[2:, 1:-1])
        ) / denom
        phi = xp.where(cathode, voltage_v, phi)
        phi = xp.where(anode, 0.0, phi)
        phi[0, :] = phi[1, :]
        phi[-1, :] = phi[-2, :]
        phi[:, 0] = phi[:, 1]
        phi[:, -1] = phi[:, -2]
        phi = xp.where(fixed, xp.where(cathode, voltage_v, 0.0), phi)
    return phi
