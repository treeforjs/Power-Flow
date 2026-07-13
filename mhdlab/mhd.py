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
        self.temperature_floor_k = float(config.get("temperature_floor_k", 1.0))
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
        jz = self._current_density(current_a)
        az = solve_poisson(
            source=-MU0 * jz,
            dx=self.raster.dx,
            dy=self.raster.dy,
            iterations=self.poisson_iterations,
        )
        bx = np.gradient(az, self.raster.dy, axis=0)
        by = -np.gradient(az, self.raster.dx, axis=1)
        b2 = bx * bx + by * by
        mag_pressure = b2 / (2.0 * MU0)

        heat = self.material.electrical_resistivity_ohm_m * jz * jz
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
        jz = xp.asarray(self._current_density(current_a))
        temp0 = xp.asarray(state.temperature_k)
        density0 = xp.asarray(state.density_kg_m3)
        displacement0 = xp.asarray(state.surface_displacement_m)

        az = solve_poisson_xp(
            source=-MU0 * jz,
            dx=self.raster.dx,
            dy=self.raster.dy,
            iterations=self.poisson_iterations,
            xp=xp,
        )
        bx = xp.gradient(az, self.raster.dy, axis=0)
        by = -xp.gradient(az, self.raster.dx, axis=1)
        b2 = bx * bx + by * by
        mag_pressure = b2 / (2.0 * MU0)

        heat = self.material.electrical_resistivity_ohm_m * jz * jz
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

    def _current_density(self, current_a: float) -> np.ndarray:
        jz = np.zeros(self.raster.shape, dtype=float)
        cell_area = self.raster.dx * self.raster.dy
        cathode_area = max(float(self.cathode_mask.sum()) * cell_area * self.depth_m, 1e-30)
        anode_area = max(float(self.anode_mask.sum()) * cell_area * self.depth_m, 1e-30)
        jz[self.cathode_mask] = current_a / cathode_area
        jz[self.anode_mask] = -current_a / anode_area
        return jz


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
