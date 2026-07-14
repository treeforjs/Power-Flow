"""Reduced 2D MHD and electrostatic field solver."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .backend import ArrayBackend
from .constants import MU0
from .conductivity_tables import ConductivityTable
from .geometry import Raster, boundary_mask
from .materials import Material


@dataclass
class MHDState:
    temperature_k: np.ndarray
    specific_enthalpy_j_kg: np.ndarray
    density_kg_m3: np.ndarray
    pressure_pa: np.ndarray
    magnetic_pressure_pa: np.ndarray
    az_wb_m: np.ndarray
    bx_t: np.ndarray
    by_t: np.ndarray
    jz_a_m2: np.ndarray
    conductivity_s_m: np.ndarray
    joule_heating_w_m3: np.ndarray
    phi_v: np.ndarray
    ex_v_m: np.ndarray
    ey_v_m: np.ndarray
    surface_displacement_m: np.ndarray


@dataclass
class CurrentClusterMesh:
    labels: np.ndarray
    active_mask: np.ndarray
    counts: np.ndarray
    area_m2: np.ndarray
    x_m: np.ndarray
    y_m: np.ndarray
    region: np.ndarray
    green_wb_m_per_a_m2: np.ndarray
    reference_length_m: float
    min_distance_m: float


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
        self.induction_solver = str(config.get("induction_solver", "open_boundary_impedance")).lower()
        self.induction_max_unknowns = int(config.get("induction_max_unknowns", 640))
        self.induction_vector_potential_chunk = int(config.get("induction_vector_potential_chunk", 8192))
        self.enforce_unidirectional_region_current = bool(config.get("enforce_unidirectional_region_current", True))
        self.current_density_limit_factor = float(config.get("current_density_limit_factor", 25.0))
        self.temperature_floor_k = float(config.get("temperature_floor_k", 1.0))
        self.conductivity_config = dict(config.get("conductivity", {"model": "temperature"}))
        self.conductivity_table = self._load_conductivity_table(self.conductivity_config)
        self.material_mask = (
            raster.mask_by_kind("cathode")
            | raster.mask_by_kind("anode")
            | raster.mask_by_kind("material")
        )
        self.vacuum_mask = raster.mask_by_kind("vacuum") | ~self.material_mask
        self.cathode_mask = raster.mask_by_kind("cathode")
        self.anode_mask = raster.mask_by_kind("anode")
        self.surface_mask = boundary_mask(self.material_mask, self.vacuum_mask)
        self.current_cluster_mesh = self._build_current_cluster_mesh(config)

    def initial_state(self) -> MHDState:
        shape = self.raster.shape
        temp = np.full(shape, self.material.initial_temperature_k, dtype=float)
        enthalpy = self.material.specific_enthalpy_from_temperature(temp)
        density = np.where(self.material_mask, self.material.density_kg_m3, 0.0)
        zero = np.zeros(shape, dtype=float)
        return MHDState(
            temperature_k=temp,
            specific_enthalpy_j_kg=enthalpy,
            density_kg_m3=density,
            pressure_pa=zero.copy(),
            magnetic_pressure_pa=zero.copy(),
            az_wb_m=zero.copy(),
            bx_t=zero.copy(),
            by_t=zero.copy(),
            jz_a_m2=zero.copy(),
            conductivity_s_m=zero.copy(),
            joule_heating_w_m3=zero.copy(),
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
        density_for_energy = np.maximum(state.density_kg_m3, 1.0)
        enthalpy = state.specific_enthalpy_j_kg + dt_s * (
            heat / density_for_energy
            + self.material.thermal_conductivity_w_m_k * lap_t / density_for_energy
        )
        temp = self.material.temperature_from_specific_enthalpy(enthalpy)
        temp = np.where(self.material_mask, np.maximum(temp, self.temperature_floor_k), state.temperature_k)
        enthalpy = np.where(self.material_mask, enthalpy, state.specific_enthalpy_j_kg)

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
            specific_enthalpy_j_kg=enthalpy,
            density_kg_m3=density,
            pressure_pa=pressure,
            magnetic_pressure_pa=mag_pressure,
            az_wb_m=az,
            bx_t=bx,
            by_t=by,
            jz_a_m2=jz,
            conductivity_s_m=conductivity,
            joule_heating_w_m3=heat,
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
        enthalpy0 = xp.asarray(state.specific_enthalpy_j_kg)
        density0 = xp.asarray(state.density_kg_m3)
        displacement0 = xp.asarray(state.surface_displacement_m)

        az = xp.asarray(az_np)
        bx = xp.gradient(az, self.raster.dy, axis=0)
        by = -xp.gradient(az, self.raster.dx, axis=1)
        b2 = bx * bx + by * by
        mag_pressure = b2 / (2.0 * MU0)

        heat = jz * jz / xp.maximum(conductivity, 1.0e-30)
        lap_t = laplacian_xp(temp0, self.raster.dx, self.raster.dy, xp)
        density_for_energy = xp.maximum(density0, 1.0)
        enthalpy = enthalpy0 + dt_s * (
            heat / density_for_energy
            + self.material.thermal_conductivity_w_m_k * lap_t / density_for_energy
        )
        temp_np = self.material.temperature_from_specific_enthalpy(self.backend.asnumpy(enthalpy))
        temp = xp.asarray(temp_np)
        temp = xp.where(material_mask, xp.maximum(temp, self.temperature_floor_k), temp0)
        enthalpy = xp.where(material_mask, enthalpy, enthalpy0)

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
            specific_enthalpy_j_kg=self.backend.asnumpy(enthalpy),
            density_kg_m3=self.backend.asnumpy(density),
            pressure_pa=self.backend.asnumpy(pressure),
            magnetic_pressure_pa=self.backend.asnumpy(mag_pressure),
            az_wb_m=self.backend.asnumpy(az),
            bx_t=self.backend.asnumpy(bx),
            by_t=self.backend.asnumpy(by),
            jz_a_m2=self.backend.asnumpy(jz),
            conductivity_s_m=self.backend.asnumpy(conductivity),
            joule_heating_w_m3=self.backend.asnumpy(heat),
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
        if self.induction_solver in {"open", "open_boundary", "open_boundary_impedance", "impedance", "impedance_matrix"}:
            return self._solve_open_boundary_impedance_current(state, current_a, dt_s)
        return self._solve_local_poisson_inductive_current(state, current_a, dt_s)

    def _solve_local_poisson_inductive_current(
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
            self._enforce_region_current_constraints(jz, self.cathode_mask, float(current_a), cell_area)
            self._enforce_region_current_constraints(jz, self.anode_mask, -float(current_a), cell_area)
            next_az = solve_poisson(
                source=-MU0 * jz,
                dx=self.raster.dx,
                dy=self.raster.dy,
                iterations=self.poisson_iterations,
            )
            az = relaxation * next_az + (1.0 - relaxation) * az
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
        self._enforce_region_current_constraints(jz, self.cathode_mask, float(current_a), cell_area)
        self._enforce_region_current_constraints(jz, self.anode_mask, -float(current_a), cell_area)
        return jz, az

    def _solve_open_boundary_impedance_current(
        self,
        state: MHDState,
        current_a: float,
        dt_s: float,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Solve the driven conductor current with an open-boundary Green matrix.

        The cluster solve is the discrete form of

            J/sigma + (A^{n+1} - A^n)/dt = E_region

        with one scalar drive field per electrode and integral current
        constraints.  The Green matrix maps every current cluster to vector
        potential with the 2D open-boundary logarithmic kernel, so the x-profile
        comes from mutual/self inductance instead of a nearby box Dirichlet
        boundary.
        """
        mesh = self.current_cluster_mesh
        if mesh is None or mesh.area_m2.size == 0:
            zero = np.zeros(self.raster.shape, dtype=float)
            return zero.copy(), zero.copy()

        dt_s = max(float(dt_s), 1.0e-30)
        conductivity = self._electrical_conductivity(state.temperature_k, state.density_kg_m3, state.jz_a_m2)
        sigma_cluster = self._cluster_harmonic_conductivity(conductivity, mesh)
        previous_az_cluster = self._cluster_average(np.asarray(state.az_wb_m, dtype=float), mesh)
        j_cluster = self._solve_cluster_currents(
            mesh,
            sigma_cluster,
            previous_az_cluster,
            cathode_current_a=float(current_a),
            anode_current_a=-float(current_a),
            dt_s=dt_s,
        )
        self._enforce_cluster_current_constraints(j_cluster, mesh, region=1, target_current_a=float(current_a))
        self._enforce_cluster_current_constraints(j_cluster, mesh, region=-1, target_current_a=-float(current_a))

        jz = np.zeros(self.raster.shape, dtype=float)
        active_labels = mesh.labels[mesh.active_mask]
        jz[mesh.active_mask] = j_cluster[active_labels]
        az = self._vector_potential_from_clusters(mesh, j_cluster)
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
        elif model in {"temperature", "knoepfel", "knoepfel_like", "ec_knoepfel"}:
            beta_cv = float(cfg.get("betacv_1_k", cfg.get("temperature_coefficient_1_k", 9.4e-4)))
            t0 = float(cfg.get("reference_temperature_k", self.material.initial_temperature_k))
            temperature_cutoff_k = float(cfg.get("temperature_cutoff_k", 100.0))
            density_exponent = float(cfg.get("density_exponent", cfg.get("alpha", 0.0)))
            effective_temperature = np.maximum(temperature_k, temperature_cutoff_k)
            denom = 1.0 + beta_cv * np.maximum(effective_temperature - t0, 0.0)
            sigma = sigma0 * np.power(np.maximum(rho_ratio, 0.0), density_exponent) / np.maximum(denom, 1.0e-30)
        elif model in {"current_dependent_resistivity", "anomalous_resistivity"}:
            eta0 = float(cfg.get("base_resistivity_ohm_m", self.material.electrical_resistivity_ohm_m))
            j0 = max(float(cfg.get("current_density_scale_a_m2", 1.0e12)), 1.0e-30)
            exponent = float(cfg.get("exponent", 2.0))
            eta = eta0 * (1.0 + np.power(np.abs(jz_a_m2) / j0, exponent))
            sigma = 1.0 / np.maximum(eta, 1.0e-30)
        elif model in {"table", "tabular", "ethos_table", "conductivity_table"}:
            if self.conductivity_table is None:
                raise ValueError("conductivity model 'table' requires conductivity.file")
            sigma = self.conductivity_table.interpolate(
                density_kg_m3,
                temperature_k,
                ensemble_index=cfg.get("ensemble_index"),
                statistic=cfg.get("ensemble_statistic"),
                log_interpolation=bool(cfg.get("log_interpolation", True)),
            )
            sigma *= float(cfg.get("multiplier", 1.0))
        else:
            raise ValueError(f"unsupported conductivity model: {model}")

        minimum = float(cfg.get("minimum_s_m", 1.0e3))
        sigma = np.where(self.material_mask, np.maximum(sigma, minimum), 0.0)
        return sigma

    @staticmethod
    def _load_conductivity_table(cfg: dict) -> ConductivityTable | None:
        model = str(cfg.get("model", "constant")).lower()
        if model not in {"table", "tabular", "ethos_table", "conductivity_table"}:
            return None
        path = cfg.get("file") or cfg.get("table") or cfg.get("path")
        if not path:
            raise ValueError("conductivity table model requires a file/table/path entry")
        return ConductivityTable.from_file(path)

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

    def _enforce_region_current_constraints(
        self,
        jz: np.ndarray,
        mask: np.ndarray,
        target_current_a: float,
        cell_area: float,
    ) -> None:
        if not self.enforce_unidirectional_region_current or not mask.any():
            return
        sign = 1.0 if target_current_a >= 0.0 else -1.0
        local = sign * jz[mask]
        local = np.maximum(local, 0.0)
        mean_abs = abs(target_current_a) / max(float(mask.sum()) * cell_area, 1.0e-30)
        limit_factor = max(self.current_density_limit_factor, 1.0)
        local = np.minimum(local, limit_factor * mean_abs)
        total = float(local.sum()) * cell_area
        if total <= 0.0:
            local = np.full(mask.sum(), mean_abs, dtype=float)
            total = float(local.sum()) * cell_area
        local *= abs(target_current_a) / max(total, 1.0e-30)
        jz[mask] = sign * local

    def _build_current_cluster_mesh(self, config: dict) -> CurrentClusterMesh | None:
        active = self.cathode_mask | self.anode_mask
        if not active.any():
            return None
        labels = np.full(self.raster.shape, -1, dtype=int)
        max_unknowns = max(int(config.get("induction_max_unknowns", self.induction_max_unknowns)), 8)
        region_masks = [(1, self.cathode_mask), (-1, self.anode_mask)]
        active_regions = [(region, mask) for region, mask in region_masks if mask.any()]
        target_per_region = max(4, int(np.ceil(max_unknowns / max(len(active_regions), 1))))

        key_blocks = []
        y_blocks = []
        x_blocks = []
        for region, mask in active_regions:
            yy_region, xx_region = np.nonzero(mask)
            x_cells = self.raster.x[xx_region]
            y_cells = self.raster.y[yy_region]
            xmin = float(x_cells.min() - 0.5 * self.raster.dx)
            xmax = float(x_cells.max() + 0.5 * self.raster.dx)
            ymin = float(y_cells.min() - 0.5 * self.raster.dy)
            ymax = float(y_cells.max() + 0.5 * self.raster.dy)
            width = max(xmax - xmin, self.raster.dx)
            height = max(ymax - ymin, self.raster.dy)
            aspect = max(width / height, 1.0)
            ny_bins = max(1, int(round(np.sqrt(target_per_region / aspect))))
            nx_bins = max(1, int(np.ceil(target_per_region / ny_bins)))
            nx_bins = min(nx_bins, np.unique(xx_region).size)
            ny_bins = min(ny_bins, np.unique(yy_region).size)
            bx = np.floor((x_cells - xmin) / width * nx_bins).astype(int)
            by = np.floor((y_cells - ymin) / height * ny_bins).astype(int)
            bx = np.clip(bx, 0, nx_bins - 1)
            by = np.clip(by, 0, ny_bins - 1)
            key_blocks.append(np.column_stack((np.full(xx_region.size, region), by, bx)))
            y_blocks.append(yy_region)
            x_blocks.append(xx_region)

        keys = np.vstack(key_blocks)
        yy = np.concatenate(y_blocks)
        xx = np.concatenate(x_blocks)
        unique, inverse = np.unique(keys, axis=0, return_inverse=True)
        labels[yy, xx] = inverse

        n = unique.shape[0]
        counts = np.bincount(inverse, minlength=n).astype(float)
        cell_area = self.raster.dx * self.raster.dy
        area = counts * cell_area
        x_cells = self.raster.x[xx]
        y_cells = self.raster.y[yy]
        x = np.bincount(inverse, weights=x_cells, minlength=n) / np.maximum(counts, 1.0)
        y = np.bincount(inverse, weights=y_cells, minlength=n) / np.maximum(counts, 1.0)
        region = unique[:, 0].astype(int)

        width = max(self.raster.x_max - self.raster.x_min, self.raster.dx)
        height = max(self.raster.y_max - self.raster.y_min, self.raster.dy)
        reference_length = float(config.get("magnetic_reference_length_m", max(width, height)))
        min_distance = float(config.get("magnetic_min_distance_m", 0.25 * min(self.raster.dx, self.raster.dy)))
        green = open_boundary_green_matrix(
            x,
            y,
            area,
            reference_length_m=reference_length,
            min_distance_m=min_distance,
        )
        return CurrentClusterMesh(
            labels=labels,
            active_mask=active,
            counts=counts,
            area_m2=area,
            x_m=x,
            y_m=y,
            region=region,
            green_wb_m_per_a_m2=green,
            reference_length_m=reference_length,
            min_distance_m=min_distance,
        )

    @staticmethod
    def _cluster_average(values: np.ndarray, mesh: CurrentClusterMesh) -> np.ndarray:
        labels = mesh.labels[mesh.active_mask]
        sums = np.bincount(labels, weights=values[mesh.active_mask], minlength=mesh.area_m2.size)
        return sums / np.maximum(mesh.counts, 1.0)

    def _cluster_harmonic_conductivity(self, conductivity: np.ndarray, mesh: CurrentClusterMesh) -> np.ndarray:
        labels = mesh.labels[mesh.active_mask]
        sigma = np.maximum(conductivity[mesh.active_mask], 1.0e-30)
        resistivity_sum = np.bincount(labels, weights=1.0 / sigma, minlength=mesh.area_m2.size)
        mean_resistivity = resistivity_sum / np.maximum(mesh.counts, 1.0)
        return 1.0 / np.maximum(mean_resistivity, 1.0e-30)

    @staticmethod
    def _solve_cluster_currents(
        mesh: CurrentClusterMesh,
        sigma_cluster: np.ndarray,
        previous_az_cluster: np.ndarray,
        cathode_current_a: float,
        anode_current_a: float,
        dt_s: float,
    ) -> np.ndarray:
        n = mesh.area_m2.size
        system = np.zeros((n + 2, n + 2), dtype=float)
        rhs = np.zeros(n + 2, dtype=float)
        system[:n, :n] = mesh.green_wb_m_per_a_m2 / dt_s
        diag = np.diag_indices(n)
        system[diag] += 1.0 / np.maximum(sigma_cluster, 1.0e-30)
        system[:n, n] = -((mesh.region == 1).astype(float))
        system[:n, n + 1] = -((mesh.region == -1).astype(float))
        rhs[:n] = previous_az_cluster / dt_s
        cathode = mesh.region == 1
        anode = mesh.region == -1
        system[n, :n] = np.where(cathode, mesh.area_m2, 0.0)
        system[n + 1, :n] = np.where(anode, mesh.area_m2, 0.0)
        rhs[n] = cathode_current_a
        rhs[n + 1] = anode_current_a
        try:
            solution = np.linalg.solve(system, rhs)
        except np.linalg.LinAlgError:
            solution = np.linalg.lstsq(system, rhs, rcond=None)[0]
        return solution[:n]

    def _enforce_cluster_current_constraints(
        self,
        j_cluster: np.ndarray,
        mesh: CurrentClusterMesh,
        region: int,
        target_current_a: float,
    ) -> None:
        if not self.enforce_unidirectional_region_current:
            return
        mask = mesh.region == region
        if not mask.any():
            return
        sign = 1.0 if target_current_a >= 0.0 else -1.0
        local = sign * j_cluster[mask]
        local = np.maximum(local, 0.0)
        mean_abs = abs(target_current_a) / max(float(mesh.area_m2[mask].sum()), 1.0e-30)
        local = np.minimum(local, max(self.current_density_limit_factor, 1.0) * mean_abs)
        total = float((local * mesh.area_m2[mask]).sum())
        if total <= 0.0:
            local = np.full(mask.sum(), mean_abs, dtype=float)
            total = float((local * mesh.area_m2[mask]).sum())
        local *= abs(target_current_a) / max(total, 1.0e-30)
        j_cluster[mask] = sign * local

    def _vector_potential_from_clusters(self, mesh: CurrentClusterMesh, j_cluster: np.ndarray) -> np.ndarray:
        current = j_cluster * mesh.area_m2
        if not np.any(current):
            return np.zeros(self.raster.shape, dtype=float)
        points_x = self.raster.xx.ravel()
        points_y = self.raster.yy.ravel()
        out = np.zeros(points_x.size, dtype=float)
        chunk = max(int(self.induction_vector_potential_chunk), 1)
        coefficient = -MU0 / (2.0 * np.pi)
        for start in range(0, points_x.size, chunk):
            stop = min(start + chunk, points_x.size)
            dx = points_x[start:stop, None] - mesh.x_m[None, :]
            dy = points_y[start:stop, None] - mesh.y_m[None, :]
            r = np.sqrt(dx * dx + dy * dy)
            r = np.maximum(r, mesh.min_distance_m)
            out[start:stop] = coefficient * (np.log(r / mesh.reference_length_m) @ current)
        return out.reshape(self.raster.shape)


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


def open_boundary_green_matrix(
    x_m: np.ndarray,
    y_m: np.ndarray,
    area_m2: np.ndarray,
    reference_length_m: float,
    min_distance_m: float,
) -> np.ndarray:
    dx = x_m[:, None] - x_m[None, :]
    dy = y_m[:, None] - y_m[None, :]
    r = np.sqrt(dx * dx + dy * dy)
    self_radius = np.sqrt(np.maximum(area_m2, 1.0e-300) / np.pi)
    r = np.where(np.eye(x_m.size, dtype=bool), self_radius[None, :], r)
    r = np.maximum(r, max(min_distance_m, 1.0e-30))
    reference = max(float(reference_length_m), max(min_distance_m, 1.0e-30))
    return (-MU0 / (2.0 * np.pi)) * np.log(r / reference) * area_m2[None, :]


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
