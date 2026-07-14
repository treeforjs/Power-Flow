"""High-level simulation orchestration."""

from __future__ import annotations

import json
import shutil
from collections.abc import Callable
from datetime import datetime
from pathlib import Path

import numpy as np

from .bolsig import BolsigRunner, approximate_bolsig_table, supplement_missing_rate_coefficients
from .chemistry import Mechanism, line_emissivity, local_reaction_rates
from .config import resolve_path
from .constants import KB, QE
from .cr_data import CRDataLibrary
from .cross_sections import CrossSectionLibrary
from .diagnostics import HDF5DiagnosticWriter, load_hdf5_arrays
from .desorption import TemkinDesorption
from .geometry import Geometry
from .materials import material_from_config
from .mhd import ReducedMHDSolver
from .neutrals import KineticNeutralSolver, VelocityGrid
from .spectra import CCDCalibration, integrate_los, synthesize_spectrum
from .traces import drive_from_parametric_current, drive_from_trace, load_trace_csv


def run_from_config(config: dict, progress: Callable[[str], None] | None = None) -> dict:
    root = resolve_path(config, config.get("project_root", ".")) or Path.cwd()
    out_root = resolve_path(config, config.get("output_root", str(root / "runs"))) or (root / "runs")
    run_dir = out_root / datetime.now().strftime("run_%Y%m%d_%H%M%S")
    run_dir.mkdir(parents=True, exist_ok=True)
    _progress(progress, f"created run directory: {run_dir}")

    config_path = config.get("_config_path")
    if config_path:
        shutil.copy2(config_path, run_dir / "config_used.yaml")

    geometry = Geometry.from_json(resolve_path(config, config["geometry"]["file"]))
    grid_cfg = config.get("grid", {})
    nx, ny, grid_note = _grid_shape_from_config(geometry, grid_cfg)
    raster = geometry.rasterize(nx=nx, ny=ny)
    _progress(
        progress,
        f"rasterized geometry on {raster.nx} x {raster.ny} cells "
        f"(dx={raster.dx * 1.0e6:.1f} um, dy={raster.dy * 1.0e6:.1f} um; {grid_note})",
    )
    material = material_from_config(config.get("material"))
    mhd_config = _resolve_mhd_config(config)
    mhd_solver = ReducedMHDSolver(raster, material, mhd_config)
    mhd_state = mhd_solver.initial_state()
    _progress(progress, f"MHD array backend: {mhd_solver.backend.name}")
    _progress(progress, _conductivity_status(mhd_config))

    mechanism = Mechanism.from_yaml(resolve_path(config, config["chemistry"]["mechanism"]))
    species = mechanism.species or ["H2O", "H", "H2", "O", "O2", "OH"]
    vcfg = config.get("velocity_grid", {})
    velocity_grid = VelocityGrid.polar(
        max_speed_m_s=float(vcfg.get("max_speed_m_s", 8000.0)),
        n_speed=int(vcfg.get("n_speed", 4)),
        n_angle=int(vcfg.get("n_angle", 12)),
    )
    neutral_solver = KineticNeutralSolver(raster, species, velocity_grid, config.get("neutrals", {}))
    _progress(progress, f"neutral kinetics backend: {neutral_solver.backend.name}")
    _progress(progress, f"estimated neutral velocity-space working set: {neutral_solver.estimated_working_set_gb():.2f} GiB")
    neutral_state = neutral_solver.initial_state()
    _progress(progress, f"initialized {len(species)} neutral species on {velocity_grid.vx.size} velocity ordinates")

    desorption = TemkinDesorption(**config.get("desorption", {}))
    coverage = desorption.initial_coverage(raster.shape)

    bolsig_table = _load_bolsig_table(config, run_dir, progress)
    allow_provisional = bool(config.get("chemistry", {}).get("allow_provisional_rates", False))
    base_rates = local_reaction_rates(mechanism, bolsig_table.rate_coefficients, allow_provisional)
    cross_section_library = _load_cross_sections(config)
    cr_data_library = CRDataLibrary.from_config(config, resolve_path)
    incident_energies_ev = {
        str(k): float(v)
        for k, v in config.get("cross_sections", {}).get("incident_energies_ev", {}).items()
    }

    t_end = float(config["time"]["end_s"])
    dt = float(config["time"]["dt_s"])
    if "sample_interval_s" in config["time"]:
        sample_every = max(1, int(round(float(config["time"]["sample_interval_s"]) / dt)))
    else:
        sample_every = int(config["time"].get("sample_every", 1))
    times = np.arange(0.0, t_end + 0.5 * dt, dt)
    drive = _load_drive_profile(config, times)
    rl_fit = drive.metadata
    _progress(progress, f"loaded {rl_fit.get('mode', 'drive')} drive with peak current {rl_fit.get('peak_current_a', 0.0) / 1.0e3:.1f} kA")

    diag_cfg = config.get("diagnostics", {})
    output_dtype = np.dtype(diag_cfg.get("output_dtype", "float32"))
    hdf5_path = run_dir / str(diag_cfg.get("hdf5_file", "fields.h5"))
    expected_samples = _expected_sample_count(len(times), sample_every)
    keep_samples = bool(diag_cfg.get("keep_samples_in_memory", False))
    samples = [] if keep_samples else None
    last_sample = None
    progress_every_samples = max(1, int(diag_cfg.get("progress_every_samples", 10)))
    _progress(
        progress,
        f"starting timestep loop: {len(times)} steps, dt={dt * 1.0e9:.3f} ns, "
        f"{expected_samples} saved samples to {hdf5_path.name}",
    )

    writer_metadata = {
        "nx": raster.nx,
        "ny": raster.ny,
        "dx_m": raster.dx,
        "dy_m": raster.dy,
        "x_min_m": raster.x_min,
        "x_max_m": raster.x_max,
        "y_min_m": raster.y_min,
        "y_max_m": raster.y_max,
        "sample_every_steps": sample_every,
        "expected_samples": expected_samples,
    }
    with HDF5DiagnosticWriter(
        hdf5_path,
        raster.shape,
        species,
        dtype=output_dtype,
        compression=diag_cfg.get("hdf5_compression", "gzip"),
        compression_level=diag_cfg.get("hdf5_compression_level", 4),
        metadata=writer_metadata,
    ) as writer:
        _write_diagnostics_manifest(run_dir, hdf5_path, 0, expected_samples, None, complete=False)
        for step, t in enumerate(times):
            voltage = drive.sample_voltage(t)
            current = drive.sample_current(t)
            mhd_state = mhd_solver.step(mhd_state, current_a=current, voltage_v=voltage, dt_s=dt)
            pre_collision_neutral = neutral_solver.asnumpy(neutral_state.total_neutral_density())
            electron_mean_energy_ev = _electron_mean_energy_map(
                config,
                bolsig_table,
                mhd_state,
                mhd_solver.surface_mask,
                mhd_solver.vacuum_mask,
                pre_collision_neutral,
            )
            source_rate, coverage = desorption.step(
                mhd_state.temperature_k,
                coverage,
                mhd_solver.surface_mask,
                dt_s=dt,
            )
            neutral_state = neutral_solver.step(
                neutral_state,
                dt_s=dt,
                surface_source_m2_s=source_rate,
                surface_temperature_k=mhd_state.temperature_k,
                source_species=desorption.species,
                surface_mask=mhd_solver.surface_mask,
                reaction_rates=base_rates,
                cross_sections=cross_section_library,
                incident_energies_ev=incident_energies_ev,
                electron_energy_ev=electron_mean_energy_ev,
            )
            if step % sample_every == 0 or step == len(times) - 1:
                moments = neutral_solver.moments(neutral_state)
                total_neutral = neutral_solver.asnumpy(neutral_state.total_neutral_density())
                en_td = _reduced_field_td(mhd_state.ex_v_m, mhd_state.ey_v_m, total_neutral)
                sample = {
                    "time_s": t,
                    "temperature_k": _sample_array(mhd_state.temperature_k, output_dtype),
                    "specific_enthalpy_j_kg": _sample_array(mhd_state.specific_enthalpy_j_kg, output_dtype),
                    "pressure_pa": _sample_array(mhd_state.pressure_pa, output_dtype),
                    "bx_t": _sample_array(mhd_state.bx_t, output_dtype),
                    "by_t": _sample_array(mhd_state.by_t, output_dtype),
                    "jz_a_m2": _sample_array(mhd_state.jz_a_m2, output_dtype),
                    "conductivity_s_m": _sample_array(mhd_state.conductivity_s_m, output_dtype),
                    "joule_heating_w_m3": _sample_array(mhd_state.joule_heating_w_m3, output_dtype),
                    "surface_displacement_m": _sample_array(mhd_state.surface_displacement_m, output_dtype),
                    "electron_density_m3": _sample_array(neutral_solver.asnumpy(neutral_state.electron_density_m3), output_dtype),
                    "electron_mean_energy_ev": _sample_array(electron_mean_energy_ev, output_dtype),
                    "total_neutral_density_m3": _sample_array(total_neutral, output_dtype),
                    "en_td": _sample_array(en_td, output_dtype),
                    "species_density_m3": {
                        sp: _sample_array(neutral_solver.asnumpy(data["density_m3"]), output_dtype)
                        for sp, data in moments.items()
                    },
                }
                writer.append(sample)
                last_sample = sample
                if keep_samples and samples is not None:
                    samples.append(sample)
                _write_diagnostics_manifest(
                    run_dir,
                    hdf5_path,
                    writer.sample_count,
                    expected_samples,
                    float(t),
                    complete=False,
                )
                if writer.sample_count == 1 or writer.sample_count % progress_every_samples == 0 or step == len(times) - 1:
                    _progress(
                        progress,
                        f"sample {writer.sample_count}/{expected_samples}: "
                        f"t={t * 1.0e9:.1f} ns, I={current / 1.0e3:.1f} kA, "
                        f"max T={float(np.nanmax(mhd_state.temperature_k)):.1f} K, "
                        f"max |Jz|={float(np.nanmax(np.abs(mhd_state.jz_a_m2))):.3e} A/m^2, "
                        f"max neutral={float(np.nanmax(total_neutral)):.3e} m^-3",
                    )

    if last_sample is None:
        raise RuntimeError("no diagnostic samples were written")
    _write_diagnostics_manifest(run_dir, hdf5_path, expected_samples, expected_samples, float(last_sample["time_s"]), complete=True)
    if bool(diag_cfg.get("write_npz_compat", False)):
        _progress(progress, "writing fields.npz compatibility archive")
        _save_hdf5_as_npz(hdf5_path, run_dir / "fields.npz")
    _progress(progress, "writing spectra, summary, and quick-look plots")
    spectra_summary = _write_spectra_products(config, run_dir, raster, mechanism, last_sample)
    _write_summary(
        run_dir,
        rl_fit,
        bolsig_table,
        base_rates,
        spectra_summary,
        cross_section_library,
        cr_data_library,
    )
    _write_plots(run_dir, last_sample, raster)
    _progress(progress, f"complete: {run_dir}")
    return {
        "run_dir": str(run_dir),
        "fields_hdf5": str(hdf5_path),
        "sample_count": expected_samples,
        "samples": samples if samples is not None else [],
        "rl_fit": rl_fit,
        "spectra": spectra_summary,
    }


def _load_drive_profile(config: dict, times: np.ndarray):
    dcfg = config.get("drive", {})
    mode = str(dcfg.get("mode", "fit_rl")).lower()
    if mode == "parametric_current":
        return drive_from_parametric_current(times, dcfg)
    trace_path = resolve_path(config, dcfg["trace_csv"])
    trace = load_trace_csv(trace_path)
    return drive_from_trace(
        trace,
        voltage_column=dcfg.get("voltage_column", "voltage_v"),
        current_column=dcfg.get("current_column", "feed_current_a"),
        mode=mode,
    )


def _sample_array(values: np.ndarray, dtype: np.dtype) -> np.ndarray:
    return np.asarray(values, dtype=dtype).copy()


def _grid_shape_from_config(geometry: Geometry, grid_cfg: dict) -> tuple[int, int, str]:
    xmin, xmax, ymin, ymax = geometry.bounds
    width = xmax - xmin
    height = ymax - ymin
    nx = grid_cfg.get("nx")
    ny = grid_cfg.get("ny")
    note = "explicit grid"
    target_cell_m = None

    if "target_cell_size_m" in grid_cfg:
        target_cell_m = float(grid_cfg["target_cell_size_m"])
    elif "target_cell_size_um" in grid_cfg:
        target_cell_m = float(grid_cfg["target_cell_size_um"]) * 1.0e-6

    min_cells = grid_cfg.get("min_cells_per_material_thickness", grid_cfg.get("min_cells_per_foil_thickness"))
    material_scales = geometry.material_length_scales()
    if min_cells is not None and material_scales:
        material_cell = min(material_scales) / max(float(min_cells), 1.0)
        target_cell_m = material_cell if target_cell_m is None else min(target_cell_m, material_cell)

    if target_cell_m is not None:
        if target_cell_m <= 0.0:
            raise ValueError("grid target cell size must be positive")
        nx = int(np.ceil(width / target_cell_m))
        ny = int(np.ceil(height / target_cell_m))
        note = f"target cell <= {target_cell_m * 1.0e6:.2f} um"

    nx = int(nx if nx is not None else 64)
    ny = int(ny if ny is not None else 48)
    max_cells = grid_cfg.get("max_cells")
    if max_cells is not None and nx * ny > int(max_cells):
        scale = np.sqrt((nx * ny) / max(int(max_cells), 1))
        nx = max(4, int(np.floor(nx / scale)))
        ny = max(4, int(np.floor(ny / scale)))
        note += f", capped at {int(max_cells)} cells"
    return nx, ny, note


def _resolve_mhd_config(config: dict) -> dict:
    mhd_config = dict(config.get("mhd", {}))
    conductivity = dict(mhd_config.get("conductivity", {}))
    model = str(conductivity.get("model", "")).lower()
    if model in {"table", "tabular", "ethos_table", "conductivity_table"}:
        for key in ("file", "table", "path"):
            if key in conductivity:
                conductivity[key] = str(resolve_path(config, conductivity[key]))
                break
    if conductivity:
        mhd_config["conductivity"] = conductivity
    return mhd_config


def _conductivity_status(mhd_config: dict) -> str:
    cfg = mhd_config.get("conductivity", {})
    model = str(cfg.get("model", "constant"))
    if model.lower() in {"table", "tabular", "ethos_table", "conductivity_table"}:
        source = cfg.get("file") or cfg.get("table") or cfg.get("path")
        ensemble = cfg.get("ensemble_index", cfg.get("ensemble_statistic", 0))
        return f"using tabular conductivity model from {source} (ensemble/stat={ensemble})"
    return f"using conductivity model: {model}"


def _load_bolsig_table(config: dict, run_dir: Path, progress: Callable[[str], None] | None = None):
    bcfg = config.get("bolsig", {})
    if not bcfg.get("enabled", True):
        _progress(progress, "BOLSIG disabled; using approximate electron-kinetics table")
        return approximate_bolsig_table()
    bolsig_dir = resolve_path(config, bcfg.get("path", "third_party/bolsigplus"))
    root = resolve_path(config, config.get("project_root", ".")) or Path.cwd()
    cache_dir = resolve_path(config, bcfg.get("cache_dir")) if bcfg.get("cache_dir") else root / ".cache" / "bolsigplus"
    runner = BolsigRunner(bolsig_dir, cache_dir)
    _progress(progress, f"loading BOLSIG table from shared cache: {cache_dir}")
    try:
        table = runner.run_table(
            species=list(bcfg.get("species", ["H2", "O2"])),
            fractions=list(bcfg.get("fractions", [0.5, 0.5])),
            en_min_td=float(bcfg.get("en_min_td", 0.1)),
            en_max_td=float(bcfg.get("en_max_td", 1000.0)),
            count=int(bcfg.get("count", 20)),
            gas_temperature_k=float(bcfg.get("gas_temperature_k", 300.0)),
            gas_density_m3=float(bcfg.get("gas_density_m3", 3.295e22)),
            timeout_s=float(bcfg.get("timeout_s", 20.0)),
        )
        if bcfg.get("supplement_missing_provisional_rates", True):
            table = supplement_missing_rate_coefficients(table, approximate_bolsig_table())
        if table.warning:
            _progress(progress, f"BOLSIG warning: {table.warning}")
        else:
            _progress(progress, f"BOLSIG ready: {table.reduced_field_td.size} E/N points from {table.source_file}")
        return table
    except Exception as exc:
        (run_dir / "bolsig_warning.txt").write_text(f"Falling back to approximate table: {exc}\n", encoding="utf-8")
        _progress(progress, f"BOLSIG failed; using approximate table ({exc})")
        return approximate_bolsig_table()


def _load_cross_sections(config: dict) -> CrossSectionLibrary | None:
    xcfg = config.get("cross_sections", {})
    if not xcfg.get("enabled", False):
        return None
    manifest = resolve_path(config, xcfg.get("manifest"))
    if manifest is None:
        return None
    return CrossSectionLibrary.from_manifest(manifest)


def _reduced_field_td(ex: np.ndarray, ey: np.ndarray, neutral_density: np.ndarray) -> np.ndarray:
    e_mag = np.sqrt(ex * ex + ey * ey)
    n = np.maximum(neutral_density, 1.0e10)
    return e_mag / n / 1.0e-21


def _electron_mean_energy_map(
    config: dict,
    bolsig_table,
    mhd_state,
    surface_mask: np.ndarray,
    vacuum_mask: np.ndarray,
    neutral_density: np.ndarray,
) -> np.ndarray:
    xcfg = config.get("cross_sections", {})
    ecfg = xcfg.get("electron_energy", {})
    mode = str(ecfg.get("source", "bolsig_with_thermionic_floor")).lower()
    en_td = _reduced_field_td(mhd_state.ex_v_m, mhd_state.ey_v_m, neutral_density)
    bolsig_energy = bolsig_table.interp_mean_energy(en_td)
    if mode in {"bolsig", "eedf"}:
        return bolsig_energy
    thermionic = _thermionic_electron_energy_ev(
        mhd_state.temperature_k,
        surface_mask,
        vacuum_mask,
        flux_mean_factor=float(ecfg.get("thermionic_flux_mean_factor", 2.0)),
        min_temperature_k=float(ecfg.get("thermionic_min_temperature_k", 0.0)),
    )
    if mode in {"thermionic", "thermionic_emission"}:
        return thermionic
    if mode in {"bolsig_with_thermionic_floor", "auto", "max"}:
        return np.maximum(bolsig_energy, thermionic)
    raise ValueError(f"unsupported cross_sections.electron_energy.source: {mode}")


def _thermionic_electron_energy_ev(
    temperature_k: np.ndarray,
    surface_mask: np.ndarray,
    vacuum_mask: np.ndarray,
    flux_mean_factor: float = 2.0,
    min_temperature_k: float = 0.0,
) -> np.ndarray:
    surface_energy = flux_mean_factor * KB * np.maximum(temperature_k, min_temperature_k) / QE
    total = np.zeros_like(temperature_k, dtype=float)
    count = np.zeros_like(temperature_k, dtype=float)
    shifts = (
        (slice(0, -1), slice(None), slice(1, None), slice(None)),
        (slice(1, None), slice(None), slice(0, -1), slice(None)),
        (slice(None), slice(0, -1), slice(None), slice(1, None)),
        (slice(None), slice(1, None), slice(None), slice(0, -1)),
    )
    for dst_y, dst_x, src_y, src_x in shifts:
        adjacent = surface_mask[src_y, src_x] & vacuum_mask[dst_y, dst_x]
        total[dst_y, dst_x] += np.where(adjacent, surface_energy[src_y, src_x], 0.0)
        count[dst_y, dst_x] += adjacent.astype(float)
    return np.divide(total, count, out=np.zeros_like(total), where=count > 0.0)


def _save_samples(path: Path, samples: list[dict]) -> None:
    arrays = {"time_s": np.asarray([s["time_s"] for s in samples])}
    for key in (
        "temperature_k",
        "specific_enthalpy_j_kg",
        "pressure_pa",
        "bx_t",
        "by_t",
        "jz_a_m2",
        "conductivity_s_m",
        "joule_heating_w_m3",
        "surface_displacement_m",
        "electron_density_m3",
        "electron_mean_energy_ev",
        "total_neutral_density_m3",
        "en_td",
    ):
        arrays[key] = np.asarray([s[key] for s in samples])
    species = sorted(samples[-1]["species_density_m3"])
    for sp in species:
        arrays[f"species_{sp}_density_m3"] = np.asarray([s["species_density_m3"][sp] for s in samples])
    np.savez_compressed(path, **arrays)


def _write_spectra_products(config: dict, run_dir: Path, raster, mechanism: Mechanism, last_sample: dict) -> dict:
    scfg = config.get("spectra", {})
    if not scfg.get("enabled", True):
        return {"enabled": False}
    cal = CCDCalibration(**scfg.get("calibration", {}))
    n_pixels = int(scfg.get("n_pixels", 1024))
    wavelength = cal.wavelength_axis(n_pixels)
    emiss_maps = line_emissivity(
        mechanism.spectral_lines,
        last_sample["electron_density_m3"],
        last_sample["species_density_m3"],
        excitation_scale=float(scfg.get("excitation_scale", 1.0e-21)),
    )
    strengths = {}
    for line in mechanism.spectral_lines:
        if geometry_los := scfg.get("los"):
            strengths[line.name] = sum(integrate_los(emiss_maps[line.name], raster, los) for los in geometry_los)
        else:
            strengths[line.name] = float(emiss_maps[line.name].sum() * raster.dx * raster.dy)
    synthetic = synthesize_spectrum(wavelength, mechanism.spectral_lines, strengths, cal)
    np.savetxt(
        run_dir / "synthetic_spectrum.csv",
        np.column_stack([wavelength, synthetic]),
        delimiter=",",
        header="wavelength_nm,intensity",
        comments="",
    )
    return {"enabled": True, "lines": strengths, "spectrum_csv": "synthetic_spectrum.csv"}


def _write_summary(
    run_dir: Path,
    rl_fit: dict,
    bolsig_table,
    rates: dict,
    spectra_summary: dict,
    cross_section_library: CrossSectionLibrary | None,
    cr_data_library: CRDataLibrary | None,
) -> None:
    summary = {
        "rl_fit": rl_fit,
        "bolsig_source": bolsig_table.source_file,
        "bolsig": {
            "source": bolsig_table.source_file,
            "log_file": bolsig_table.log_file,
            "return_code": bolsig_table.return_code,
            "warning": bolsig_table.warning,
            "points": int(bolsig_table.reduced_field_td.size),
            "rate_coefficient_count": len(bolsig_table.rate_coefficients),
            "rate_coefficient_names": sorted(bolsig_table.rate_coefficients),
            "mean_energy_ev_min": float(np.nanmin(bolsig_table.mean_energy_ev)),
            "mean_energy_ev_max": float(np.nanmax(bolsig_table.mean_energy_ev)),
        },
        "reaction_rates": {k: float(v) for k, v in rates.items()},
        "cross_sections": cross_section_library.summary() if cross_section_library else {},
        "cr_model": cr_data_library.summary() if cr_data_library else {},
        "spectra": spectra_summary,
    }
    (run_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


def _write_plots(run_dir: Path, sample: dict, raster) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fields = [
        ("temperature_k", "Temperature K"),
        ("total_neutral_density_m3", "Neutral density m^-3"),
        ("electron_density_m3", "Electron density m^-3"),
        ("surface_displacement_m", "Surface displacement m"),
        ("jz_a_m2", "Current density Jz A/m^2"),
        ("joule_heating_w_m3", "Joule heating W/m^3"),
        ("conductivity_s_m", "Electrical conductivity S/m"),
        ("bx_t", "Magnetic field Bx T"),
        ("by_t", "Magnetic field By T"),
    ]
    derived = {
        "bmag_t": np.sqrt(sample["bx_t"] * sample["bx_t"] + sample["by_t"] * sample["by_t"]),
    }
    extent_mm = [raster.x_min * 1.0e3, raster.x_max * 1.0e3, raster.y_min * 1.0e3, raster.y_max * 1.0e3]
    for key, title in fields:
        fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
        image = np.log10(np.maximum(sample[key], 1.0e-30)) if key in {"joule_heating_w_m3", "conductivity_s_m"} else sample[key]
        im = ax.imshow(image, origin="lower", aspect="auto", extent=extent_mm)
        ax.set_title(title)
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        fig.colorbar(im, ax=ax)
        fig.savefig(run_dir / f"{key}.png", dpi=160)
        plt.close(fig)
    fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
    im = ax.imshow(derived["bmag_t"], origin="lower", aspect="auto", extent=extent_mm)
    ax.set_title("Magnetic field |B| T")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    fig.colorbar(im, ax=ax)
    fig.savefig(run_dir / "bmag_t.png", dpi=160)
    plt.close(fig)


def _expected_sample_count(n_steps: int, sample_every: int) -> int:
    count = sum(1 for step in range(n_steps) if step % sample_every == 0 or step == n_steps - 1)
    return max(count, 1)


def _write_diagnostics_manifest(
    run_dir: Path,
    hdf5_path: Path,
    completed_samples: int,
    expected_samples: int,
    last_time_s: float | None,
    complete: bool,
) -> None:
    manifest = {
        "format": "hdf5",
        "fields_hdf5": hdf5_path.name,
        "completed_samples": int(completed_samples),
        "expected_samples": int(expected_samples),
        "last_time_s": last_time_s,
        "complete": bool(complete),
    }
    (run_dir / "diagnostics_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def _save_hdf5_as_npz(hdf5_path: Path, npz_path: Path) -> None:
    arrays = load_hdf5_arrays(hdf5_path)
    np.savez_compressed(npz_path, **arrays)


def _progress(progress: Callable[[str], None] | None, message: str) -> None:
    if progress:
        progress(message)
