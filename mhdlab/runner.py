"""High-level simulation orchestration."""

from __future__ import annotations

import json
import shutil
from datetime import datetime
from pathlib import Path

import numpy as np

from .bolsig import BolsigRunner, approximate_bolsig_table, supplement_missing_rate_coefficients
from .chemistry import Mechanism, line_emissivity, local_reaction_rates
from .config import resolve_path
from .cr_data import CRDataLibrary
from .cross_sections import CrossSectionLibrary
from .desorption import TemkinDesorption
from .geometry import Geometry
from .materials import material_from_config
from .mhd import ReducedMHDSolver
from .neutrals import KineticNeutralSolver, VelocityGrid
from .spectra import CCDCalibration, integrate_los, synthesize_spectrum
from .traces import fit_effective_rl, load_trace_csv, synthesize_current_from_rl


def run_from_config(config: dict) -> dict:
    root = resolve_path(config, config.get("project_root", ".")) or Path.cwd()
    out_root = resolve_path(config, config.get("output_root", str(root / "runs"))) or (root / "runs")
    run_dir = out_root / datetime.now().strftime("run_%Y%m%d_%H%M%S")
    run_dir.mkdir(parents=True, exist_ok=True)

    config_path = config.get("_config_path")
    if config_path:
        shutil.copy2(config_path, run_dir / "config_used.yaml")

    geometry = Geometry.from_json(resolve_path(config, config["geometry"]["file"]))
    grid_cfg = config.get("grid", {})
    raster = geometry.rasterize(nx=int(grid_cfg.get("nx", 64)), ny=int(grid_cfg.get("ny", 48)))
    material = material_from_config(config.get("material"))
    mhd_solver = ReducedMHDSolver(raster, material, config.get("mhd", {}))
    mhd_state = mhd_solver.initial_state()

    trace = load_trace_csv(resolve_path(config, config["drive"]["trace_csv"]))
    voltage_name = config["drive"].get("voltage_column", "voltage_v")
    current_name = config["drive"].get("current_column", "feed_current_a")
    rl_fit = fit_effective_rl(trace.time_s, trace.values[voltage_name], trace.values[current_name])
    drive_current = trace.values[current_name]
    if config["drive"].get("mode", "fit_rl") == "fit_rl":
        drive_current = synthesize_current_from_rl(
            trace.time_s,
            trace.values[voltage_name],
            rl_fit["resistance_ohm"],
            rl_fit["inductance_h"],
        )

    mechanism = Mechanism.from_yaml(resolve_path(config, config["chemistry"]["mechanism"]))
    species = mechanism.species or ["H2O", "H", "H2", "O", "O2", "OH"]
    vcfg = config.get("velocity_grid", {})
    velocity_grid = VelocityGrid.polar(
        max_speed_m_s=float(vcfg.get("max_speed_m_s", 8000.0)),
        n_speed=int(vcfg.get("n_speed", 4)),
        n_angle=int(vcfg.get("n_angle", 12)),
    )
    neutral_solver = KineticNeutralSolver(raster, species, velocity_grid, config.get("neutrals", {}))
    neutral_state = neutral_solver.initial_state()

    desorption = TemkinDesorption(**config.get("desorption", {}))
    coverage = desorption.initial_coverage(raster.shape)

    bolsig_table = _load_bolsig_table(config, run_dir)
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
    sample_every = int(config["time"].get("sample_every", 1))
    times = np.arange(0.0, t_end + 0.5 * dt, dt)
    samples = []

    for step, t in enumerate(times):
        voltage = trace.sample(voltage_name, t)
        current = float(np.interp(t, trace.time_s, drive_current))
        mhd_state = mhd_solver.step(mhd_state, current_a=current, voltage_v=voltage, dt_s=dt)
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
        )
        if step % sample_every == 0 or step == len(times) - 1:
            moments = neutral_solver.moments(neutral_state)
            total_neutral = neutral_state.total_neutral_density()
            en_td = _reduced_field_td(mhd_state.ex_v_m, mhd_state.ey_v_m, total_neutral)
            samples.append(
                {
                    "time_s": t,
                    "temperature_k": mhd_state.temperature_k.copy(),
                    "pressure_pa": mhd_state.pressure_pa.copy(),
                    "bx_t": mhd_state.bx_t.copy(),
                    "by_t": mhd_state.by_t.copy(),
                    "jz_a_m2": mhd_state.jz_a_m2.copy(),
                    "surface_displacement_m": mhd_state.surface_displacement_m.copy(),
                    "electron_density_m3": neutral_state.electron_density_m3.copy(),
                    "total_neutral_density_m3": total_neutral.copy(),
                    "en_td": en_td.copy(),
                    "species_density_m3": {sp: data["density_m3"].copy() for sp, data in moments.items()},
                }
            )

    result_path = run_dir / "fields.npz"
    _save_samples(result_path, samples)
    spectra_summary = _write_spectra_products(config, run_dir, raster, mechanism, samples[-1])
    _write_summary(
        run_dir,
        rl_fit,
        bolsig_table,
        base_rates,
        spectra_summary,
        cross_section_library,
        cr_data_library,
    )
    _write_plots(run_dir, samples[-1])
    return {"run_dir": str(run_dir), "samples": samples, "rl_fit": rl_fit, "spectra": spectra_summary}


def _load_bolsig_table(config: dict, run_dir: Path):
    bcfg = config.get("bolsig", {})
    if not bcfg.get("enabled", True):
        return approximate_bolsig_table()
    bolsig_dir = resolve_path(config, bcfg.get("path", "third_party/bolsigplus"))
    runner = BolsigRunner(bolsig_dir, run_dir / "bolsig_cache")
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
        return table
    except Exception as exc:
        (run_dir / "bolsig_warning.txt").write_text(f"Falling back to approximate table: {exc}\n", encoding="utf-8")
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


def _save_samples(path: Path, samples: list[dict]) -> None:
    arrays = {"time_s": np.asarray([s["time_s"] for s in samples])}
    for key in ("temperature_k", "pressure_pa", "bx_t", "by_t", "jz_a_m2", "surface_displacement_m", "electron_density_m3", "total_neutral_density_m3", "en_td"):
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


def _write_plots(run_dir: Path, sample: dict) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fields = [
        ("temperature_k", "Temperature K"),
        ("total_neutral_density_m3", "Neutral density m^-3"),
        ("electron_density_m3", "Electron density m^-3"),
        ("surface_displacement_m", "Surface displacement m"),
    ]
    for key, title in fields:
        fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
        im = ax.imshow(sample[key], origin="lower", aspect="auto")
        ax.set_title(title)
        fig.colorbar(im, ax=ax)
        fig.savefig(run_dir / f"{key}.png", dpi=160)
        plt.close(fig)
