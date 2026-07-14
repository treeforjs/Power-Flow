"""Build starter electrical-conductivity tables for MHD runs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from mhdlab.materials import material_from_config


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate an ETHOS-style starter conductivity table.")
    parser.add_argument("--output", required=True, help="Output .csv, .npz, .h5, or .hdf5 table")
    parser.add_argument("--material", default="SS304", help="Material preset")
    parser.add_argument("--density-min", type=float, default=4000.0, help="Minimum density in kg/m^3")
    parser.add_argument("--density-max", type=float, default=12000.0, help="Maximum density in kg/m^3")
    parser.add_argument("--n-density", type=int, default=9, help="Number of density grid points")
    parser.add_argument("--temperature-min", type=float, default=100.0, help="Minimum temperature in K")
    parser.add_argument("--temperature-max", type=float, default=20000.0, help="Maximum temperature in K")
    parser.add_argument("--n-temperature", type=int, default=80, help="Number of temperature grid points")
    parser.add_argument("--density-exponent", type=float, default=0.0, help="Knoepfel density exponent")
    parser.add_argument("--betacv-1-k", type=float, default=9.4e-4, help="Knoepfel temperature coefficient")
    parser.add_argument("--sigma-min", type=float, default=1.0e5, help="Conductivity floor in S/m")
    parser.add_argument("--ensemble-count", type=int, default=1, help="Number of uncertainty table members")
    parser.add_argument("--uncertainty-fraction", type=float, default=0.0, help="One-sigma lognormal multiplier uncertainty")
    parser.add_argument("--seed", type=int, default=1234, help="Random seed for ensemble generation")
    args = parser.parse_args()

    material = material_from_config({"preset": args.material})
    density = np.geomspace(args.density_min, args.density_max, args.n_density)
    temperature = np.geomspace(args.temperature_min, args.temperature_max, args.n_temperature)
    base = knoepfel_table(
        density,
        temperature,
        sigma0=1.0 / material.electrical_resistivity_ohm_m,
        rho0=material.density_kg_m3,
        t0=material.initial_temperature_k,
        density_exponent=args.density_exponent,
        betacv_1_k=args.betacv_1_k,
        sigma_min=args.sigma_min,
    )
    values = build_ensemble(base, args.ensemble_count, args.uncertainty_fraction, args.seed)
    write_table(Path(args.output), density, temperature, values)
    print(args.output)
    return 0


def knoepfel_table(
    density: np.ndarray,
    temperature: np.ndarray,
    sigma0: float,
    rho0: float,
    t0: float,
    density_exponent: float,
    betacv_1_k: float,
    sigma_min: float,
) -> np.ndarray:
    rho_ratio = np.maximum(density[:, None] / max(rho0, 1.0e-30), 0.0)
    denom = 1.0 + betacv_1_k * np.maximum(temperature[None, :] - t0, 0.0)
    sigma = sigma0 * np.power(rho_ratio, density_exponent) / np.maximum(denom, 1.0e-30)
    return np.maximum(sigma, sigma_min)


def build_ensemble(base: np.ndarray, count: int, uncertainty_fraction: float, seed: int) -> np.ndarray:
    count = max(int(count), 1)
    uncertainty_fraction = max(float(uncertainty_fraction), 0.0)
    if count == 1:
        return base
    rng = np.random.default_rng(seed)
    sigma_ln = np.log1p(uncertainty_fraction)
    perturb = rng.normal(loc=0.0, scale=sigma_ln, size=(count, *base.shape))
    perturb[0, :, :] = 0.0
    return base[None, :, :] * np.exp(perturb)


def write_table(path: Path, density: np.ndarray, temperature: np.ndarray, conductivity: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    suffix = path.suffix.lower()
    if suffix == ".csv":
        write_csv(path, density, temperature, conductivity)
    elif suffix == ".npz":
        np.savez_compressed(path, density_kg_m3=density, temperature_k=temperature, conductivity_s_m=conductivity)
    elif suffix in {".h5", ".hdf5"}:
        try:
            import h5py
        except ImportError as exc:
            raise RuntimeError("h5py is required for HDF5 conductivity table output") from exc
        with h5py.File(path, "w") as h5:
            h5.create_dataset("density_kg_m3", data=density)
            h5.create_dataset("temperature_k", data=temperature)
            h5.create_dataset("conductivity_s_m", data=conductivity)
            h5.attrs["schema"] = "mhdlab-conductivity-table-v1"
            h5.attrs["interpolation"] = "log-density/log-temperature/log-conductivity"
    else:
        raise ValueError(f"unsupported output format: {path}")


def write_csv(path: Path, density: np.ndarray, temperature: np.ndarray, conductivity: np.ndarray) -> None:
    values = conductivity[None, :, :] if conductivity.ndim == 2 else conductivity
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["ensemble_index", "density_kg_m3", "temperature_k", "conductivity_s_m"])
        for member in range(values.shape[0]):
            for i, rho in enumerate(density):
                for j, temp in enumerate(temperature):
                    writer.writerow([member, f"{rho:.12g}", f"{temp:.12g}", f"{values[member, i, j]:.12g}"])


if __name__ == "__main__":
    raise SystemExit(main())
