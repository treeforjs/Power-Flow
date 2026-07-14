"""Post-process a Power Flow run directory into quick-look plots."""

from __future__ import annotations

import argparse
import io
from pathlib import Path

import imageio.v2 as imageio
import numpy as np


FIELD_LABELS = {
    "temperature_k": "Temperature (K)",
    "total_neutral_density_m3": "Total neutral density (m^-3)",
    "electron_density_m3": "Electron density (m^-3)",
    "surface_displacement_m": "Surface displacement (m)",
    "pressure_pa": "Pressure (Pa)",
    "bx_t": "Bx (T)",
    "by_t": "By (T)",
    "jz_a_m2": "Jz (A/m^2)",
    "en_td": "E/N (Td)",
}

LOG_FIELDS = {"total_neutral_density_m3", "electron_density_m3", "pressure_pa", "en_td"}


def main() -> int:
    parser = argparse.ArgumentParser(description="Create quick-look plots from a run directory.")
    parser.add_argument("run_dir", help="Run directory containing fields.npz")
    parser.add_argument("--output", default=None, help="Output directory; defaults to run_dir/visualization")
    parser.add_argument("--time-ns", type=float, default=None, help="Also write an overview nearest this time")
    parser.add_argument("--gif", action="append", default=[], help="Field key to animate; repeatable")
    parser.add_argument("--all-gifs", action="store_true", help="Animate the standard overview fields")
    parser.add_argument("--fps", type=float, default=8.0, help="GIF frame rate")
    args = parser.parse_args()

    run_dir = Path(args.run_dir).resolve()
    output_dir = Path(args.output).resolve() if args.output else run_dir / "visualization"
    output_dir.mkdir(parents=True, exist_ok=True)

    fields_path = run_dir / "fields.npz"
    if not fields_path.exists():
        raise FileNotFoundError(fields_path)

    with np.load(fields_path) as data:
        arrays = {key: data[key] for key in data.files}

    write_overview(arrays, output_dir / "overview_last.png", index=-1)
    if "total_neutral_density_m3" in arrays:
        neutral = arrays["total_neutral_density_m3"].reshape(arrays["total_neutral_density_m3"].shape[0], -1)
        peak_idx = int(np.argmax(neutral.max(axis=1)))
        write_overview(arrays, output_dir / "overview_peak_neutral.png", index=peak_idx)
    if args.time_ns is not None:
        idx = int(np.argmin(np.abs(arrays["time_s"] * 1.0e9 - args.time_ns)))
        write_overview(arrays, output_dir / f"overview_{arrays['time_s'][idx] * 1.0e9:.1f}ns.png", index=idx)
    write_time_traces(arrays, output_dir / "time_traces.png")

    gif_fields = list(args.gif)
    if args.all_gifs:
        gif_fields.extend(["temperature_k", "total_neutral_density_m3", "electron_density_m3", "surface_displacement_m"])
    for key in dict.fromkeys(gif_fields):
        if key not in arrays:
            raise KeyError(f"field not found in fields.npz: {key}")
        write_gif(arrays, key, output_dir / f"{key}.gif", fps=args.fps)

    print(output_dir)
    return 0


def write_overview(arrays: dict[str, np.ndarray], path: Path, index: int) -> None:
    import matplotlib.pyplot as plt

    keys = ["temperature_k", "total_neutral_density_m3", "electron_density_m3", "surface_displacement_m"]
    time_ns = float(arrays["time_s"][index] * 1.0e9)
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    for ax, key in zip(axes.ravel(), keys):
        image = display_array(arrays[key][index], key)
        im = ax.imshow(image, origin="lower", aspect="auto")
        ax.set_title(f"{FIELD_LABELS.get(key, key)} at {time_ns:.1f} ns")
        ax.set_xlabel("x cell")
        ax.set_ylabel("y cell")
        cbar = fig.colorbar(im, ax=ax)
        if key in LOG_FIELDS:
            cbar.set_label("log10")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_time_traces(arrays: dict[str, np.ndarray], path: Path) -> None:
    import matplotlib.pyplot as plt

    time_ns = arrays["time_s"] * 1.0e9
    traces = [
        ("temperature_k", "max", "Max temperature (K)", 1.0),
        ("total_neutral_density_m3", "sum", "Total neutral inventory (cell-weighted)", 1.0),
        ("electron_density_m3", "max", "Max electron density (m^-3)", 1.0),
        ("surface_displacement_m", "max", "Max surface displacement (um)", 1.0e6),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    for ax, (key, reducer, label, scale) in zip(axes.ravel(), traces):
        values = arrays[key]
        if reducer == "sum":
            y = values.reshape(values.shape[0], -1).sum(axis=1) * scale
        else:
            y = values.reshape(values.shape[0], -1).max(axis=1) * scale
        ax.plot(time_ns, y, lw=2)
        ax.set_xlabel("time (ns)")
        ax.set_ylabel(label)
        ax.grid(True, alpha=0.25)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_gif(arrays: dict[str, np.ndarray], key: str, path: Path, fps: float) -> None:
    import matplotlib.pyplot as plt

    values = arrays[key]
    time_ns = arrays["time_s"] * 1.0e9
    display_values = display_array(values, key)
    finite = display_values[np.isfinite(display_values)]
    vmin, vmax = np.percentile(finite, [1.0, 99.0]) if finite.size else (0.0, 1.0)
    if np.isclose(vmin, vmax):
        vmax = vmin + 1.0

    frames = []
    for idx in range(values.shape[0]):
        fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
        im = ax.imshow(display_values[idx], origin="lower", aspect="auto", vmin=vmin, vmax=vmax)
        ax.set_title(f"{FIELD_LABELS.get(key, key)}  t={time_ns[idx]:.1f} ns")
        ax.set_xlabel("x cell")
        ax.set_ylabel("y cell")
        fig.colorbar(im, ax=ax)
        frames.append(figure_to_rgb(fig))
        plt.close(fig)
    imageio.mimsave(path, frames, duration=1.0 / max(fps, 1.0e-6))


def display_array(values: np.ndarray, key: str) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if key in LOG_FIELDS:
        positive = arr[arr > 0.0]
        floor = max(float(np.nanmin(positive)) * 0.1, 1.0e-300) if positive.size else 1.0e-300
        return np.log10(np.maximum(arr, floor))
    return arr


def figure_to_rgb(fig) -> np.ndarray:
    fig.canvas.draw()
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=130)
    buffer.seek(0)
    return imageio.imread(buffer)


if __name__ == "__main__":
    raise SystemExit(main())
