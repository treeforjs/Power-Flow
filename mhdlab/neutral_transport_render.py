"""Render saved neutral-transport maps without evolving or altering populations."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import platform

import numpy as np


SPECIES = ("H", "H2", "H2O", "OH")
REQUIRED_ARRAYS = ("time_s", "x_edges_m", "y_edges_m", "species",
                   "number_density_m3", "foil_bounds_m")


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _finite_array(value, name, ndim):
    raw = np.asarray(value)
    if raw.ndim != ndim or not raw.size or raw.dtype.kind not in "iuf":
        raise ValueError(f"{name} must be a nonempty real rank-{ndim} array")
    value = np.asarray(raw, dtype=float)
    if not np.isfinite(value).all():
        raise ValueError(f"{name} must contain finite values")
    return value


def _text(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be an explicit nonempty string")
    return value


def _metadata(path):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError(f"duplicate JSON key: {key}")
            result[key] = value
        return result

    def constant(value):
        raise ValueError(f"nonfinite JSON constant: {value}")

    value = json.loads(Path(path).read_bytes(), object_pairs_hook=pairs, parse_constant=constant)
    if value.get("status") not in {"manufactured", "provisional", "owner_supplied"}:
        raise ValueError("status must be manufactured, provisional or owner_supplied")
    for name in ("time_registration", "model_description"):
        _text(value.get(name), name)
    length = value.get("invariant_current_length_m")
    if isinstance(length, bool) or not isinstance(length, (int, float)) or not math.isfinite(length) or length <= 0:
        raise ValueError("invariant_current_length_m must be explicit, finite and positive")
    assumptions = value.get("assumptions")
    if not isinstance(assumptions, list) or not assumptions:
        raise ValueError("assumptions must be a nonempty list of explicit strings")
    for assumption in assumptions:
        _text(assumption, "assumption")
    return value


def load_transport_maps(path, metadata_path):
    """Validate the shared saved-map contract; no interpolation or repair."""
    metadata = _metadata(metadata_path)
    with np.load(path, allow_pickle=False) as archive:
        missing = set(REQUIRED_ARRAYS) - set(archive.files)
        if missing:
            raise ValueError(f"missing transport arrays: {sorted(missing)}")
        arrays = {name: np.array(archive[name], copy=True) for name in REQUIRED_ARRAYS}
    species = arrays["species"]
    if species.ndim != 1 or species.dtype.kind not in "US" or tuple(species.astype(str)) != SPECIES:
        raise ValueError("species must be the exact ordered string vector H,H2,H2O,OH")
    time = _finite_array(arrays["time_s"], "time_s", 1)
    x = _finite_array(arrays["x_edges_m"], "x_edges_m", 1)
    y = _finite_array(arrays["y_edges_m"], "y_edges_m", 1)
    if len(time) < 2 or np.any(np.diff(time) <= 0):
        raise ValueError("time_s must have at least two strictly increasing samples")
    for name, edges in (("x_edges_m", x), ("y_edges_m", y)):
        if len(edges) < 2 or np.any(np.diff(edges) <= 0) or not math.isfinite(float(edges[-1]-edges[0])):
            raise ValueError(f"{name} must be a finite strictly increasing partition")
    density = _finite_array(arrays["number_density_m3"], "number_density_m3", 4)
    if density.shape != (len(time), 4, len(y)-1, len(x)-1) or np.any(density < 0):
        raise ValueError("density must be nonnegative with shape(time,4,normal_bin,width_bin)")
    foil = _finite_array(arrays["foil_bounds_m"], "foil_bounds_m", 1)
    if foil.shape != (4,) or foil[1] <= foil[0] or foil[3] <= foil[2]:
        raise ValueError("foil_bounds_m requires xmin<xmax and ymin<ymax")
    if foil[0] < x[0] or foil[1] > x[-1] or foil[2] < y[0] or foil[3] > y[-1]:
        raise ValueError("the full foil rectangle must lie inside the displayed map")
    arrays.update(time_s=time, x_edges_m=x, y_edges_m=y, number_density_m3=density, foil_bounds_m=foil)
    return arrays, metadata


def log_display_limits(density, color_span_decades=None):
    """Fix one global log scale per species; zero never acquires a density floor."""
    result = []
    for species_index, name in enumerate(SPECIES):
        positive = density[:, species_index][density[:, species_index] > 0]
        if not positive.size:
            result.append({"species": name, "positive_min_m3": None, "positive_max_m3": None,
                           "log10_color_limits": None, "normalization": "no_positive_density"})
            continue
        low, high = float(np.min(positive)), float(np.max(positive))
        lo_log, hi_log = math.log10(low), math.log10(high)
        constant = lo_log == hi_log
        display_low = max(lo_log, hi_log-color_span_decades) if color_span_decades is not None else lo_log
        result.append({"species": name, "positive_min_m3": low, "positive_max_m3": high,
                       "log10_color_limits": [lo_log-.5, hi_log+.5] if constant else [display_low, hi_log],
                       "normalization": "fixed_per_species_log10",
                       "positive_bins_below_color_range": int(np.count_nonzero(np.log10(positive) < display_low)),
                       "single_positive_value_color_axis_expanded_half_dex": constant})
    return result


def render_transport(npz_path, metadata_path, output_dir, *, fps=10., dpi=110, color_span_decades=None,
                     view_bounds_m=None):
    """Write a four-panel animated GIF plus actual decoded-GIF previews.

    Saved density values and times are left untouched. Every source snapshot is
    rendered once with equal frame dwell; physical timestamps remain explicit.
    """
    from PIL import Image, __version__ as pillow_version
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.patches import Rectangle
    from matplotlib.ticker import FuncFormatter, MaxNLocator

    npz_path, metadata_path, output = Path(npz_path).resolve(), Path(metadata_path).resolve(), Path(output_dir).resolve()
    if output.exists():
        raise ValueError("use a new render output directory")
    if isinstance(fps, bool) or not isinstance(fps, (int, float)) or not math.isfinite(fps) or not 0 < fps <= 100:
        raise ValueError("fps must be finite and in (0,100]")
    if type(dpi) is not int or not 60 <= dpi <= 240:
        raise ValueError("dpi must be an integer from 60 to 240")
    if color_span_decades is not None and (isinstance(color_span_decades, bool)
            or not isinstance(color_span_decades, (int, float)) or not math.isfinite(color_span_decades)
            or not 1 <= color_span_decades <= 300):
        raise ValueError("color_span_decades must be finite and in [1,300]")
    input_hashes = {"npz": _sha(npz_path), "metadata": _sha(metadata_path)}
    arrays, metadata = load_transport_maps(npz_path, metadata_path)
    full_bounds = np.array([arrays["x_edges_m"][0], arrays["x_edges_m"][-1],
                            arrays["y_edges_m"][0], arrays["y_edges_m"][-1]])
    view = full_bounds.copy() if view_bounds_m is None else _finite_array(view_bounds_m, "view_bounds_m", 1)
    if (view.shape != (4,) or view[0] >= view[1] or view[2] >= view[3]
            or view[0] < full_bounds[0] or view[1] > full_bounds[1]
            or view[2] < full_bounds[2] or view[3] > full_bounds[3]):
        raise ValueError("view_bounds_m must be ordered xmin,xmax,ymin,ymax within the full saved map")
    cropped_view = not np.array_equal(view, full_bounds)
    view_mm = view * 1000
    normal_magnification = metadata.get("normal_axis_magnification", 1.0)
    if (isinstance(normal_magnification, bool) or not isinstance(normal_magnification, (int, float))
            or not math.isfinite(normal_magnification) or not 1 <= normal_magnification <= 20):
        raise ValueError("normal_axis_magnification must be a declared finite scale in [1,20]")
    density = arrays["number_density_m3"]
    times = arrays["time_s"]
    areas = np.outer(np.diff(arrays["y_edges_m"]), np.diff(arrays["x_edges_m"]))
    with np.errstate(over="ignore", invalid="ignore"):
        counts = np.einsum("tsij,ij->ts", density, areas) * metadata["invariant_current_length_m"]
    if not np.isfinite(counts).all() or not np.isfinite(np.sum(counts, axis=1)).all():
        raise ValueError("integrated map counts exceed the finite numerical domain")
    x, y, foil = (arrays[name] * 1000 for name in ("x_edges_m", "y_edges_m", "foil_bounds_m"))
    limits = log_display_limits(density, color_span_decades)
    # Preserve SI arrays/count integration; convert display logarithms by -6 dex
    # so even subnormal positive native densities remain visible.
    for scale in limits:
        for key in ("positive_min", "positive_max"):
            value = scale.pop(key + "_m3")
            scale[key + "_cm3"] = None if value is None else value / 1e6
        if scale["log10_color_limits"] is not None:
            scale["log10_color_limits"] = [v-6 for v in scale["log10_color_limits"]]
    cmap = plt.get_cmap("YlGnBu").with_extremes(bad="white")
    fig, axes = plt.subplots(2, 2, figsize=(10.4, 8.6), layout="constrained")
    fig.set_dpi(dpi)
    title = metadata.get("case_title", "Neutral desorption and vacuum transport")
    status = metadata["status"].upper().replace("_", " ")
    view_description = "cropped physical view" if cropped_view else "foil width–normal plane"
    fig.suptitle(f"{title}\n{status} — number density in the {view_description}", fontsize=13)
    # A fixed palette prevents frame-dependent GIF quantization from changing colors.
    palette_image = Image.new("P", (1, 1))
    colors = [tuple(round(255 * channel) for channel in cmap(i / 223)[:3]) for i in range(224)]
    colors.extend([(round(255*i/31),) * 3 for i in range(32)])
    palette_image.putpalette([channel for rgb in colors for channel in rgb])
    artists = []
    for species_index, (ax, scale) in enumerate(zip(axes.flat, limits)):
        color_limits = scale["log10_color_limits"]
        field = np.ma.array(np.zeros(density.shape[2:], dtype=float), mask=True)
        normalizer = Normalize(*color_limits) if color_limits is not None else Normalize(0, 1)
        artist = ax.pcolormesh(x, y, field, shading="flat", cmap=cmap, norm=normalizer, rasterized=True)
        artists.append(artist)
        ax.add_patch(Rectangle((foil[0], foil[2]), foil[1]-foil[0], foil[3]-foil[2],
                               facecolor="0.35", edgecolor="black", linewidth=.8, zorder=5))
        ax.set(xlabel="Foil width x [mm]", ylabel="Thickness / vacuum normal y [mm]",
               xlim=tuple(view_mm[:2]), ylim=tuple(view_mm[2:]), title=SPECIES[species_index])
        ax.set_aspect(normal_magnification, adjustable="box")
        ax.xaxis.set_major_locator(MaxNLocator(5))
        ax.yaxis.set_major_locator(MaxNLocator(5))
        if color_limits is not None:
            colorbar = fig.colorbar(artist, ax=ax, fraction=.046, pad=.035,
                                   extend="min" if scale["positive_bins_below_color_range"] else "neither")
            colorbar.set_label(r"Number density [cm$^{-3}$]", fontsize=9)
            colorbar.locator = MaxNLocator(5)
            colorbar.formatter = FuncFormatter(lambda value, _: r"$10^{%g}$" % value)
            colorbar.update_ticks()
        else:
            ax.text(.5, .94, "No positive density in this species", transform=ax.transAxes,
                    ha="center", va="top", fontsize=9, bbox={"facecolor": "white", "edgecolor": "0.7"})
    maximum_time = float(np.max(np.abs(times)))
    unit, multiplier = ("ns", 1e9) if maximum_time < 1e-6 else ("µs", 1e6) if maximum_time < 1e-3 else ("ms", 1e3) if maximum_time < 1 else ("s", 1.)
    footer = fig.supxlabel("", fontsize=10)
    frames = []
    for frame, stamp in enumerate(times):
        for species_index, artist in enumerate(artists):
            values = density[frame, species_index]
            positive = values > 0
            logged = np.zeros_like(values)
            np.log10(values, out=logged, where=positive)
            logged[positive] -= 6
            display = np.ma.array(logged, mask=~positive)
            artist.set_array(display)
        aspect_note = ("Equal length scale" if normal_magnification == 1
                       else f"Normal axis enlarged {normal_magnification:g}×")
        color_note = ("No time interpolation" if color_span_decades is None
                      else f"Color span ≤{color_span_decades:g} decades; lower positive values pale")
        footer.set_text(f"t = {stamp*multiplier:.5g} {unit}    |    frame {frame+1}/{len(times)}\n"
                        f"Gray: solid foil   •   White: zero density   •   {aspect_note}\n{color_note}")
        fig.canvas.draw()
        rgb = np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()
        image = Image.fromarray(rgb).quantize(palette=palette_image, dither=Image.Dither.NONE)
        frames.append(image)
    plt.close(fig)
    dwell_ms = max(10, int(round(100 / fps)) * 10)
    output.mkdir(parents=True)
    gif_path = output / "neutral_transport.gif"
    frames[0].save(gif_path, save_all=True, append_images=frames[1:], duration=dwell_ms,
                   loop=0, disposal=2, optimize=False)
    peak_index = int(np.argmax(np.sum(counts, axis=1)))
    with Image.open(gif_path) as decoded:
        frame_count = decoded.n_frames
        if frame_count != len(times):
            raise ValueError("encoded GIF did not retain every source snapshot")
        dimensions = list(decoded.size)
        durations = []
        for index in range(frame_count):
            decoded.seek(index)
            durations.append(decoded.info.get("duration", 0))
            if index in {0, peak_index, frame_count-1}:
                for filename, target in (("preview_first.png", 0), ("preview_peak.png", peak_index),
                                          ("preview_last.png", frame_count-1)):
                    if index == target:
                        decoded.convert("RGB").save(output / filename)
    if input_hashes != {"npz": _sha(npz_path), "metadata": _sha(metadata_path)}:
        raise ValueError("source maps changed during rendering")
    report = {"schema": "mhdlab-neutral-transport-gif-v1", "status": metadata["status"],
        "input_npz": str(npz_path), "input_metadata": str(metadata_path), "input_sha256": input_hashes,
        "source_metadata": metadata, "frame_count": frame_count, "dimensions_pixels": dimensions,
        "requested_fps": fps, "frame_duration_ms": durations, "total_animation_duration_ms": sum(durations),
        "physical_time_s": times.tolist(), "display_color_scales": limits,
        "palette": "fixed_224_YlGnBu_colors_plus32_greys; no dithering",
        "display_density_units": "cm^-3", "native_archive_density_units": "m^-3",
        "display_conversion_factor_from_m3": 1e-6,
        "density_values_changed": False, "zeros_replaced_by_floor": False,
        "time_interpolation": False,
        "geometry_aspect": "equal_millimetre_units" if normal_magnification == 1 else "declared_normal_axis_magnification",
        "normal_axis_magnification": normal_magnification,
        "full_map_bounds_m": full_bounds.tolist(), "display_view_bounds_m": view.tolist(),
        "view_crops_axes_only": cropped_view,
        "requested_color_span_decades": color_span_decades,
        "positive_under_range_display": "lowest colormap color; values retained and included in particle totals",
        "foil_display": "solid rectangle overlay; source density arrays unchanged",
        "peak_number_frame_index": peak_index, "integrated_plot_particle_numbers": counts.tolist(),
        "integrated_particle_number_scope": "whole saved map, including any area outside the display view",
        "physical_model_qualified_by_renderer": False,
        "runtime": {"python": platform.python_version(), "numpy": np.__version__,
                    "matplotlib": matplotlib.__version__, "pillow": pillow_version},
        "renderer_source_sha256": _sha(__file__),
        "output_sha256": {path.name: _sha(path) for path in sorted(output.iterdir()) if path.is_file()}}
    (output / "render_report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return report
