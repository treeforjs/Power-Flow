# Saved neutral-transport GIF

`tools/render_neutral_transport_gif.py` renders saved neutral number-density maps
without changing populations, evolving transport, adding chemistry or assigning
physical validity. It displays H, H2, H2O and OH in the foil-width / thickness-
normal plane, including surrounding vacuum. All coordinates remain physical
millimetres. The default aspect ratio is equal; optional sidecar
`normal_axis_magnification` in [1,20] enlarges the normal display direction and
is printed in the figure footer and render report. It changes no coordinates,
density values or timestamps. The solid foil is a gray overlay.

```powershell
python tools/render_neutral_transport_gif.py --input runs/case/transport_maps.npz --metadata runs/case/metadata.json --output-dir runs/case/new_visualization
```

The output directory must be new. It contains `neutral_transport.gif`, previews
decoded from the actual GIF's first, maximum in-domain number and last frames,
and `render_report.json`. The report records input/source/output hashes, original
timestamps, frame count and durations, fixed color scales and source metadata.

Required NPZ arrays:

| Name | Meaning |
| --- | --- |
| `time_s` | At least two finite strictly increasing physical snapshot times |
| `x_edges_m` | Increasing foil-width grid edges |
| `y_edges_m` | Increasing thickness/vacuum-normal grid edges |
| `species` | Exact ordered string array `H,H2,H2O,OH` |
| `number_density_m3` | Finite nonnegative array `(time,4,normal_bin,width_bin)` |
| `foil_bounds_m` | Solid rectangle `[xmin,xmax,ymin,ymax]` inside displayed support |

The JSON sidecar requires `status` (`provisional`, `manufactured` or
`owner_supplied`), `time_registration`, `model_description`, positive
`invariant_current_length_m`, and a nonempty string-list `assumptions`. Other
source/geometry/chemistry/ledger records are retained unchanged. Optional
`case_title` supplies the figure title. A rendering status is not scientific
qualification.

Each species has one fixed log10 color range over every input snapshot. Exact
zeros are masked white; there is no density floor or artificial population.
An entirely empty species receives an explicit empty-panel annotation without a
color scale. A single positive value gets a documented half-decade display-axis
bracket; the density array remains unchanged. The GIF uses a fixed palette across
frames. Every input time is rendered once with equal frame dwell and its own
physical timestamp; no unobserved time or transport is interpolated. GIF dwell
times are quantized to its 10-ms format precision and recorded.

Optional `--color-span-decades 8` restricts each fixed color scale to at most
eight decades below that species' global maximum. Smaller positive densities
remain visible in the lowest color, with a colorbar under-range extension;
they are not set to zero, floored, or removed from particle totals. The footer
declares the span, and the report preserves actual extrema and counts of
positive bins below the color range. This avoids negligible depleted-pool tails
spanning hundreds of decades and washing out the resolved plume contrast.

Optional `--view-bounds-mm XMIN XMAX YMIN YMAX` crops only the displayed axes
to a finite, ordered physical rectangle inside the saved map. The API argument
is `view_bounds_m` in metres. This supports an equal-scale edge close-up while
retaining the full density arrays, original times, and global per-species color
scales. It does not resample, magnify one coordinate relative to another, change
velocities or evolve additional flight. The title labels a cropped view, and
the report records both full-map and displayed bounds. Its integrated particle
numbers always cover the whole saved map, including any part outside the view.

The renderer verifies that the encoded GIF retains every input snapshot and
provides integrated number diagnostics using the declared invariant current
length. These are diagnostics of the supplied arrays, not source/chemistry or
conservation verification. A successful render does not establish an MFEM foil
calculation, water fragment yields, a neutral temperature, or agreement with the
dissertation.
