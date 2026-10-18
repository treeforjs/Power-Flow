# FLASH surface-temperature interchange

`mhdlab.flash_surface` reads and writes a versioned HDF5 contract for a
**stationary planar single face**. `mhdlab.flash_extract` now supplies a configured
exporter for supported structured FLASH HDF5 scalar layouts. Neither module
launches FLASH or invokes the legacy reduced runner.

Use one file per emitting face. The positive x coordinate runs across the full
foil width. y is the gap normal and z the current direction. The file describes
surface temperature averaged over each x bin, per unit length in z. Distinct face
identifiers prevent anode and cathode from being silently combined. The contract
does not encode a face normal/location, moving surface, velocity, electromagnetic
fields, or a volume state; those require additional interfaces before transport or
fully coupled discharge predictions are possible.

## Layout

| HDF5 item | Shape/value | Meaning |
| --- | --- | --- |
| `/time_s` | `(nt,)`, units attribute `s` | At least two finite nonnegative times, strictly increasing |
| `/x_edges_m` | `(nx+1,)`, units attribute `m` | Finite strictly increasing full-width bin edges |
| `/temperature_k` | `(nt,nx)`, units attribute `K` | Positive finite material surface/lattice temperature |
| root `schema` | `mhdlab-flash-surface-v2` | Exact schema identifier; v1 lacks required temperature-association metadata |
| root `coordinate_system` | `x=foil_width;y=gap_normal;z=current` | Exact convention |
| root `surface_model` | `stationary_planar_single_face` | Exact geometry restriction |
| root `surface_region` | nonempty string | Explicit face identifier, e.g. cathode gap-facing face |
| root `source_kind` | `flash` or `synthetic` | Source classification; synthetic fixtures stay labeled |
| root `flash_version` | nonempty for `flash` | Installed version used for the actual run |
| root `run_id` | nonempty string | Identifies the external calculation |
| root `source_reference` | nonempty string | Recoverable run manifest/path/hash reference |
| root `temperature_location` | `surface_boundary`, `material_sample_plane`, or `boundary_temperature_proxy` | Explicit relationship to the physical boundary |
| root `material_sample_depth_m` | finite nonnegative number | Positive for a sampled material plane or boundary proxy; zero for a boundary value |
| root `material_sample_depth_units` | `m` | Exact depth units |
| root `boundary_temperature_reference` | string; nonempty for a proxy | Evidence supporting the boundary-temperature approximation |
| root `current_averaging_interval_m_json` | JSON `null` or increasing pair | Physical interval averaged along current; null for an invariant-direction history |
| root `current_uniformity_reference` | string | Required before using a current-averaged temperature in nonlinear release |

Choose and record the simulation time origin in the source run manifest, with
nonnegative stored simulation times. Transform diagnostic trigger offsets explicitly;
do not silently align each spectrum to an arbitrary temperature feature.

The `FlashSurfaceHistory` constructor validates arrays/provenance;
`write_surface_history(path, history)` writes the canonical format and
`load_surface_history(path)` validates it on read. Dataset units are mandatory.
There is no guessed CGS conversion, guessed field name, or temperature-floor repair.
Metadata classifies provenance; it cannot authenticate an external simulation.

## Inspection

From the repository root with normal project dependencies installed:

```bash
python tools/inspect_flash_surface.py /path/to/cathode_surface.h5
python tools/inspect_flash_surface.py /path/to/cathode_surface.h5 --time-s 4e-8
```

The JSON report explicitly includes the source classification. It reports left
and right edge temperatures and a center temperature, with exact bin/region overlap
weights even when bins have unequal widths. Defaults are the outer 10% on each
side and the centered 20%; options can change these fractions but cannot overlap
the regions. Temperature is piecewise constant spatially within a bin, and linear
interpolation is available in time with extrapolation rejected.

These are **thermal metrics**, not predicted desorption ratios. The ratio of
mean temperatures is not a ratio of release rates. Width binning and temporal
sampling must resolve temperature variations before a nonlinear desorption law
is evaluated; applying that law to an underresolved bin average can bias release.

## Configured raw-output exporter

```bash
python tools/export_flash_surface.py --demo
python tools/export_flash_surface.py --config /path/to/extraction.json \
  --output-dir runs/flash_surface_export
```

The demonstration creates manufactured leaf blocks and uses the same extraction
path as external data. It does not represent a FLASH run. The
[`flash_extract_example.json`](../configs/flash_extract_example.json) template
requires explicit configuration before execution.

The configuration specifies dataset paths, storage/bounding-box axes, physical
axis mapping and signs, unit conversions, leaf/parent values, material purity,
times, width bins, material-side sampling depth and current-direction interval.
The exporter handles supported Cartesian cell-centred structured layouts with
no guard cells. It rejects unsupported layouts, vacuum/mixed selected cells,
missing coverage and overlaps. Tolerance does not excuse unresolved topology.
Raw files, configuration and final output are hashed in the extraction manifest.

The result remains a `material_sample_plane` diagnostic. A positive-depth cell
sample is not automatically a boundary value. The exporter does not offer a
silent conversion; the downstream source workflow requires an explicit justified
proxy decision. Three-dimensional current-direction averaging is also preserved.
Evaluate local nonlinear release on resolved patches before area averaging when
that direction is nonuniform.

## Workstation and physical verification still required

The workstation configuration must identify the physical face in the installed
FLASH version, supply the appropriate material/lattice temperature rather than
an electron temperature, verify conversions and field conventions, and preserve
the face's true width coverage. Its source
manifest should record the FLASH version, configuration, material tables, drive,
time origin, unit conversions, source output checksums and convergence evidence.
Verify the surface against raw local visualization and a known-temperature benchmark.

For shots after material expansion becomes appreciable, add and verify a moving-face
contract with coordinates, normals, areas and velocities. A stationary early-time
temperature export cannot be extended to those gates merely by adding more samples.
The external temperature coupling is useful preparation; the full neutral, reaction,
rovibrational and detector models remain to be implemented and validated.
