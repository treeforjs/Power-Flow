# FLASH surface-temperature interchange

`mhdlab.flash_surface` reads and writes a deliberately small, versioned HDF5
contract for a **stationary planar single face**. This is the first coupling
boundary for independently run FLASH heating calculations. It neither reads raw
FLASH AMR output nor invokes the existing reduced MHD/neutral/CR runner.

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
| root `schema` | `mhdlab-flash-surface-v1` | Exact schema identifier |
| root `coordinate_system` | `x=foil_width;y=gap_normal;z=current` | Exact convention |
| root `surface_model` | `stationary_planar_single_face` | Exact geometry restriction |
| root `surface_region` | nonempty string | Explicit face identifier, e.g. cathode gap-facing face |
| root `source_kind` | `flash` or `synthetic` | Source classification; synthetic fixtures stay labeled |
| root `flash_version` | nonempty for `flash` | Installed version used for the actual run |
| root `run_id` | nonempty string | Identifies the external calculation |
| root `source_reference` | nonempty string | Recoverable run manifest/path/hash reference |

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

## Exporter work still required

The workstation exporter must identify the physical face in the installed FLASH
version, convert its temperature units and field conventions explicitly, reject
vacuum contamination, and preserve the face's true width coverage. Its source
manifest should record the FLASH version, configuration, material tables, drive,
time origin, unit conversions, source output checksums and convergence evidence.
Verify the surface against raw local visualization and a known-temperature benchmark.

For shots after material expansion becomes appreciable, add and verify a moving-face
contract with coordinates, normals, areas and velocities. A stationary early-time
temperature export cannot be extended to those gates merely by adding more samples.
The external temperature coupling is useful preparation; the full neutral, reaction,
rovibrational and detector models remain to be implemented and validated.
