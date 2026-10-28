# MFEM Joule surface-temperature handoff

`mhdlab.mfem_surface` accepts a deliberately explicit SI sample table from a
native MFEM surface evaluator. It produces the existing resolved planar surface
HDF5 contract used by conservative fixed-inventory and Temkin release. MFEM data
retain `source_kind=mfem`, `heating_backend=mfem_joule`, the MFEM version and full
source revision, and original source schema `mhdlab-mfem-surface-v1`. No MFEM
calculation is labeled as FLASH, and historical FLASH readers remain unchanged.

The stock `joule.cpp` rod/coil examples are not planar foil calculations. Their
mesh attributes, example material constants, drive, temperature zero and units
need independent configuration. Exporting stock native grid functions directly
through this interface is unsupported. It neither parses a VisIt collection nor
samples arbitrary MFEM meshes; the native producer must supply the table and its
hash-bound run, mesh, field and quadrature evidence.

## Required table and configuration

[mfem_surface_example.json](../configs/mfem_surface_example.json) is a
non-executable template. Replace its model, geometry, clock and sample-table
fields and explicitly set `configuration_status=configured` for an MFEM source.
The alternative manufactured status requires `source_kind=synthetic`. Placeholder
mapping and zero-convention references reject before export.

The UTF-8 CSV columns, in order, are:

```text
time_s,z_index,x_index,temperature_value
```

Times already use seconds, the supplied face partitions use metres, and
temperature values use kelvin or kelvin increments under an explicitly declared
zero convention. No guessed conversion from a toy model's dimensionless numbers
is performed. Every configured time/current/width patch must appear exactly once.
Duplicate, missing, out-of-partition, nonfinite or unknown-time rows reject.
Rows may arrive in any order; integer patch IDs supply their exact association.
No regridding, temporal interpolation or spatial temperature averaging occurs.

For `value_basis=absolute_kelvin`, `initial_temperature_k` must be zero and the
table values remain unchanged. For `temperature_increment_kelvin`, the caller
must supply a positive, uniform initial material temperature and a justified zero
convention. That offset is added once, and its value and reference remain in the
surface metadata. This accommodates a qualified delta-temperature formulation;
it does not establish that a stock `T=0` toy field has physical kelvin units.
All resulting absolute temperatures must be positive and finite.

## Material side, face and unresolved gradients

The legacy broad face lies at a declared y coordinate with an outward normal sign,
material boundary attribute and selection/complete-width evidence. The surface
uses x across foil width and z along current. Both x and z partitions retain
resolved local patches. A whole-current averaged temperature cannot be inflated
into resolved current patches. Distinct emitting faces require distinct files.

Thickness-side faces instead declare
`coordinate_system="u=foil_thickness;v=foil_width;z=current"`,
`physical_in_plane_axis="y"` and `physical_normal_axis="x"`. Their local
`x_edges_m` partition physical y, while `normal_position_m` is physical x and
the outward sign means left/right. The shared validator checks this frame
through direct construction, HDF reload and release. Packet coordinates remain
explicitly face-local until `surface_coordinates.to_physical_coordinates`
permutes both position and velocity into physical xyz, preserving speed and
kinetic energy. The [four-face workflow](neutral_perimeter_release.md) applies
that mapping and accounts for the smaller side areas without an extra factor.

MFEM Joule temperature lies in a discontinuous L2 finite-element space. A true
surface export must evaluate its **one-sided material trace** on the selected
boundary face and integrate that trace over each physical x/z patch. Neighboring
air/vacuum element values and arbitrary volume-element means cannot silently
become the material boundary temperature. The configuration identifies this
evaluation as `one_sided_material_boundary_trace_patch_average`.

A positive-depth material sample instead declares
`resolved_material_sample_plane_patch_average` and `location=material_sample_plane`.
It may be inspected, but existing release functions reject it until a separately
justified boundary proxy is recorded. The proxy preserves its positive sampling
depth and evidence reference. Neither this depth, foil thickness nor the thermal
skin depth supplies an optical path or a neutral temperature.

The generated `PlanarSurfaceHistory` retains all metadata through HDF5 and
desorption, including source sample/configuration hashes. Surface area enters
the existing particle/energy ledgers once from `dx*dz`; no extra current length,
depth, gate duration or density multiplies the temperature or particle count.
Refine local x/z/time sampling before applying nonlinear release to patch means.
The shared planar validator applies the MFEM trace, zero-basis, face, source
identity and solver-scope checks during direct construction, HDF5 reload and
the internal Temkin batch path. A saved metadata edit cannot turn a volume mean
into an admitted material-boundary trace. Metadata and hashes remain declarations
whose physical origin must be reviewed independently.

## Export and manufactured verification

```powershell
python tools/export_mfem_surface.py --config path/to/configuration.json --output-dir runs/new_mfem_surface
python tools/export_mfem_surface.py --demo --output-dir runs/new_mfem_surface_demo
```

Configured output contains `surface.h5` and `extraction_manifest.json`. The demo
also writes its authored input table/configuration under `manufactured_inputs/`,
with the exported history under `export/`. It performs no heating solve and is
explicitly synthetic. Export destinations must be fresh; source files and prior
outputs are preserved.

Load the file with `mhdlab.planar_release_io.load_planar_surface_history` and pass
it to the existing `evolve_planar_release` or `evolve_temkin_planar_release` APIs.
The existing `tools/run_temkin_planar.py --surface ... --config ...` accepts this
canonical planar file with explicit surface parameters. No new desorption or
heating solver is introduced.

The source model declares conductivity-temperature feedback, Hall induction,
material motion and qualification scope separately. The stock one-way Joule
example does not gain those mechanisms through this adapter. Moving material
requires a different surface contract. Required physical foil drive, material,
energy/current closure, temperature feedback, benchmark/convergence and surface
evaluation evidence remains external to a successful table conversion.
