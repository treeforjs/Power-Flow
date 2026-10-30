# Matched native thermal input and mechanical-support audit

The expansion benchmark uses only
`runs/quick_heating_256x32_domain40mm_height25mm_20261029`, pinned by the external
SHA256 of `source_run.json`:
`df79133957684abbe5678b10b0c4b7625cbb0c8b7e17d31712410decd5efca05`.
The loader verifies consumed files against this manifest and each native field
against `em_field_snapshots.json`. The manifest contains an obsolete self-entry
under `outputs_sha256/source_run.json`; this self-entry is not a valid hash of the
final JSON and is not used. The external hash above is the identity anchor.

The original run source snapshot, particularly
`source_snapshot/native/mfem_edge_source/edge_driver.cpp`, determines the data
interpretation. `Thermal` initializes its scalar FE unknown to zero (line 126).
`Inspect` adds 300 K to the upper volume cell average (line 218). It emits volume
rows only for `conductor==0` (line 214). Surface writers emit increments for both
electrodes (lines 291–297), and their metadata explicitly declares the 300 K
offset. Consequently, the adapter adds 300 K to face tables once and adds nothing
to `foil_snapshots.csv`.

There are 67 bulk/field snapshots from 0.749788055492 ns through 60 ns. Every
snapshot has 256 width cells and 32 thickness cells. Native element order is
not Cartesian row order; the loader indexes cell coordinates and checks complete,
unique coverage. It returns arrays ordered `(time, ascending y, ascending x)`.
The upper foil occupies x = ±2.5 mm, y = 2.5…2.7 mm. The lower foil occupies
x = ±2.5 mm, y = −2.7…−2.5 mm. The modeled thermal/EM fields are invariant along
the 25 mm axial height; there is no measured/resolved axial thermal gradient.

All eight actual surface tables are loaded. **The lower bulk field was not
exported.** The default loader leaves it unavailable. Passing
`lower_volume="mirror_with_surface_validation"` explicitly selects a benchmark
reflection of the upper bulk. It checks upper top against lower bottom, upper
bottom against lower top, and left/right faces with reversed thickness ordering.
On the specified matched run the maximum discrepancy across all recorded faces
and times is 6.87696e−7 K, below the declared 1e−4 K tolerance. This supports the
symmetry scenario; it does not turn the missing lower volume into independently
exported data.

The full history is unsuitable for unchecked solid elasticity. Its maximum
upper cell-average temperature is 12,269.544 K, and its maximum directly exported
surface-patch temperature is 13,300.599 K. Neither is clipped by the loader.
For a 293.15–773.15 K material interval, the first invalid bulk snapshot is
19.265479288 ns and the preceding bulk snapshot is 18.339694727 ns. Direct face
temperatures exceed the upper limit earlier, at 17.622896529 ns; the preceding
face snapshot is 17.413910165 ns. Thus the last common bulk/field snapshot before
the observed face violation is **17.413910165 ns**. Later snapshots must be
rejected by this solid model. Even accepted cell/patch averages do not prove that
unexported nodal or subcell maxima stay within the interval. The source remains
a constant-property heating scenario without latent heat or phase change.

`MatchedThermalHistory.admissibility(minimum_k, maximum_k)` records both bulk and
surface violations and the last common admissible bulk time. This is a material
input guard, not a statement of mechanical, mesh or experimental convergence.
The detailed consumed-file hashes and numerical source audit are in the report
directory `outputs/ss304-thermal-expansion-20261029/source-audit.json`.

## What the models actually support mechanically

No foil displacement, traction, clamp stiffness, weld compliance, residual
strain, mechanical contact or support condition is specified in this native
heating source. `Thermal` assembles scalar heat capacity and conductivity with
natural insulated boundaries (lines 118–138). `EMGrid` constrains magnetic
potential on the outer electromagnetic boundary (lines 70–79); these are not
displacement constraints. `examples/geometry/mykonos_foil_gap_5mm.json` supplies
stationary material polygons and a line of sight, with no mechanical conditions.

The finite-axial PIC example in `raikou-mfem/examples/mfem_heating_history_pic.cpp`
describes an ideal electrical short at z = −12.5 mm, feed at z = +12.5 mm and
observation at z = 0 (lines 126–143). It imposes the axial voltage profile and
particle absorption/escape. The 5×5 mm short footprint supplies electrical and
particle geometry; it does not identify a clamp, its motion, a finite-stiffness
joint, or a structural short mesh. The preserved run description likewise does
not provide a measured mechanical fixture model.

Accordingly, free, fixed-end and fully constrained calculations must remain
explicit benchmarks. The free case requires a rigid-motion origin/gauge; it
does not fix the experiment's absolute gap drift. An affine foil-only geometry
interface must not pretend to move the short or feed consistently. Production
coupling requires an explicit support/short model and an admissible constitutive
model, plus compatible volume/surface/field transfers and their conservation
checks. Static equilibrium also cannot substitute for transient elastic-wave
response during a nanosecond heating pulse.

## Interface and bounded verification

```python
from mhdlab.thermal_expansion_source import load_matched_history

h = load_matched_history(run_directory)  # Exact preferred quick-run hash by default.
guard = h.admissibility(293.15, 773.15)
# h.upper_temperature_k: absolute K, native cell averages.
# h.lower_temperature_k: None unless explicit symmetry reconstruction requested.
# h.native_fields: verified Az and dAz/dt paths, hashes and identical native times.
```

No mesh interpolation, temporal interpolation, geometric deformation or EM solve
occurs in this reader. Temperature is not relabeled as conserved thermal energy.
The companion geometry/remap module owns those explicit representations and
their separate guards. The 19 source-adapter tests exercise offsets, native
ordering, lower-field labeling/reflection, all-file integrity, exact-time
matching, incomplete/duplicate cells, geometry, hot-cell guards and asymmetry
rejection. They pass in the frozen userspace with CPU 7, nice 15, one thread and
the shared `testing-core7.lock`.
