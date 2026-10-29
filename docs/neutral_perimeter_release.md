# Four-face neutral release in the early VUV window

The perimeter case adds the left and right material sides to the two broad
emitting faces. It retains the same shot 17133 current, 60 ns stop time,
constant-property MFEM heating and illustrative independent surface reservoirs.
It also removes the previous display's fourfold enlargement of the normal
direction. No launch speed is increased to change the appearance.

## Why the previous movie looked directional

The previous source emitted only from y-normal broad faces. Its x velocity
was the tangential component of the thermal flux law. Under that law the
normal second moment is twice either tangential second moment. This gives a
physical RMS normal/tangential speed ratio of sqrt(2). The old display enlarged
y by four, making that ratio appear approximately 5.66 rather than 1.41.

An audit conditioned on each old packet's actual birth temperature and weight
found normalized second moments approximately (1,2,1). At the 60 ns frame:

| Species | RMS x displacement since birth | RMS y displacement since birth |
| --- | --- | --- |
| H | 45.59 micrometres | 64.50 micrometres |
| H2 | 32.15 micrometres | 45.61 micrometres |
| H2O | 10.77 micrometres | 15.24 micrometres |
| OH | 11.08 micrometres | 15.71 micrometres |

The weighted mean packet age was only 11.06 ns; particles had not all flown
for the full current-pulse duration. These numbers describe the earlier
two-broad-face control, not the new perimeter source. The thermal flux law is
also described in the [COMSOL molecular-flow derivation](https://doc.comsol.com/6.3/doc/com.comsol.help.molec/molec_ug_fluidflow_rarefied.6.36.html).

## Actual material boundaries

`native/mfem_edge_source/edge_driver.cpp` evaluates each side from the thermal
finite-element material boundary using separate quadrature. The side traces
are not estimated from the broad faces or from volume-cell averages.

| Face | Physical normal | In-plane patch direction | Area over 10 mm current length |
| --- | --- | --- | --- |
| Top/back | +y at 2.7 mm | x across 5 mm width | 5e-5 m2 |
| Bottom/gap | -y at 2.5 mm | x across 5 mm width | 5e-5 m2 |
| Left | -x at -2.5 mm | y through 0.2 mm thickness | 2e-6 m2 |
| Right | +x at +2.5 mm | y through 0.2 mm thickness | 2e-6 m2 |

The sides together have 4% of the broad pair's area. The case assigns equal
per-area initial coverage and rate laws on every face; this is an explicit
surface-condition assumption. Side inventories use their actual smaller areas.
Different surface preparation or a nonthermal desorption law would require
separate inputs. No side source multiplier is fitted.

For side-face sampling, local (u,v,z) means physical (y,x,z). The existing
planar law samples the signed local v-normal. The explicit orthogonal axis
permutation maps positions and velocities into physical xyz, preserving speed,
kinetic energy and the left/right outward signs. Metadata retains both the
source-local frame and the physical packet frame through saved histories.

The new source run is `runs/mfem_edge_current_source_four_faces`. Its broad
face tables, diagnostics, cell snapshots and schedule are byte-identical to
the earlier v3 run, as recorded in `four_face_verification.json`. Side traces
reach about 12074 K, further emphasizing the fixed-property stationary-solid
model's limited physical applicability near the corners. The published ALEGRA
comparison and its existing qualifications therefore remain applicable.

## Numerical representation

Hot side pools can become depleted to positive binary64 subnormal tails.
Multiplying those per-area tails by the small side-patch area can underflow a
packet's represented number. The original law history remains unmodified.
`prepare_sampling_interval` omits only components below the case's explicit
1e-250-particle representation cutoff when constructing a validated sampling
view. It reuses the same frozen-temperature law and checks that all represented
release and energy components are bit-identical. Omitted component counts,
exact-decimal totals and conservative number bounds are recorded separately in
each face's `packet_numerical_omissions`; no density floor is applied.

## Reproduce the source and views

```powershell
python tools/run_mfem_edge_source.py --wsl-distro Ubuntu-24.04 --output-dir runs/new_four_face_heating
python tools/run_neutral_desorption_map.py --top-surface runs/new_four_face_heating/top_surface.json --bottom-surface runs/new_four_face_heating/bottom_surface.json --left-surface runs/new_four_face_heating/left_surface.json --right-surface runs/new_four_face_heating/right_surface.json --case configs/neutral_desorption_map_17133_perimeter.json --output-dir runs/new_perimeter_maps
python tools/render_neutral_transport_gif.py --input runs/new_perimeter_maps/neutral_maps.npz --metadata runs/new_perimeter_maps/metadata.json --output-dir runs/new_perimeter_overview --fps 12 --dpi 110 --color-span-decades 8
python tools/render_neutral_transport_gif.py --input runs/new_perimeter_maps/neutral_maps.npz --metadata runs/new_perimeter_maps/metadata.json --output-dir runs/new_perimeter_edge --fps 12 --dpi 110 --color-span-decades 8 --view-bounds-mm 2.25 2.9 2.2 3.0
```

Both views use equal physical x/y length scales and the same saved maps, times
and fixed color scales. The edge view crops only the axes. Collisionless
transport, independent assumed species pools, passive removal-energy accounting
and the absence of a calibrated chemistry or spectral fit remain unchanged.

## Current observation grid and physics extensions

The completed Arch source now has a65-frame four-face map with8,824,000
weighted packets. The original208×144 observation grid changed width bins
from62.5um over the foil to6.25um beside it. Refining all bins fourfold retained
that jump. A uniform3.125um grid (1856×320) removes this observation-grid
seam while retaining the same packets, physical bounds and frame times.
Coarsening its cell inventories recovers the original density to6.71e-14
relative, under a2e-12 numerical comparison threshold. Fine H/H2 cells remain
limited by weighted packet statistics; no smoothing was applied.

The residual-gas request uses the50 nominal mass channels of dissertation
Figure2.3. The450 encoded PDF rectangles are graphics provenance, not450
independent observed bars. The reported pressure is approximately5e-5Torr;
it is not imposed as an exact normalization. Background gas interactions,
plume self-collisions and a calibrated RGA parent composition are not present
in the current GIFs. A conditional collision-size screening has been computed;
it does not enable a collision operator or replace energy-dependent data.

Water fragmentation work now explicitly uses Raikou. A sourced Table10 OH-yield
screening advanced99 prescribed monoenergetic-electron cases through60ns with
H/O/charge checks and an analytic single-channel comparison. It does not close
the full19-species water model: after binding this OH coefficient, the full
mechanism still refuses missing electron/heavy rates at zero accepted time.
The water-monolayer count remains adjustable; no RGA-derived surface coverage
or independent fragment pools are assigned to this parent-water source.

The requested surface-temperature-dependent electron current is now evaluated
with the unrestricted Richardson-Dushman law, with explicit work function and
Richardson coefficient. `mhdlab/thermionic_supply.py` is an analytic boundary
supply calculation, not a second electron-kinetics solver. Its outputs have no
space-charge cap, field enhancement or collected-current claim. The next
physical closure must normalize Raikou's incoming electron distribution to
this supply, solve the transient Poisson/sheath/circuit field, and obtain the
EEDF from that field and collisions. The axial conductor voltage gradients
from the heating solve are not a cross-gap voltage. Surface phase changes,
work-function evolution, magnetic insulation and full reaction energy closure
remain required physical qualifications.
