# Early-time current-driven neutral transport

The GIF shows **H, H2, H2O and OH** desorption and collisionless flight in the
foil width/thickness-normal plane from **0.749788 ns to 60 ns**. The source uses
the recovered dissertation **shot 17133, Fig. C.19 current trace**. The run stays
within that trace's retained support and is capped at 60 ns for the early VUV
window. It supersedes the earlier uniform-voltage, metal-only demonstration.

Open `neutral_desorption_H_H2_H2O_OH.gif` and `mfem_alegra_comparison.png`.
The three previews were decoded from the actual GIF at its first, middle and
last frames. All 65 supplied times are retained; loop duration is 5.2 seconds.

## What the calculation establishes

The custom MFEM 2D A_z magnetoquasistatic driver solves a finite-width foil,
oppositely driven return conductor and exactly nonconducting exterior. Both
integrated conductor currents are constrained to the C.19 drive. Actual J^2/sigma
quadrature supplies separate conductor-only thermal finite-element solves.
Edge concentration follows from the fields; no edge temperature multiplier is
imposed. This is an MFEM extension, not the unchanged stock Joule kernel or a
full MHD calculation.

At 60 ns, the broad edge-band/central-band mean Joule ratio is approximately
7.806. Edge-band mean temperature is 2625 K; central-band mean is 570 K. These
bands average through the thickness and differ from the nearest edge-cell
lineouts used in the ALEGRA comparison.

Doubling the conductor mesh from 64x8 to 128x16, and doubling the other mesh
counts, preserves edge concentration. The native band Joule ratio changes by
-3.03%; edge-band mean T changes by -68.4 K. Expanding the field half-domain
from 10 to 20 mm changes that ratio by -0.204%. The separate sensitivity report
also quantifies surface/interior lineout changes. These checks do not establish
convergence; corner maxima increase under refinement.

The comparison uses the dissertation's published ALEGRA displays in Figs.
4.28-4.30, with exact 40/60 ns model outputs. At the first retained inward depth
(approximately 7.13 micrometres):

| Time | MFEM nearest edge / ALEGRA edge | MFEM middle / ALEGRA middle |
| --- | --- | --- |
| 40 ns | 3354 / 1212 K | 471 / 540 K |
| 60 ns | 9089 / 3225 K | 862 / 1173 K |

The interior is closer to the published temperature profiles, while the nearest
resolved edge cells overshoot substantially. The original ALEGRA edge probe
offset, input waveform and clock registration are unavailable. A separately
plotted outer-band position sensitivity is not statistical uncertainty. The
Fig. 4.30 eV-to-K interpretation remains conditional. Published readout envelopes
represent glyph recovery bounds, not experimental temperature uncertainty.

## Explicit model assumptions

- Upper foil x=-2.5 to +2.5 mm, y=+2.5 to +2.7 mm; thickness 200 micrometres.
  The identical return foil is at y=-2.7 to -2.5 mm. The invariant current
  length is 10 mm, used once in power, particle-number and volume accounting.
- At the first supported C.19 knot, I=6012.071 A, initialize uniform DC current
  and its static magnetic field, with T=300 K. Prehistory is unknown. Published
  times are unshifted; comparison with the ALEGRA clock is an unverified scenario.
  I(40 ns)=308951.866 A and I(60 ns)=500734.776 A.
- Constant illustrative sigma=3e6 S/m, Cv=2.8e6 J/(m3 K), k=16 W/(m K), mu=mu0;
  no temperature feedback, phase change, material motion or Hall induction.
  Face patches reach approximately 9811 K at 60 ns, and cell-average temperatures
  already exceed 3300 K at 40 ns. The stationary solid/property assumptions fail
  near the hottest corners. These results are a provisional numerical scenario,
  not a qualified stainless-steel prediction or reproduction of ALEGRA.
- Four independent assumed adsorbate reservoirs emit from the two broad faces;
  side edges do not emit. Initial per-face coverages are H=2e15, H2=2e15,
  H2O=1.4e16, OH=2e15 per m2. Common first-order attempt rate 1e12/s and barrier
  0.8 eV are illustrative. There is no inferred water fragmentation or reaction
  branching and no calibration of absolute yields.
- Existing finite-inventory release and exact conditional birth-time sampling
  supply 3,432,000 weighted packets. Launches use flux-weighted thermal normal
  velocities and Gaussian tangential velocities. All three velocity components
  are retained; current-direction transport is invariant/periodic.
- Flight is collisionless and force-free; map exits and solid entries absorb.
  Removal energy is a passive diagnostic. There is no plasma, CR/spectral fit,
  neutral density/temperature inference or radiation feedback in this GIF.

## Display and verification

Coordinates are physical millimetres, with the normal direction explicitly
enlarged 4x. Fixed per-species log color scales span at most eight decades below
the global maximum. Smaller positive values remain in the lowest color and in
particle totals; only exact zeros are white. No density floor, smoothing,
rescaling of populations, or time interpolation is applied.

The actual encoded GIF was decoded and its first, middle and last frames were
visually inspected. The particle ledger's largest relative number residual is
2.16e-17 of initial inventory; histogram residual is 5.68e-15. Current constraint
error is below 2.48e-8 A; maximum EM and thermal increment-balance residuals are
6.14e-12 J and 8.00e-14 J. Expected and sampled launch energies are separate.
Seven renderer tests pass, including an extreme under-range-tail test that
checks unchanged density values and particle totals. Earlier source/transport
checks and comparison contract tests are documented in the project.

`inputs` contains the actual heating faces, current schedule, diagnostics,
snapshots, transport arrays, ledgers, comparison and sensitivity records. Large
packet and release-history files remain in the original run directory and can
be regenerated with the saved seed. `manifest.json` binds every packaged file
except itself by SHA256. Original run records retain their original paths and
hashes. The replay runner includes a post-run guard rejecting T0 other than
300 K; the saved native case already used 300 K and its fields are unchanged.

## Replay

Use Python >=3.11 with numpy, scipy, h5py, matplotlib, pillow, PyYAML and imageio.
From this package's `replay` directory, using fresh output directories:

```powershell
python tools/render_neutral_transport_gif.py --input ../inputs/transport/neutral_maps.npz --metadata ../inputs/transport/metadata.json --output-dir runs/replayed_gif --fps 12 --dpi 110 --color-span-decades 8
python tools/compare_mfem_alegra.py --snapshots ../inputs/heating/foil_snapshots.csv --diagnostics ../inputs/heating/diagnostics.csv --output-dir runs/replayed_comparison
python tools/run_neutral_desorption_map.py --top-surface ../inputs/heating/top_surface.json --bottom-surface ../inputs/heating/bottom_surface.json --case ../inputs/transport/case.json --output-dir runs/replayed_transport
```

Fresh heating additionally needs installed MFEM 4.10 with MPI/Hypre/Metis,
CMake/Ninja and the existing Ubuntu-24.04 WSL environment. The owned native
driver is included; MFEM itself is not vendored:

```powershell
python tools/run_mfem_edge_source.py --wsl-distro Ubuntu-24.04 --output-dir runs/recomputed_heating
```
