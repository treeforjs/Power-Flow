# Conditional comparison with the dissertation ALEGRA displays

`tools/compare_mfem_alegra.py` compares MFEM temperature snapshots with recovered
simulation curves from dissertation Figs.4.28–4.30. These are published ALEGRA
displays, labeled “Allegra” in the dissertation, not measured temperatures or
the original simulation arrays. The comparison does not fit material constants,
repair drive amplitudes, shift clocks or claim physical reproduction.

The drive context is shot **17133**, named for the heating calculation. Its
[Fig.C.19 current](../data/dissertation/shot_17133_current.csv) is distinct from
the finite current insets in Figs.4.28 and 4.29. Their central displayed curves
differ on their literal coordinates; original simulation input and physical
clock registration are unavailable. The current panel preserves those separate
curves and can overlay the actual current enforced by the MFEM driver.

## MFEM snapshot contract

The input CSV columns are exactly:

```text
time_s,x_left_m,x_right_m,depth_left_m,depth_right_m,temperature_k,current_density_a_m2,joule_power_w_m3
```

A required JSON sidecar uses `mhdlab-mfem-foil-snapshots-v1`, explicit status,
`sampling=cell_average`, fixed SI units, CSV SHA256, and geometry identifying
the upper foil and depth inward from its gap-facing surface. Its width and
thickness must match a complete nonoverlapping rectangular cell partition.
Exact snapshots at 40 and 60 ns are required; the tool performs no temporal
interpolation. Additional drive/model/source metadata are retained in full.

```powershell
python tools/compare_mfem_alegra.py --snapshots runs/case/foil_snapshots.csv --diagnostics runs/case/diagnostics.csv --output-dir runs/new_mfem_alegra_comparison
```

The optional diagnostics table uses `time_s` and `foil_current_a` for the drive
overlay. The output directory must be new. JSON, PNG and PDF outputs preserve
input/source hashes, all retained readouts, prediction differences and the
sampling conventions below.

## Comparison quantities and sampling

Figs.4.28 and 4.29 retain 38 points each for the foil-edge and foil-middle depth
profiles at 40 and 60 ns. Their depth coordinates reference the initial 5-mm
gap face and 200-µm foil; actual moving material boundaries are not recovered.
Their temperature axes explicitly use Kelvin.

The MFEM cell averages are read as a piecewise constant reconstruction. Interior
cell ties belong to the right/deeper cell; the final boundary belongs to the
last cell. This is a declared observation approximation, not a native point
temperature or material-surface trace. No spatial extrapolation is allowed.

The precise transverse location of the published foil-edge lineout is unknown.
The tool therefore reports the mean of the nearest left/right MFEM edge cells,
the central cell(s) touching the foil middle, and a separately labeled min/max
position sensitivity across the outer 10% of model columns per side. This
sensitivity is not statistical uncertainty and does not establish identical
probe placement. The chosen cell intervals are saved explicitly.

Fig.4.30 supplies a width lineout at nominal depth 12.5 µm and time 40 ns. Its
published eV label is converted to Kelvin only under the existing conditional
interpretation `temperature_eV=k_B*T/e`. The original temperature-variable and
energy-mode ownership remain unverified. The tool checks that the committed
conditional Kelvin column matches this conversion before using it.

Plotted bands are correlated glyph readout envelopes. They are not diagnostic
noise, physical confidence intervals or independent likelihood weights.
Differences are descriptive maximum absolute departures on retained support;
there is no parameter fit or amplitude normalization. Unverified published
current-density normalization is not used to infer a physical Joule-power fit.

Physical comparison still depends on the actual material law, edge radius,
thermal/field boundaries, support geometry, drive history and clock conventions,
phase response and convergence. A provisional MFEM constant-property model can
show whether field-derived edge concentration occurs, but visual similarity
alone does not qualify the original ALEGRA model or the new calculation.

## First current-driven comparison, capped at 60 ns

The native source in `runs/mfem_edge_current_source_v3` was compared with the
published 40 and 60 ns displays. The resulting plot and complete numerical
report are in `runs/mfem_alegra_comparison_v3`. The run starts at the first
supported C.19 knot, 0.749788 ns, with a declared 300 K initial state and static
magnetic initialization at that current. Deposited heat before that point and
physical registration with the ALEGRA clock remain unknown. No later evolution
is needed for this early-window comparison.

At the first retained depth coordinate, 7.13 µm, the comparison is:

| Time | Profile selection | Published ALEGRA [K] | MFEM cell average [K] | Maximum absolute departure over all 38 retained depths [K] |
| --- | --- | ---: | ---: | ---: |
| 40 ns | nearest edge cells | 1212 | 3354 | 2225 |
| 40 ns | central cells | 540 | 471 | 93 |
| 60 ns | nearest edge cells | 3225 | 9089 | 6035 |
| 60 ns | central cells | 1173 | 862 | 404 |

At nominal depth 12.5 µm, the resolved nearest-edge/central Joule-power ratios
are 16.33 at 40 ns and 13.52 at 60 ns. These are properties of the actual native
field, not an imposed temperature or edge factor. The source driver's separate
outer-band diagnostic averages over a wider region and is a different sampling
quantity. The comparison selects the first depth cell, 0–25 µm; its cell average
does not resolve a pointwise temperature at 12.5 µm or the surface temperature.

For the conditional Fig.4.30 Kelvin interpretation, the nearest retained
width coordinate to the middle is −5.76 µm. Its temperature is 492 K in the
published curve and 471 K in MFEM. The maximum absolute departure over all 137
retained width points is 1432 K, dominated by the sharp outer cells.

The model therefore produces edge concentration, but does not reproduce the
published edge temperatures. Its sharp rectangular edge, unresolved probe
placement and fixed material properties require further physical qualification.
The 60 ns edge temperatures particularly expose the omission of phase response
and temperature feedback. The outer-band position envelope is sensitivity to a
declared sampling choice; overlap with a published curve does not establish an
edge-location match. No parameters, time shifts or edge-cell choices were fitted
to reduce these departures.

The comparison command exited successfully, and its PNG was visually inspected.
The JSON records the source snapshot, sidecar, diagnostics, dissertation data
and tool hashes. Focused tests exercise grid completeness, exact-time support,
SI/quantity association, P0 sampling and the distinction between actual Joule
contrast and an authored temperature contrast.
