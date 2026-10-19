# Signed current-density depth profiles

The two recoveries retain the foil-edge and foil-middle J_z curves in the
lower-left panels of dissertation Figs. 4.28 and 4.29. These are published
legacy “Allegra” simulation results at 40 and 60 ns, driven by Mykonos shot
17133 with a reported 200 um SS304 foil and 5 mm initial gap. They are
reconstructed raster lineouts, rather than experimental current-density data
or FLASH output.

| Figure, printed / PDF page | Data and provenance | Samples |
| --- | --- | --- |
| 4.28, 94 / 119 | [CSV](../data/dissertation/figure_4_28_current_depth.csv), [metadata](../data/dissertation/figure_4_28_current_depth.json) | 141 per profile, 40 ns |
| 4.29, 96 / 121 | [CSV](../data/dissertation/figure_4_29_current_depth.csv), [metadata](../data/dissertation/figure_4_29_current_depth.json) | 141 per profile, 60 ns |

J_z is plotted horizontally; the vertical axis is distance from the middle of
the A–K gap. Fresh tick calibration uses each panel's own axes. The signed
literal MA/cm² values remain primary. All retained values are negative; no
absolute-value replacement, sign repair or host-frame rotation is applied.
The `Jz_literal_conditional_A_per_m2` columns multiply the printed values by
10¹⁰ A/m² per MA/cm². **The original simulation normalization is unverified.**
These columns cannot establish physical current density or J²/sigma heating.

Native rows 640–780 provide 564 samples across the four curves. Their plotted
gap coordinate spans approximately 2.52232–2.66407 mm. The separately named
nominal initial depth spans 22.3201–164.068 um from the declared initial
+2.5 mm gap face. That reference does not recover a moving material boundary,
a surface value or the precise transverse location of the foil-edge lineout.

Each sample retains its native row, exact-color run and median x coordinate.
The conservative interior mask excludes visible corners, broad boundary ramps,
exterior/zero branches and legends. All 396 candidate rows per profile are
catalogued, including 255 excluded rows per profile, for 1,020 exclusions in
total. No missing region is filled by interpolation or smoothing.

Current readout bounds include the full matched color run, an additional
two pixels on each side and the largest current tick-fit residual. Position
bounds use two pixels and the corresponding coordinate-fit residual. These
are correlated display bounds, not physical uncertainty, confidence intervals
or independent likelihood errors. Adjacent raster rows are not independent
measurements or recovered simulation mesh nodes.

Reproduce the new profiles and their PNG/PDF QA figures from the hash-matched
dissertation:

```sh
python tools/recover_dissertation_depth_current.py --pdf /path/to/smtrevor_1.pdf --output-dir data/dissertation --qa-dir runs/current_depth_recovery
```

The tool checks the PDF, decoded page contents, native RGB pixels and detected
tick geometry before extraction. Sidecars bind the source, calibration,
retained/excluded rows, generator, CSV and QA/source-geometry files. Existing
temperature and electromagnetic width recoveries remain separate artifacts.

The recovered shot 17133 current display supplies same-shot context, while
current-to-simulation clock registration remains unresolved. These partial
depth profiles supply no whole-depth or full cross-section current integral,
Ampere-law closure, thermal history, conductivity or deposited energy. The
magnetic gap panels remain pending. The [finite simulation-drive insets](legacy_simulation_drive_recovery.md) are separately recovered; original inputs and physical clock registration remain unresolved.
