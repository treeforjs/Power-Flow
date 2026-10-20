# Published magnetic-field profiles across the gap

The upper-right panels of dissertation Figures 4.28 and 4.29 contain separate
foil-edge and foil-middle magnetic-field lineouts. This recovery preserves
their signed plotted values and native raster support. These are published
legacy **Allegra simulation** profiles, associated by their captions with
Mykonos shot 17133 at nominal 40 and 60 ns. They are reconstructed display
coordinates, not measured fields, original simulation arrays or FLASH output.

| Figure; printed / PDF page | Recovered files | Caption snapshot |
| --- | --- | --- |
| 4.28; 94 / 119 | [CSV](../data/dissertation/figure_4_28_magnetic_gap.csv), [metadata](../data/dissertation/figure_4_28_magnetic_gap.json) | 40 ns |
| 4.29; 96 / 121 | [CSV](../data/dissertation/figure_4_29_magnetic_gap.csv), [metadata](../data/dissertation/figure_4_29_magnetic_gap.json) | 60 ns |

Each profile retains 131 native rows, for 524 samples in total. Their central
support is approximately −1.658228 to +1.632911 mm, with no inserted zero-gap
sample or padded endpoint. All retained B_x values remain negative on both
coordinate halves. The complete catalogs account for 1,572 candidate rows,
including 1,048 exclusions.

The literal field label is **B_x [100 Tesla]** and the spatial coordinate is
distance from the middle of the A–K gap in mm. The primary data retain that
literal field ordinate. Separate explicitly conditional SI columns multiply
it by the printed factor of 100 tesla. The original electromagnetic
normalization and unit-export definitions remain unavailable. No scale repair,
sign change, magnitude operation or field-frame rotation is applied.

## Native support and readout

The [recovery tool](../tools/recover_dissertation_gap_magnetic.py) validates the
dissertation PDF, page content and decoded native RGB image identities before
using inspected axis calibrations. Each sidecar binds the selected curve
colors, tick anchors, calibration coefficients and residuals, retained exact
color runs, excluded candidate coordinates and generator hash. Source pixel
locations remain available alongside the calibrated values.

The conservative mask retains native rows 150–280 inclusive. The field-axis
tick centres are 606.5, 706.5, 806.5 and 906.5 pixels for −1, −0.5, 0 and +0.5
literal field units. The gap tick centres are 56.5, 135.5, 214.5, 293.5 and
372.5 pixels for +4, +2, 0, −2 and −4 mm. Both exact rational fits have zero
residual on those declared tick centres; this does not eliminate pixel-readout
or physical uncertainty.

Each profile has directly visible finite support. Omitted or ambiguous
parts are not bridged, smoothed, interpolated or extended. Legends, overlapping
strokes, shared branches, sharp boundary neighborhoods and unsupported
endpoints must not become inferred curve values. The complete selection
catalog records exclusions rather than silently discarding them.

Readout bounds include the complete matched run, an additional ±2 pixels
and the maximum tick-fit residual. They describe correlated raster-display
readout, not independent observations, physical uncertainty or a statistical
confidence interval. Native-panel overlays and signed-field plots accompany
the local verification evidence.

## Reproduction and independent checks

Run into new output locations with the source PDF available. NumPy, PyMuPDF,
Pillow and Matplotlib must be available for this source-recovery command:

```sh
python tools/recover_dissertation_gap_magnetic.py --pdf /path/to/smtrevor_1.pdf --output-dir runs/magnetic_gap/data --qa-dir runs/magnetic_gap/qa
```

The public CSVs and sidecars are the saved recovered inputs. Source geometry,
full selection catalogs, native inspection and QA artifacts remain associated
with the local handoff. Independent verification reads the saved artifacts
against native source pixels and separately calculated axis transformations;
it does not replay the recovery producer.

Four new invalid-input controls rejected an altered PDF, a false native RGB
identity, an erased tick mask and two ambiguous exact-color runs. One complete
independent readback checked all 524 retained rows, all 1,048 exclusions,
7,860 numeric CSV cells and 18 tick anchors. Its original tick-tail scan and
the producer-declared scan independently give the same exact calibration.
The first attempt stopped on that scan-row declaration mismatch before the
catalog and CSV checks; the recorded work-only correction changed neither
the producer nor the recovered inputs. No earlier suite or positive recovery
was replayed.

## Scientific use and limits

These coarse gap profiles add a different coordinate cut to the separately
recovered [current-density depth profiles](dissertation_current_depth_recovery.md)
and [electromagnetic width profiles](dissertation_recovery.md#recovered-electromagnetic-width-profiles).
They do not recover a complete vector field or a spatial field mesh. The
nominal 200 µm foil occupies only a few native pixels on this gap-scale panel,
so it supplies no resolved skin profile, field derivative or exact surface
field. The gap-centred plot coordinate does not identify an evolving material
boundary.

No Ampere-law closure, current reconstruction, magnetic-energy integral,
Joule heating, Hall drift or charged-particle magnetization is inferred.
Those require qualified original units, geometry and field definitions.
Caption snapshot labels also do not register these profiles to the
experimental diagnostic clock. The original drive arrays, evolving material
fractions, mesh and material closure remain unavailable.

The recovered temperature interiors and finite current insets remain separate
quantities. These magnetic profiles install no heating history, desorption
boundary temperature, gap transport or owner-solver interface. They provide
source-qualified comparison shapes while physical normalization and coupling
remain unresolved.
