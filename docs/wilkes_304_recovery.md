# Partial Wilkes 304 resistivity recovery

The recovered [CSV](../data/materials/wilkes_1997_304_resistivity.csv) contains **ten identifiable measurement markers for the 0.203 mm 304 foil**: six heating and four cooling points. The [sidecar](../data/materials/wilkes_1997_304_resistivity.json) records source hashes, calibration, marker-quality metrics, exclusions and readout envelopes. This is a partial recovery of published plotted measurements, not raw laboratory data, a fitted material law or a runtime default.

Source: K. E. Wilkes, J. P. Strizak, F. J. Weaver, J. E. Besser and D. L. Smith, *Thermophysical properties of stainless steel foils* (1997), ORNL/CP-96813, [DOI 10.2172/658177](https://doi.org/10.2172/658177), [official OSTI PDF](https://www.osti.gov/servlets/purl/658177). The validated 15-page scan has SHA-256 `5738cdc4407bf6d09da8a8f2379105ff7a32c6d77c6fbb495d5391be67d154e7`. Manuscript pages are unnumbered; page references here are one-based PDF pages. No source PDF or scanned figure is redistributed in the repository.

## What the figures actually support

| Evidence | Identification | Recovery decision |
| --- | --- | --- |
| Table I, PDF p. 4 | Specimen 1: Type 304, nominal thickness 0.305 mm. Specimen 2: Type 304, 0.203 mm. | Keep the specimens distinct. |
| Figure 1 and discussion, PDF p. 6 | Absolute electrical **resistivity** of specimen 2; ordinate is micro-ohm cm, abscissa is °C. Filled diamonds mean heating; open diamonds mean cooling. Regression and three CINDAS reference lines are also present. | Recover only ten visible, separately identifiable diamonds. No reference/regression line is digitized. |
| Figure 2, PDF p. 7 | Ordinate is **deviation percent from CINDAS**. Open diamonds identify specimen 1; open squares identify specimen 2. The legend does not distinguish heating/cooling. | No absolute resistivity or branch labels are reconstructed from this figure. The 0.305 mm foil's absolute curve remains unresolved. |

The source describes heating from room temperature to about 500 °C, then cooling. Retained Figure 1 coordinates span approximately **25.67–505.01 °C** (298.82–778.16 K). The last point is its plotted coordinate; no point is snapped to a nominal 500 °C setpoint. This coordinate span is not continuous curve support. The retained cooling points stop near 102.30 °C; they do not constitute a recovered full cooling curve or hysteresis loop.

Three merged marker groups near 200, 300 and 400 °C are masked. Their hidden point count and branch assignments remain unknown. Visible white faces inside some of these clusters do not justify inventing the centres of their partly covered neighbours. Figure 2's identifiable specimen symbols cannot supply absolute resistivity without a separately justified reference reconstruction, which is outside this recovery.

## Native-pixel method and quality

The [repeatable tool](../tools/recover_wilkes_304_resistivity.py) checks the entire PDF, page-content streams and native decoded monochrome pixels for PDF pp. 4, 6 and 7. Figure 1 is a **2528 × 3322**, one-bit source scan, decoded to 0/255 grayscale pixels at image xref 32. The Figure 1 native-pixel SHA-256 is `81c0100d553ada5509f26a78d3fa4a9231afb430ff7fd2d39558b0fb7f2a3bf5`. Native pixels are used directly, avoiding renderer-dependent page resampling.

Fourteen tick anchors are identified in fixed, inspected regions: seven temperature ticks from 0 to 600 °C and seven resistivity ticks from 60 to 120 micro-ohm cm. A two-dimensional affine least-squares transform accounts for the slight scan skew. The origin anchor is shared by both axes and included in both sets. Maximum absolute calibration residuals are about **0.464 °C** and **0.0927 micro-ohm cm**; these are axis-fit residuals, not source measurement errors.

Each retained marker has a fixed 33 × 33 native-pixel region, selected after visual inspection. The source legend determines the candidate branch; pixel geometry must then support that classification:

- Cooling uses the centroid of one enclosed white diamond face. The retained faces contain **48–56 pixels**; the local black-stroke depth is 4–6 pixels. Components touching the region boundary are excluded from the enclosed-face test.
- Heating requires a distinct filled diamond with no enclosed white face. The centre is the centroid of the maximum taxicab-distance plateau inside the inspected centre neighbourhood. All six retained markers have an **11-pixel** black-core depth, with plateau sizes **1–4 pixels**. Thin reference/regression strokes do not meet this core criterion.
- Three explicit overlap rectangles are excluded. The tool does not infer any points inside them, recover a regression polynomial or resample a drawn line. Per-marker regions, pixel hashes, centre coordinates, component counts and plateau spans are in the sidecar.

The declared centre-readout allowance is ±4 native pixels per coordinate, enlarged by half the filled-core plateau span where applicable. Its rectangular bounds propagate through the absolute affine coefficients; the maximum axis-fit residual is then added. These bounds describe display readout and centre selection. They are **not** confidence intervals, physical error bars, repeatability estimates or between-heat variability. The many stored decimal digits preserve the reconstruction calculation; they do not imply laboratory precision.

The CSV retains source °C and micro-ohm cm values alongside SI conversions: `T_K=T_degC+273.15` and `1 micro-ohm cm=1e-8 ohm m`. It contains resistivity, not its reciprocal conductivity.

## Physical scope

The report supplies no batch chemical analysis for these foil specimens. It says the rolled foils are most likely cold worked, rather than establishing a measured treatment state. These are distinct source specimens, not the dissertation's particular foil heat.

The reported 1.5–2% offset of specimen 2 from the CINDAS annealed reference is a **comparison**, not its uncertainty. The report discusses a 0.0013 mm thickness measurement resolution and effects of roughness/stock condition. Its quoted 0.4–2.5% resolution-effect range concerns the thickest and thinnest specimens across several grades; it is not adopted as a pointwise physical error for specimen 2. No Flynn error estimate is transferred to this dataset.

No thermal-conductivity curve, phase/enthalpy property, high-field correction, nanosecond calibration or absolute curve for specimen 1 is produced. These sparse markers cannot establish a width-dependent foil heating/desorption result by themselves.

## Reproduction and checks

```sh
.venv/bin/python tools/recover_wilkes_304_resistivity.py \
  --pdf /workspace/research_sources/ss304_properties/wilkes_1997_ornl_cp96813.pdf \
  --qa-overlay /workspace/research_work/mitl-upstream-20261018T1858/wilkes_304_resistivity_qa.png
.venv/bin/python -m pytest -q tests/test_wilkes_304_recovery.py
```

The optional integration test skips when the external source PDF or optional PDF-extraction dependency is unavailable. Parsing, morphology, axis, unit, exclusion and committed-artifact tests do not require that external PDF. Source inspection crops, a labelled native-pixel overlay and exact command evidence stay outside Git in `/workspace/research_work/mitl-upstream-20261018T1858/`; the sidecar records the overlay hash.
