# Dissertation data reconstruction

Trevor Johannes Smith's 2025 dissertation contains gated VUV emission spectra, spatial lineouts across the anode–cathode gap, electrical traces and selected simulation profiles that can support a publication reconstruction. These are published observations and derived estimates. Original detector data, original electrical waveforms and complete per-shot calibration have not been recovered.

The source is [Characterizing Desorbed Surface Contaminants Leading to Plasma Formation in Power Feeds of Pulsed Power Accelerators](https://deepblue.lib.umich.edu/items/9ccfa886-0b19-4bd8-b3b2-bb685b9ddefc). The [current repository PDF endpoint](https://backend.production.deepblue-documents.lib.umich.edu/server/api/core/bitstreams/24e682a8-896a-410d-85bb-ca4654a9db02/content) supplied the validated 203-page PDF. Its SHA-256 is `1d7b1f9c191a37e94f1cc3634587750f1af35b07eb13e0b94c60f3da47495ed7`. For its Arabic numbered pages, add 25 to the printed page to get the one-based PDF page. The machine-readable [recovery manifest](../data/dissertation/recovery_manifest.yaml) records source provenance, 44 figure/table targets, published table entries, partial curve recoveries and unresolved metadata. The PDF and extracted images remain outside this repository.

## What the VUV diagnostic measures

The Mykonos target has a 5 mm gap, 5 mm electrode width and 25 mm height (Fig. 4.2, printed p. 54; PDF p. 79). The support uses Ultem 1010 and epoxy (p. 55). The VUV diagnostic records **emission**, using a McPherson 234/302 spectrometer and Photek iCMOS 160 camera over approximately 115–165 nm with a 3 ns gate. The reported 0.03 nm/pixel describes spectral sampling; a calibrated instrumental line-spread function is still required for spectral inference (pp. 78–79).

Each shot supplies one gated spectral frame with a spatial axis across the anode–cathode gap. The other spatial direction, across the foil width, is integrated along the line of sight. The dissertation assumes an emitting depth of 5 mm in its intensity conversion (Eq. 4.8, p. 92). Therefore the images do not directly measure the transverse desorption profile. They can constrain spectra projected from a width-dependent source, subject to degeneracies in excitation, emitting depth and source inventory.

The directly identified VUV features are atomic H Lyman-alpha and molecular H2 Lyman/Werner emission bands (pp. 78–83). Figure 2.3 (p. 19) provides a background residual-gas spectrum consistent with water and hydrogen being important constituents. It does not measure the water parent inventory of an individual discharge. Visible spectroscopy supports oxygen-bearing species, including O II near 464 nm, and later metal emission (pp. 78 and 100). Section 2.4 also discusses hydrogen trapped in SS304 bulk (p. 17). An H2O breakdown network is a physical hypothesis to constrain alongside these alternatives; the published VUV data do not uniquely establish a water-only parent population.

## Reconstruction targets

| Evidence | Printed pages | Reconstruction use |
| --- | --- | --- |
| Fig. 4.2 and Figs. 4.3, 4.5 | 54, 56, 58 | Target dimensions, campaign sight lines and collection regions |
| Figs. 4.17–4.19 and accompanying text | 79–82 | Instrument assembly, D2 lamp calibration and camera gain response |
| Figs. 4.20–4.25 | 84–89 | H/H2 gap profiles, displayed uncertainty bands and three wavelength lineouts |
| Table 4.1 | 89 | Published by-eye Yacora estimates, transcribed below and in the manifest |
| Figs. A.1–A.17 | 102–118 | Seventeen VUV shots, associated current/resistive-voltage traces and gate markers |
| Figs. 4.10–4.12, 4.14–4.16 | 66, 68, 70, 73, 75, 77 | Visible spectral/timing context with campaign-specific collection regions |
| Figs. 4.27–4.30 | 93–98 | Existing simulation geometry; temperatures and Fig. 4.30 EM width profiles partially recovered; other field/current panels pending |
| Fig. C.11 | 157 | Shot 15210 encoded current polyline and VUV gate geometry recovered; voltage/SVS curves pending |

Most published spectral figures are embedded raster images. Digitization can recover displayed curves with finite plot resolution, but the pseudocolor image cannot be treated as original detector counts. Saturation, clipping, interpolation, averaging and plotting transformations may be irreversible. The first partial recovery is [Fig. 4.20's central gap profiles](../data/dissertation/figure_4_20_profiles.csv): 106 correlated raster samples, 53 each for H Lyman-alpha and the H2 161 nm Lyman band. The [provenance record](../data/dissertation/figure_4_20_profiles.json) includes source hashes, pixel coordinates, axis calibration, and pixel-reading uncertainty. These are plotted transition densities, not total ground-state neutral densities. Physical error bands and wavelength lineouts remain pending. The heating and current recoveries below are separately labeled. All other specified quantities retain their pending status; numerical table transcription is marked separately.

To reproduce this recovery with the exact source PDF, install the optional recovery dependency with `python -m pip install 'PyMuPDF>=1.24'` alongside the project's normal requirements, then run:

```bash
python tools/recover_dissertation_profiles.py --pdf /path/to/smtrevor_1.pdf \
  --qa-plot runs/recovery/figure_4_20_qa.png
```

The script rejects a different PDF or changed embedded image before using its manually inspected calibration. Digitization estimates are approximately +/-0.0262 dex horizontally and +/-0.01825 mm vertically; these do not represent experimental confidence intervals. The QA overlay checks the recovered samples against the solid central curves.

The seventeen Appendix A shot numbers are 15156, 15166, 15169, 15174, 15179, 15182, 15187, 15191, 15194, 15210, 15215, 15224, 15750, 15765, 15778, 15795 and 15802. Spectra at different gates on different shots require a shot-repeatability model before they can be interpreted as a temporal sequence.

## Recovered legacy simulation temperatures

The thermal benchmarks are published “Allegra” simulation profiles driven by **shot 17133**. They use a 200 um SS304 foil and 5 mm gap, and omit the dielectric support (section 4.5, printed p. 93; PDF p. 118). They are not measured temperatures or FLASH output. The separately recovered shot 15210 current below does not supply their drive.

| Figure and source page | Recovered artifact | Recovered quantity |
| --- | --- | --- |
| 4.28, printed 94 / PDF 119 | [CSV](../data/dissertation/figure_4_28_temperature.csv), [metadata](../data/dissertation/figure_4_28_temperature.json) | 38 interior temperature points each at the foil-edge and foil-middle labels, 40 ns |
| 4.29, printed 96 / PDF 121 | [CSV](../data/dissertation/figure_4_29_temperature.csv), [metadata](../data/dissertation/figure_4_29_temperature.json) | 38 interior temperature points each at the foil-edge and foil-middle labels, 60 ns |
| 4.30, printed 98 / PDF 123 | [CSV](../data/dissertation/figure_4_30_heating.csv), [metadata](../data/dissertation/figure_4_30_heating.json) | 137 width-dependent temperature points at 12.5 um depth, 40 ns |

Figures 4.28 and 4.29 explicitly label their temperature axes **K**. Their vertical coordinate is distance from the middle of the AK gap in mm. The reported initial geometry places the cathode gap face at +2.5 mm and its back face at +2.7 mm. Each CSV retains the absolute plotted coordinate and a separately named **nominal depth from the initial gap face**, covering approximately 7.13–194.44 um. The original mesh, material fractions and evolving material boundaries are unavailable. Boundary ramps, exterior zero segments and the legend region were excluded, so neither endpoint is a recovered surface value. The foil-edge label is tied to the schematic in Fig. 4.27; its exact transverse offset from the geometrical edge is unspecified.

At the retained point approximately 7.13 um from the initial face, the recovered edge/middle values are about 1212/540 K at 40 ns and 3225/1173 K at 60 ns. These are interior plot samples, not desorption boundary temperatures. Estimated coordinate readout envelopes are +/-22.4 K and +/-2.49 um, including the maximum tick-fit residual. They do not include physical uncertainty. The dashed “SS304 Melt” annotation maps to approximately 1651 K on these axes; it is recorded as a figure annotation rather than adopted as an independently sourced material property.

Figure 4.30 instead labels temperature in **eV**. The recovery interprets that label literally as `k_B T` in electronvolts; the original simulation temperature-variable definition and energy-mode ownership have not been independently recovered. Under that interpretation it converts to kelvin with `T_K = T_eV * e / k_B`, approximately 11604.518 K/eV, without a Celsius offset. The horizontal coordinate runs across the foil width, not the gap. Retained samples cover approximately -2.475 to +2.471 mm from the foil middle; the nearest-middle value is about 492 K. The largest retained near-edge value is about 2101 K, while the exact edge jumps and maxima remain excluded. Estimated readout envelopes are +/-29.8 K and +/-14.7 um laterally. These are correlated raster samples, not independent observations, surface values or a temperature history.

Reproduce the temperature recovery with:

```bash
python tools/recover_dissertation_depth_heating.py --pdf /path/to/smtrevor_1.pdf \
  --qa-dir runs/recovery
python tools/recover_dissertation_heating.py --pdf /path/to/smtrevor_1.pdf \
  --qa-plot runs/recovery/figure_4_30_heating_qa.png
```

Both tools check the source PDF and decoded native image hashes before applying fixed inspected calibrations. The sidecars retain pixel coordinates, units, exclusions, source pages and generator identity. No smoothing, interpolation or extrapolation is performed. The magnetic-field/current-density depth panels in Figs. 4.28 and 4.29, and the 17133 drive, remain unrecovered. The separate partial Fig. 4.30 electromagnetic recovery is described below.

## Recovered electromagnetic width profiles

The separate [Fig. 4.30 electromagnetic CSV](../data/dissertation/figure_4_30_electromagnetic.csv) and [metadata](../data/dissertation/figure_4_30_electromagnetic.json) recover **131 signed samples each** for the purple B_x and green J_z curves. They share the 40 ns, 12.5 um subsurface, shot 17133 legacy simulation context above. Retained points span approximately -2.438 to +2.434 mm across the foil width. Sharp boundary jumps, outermost edge regions, narrow curve-crossing neighborhoods and legend glyphs are excluded. No interpolation fills those gaps, and extrema describe retained samples rather than exact edge peaks. Fig. 4.30 remains partially recovered overall.

The primary recovered quantity is the signed plotted **Signal Height**. The source legends are preserved as **B_x [100 Tesla]** and **J_z [MA/cm²]**. Separate CSV columns make only the declared literal conversions: B_x ordinate times 100 gives tesla, and J_z ordinate times `10^10` gives A/m². Section 4.5 and the figure caption provide no explicit normalization recipe, and the original simulation-unit export is unavailable. These conversions remain unverified; no alternative area denominator, factor correction, sign flip or host-frame rotation is applied. The converted J_z values therefore **cannot support a quantitative J²/sigma heating benchmark until their units and normalization are confirmed**.

The readout envelopes include the visible exact-color line-run extent, an estimated +/-2 pixels and maximum tick-fit residual. They describe digitization rather than physical uncertainty. A width lineout at one depth does not determine the full current-density cross-section: no total-current integral, Ampere-law agreement, magnetic-divergence check or Joule-power closure is claimed. The recovered 15210 current is a different shot and supplies no normalization for these curves.

```bash
python tools/recover_dissertation_em_profiles.py --pdf /path/to/smtrevor_1.pdf \
  --qa-plot runs/recovery/figure_4_30_electromagnetic_qa.png
```

The tool validates the PDF and native pixel hashes before applying its fixed inspected calibration. Its output retains the raw ordinates, literal-legend SI interpretation, signed values, accepted source pixels and excluded crossing columns. The temperature artifacts remain separate.

## Recovered shot 15210 current

The [current CSV](../data/dissertation/shot_15210_current.csv) and [provenance record](../data/dissertation/shot_15210_current.json) recover the red “Average Current [MA]” curve from **Fig. C.11, printed p. 157 / PDF p. 182**. The current ramp in the same-shot Fig. A.10 clips near 0.24 MA, so it cannot supply the complete current history. C.11 contains native vector geometry, allowing 118 original encoded vertices inside the displayed axes to be recovered. Twenty vertices behind the later white legend are retained and marked `occluded_by_legend`; they are encoded source geometry rather than values interpolated across a hidden region.

The CSV provides `time_s` and `current_A`, alongside the original plotted ns/MA coordinates. Its finite time support is approximately 2.250–198.250 ns; the largest retained current is approximately 0.848 MA at 136.75 ns. The encoded endpoints outside the displayed axes were excluded. No zero-time or 200 ns endpoint, pre-pulse zero, baseline correction, smoothing or extrapolation was inserted. Small early negative currents, down to approximately -259 A, are preserved. These are published polyline vertices with unknown original filtering or plot simplification, not recovered oscilloscope samples.

The sidecar's display envelope is approximately +/-0.327 ns and +/-2000 A, derived from half the published line width plus maximum tick-fit residual. It is not diagnostic calibration uncertainty or physical current uncertainty. The plotted VUV rectangle spans approximately 52.055–55.055 ns. Table 4.1 reports acquisition time 55 ns and “After” breakdown; neither that entry nor the rectangle establishes the actual breakdown time, electronic gate-edge definition or trigger origin. Those remain unresolved.

Table 4.1, Fig. 4.23, Fig. A.10 and the C.11 plot title report 200 um SS304. The B.9 and C.11 captions report 500 um. The recovery records that conflict; current extraction cannot determine thickness. Shot 15210 is independent of the 17133-driven thermal profiles above.

```bash
python tools/recover_dissertation_drive.py --pdf /path/to/smtrevor_1.pdf \
  --qa-plot runs/recovery/shot_15210_current_qa.png
```

For downstream current queries, use the explicit [current loader](../mhdlab/drive_trace.py), which returns `DriveTrace`:

```python
from mhdlab.drive_trace import load_current_csv

trace = load_current_csv(
    "data/dissertation/shot_15210_current.csv",
    time_column="time_s",
    current_column="current_A",
    time_scale_to_s=1.0,
    current_scale_to_a=1.0,
    diagnostic_time_offset_s=0.0,
    source_kind="reconstructed",
    shot_id="Mykonos 15210",
    source_reference="Smith dissertation (2025), printed p. 157 / PDF p. 182",
    figure_reference="Figure C.11, Average Current [MA]",
)
```

The explicit zero offset preserves the published coordinate; it does not establish alignment to a physical trigger. `DriveTrace` preserves signed current, records the CSV hash, interpolates linearly only within the supplied support and rejects extrapolation. Its slopes describe that interpolation, not an independently measured `dI/dt`. The loader reads only the selected columns; consumers must separately inspect the sidecar's display bounds, legend flags and unresolved timing. Loading a trace does not install a FLASH boundary condition, infer a circuit voltage or restore the missing 17133 drive.

## Published estimates from Table 4.1

These values are author-derived Yacora estimates copied from printed p. 89 (PDF p. 114). They are not raw measured H/H2 densities or independently validated model inputs. The fits were assessed by eye (p. 90), and quantitative fit covariance is unavailable.

| Shot | Foil thickness in table, um | Gate relative to breakdown | Reported time, ns | Temperature, K | H2, cm^-3 | H, cm^-3 |
| --- | ---: | --- | ---: | ---: | ---: | ---: |
| 15750 | 200 | Before | 30 | 5000 | 1e14 | 1e15 |
| 15210 | 200 | After | 55 | 6600 | 1e16 | 1e15 |
| 15765 | 400 | Before | 29 | 5800 | 2.5e13 | 2e12 |
| 15795 | 400 | After | 180 | 8600 | 5e13 | 4e14 |
| 15778 | 500 | Before | 24 | 6000 | 6e13 | 1e12 |
| 15802 | 500 | After | 350 | 5000 | 5e13 | 6e14 |

Breakdown is defined through resistive voltage collapse, supported by visible onset (pp. 59–61). The approximately 50 ns value in the conclusions should not replace each shot's electrical trace and timing uncertainty.

## Existing modeling and source conflicts

The dissertation already reports preferential edge heating in its 2D “Allegra” simulations (section 4.5, pp. 93–98). These used a 200 um foil, a 5 mm gap and the current from a separate shot, 17133. Dielectric supports were omitted. Figures 4.28 and 4.29 show 40 ns and 60 ns profiles; Fig. 4.30 shows the width dependence at 12.5 um depth at 40 ns. The edge-heated region is described as hundreds of micrometres wide. Establishing hotter foil edges alone would repeat an existing result; the publication opportunity is quantifying its effect on the measured emission through a validated source and diagnostic reconstruction.

Several source inconsistencies must remain visible in any reconstruction:

- **Shot 15778 thickness:** Fig. 4.22, the Fig. 4.7 title, Table 4.1 and the embedded Fig. A.15 plot title say 500 um. The Fig. A.15 caption says 400 um. The manifest retains both.
- **Local density versus average density:** p. 97 says a shorter emitting depth raises the inferred local edge density by about an order of magnitude. The conclusion on p. 99 instead calls the inferred densities overestimated because the source is nonuniform. A forward model must distinguish column density, width-averaged volume density and local edge density before deciding the direction of any correction.
- **Density ranges:** the abstract's approximate H density ranges (printed p. xxiii; PDF p. 24) do not include the 1e15 cm^-3 entries for shots 15750 and 15210 in Table 4.1. Keep the explicit table values and the discrepancy.
- **Simulation temperature wording:** p. 93 calls the foil-middle surface temperature 1600 C at 60 ns, while p. 95 says that location remains well below melt. The recovered interior profiles use the explicit Kelvin axes; they do not recover actual surface locations or resolve that wording conflict.
- **Repeated Appendix B and C captions:** for example, shot 15174 is 100 um in Fig. 4.12 and Fig. A.4 but 500 um in Fig. B.3. For shot 15210, Table 4.1 and the C.11 plot title report 200 um, while B.9 and C.11 captions report 500 um. Verify material and thickness shot by shot.

The original CR comparison also has explicit limits. The public Yacora output used lacked atomic angular-momentum sublevels and molecular rotational populations, and the dissertation describes assigning whole-level populations to individual transitions, limiting the result to general agreement (pp. 91–92). Its statement that every H atom in n=2 must emit Lyman-alpha overlooks the separate 2s and 2p populations and two-photon/collisional pathways. A state-resolved model must use independently sourced population and transition rates rather than adopt that statement as a closure.

## Recovery limits

Publication reconstruction can recover geometry, shot references, approximate displayed curves, table estimates, visible feature identification and qualitative MHD comparisons. It cannot restore the original camera frames, dark/flat corrections, exact gain settings, response functions, scope waveforms, trigger/cable-delay calibration, Yacora decks, or Allegra mesh/material/output files from the figures alone. These are recorded as missing evidence in the manifest. Any digitized observation should retain its source figure, axis mapping, extraction method and reconstruction uncertainty, and be labeled as reconstructed published data.
