# Dissertation data reconstruction

Trevor Johannes Smith's 2025 dissertation contains gated VUV emission spectra, spatial lineouts across the anode–cathode gap, electrical traces and selected simulation profiles that can support a publication reconstruction. These are published observations and derived estimates. Original detector data, original electrical waveforms and complete per-shot calibration have not been recovered.

The source is [Characterizing Desorbed Surface Contaminants Leading to Plasma Formation in Power Feeds of Pulsed Power Accelerators](https://deepblue.lib.umich.edu/items/9ccfa886-0b19-4bd8-b3b2-bb685b9ddefc). The [current repository PDF endpoint](https://backend.production.deepblue-documents.lib.umich.edu/server/api/core/bitstreams/24e682a8-896a-410d-85bb-ca4654a9db02/content) supplied the validated 203-page PDF. Its SHA-256 is `1d7b1f9c191a37e94f1cc3634587750f1af35b07eb13e0b94c60f3da47495ed7`. For its Arabic numbered pages, add 25 to the printed page to get the one-based PDF page. The machine-readable [recovery manifest](../data/dissertation/recovery_manifest.yaml) records source provenance, 43 figure/table targets, published table entries and unresolved metadata. The PDF and extracted images remain outside this repository.

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
| Figs. 4.27–4.30 | 93–98 | Existing simulation geometry and edge/middle field, current and temperature profiles |

Most published spectral figures are embedded raster images. Digitization can recover displayed curves with finite plot resolution, but the pseudocolor image cannot be treated as original detector counts. Saturation, clipping, interpolation, averaging and plotting transformations may be irreversible. The first partial recovery is [Fig. 4.20's central gap profiles](../data/dissertation/figure_4_20_profiles.csv): 106 correlated raster samples, 53 each for H Lyman-alpha and the H2 161 nm Lyman band. The [provenance record](../data/dissertation/figure_4_20_profiles.json) includes source hashes, pixel coordinates, axis calibration, and pixel-reading uncertainty. These are plotted transition densities, not total ground-state neutral densities. Physical error bands and wavelength lineouts remain pending. Other curves retain `not_digitized` status; numerical table transcription is marked separately.

To reproduce this recovery with the exact source PDF, install the optional recovery dependency with `python -m pip install 'PyMuPDF>=1.24'` alongside the project's normal requirements, then run:

```bash
python tools/recover_dissertation_profiles.py --pdf /path/to/smtrevor_1.pdf \
  --qa-plot runs/recovery/figure_4_20_qa.png
```

The script rejects a different PDF or changed embedded image before using its manually inspected calibration. Digitization estimates are approximately +/-0.0262 dex horizontally and +/-0.01825 mm vertically; these do not represent experimental confidence intervals. The QA overlay checks the recovered samples against the solid central curves.

The seventeen Appendix A shot numbers are 15156, 15166, 15169, 15174, 15179, 15182, 15187, 15191, 15194, 15210, 15215, 15224, 15750, 15765, 15778, 15795 and 15802. Spectra at different gates on different shots require a shot-repeatability model before they can be interpreted as a temporal sequence.

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
- **Simulation temperature wording:** p. 93 calls the foil-middle surface temperature 1600 C at 60 ns, while p. 95 says that location remains well below melt. Numerical comparison needs the actual plotted profile and its position rather than either sentence alone.
- **Repeated Appendix B captions:** for example, shot 15174 is 100 um in Fig. 4.12 and Fig. A.4 but 500 um in Fig. B.3. Several Appendix B captions repeat 500 um SS304 across unlike shots. Verify material and thickness shot by shot.

The original CR comparison also has explicit limits. The public Yacora output used lacked atomic angular-momentum sublevels and molecular rotational populations, and the dissertation describes assigning whole-level populations to individual transitions, limiting the result to general agreement (pp. 91–92). Its statement that every H atom in n=2 must emit Lyman-alpha overlooks the separate 2s and 2p populations and two-photon/collisional pathways. A state-resolved model must use independently sourced population and transition rates rather than adopt that statement as a closure.

## Recovery limits

Publication reconstruction can recover geometry, shot references, approximate displayed curves, table estimates, visible feature identification and qualitative MHD comparisons. It cannot restore the original camera frames, dark/flat corrections, exact gain settings, response functions, scope waveforms, trigger/cable-delay calibration, Yacora decks, or Allegra mesh/material/output files from the figures alone. These are recorded as missing evidence in the manifest. Any digitized observation should retain its source figure, axis mapping, extraction method and reconstruction uncertainty, and be labeled as reconstructed published data.
