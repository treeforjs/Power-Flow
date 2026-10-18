# Observed spectrum recovery and VUV eligibility

No detector-count, arbitrary-intensity or unambiguously radiance-valued **VUV**
lineout was identified in the bounded candidate check below. The published
**visible** SVS lineout in Fig. 4.15 provides a separate eligible observation:
[136 partial display samples](../data/dissertation/observed_spectrum_4_15.csv)
with a [provenance sidecar](../data/dissertation/observed_spectrum_4_15.json).
Its literal quantity remains `Intensity [arb.]`; it supplies no absolute
VUV amplitude or raw detector counts. No provisional Fig. 4.21 artifact was
created: those experimental-origin curves have a derived density ordinate.

The source is Trevor Johannes Smith's 2025 dissertation, [official DeepBlue
PDF](https://backend.production.deepblue-documents.lib.umich.edu/server/api/core/bitstreams/24e682a8-896a-410d-85bb-ca4654a9db02/content),
SHA256 `1d7b1f9c191a37e94f1cc3634587750f1af35b07eb13e0b94c60f3da47495ed7`.
It contains 203 pages; Arabic printed page plus 25 gives the one-based PDF page.
The existing [reconstruction guide](dissertation_recovery.md) and
[manifest](../data/dissertation/recovery_manifest.yaml) already distinguish
published display data from unavailable raw frames. This check does not change
their pending recovery statuses.

## Inspected VUV candidates

| Candidate | Printed / PDF page; native image | Caption and literal axis evidence | Outcome |
| --- | --- | --- | --- |
| Fig. 4.20, Mykonos 15750, SS304 200 µm | 84 / 109; xref 2707, 1040×692 RGB | Caption identifies three experimental spectra; wavelength lineout ordinate is `Transition Density [Emitters/cm³]` | Experimental-origin derived display, not the requested observation quantity |
| Fig. 4.21, Mykonos 15765, SS304 400 µm | 85 / 110; xref 2721, 1040×692 RGB | Caption calls them “three experimental spectra”; lower-right lineouts use `Transition Density [Emitters/cm³]` | Do not reinterpret as counts, arbitrary intensity or radiance |
| Fig. 4.26, Mykonos 15802 | 91 / 116; xref 2807, 1600×1200 RGB | Top/Middle/Bottom measured-origin lineouts share `Transition Density [Emitters/cm³]` with model curves; caption says simulations are normalized to the 161-nm height of the 6000-K simulation | Not an alternative detector observable; model-referenced comparison remains separate |
| Fig. A.14, same Mykonos 15765 | 115 / 140; xref 3033, 1024×512 RGB | Caption identifies top, middle and bottom spectral lineouts; wavelength lineout ordinate remains `Transition Density [Emitters/cm³]` | The same-shot appendix does not recover a pre-conversion signal |

All four candidates are embedded rasters, not vector spectral paths. Their
pseudocolor panels cannot be inverted into raw camera counts. On printed p. 83
/ PDF p. 108, the narrative calls the main figures “intensity corrected photon
densities,” while the plotted axes call the lineout quantity transition density.
That terminology difference is unresolved. It does not justify relabelling a
published ordinate or applying a new amplitude factor.

Fig. 4.26 also has a source inconsistency: its visible model legend names
5000, 6000 and 7000 K, while the preceding discussion on p. 90 describes changes
at 4000 and 6000 K around an estimated 5000 K. Its caption specifies normalization
to the 6000-K result. No model-curve reassignment is made here.

## Bounded check of other figures

The figure list and relevant Chapter 3 candidates were checked to avoid missing
an observed VUV trace with a simpler ordinate:

- Fig. 3.10, printed 46 / PDF 71, is observed XUV **imaging** from MAIZE 2317,
  with a current/timing trace. It is not a wavelength-resolved spectrum.
- Fig. 3.12, printed 48 / PDF 73, contains an observed visible spectrum from
  MAIZE 2641, 50 µm Al, with an average lineout labelled `Intensity (arbitrary)`.
  Its approximately 400–510 nm support is visible light, not the requested VUV
  observable or the Mykonos SS304 experiment.
- Figs. 3.13 and 3.14, printed 50–51 / PDF 75–76, explicitly describe
  “PrismSPECT simulation results.” Their VUV portions are model outputs.
- Figs. 4.18–4.19 are lamp/camera calibration displays, not experimental MITL
  spectra; their [separate recovery](absolute_spectral_calibration.md) retains
  the unresolved literal unit and response definitions.

This is a bounded eligibility result for the inspected candidates, not a claim
that every dissertation image or an unavailable original file has been searched.
No model trace, spectral pseudocolor inversion or derived transition-density
lineout is substituted for the missing requested observation.

## Required VUV observation record

A reduction owner must identify the actual observation quantity and its mapping
to the displayed curves: original counts or a recoverable radiance/intensity
trace; spectral and spatial averaging support; gate/gain and background
conventions; and the density/volume or transition conversion already applied.
The original reduction record or a defensible inverse mapping is required before
the published densities can be used as a different observable. No covariance,
absolute coverage or monolayer estimate follows from display digitization.

## Recovered visible observation

Fig. 4.15, printed 75 / PDF 100, is an observed streaked visible spectrum from
Mykonos 15765 with a 400 µm SS304 foil. Its native image is xref 2554,
1040×1320 RGB. The middle linear panel shows measured SVS lineouts at 65, 95
and 330 ns. Only the green `SVS Lineout 3: 330 ns` trace is selected. The
lower semilog panel and the upper pseudocolor panel are not additional samples.
No species assignments or new physical parameters are extracted.

The [recovery tool](../tools/recover_dissertation_observed_spectrum.py) requires
the exact PDF hash above, the page-content SHA256
`064dfb491db122b454c86a3ee2a09e74fe40b3bdeb1c4928247969c180ae5442`, and the
native RGB-pixel SHA256
`27a9beba2c5aa9edb0de720ff0f50866157eb845b12c1318864279f1b16d5276`.
The native image was visually inspected before selecting calibration anchors;
a source-local overlay was also inspected after recovery. No source image is
distributed with the CSV or tool.

Nine printed wavelength ticks, 350–750 nm, and the three printed intensity
ticks, 50, 100 and 150, calibrate the native pixel axes. Fixed black-tick
regions and their detected centres are recorded in the sidecar. An exact
rational least-squares calculation maps these tick centres to the printed
numbers, then rounds the two coefficients once to float64. The largest
intensity tick residual is about 0.202 arbitrary units. This is a common axis
readout limitation, not independent noise assigned to each sample.

The tool considers 156 native columns, x=220 through 840 at a stride of four
pixels, in y=516 through 788. It uses only literal green RGB `(60,179,113)`;
antialiased colors and other curves are excluded. The legend rectangle is
masked. A column is retained only if it contains one contiguous green run of
3–25 pixels that does not touch the selection boundary. The median row of the
run gives the displayed centre. This yields 136 retained points spanning about
380.4–713.7 nm and 20 omitted columns with broad projected strokes. Every
omitted column and reason is preserved. Missing columns are not interpolated;
several sharp features and their peak maxima are therefore incomplete. No
line-strength integral is justified by this partial recovery.

The CSV includes wavelength and intensity **readout envelopes**. Wavelength
uses one native pixel on each side plus the maximum tick-fit residual.
Intensity uses the full selected green run padded by one pixel at each end,
plus the maximum tick-fit residual. These are glyph extents for assessing
digitization precision. They are not independent measurement errors,
confidence intervals or a covariance model. The plotted line is already
correlated display data; four-pixel column spacing is not a detector spectral
bin width or instrumental resolution.

The source label `330 ns` is retained without inferring a measured trigger or
breakdown time origin. The visible lineout exposure/averaging width, raw
counts, background, gain and covariance are unavailable in this recovery.
The same shot's Table 4.1 VUV entry, 29 ns before breakdown with the 3-ns VUV
diagnostic gate, is separate context. It is not the visible lineout interval
or a surface-emission gate. The SVS side-on gap collection region also differs
from the VUV collection. None of those records establishes a common amplitude
calibration or justifies an absolute monolayer regression from this trace.

The second-campaign description on printed p. 71 / PDF p. 96 reports a focal
point 450 µm from the cathode, a 175-µm spot radius and collection extending to
275 µm from the cathode. It reports a 70.3-mm primary lens focal length and
object distance, a 40.2-mm distance associated with the second focusing optic,
and a demagnification of 0.57. Fig. 4.13 on printed p. 72 / PDF p. 97 labels
this setup for shots 15750–15802. These reported geometry values are retained
in the sidecar. They do not define a top-hat collection function, an étendue,
a two-dimensional response kernel or the actual aperture acceptance. The
printed p. 74 discussion also says the VUV system started recording shot 15765
before the SVS system started recording. This states acquisition order, not
emission onset, and supplies no new mapping of the 330-ns lineout to breakdown.

With the original PDF held outside the repository, reproduce the public CSV
and sidecar using:

```sh
python tools/recover_dissertation_observed_spectrum.py --pdf /path/to/smtrevor_1.pdf
```

`--csv` and `--metadata` select output paths. `--qa-plot` optionally writes a
source-pixel overlay and must point outside the repository and any delivery
archive. The source-local QA hash is recorded when that option is used;
metadata without this optional QA entry still describes the same CSV.
The [focused tests](../tests/test_dissertation_observed_spectrum.py) check
pixel selection, omission, axis mapping, artifact semantics and source guards.
The exact-PDF integration check can be skipped where the external source is
unavailable; it is not replaced by a redistributed PDF fixture.
