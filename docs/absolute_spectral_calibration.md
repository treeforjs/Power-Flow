# Absolute spectral calibration: recovered evidence and unresolved scale

The dissertation's lamp and gain figures are now recoverable as encoded vector
display data. They do **not** yet establish a numerical correction to absolute
VUV intensity. The literal ordinate unit, missing detector-response definitions
and unreported relay throughput prevent selecting a justified correction factor.
No squared-millimetre substitution or count-to-photon factor is installed here.

The source is Trevor Johannes Smith's 2025 dissertation, [official DeepBlue
PDF](https://backend.production.deepblue-documents.lib.umich.edu/server/api/core/bitstreams/24e682a8-896a-410d-85bb-ca4654a9db02/content),
SHA256 `1d7b1f9c191a37e94f1cc3634587750f1af35b07eb13e0b94c60f3da47495ed7`.
Printed pages 79–83 correspond to one-based PDF pages 104–108. The narrowly
related intensity conversion is Eq. 4.8, printed p. 92 / PDF p. 117. Original
camera frames, lamp certificate and calibration files have not been recovered.

## Recovered published curves

| Source | Artifact | Content |
| --- | --- | --- |
| Fig. 4.18, printed 80 / PDF 105 | [CSV](../data/dissertation/figure_4_18_calibration.csv), [sidecar](../data/dissertation/figure_4_18_calibration.json) | 243 original encoded lamp-curve vertices; approximately 110.094–169.812 nm |
| Fig. 4.19, printed 81 / PDF 106 | [CSV](../data/dissertation/figure_4_19_gain.csv), [sidecar](../data/dissertation/figure_4_19_gain.json) | 101 encoded fit vertices, 36 measured-display marker centres and 203 polygon-boundary vertices including repeated closure |

Both figures literally label the ordinate **Intensity [µW/(sr·nm·mm)]**, with
one `mm`. The CSV quantity is therefore named `plotted_intensity`. It is not
labelled SI radiance, detector counts, a certified lamp output or corrected
experimental intensity. The Fig. 4.18 caption describes 100–200 nm, while its
displayed tick range is 110–170 nm and its retained vertices cover the narrower
range above. No missing wavelengths are extrapolated.

The lamp polyline also contains one encoded backstep: vertex 159 to 160 moves
from approximately 148.934146 to 148.839969 nm. Its original vertex order and
both ordinate values are retained and flagged in the sidecar. The recovery is
not an interpolation-ready single-valued response table. Sorting, dropping or
combining those vertices would be a separate reviewed modeling decision.

The gain file's `series` distinguishes `exponential_fit`, `measured_marker` and
`fit_error_boundary`; `record_kind` distinguishes original vertices from derived
centres of the circular glyphs. The red legend marker is excluded. Measured
centres map to approximately gain settings 20–55. The published fit extends to
approximately 0–100 and includes extrapolation. Native rounding and tick fitting
leave small departures from integer gain settings; no values are snapped to
integers or endpoints. `gain_region` uses the display envelope to identify
extrapolated portions. The teal polygon is preserved as the literal **Fit Error**
boundary; its confidence level, probability interpretation and construction are
unreported. It must not become a pointwise experimental error distribution.

The [recovery tool](../tools/recover_dissertation_calibration.py) guards the PDF,
page content streams, complete native drawing sets and selected path hashes.
It fits linear tick maps for wavelength/gain and intensity or log10 intensity.
Least-squares coefficients are evaluated with exact rational arithmetic on the
native binary tick coordinates and rounded once to float64, avoiding library
linear-algebra differences in serialized coordinates.
It retains encoded geometry without smoothing, resampling or refitting.
Seventeen-digit CSV coordinate serialization preserves extracted geometry;
it does not imply that precision in the calibration or physical observations.
Successive vertices are correlated display information. Display envelopes use
half the published stroke width or marker radius plus the largest tick-fit
residual. They describe glyph extent, not centre uncertainty or physical errors.

```sh
python tools/recover_dissertation_calibration.py --pdf /path/to/smtrevor_1.pdf \
  --output-dir data/dissertation --qa-dir /tmp/dissertation-calibration-qa
```

Recovery needs the optional PyMuPDF dependency. Optional QA overlays contain
source pixels and remain local, outside the public delivery archive; sidecars
record their hashes. Portable tests use committed facts and synthetic geometry,
with exact-PDF integration skipped when that external cache is unavailable.

## Documented acquisition and geometry

| Quantity | Source statement | Reference |
| --- | --- | --- |
| Spectrometer | McPherson 234/302; 1200 grooves/mm Al+MgF2 grating; approximately 115–165 nm band | p. 79 |
| Camera | Photek iCMOS 160, S20 intensifier and MgF2 window | p. 79 |
| Shot exposure | One frame integrated over 3 ns | p. 79 |
| Dispersion | Approximately 0.03 nm/pixel | p. 79 |
| Grating reflectance | Approximately 75–85% above 115 nm | p. 79 |
| Camera QE | 19–23% over 115–180 nm, with the stated intensifier/window | p. 79 |
| Lamp | McPherson 632 D2, affixed to the spectrometer entrance slit | pp. 79–80 |
| Lamp slit width | 50 µm | pp. 80, 82 |
| Gain measurements | Settings 20–55; gates 3, 5, 10 and 100 ns | p. 80 |
| Gain limits | Below 20 near noise; 55 approaches pixel saturation at 3 ns | p. 80 |
| Lamp étendue | 0.037 mm²·sr for the lamp–50 µm slit configuration | p. 82 |
| Mykonos relay | Two 100-mm, diameter 1-inch Al+MgF2 concave mirrors and an added MgF2 window | p. 82 |
| Relay response leads | Mirror reflectance from eSource Optics [68], window transmittance from ThorLabs [69] | p. 82 |

The gain plot does not identify which gate produced each marker, how exposure
was normalized, which spectral/spatial region was reduced, or the detector-count
definition. Its increasing ordinate is labelled with the same physical unit as
the fixed lamp curve. Thus its detector-response interpretation requires the
original procedure; the plotted ordinate cannot simply be treated as ADU or as
a gain multiplier. The original exponential coefficients, offset convention and
fit covariance are not supplied. A curve recovered beyond the measured gain
range remains an extrapolated display, not a validated camera response.

## Étendue and dimensional checks

The dissertation defines

\[
G=\int_S\int_\Omega n^2\cos\theta\,dS\,d\Omega \tag{4.4}
\]

and relates radiant flux and radiance by

\[
L_{e,\Omega}=n^2\frac{\partial\Phi_e}{\partial G}. \tag{4.5}
\]

Here vacuum has `n=1`. Ordinary spectral energy radiance has an **area**
denominator. Multiplying the figures' literal per-`mm` ordinate by the stated
`mm²·sr` étendue does not produce spectral radiant power. A missing exponent in
the label, a line-source normalization, a scale error or another undocumented
convention cannot be distinguished from these figures. The lamp certificate and
original reduction must resolve it; no particular replacement factor follows.

Equation 4.6 gives the centered, on-axis rectangular solid angle:

\[
\Omega=4\arctan\left[\frac{hw}{2d\sqrt{4d^2+h^2+w^2}}\right].
\]

Equation 4.7 prints `Ω=πr²/d²` for a circle without a stated approximation.
That is a small-angle expression. For a point on the axis of a circular aperture,
the exact solid angle is `2π[1-d/sqrt(d²+r²)]`. Finite emitting areas and projected
angular factors still require the étendue integral. The inspected section does
not supply the lamp slit height, emitter dimensions and distances needed to
recompute the reported value, so this observation does not quantify an actual
experimental correction.

The 0.037 mm²·sr value must not be assigned to the Mykonos relay. The text names
a different limiting first-mirror/window path but gives no numerical relay
étendue or per-spatial-pixel allocation. It calls the resulting calibration
approximate. The continuation on p. 83 only completes the statement about photons
within the collection volume; it adds no throughput or detector conversion.

## Photon, energy, bandwidth and exposure accounting

A dimensional audit can use the following identity **after** its quantities and
factor ownership have been established. For energy spectral radiance per nm,
the photons reaching a spatial/spectral detector bin are

\[
N_{\gamma,p}=\int_{\rm gate}\!dt\int_{\rm bin}\!d\lambda_{\rm nm}\,
L_{e,\lambda}\,G_p(\lambda)\,\mathcal T_{\rm extra}(\lambda)
\frac{\lambda_{\rm m}}{hc}.
\]

This identity is an audit convention, not a recovered dissertation calibration
equation. It requires `L_e,lambda` in W/(m²·sr·nm), `G_p` in m²·sr, time in
seconds and wavelength in metres in the photon-energy factor. Photon spectral
radiance instead already counts photons: multiplying by `lambda/(hc)` again
would be wrong. A per-nm ordinate needs integration over the actual spectral
bin; a line-integrated ordinate already contains that bandwidth. Similarly,
integrated camera counts and per-second rates differ by the effective gate
exposure. The stated 0.03 nm/pixel is approximate dispersion, not a recovered
bandwidth integration rule or instrumental line-spread function. The listed
lamp gates do not establish a count-to-rate normalization.

Camera QE converts incident photons to detected photoelectrons. Intensifier
response, camera conversion to ADU, dark/background/stray-light subtraction,
binning, saturation and effective gate shape require their own definitions.
A calibrated count-to-photon response inferred through the complete lamp setup
may already include the grating, camera QE and camera window. Multiplying those
losses again could double-count them. Conversely, the Mykonos relay adds two
reflections and a separate window absent from that lamp setup. A product
`R1(lambda)*R2(lambda)*T_added_window(lambda)` is a plausible optical model to
verify, not a recovered multiplication rule. The text provides no numerical
curves or explicit powers. Percent versus unit-fraction conventions must also be
specified before numerical use.

Equation 4.8, p. 92, is
`I_ik=A_ik*h*nu_ik*n_k*L/(4*pi)`, with emitting depth set to foil width, 5 mm.
The displayed unit includes nm^-1, while the surrounding text describes total
line intensity represented by a delta function. Dimensional line-integrated
energy radiance from that expression has no per-nm factor until a normalized
spectral profile or bin rule is supplied. The later line broadening is not a
detector calibration. This unresolved distinction prevents inferring a missing
bandwidth correction from Eq. 4.8 alone.

The dimensionally consistent separation, under the adopted emission-model
assumptions, is

\[
I_{\rm line}=\frac{A_{ik}h\nu_{ik}}{4\pi}\int n_k\,d\ell
\quad[\mathrm{W\,m^{-2}\,sr^{-1}}],\qquad
I_{\lambda}=I_{\rm line}\,\phi_{\lambda},\qquad
\int\phi_{\lambda}\,d\lambda_{\rm nm}=1.
\]

Here `phi_lambda` has units nm^-1 and can include a separately normalized
instrumental line profile. A spectral-bin mean is
`Ibar_p = integral_bin(I_lambda d_lambda_nm)/Delta_lambda_nm`; its photon
conversion applies `lambda_m/(h*c)` inside the unit-consistent bin integral.
Replacing that wavelength dependence by a bin-centre value requires an explicit
approximation, or the appropriate bin average under a stated spectral model.
This corrects the equation/unit distinction. It does not determine a new
absolute experimental amplitude factor.

## What remains unresolved for inference

The required evidence is a dimensionally identified lamp curve/certificate;
original detector response and background/gain/exposure normalization; numerical
relay geometry and wavelength throughput; spatial/spectral bin definitions;
actual shot gain/gate settings and saturation status; and line-profile/bin
normalization. The recovered curves constrain what the dissertation displayed,
but they cannot reconstruct clipped counts or uniquely determine these factors.
Absolute-response nuisance factors and their priors must remain explicit in
spectral inference until this evidence is obtained. Converting fitted source
inventory to monolayers separately requires a declared site-density convention;
an instrumental normalization cannot establish that surface quantity.
