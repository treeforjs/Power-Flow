# Manufacturer follow-up on absolute VUV calibration

Public manufacturer documents clarify several quantity definitions, but they
do not identify an absolute correction for the dissertation spectra. In
particular, the legacy McPherson lamp datasheet itself uses inconsistent
radiance/intensity terminology. The original lamp certificate, reduction files
and shot-specific camera records remain necessary. This follow-up supplements
the [dissertation calibration audit](absolute_spectral_calibration.md); it does
not change the recovered curves or install a calibration factor.

The [source manifest](../data/dissertation/calibration_primary_source_manifest.json)
records official URLs, retrieval times, actual response sizes and SHA-256
digests, HTTP/header metadata and pinpoint findings. Manufacturer PDFs, HTML,
extracted article text and pixel inspections remain in the external local cache
and are excluded from the public repository and delivery archive. These are
generic product documents, not certificates for the experimental units.

## McPherson 632 lamp: the area convention remains unresolved

The official [legacy Model 632 datasheet](https://www.mcphersoninc.com/pdf/632.pdf)
is still accessible. Page 1 offers a source “available as a calibrated standard
of spectral radiance,” with accompanying documentation and measurements every
nanometre from 115 to 400 nm. It describes an approximately 1-mm emitting
diameter, an MgF2 window and an approximately f/7 output. These nominal dimensions
do not establish the experimental source-to-slit geometry or calibration plane.

Crucially, its page-2 sample graph labels the ordinate
**“Spectral Radiant Intensity, µW / (sr nm)”**, with no area denominator. The label
is raster artwork and was visually inspected; it is not present in extracted
PDF text. The graph is not a supplied serial-number certificate, and the
datasheet does not reconcile this label with its page-1 radiance terminology.
Neither wording establishes whether the dissertation's literal per-mm ordinate
is a missing exponent, a line normalization, a conversion error or another
convention. No spectral amplitude is adopted from the generic graph.

The current [model632.html page](https://www.mcphersoninc.com/lightsources/model632.html)
now describes **Model 634**, with approximately f/6 output, and links a
[634 datasheet](https://www.mcphersoninc.com/pdf/634.pdf). This is separate family
information. It must not replace the legacy 632 geometry, spectral output or
calibration identity.

For a certificate specifying spectral radiant intensity in W/(sr·nm), collected
power uses its accepted solid angle. For spectral radiance in W/(m²·sr·nm), it
uses the projected area–solid-angle integral, or étendue. Selecting one of these
relations requires the actual certificate quantity and reference plane. The
reported `0.037 mm²·sr` does not settle that choice or identify the Mykonos relay
throughput. A nominal f-number and emitting diameter are insufficient to repair
the missing geometry.

## Spectrometer: dispersion, diffraction efficiency and throughput differ

The official [legacy 234/302 datasheet](https://www.mcphersoninc.com/pdf/234302.pdf)
gives a 200-mm focal length, f/4.5 aperture, 32-degree incidence and settable slit
height on page 1. Page 2 lists 4 nm/mm dispersion for the 1200-grooves/mm grating.
Its 0.1-nm FWHM specification was tested in scanning mode at 184.9 or 312.6 nm
with 10-µm slits; it does not supply the experimental VUV instrumental line
profile or convert the dissertation's approximate nm/pixel into a binning rule.

Page 2 also contains a theoretical grating-selection graph with a
1200-grooves/mm Al+MgF2 curve. Its ordinate literally says
**“Theor. Efficiency (percent)”**, while ticks run from 0 to 0.4. That labeling
does not unambiguously specify fraction versus percent. No numerical efficiency
table is recovered or adopted. Page 3's sample spectrum and the current product
page describe relative response, not a calibrated whole-train throughput.

Coating reflectance, useful-order diffraction efficiency and the complete
spectrometer/relay response are different quantities. A coating-reflectance
estimate cannot be substituted for diffraction efficiency. A complete empirical
lamp-to-count calibration may already contain grating and detector losses;
multiplying them again would double-count those factors. The current 234/302 URL
also advertises the newer Monarch, so current specifications must be kept
distinct from the legacy instrument and its actual grating.

## Photek: QE and gate specifications do not provide ADU normalization

[Photek AN002, revision 02, May 2021](https://www.photek.com/wp-content/uploads/2021/06/PH-AN002-ImageIntensifier-Rev02-May21.pdf)
defines photocathode QE using photons incident on the photocathode (page 2).
It treats MCP electron amplification and phosphor conversion as separate stages
(pages 3–4), with phosphor efficiency measured in photons/electron. Page 2 lists
MgF2 input-window options and explains that the window affects short-wavelength
response. Generic S20 or window curves do not identify the reference plane or
response of the historical custom assembly. QE alone is not ADU per photon.

The [current iCMOS 160 page](https://www.photek.com/icmos-camera/) distinguishes
standard 50-ns and optional high-speed 3-ns minimum gates. It also distinguishes
one gate trigger per frame from multiple gate triggers per frame. AN002 page 6
likewise describes integrating multiple excitations into one frame. Consequently,
frame duration does not independently determine the effective optical exposure.
These family specifications do not recover actual gate waveforms, gate counts,
gain, linearity, saturation, background or sensor-to-ADU conversion for any shot.

The missing evidence remains the actual lamp certificate and its measurement
geometry; the original lamp count/gain/exposure reduction; complete wavelength
throughput with factor ownership; and shot-specific acquisition settings.
Absolute response therefore remains an explicit unresolved scale in coverage
inference. Generic manufacturer curves cannot determine a unique monolayer
coverage or replace that experimental evidence.
