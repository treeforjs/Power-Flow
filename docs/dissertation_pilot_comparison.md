# Dissertation pilot observations and reconstruction scenarios

The [versioned comparison manifest](../data/dissertation/pilot_comparison_manifest.json)
binds the pilot observations, their literal units, source identities, calibration
scenarios and current admission status. Shot **15750** is the pilot; **15765** is
the first additional comparison. Shot **17133** remains a distinct heating
benchmark. Its drive is not a substitute for either VUV shot's drive.

## Available constraints

| Shot and figure | Recovered artifact | Literal quantity and support |
| --- | --- | --- |
| 15750, Fig. 4.20 | [Gap profiles](../data/dissertation/figure_4_20_profiles.csv) | 106 correlated display samples of H Lyman-alpha and H2 approximately 161-nm central gap profiles; transition density in emitters/cm^3 |
| 15750, Fig. 4.20 | [Wavelength lineouts](../data/dissertation/figure_4_20_wavelength_lineouts.csv), [sidecar](../data/dissertation/figure_4_20_wavelength_lineouts.json) | 157 retained samples: top 40, middle 38, bottom 79; transition density in emitters/cm^3 |
| 15765, Fig. 4.21 | [Wavelength lineouts](../data/dissertation/figure_4_21_wavelength_lineouts.csv), [sidecar](../data/dissertation/figure_4_21_wavelength_lineouts.json) | 183 retained samples: top 46, middle 54, bottom 83; transition density in emitters/cm^3 |
| 15765, Fig. 4.15 | [Visible lineout](../data/dissertation/observed_spectrum_4_15.csv) | 136 partial samples, literal arbitrary intensity at the displayed 330-ns label; separate acquisition and collection region |

The new VUV lineouts recover the experimental-origin curves, excluding the
Yacora overplot. They do not invert pseudocolor panels. Original raw frames are
unavailable. The derived ordinate cannot yet be mapped uniquely to counts,
radiance, a single upper-state density, or ground-state neutral density.

The dissertation is Trevor Johannes Smith's 2025
[Characterizing Desorbed Surface Contaminants Leading to Plasma Formation in
Power Feeds of Pulsed Power Accelerators](https://deepblue.lib.umich.edu/items/9ccfa886-0b19-4bd8-b3b2-bb685b9ddefc).
The exact PDF SHA256 is
`1d7b1f9c191a37e94f1cc3634587750f1af35b07eb13e0b94c60f3da47495ed7`.
Figures 4.20 and 4.21 are on printed pages 84 and 85 / one-based PDF pages 109
and 110. Their source-page and native RGB identities are independently guarded.

## Raster readout contract

The [recovery tool](../tools/recover_dissertation_vuv_lineouts.py) uses the native
1040-by-692 RGB images, with no image scaling. The lower-right wavelength axes
are calibrated against five inspected major ticks and four density ticks.
Coefficients are fitted in exact rational arithmetic and rounded once.

Only exact olive, firebrick and slate-blue curve colors are admitted. Every fourth
native column is considered; contiguous strokes of one through eight exact-color
pixels are retained. Legend-covered, separated, wide and boundary-touching runs
are omitted. The sidecars retain every omitted column and its reason. The
`support_segment_id` changes whenever an omitted sampled column interrupts a
series. It does not declare detector bins or independent resolution elements.

The envelope includes the full selected stroke, one extra native pixel, and the
maximum tick-fit residual. These are correlated glyph readout bounds, not
measurement noise or confidence intervals. Source-local overlays for both
figures were inspected; source pixels remain outside repository artifacts.
Missing segments and partially hidden peaks prevent complete line integrals.

With the exact PDF retained outside version control:

```sh
python tools/recover_dissertation_vuv_lineouts.py --pdf /path/to/smtrevor_1.pdf --figure 4.20
python tools/recover_dissertation_vuv_lineouts.py --pdf /path/to/smtrevor_1.pdf --figure 4.21
python tools/validate_pilot_comparison.py --report runs/dissertation_pilot_readiness.json
```

`--qa-plot` optionally creates a local source overlay; retain it outside public
artifacts. Validation checks 14 referenced artifact hashes. It returns success
for a valid evidence record even when physical inference is unavailable.
`--require-physical-ready` instead returns exit code 2 in that situation.

## Calibration and reduction scenarios

The manifest keeps dimensional interpretation, reduction interpretation and
wavelength response as separate categorical assumptions:

- **Literal per-mm unit:** preserves the calibration figures' single `mm`
  denominator. It does not define ordinary area radiance.
- **Area-label hypothesis:** records a possible missing exponent as an
  unverified interpretation. No numerical conversion or squared-mm substitution
  is active.
- **Literal transition-density display:** the currently usable representation
  of the VUV curves. It remains unscaled.
- **Line-integrated upper-state hypothesis:** requires a known state, transition
  rate, angular convention, reduction and original width convention for each
  datum. The broad molecular lineout does not itself provide these.
- **Spectral-density hypothesis:** requires the per-bin or per-transition
  reduction, effective bandwidth and normalized profile convention.
- **Shot relay response:** requires the separately sourced two-mirror/window
  throughput and allocation of losses already present in lamp calibration.

No unverified scenario is admissible for absolute inference. The lamp/slit
etendue, 0.037 mm^2 sr, remains separate from the unknown shot relay etendue.
The [calibration audit](absolute_spectral_calibration.md) contains the underlying
dimensional and factor-ownership evidence.

For each future admissible nominal radiance reconstruction, the sensitivity grid
has 41 logarithmic points, `10**linspace(-2,2,41)`. Its convention is
`scenario_corrected_radiance = C_abs * nominal_radiance`. This 0.01--100 range is
exploratory; it is not a calibration confidence interval. It cannot repair a
dimensional inconsistency and is not applied to the current derived-density
displays. Those observations permit only `C_abs = 1` in their literal space.

## Comparison admission and remaining evidence

An owner prediction must match the observation's actual quantity, retained
coordinates, spatial support and acquisition. The display comparison score is
the maximum departure from the recovered center divided by the appropriate
lower or upper readout-envelope half-width. A score of at most one means
readout compatibility only. It supplies no likelihood or statistical weight
from the density of digitized samples. Coordinate readout bounds remain separate
metadata; they are not a guessed detector response.

Each observation's readiness report exposes reduction-mapping, signed photon
clock and collection-kernel admission separately. The 3-ns VUV gate duration
does not establish its endpoints in the simulation clock. Table 4.1's reported
30-ns and 29-ns before-breakdown labels are not silently assigned as starts,
centers or ends. The visible 330-ns label is not registered to the VUV gate.

Physical inference additionally requires qualified MFEM heating/source provenance,
Raikou boundary-source acceptance and a resolved physical population ensemble,
state/rate coverage, optical-depth or transfer admission, and the physical
neutral-temperature interpretation. Coverage requires a sourced monolayer
convention. These requirements remain explicit; this evidence milestone makes
no numerical absolute correction or neutral-parameter claim.

Focused verification is in
[test_pilot_comparison_evidence.py](../tests/test_pilot_comparison_evidence.py).
Set `DISSERTATION_PDF` to the external exact PDF to include both source integration
checks; otherwise those two checks skip explicitly.
