# A required invariance control for the width comparison

A nonzero source correlation such as `Cov(time,(x-face_center)^2)` establishes
a difference in the emitted population. Spectral discrimination additionally
requires that the measurement response retain information about that
difference. The recovered width-integrated sightline makes this distinction
relevant to the publication comparison.

Let `S(x,xi)` be the emitted source measure per width, where `xi` contains every
other retained launch variable: time, exact state, velocity, gap-normal and
current coordinates, and any additional variables used by the forward model.
Consider a resolved source and its paired uniform-width control satisfying

```math
M(\xi)=\int_{x_-}^{x_+}S_{\rm resolved}(x,\xi)\,dx
      =\int_{x_-}^{x_+}S_{\rm uniform}(x,\xi)\,dx.
```

This is equality of the full retained-variable marginal. Equal total particle
number and mean energy alone do not imply it. In particular, replacing a
mixture of launch temperatures or states by a single averaged distribution can
break this control even when those two moments agree.

For a **fixed response linear in the emitted source**, suppose each channel
has a kernel `K_a(xi)` independent of launch width position. The source domain,
boundaries and diagnostic collection must also permit the same complete width
integration in both cases. With the integrals well defined,

```math
O_a[S]=\int d\xi\int dx\,K_a(\xi)S(x,\xi)
      =\int d\xi\,K_a(\xi)M(\xi).
```

Every such channel is identical for the two hypotheses. Arbitrarily many
channels cannot recover a source coordinate that this response marginalizes
away. The result also holds for a common fixed linear instrument response
applied afterward. It does not require the sources to have the same local
time/position covariance.

The assumptions refer to the complete source-to-measurement operator, including
transport, collection and boundaries. Uniform imposed coefficients by
themselves do not guarantee an x-independent response: finite-width losses,
side boundaries or a spatially restricted sightline can retain launch-position
dependence. Likewise, a model can be linear in a population at prescribed
conditions while its coupled source-to-signal map becomes nonlinear when those
conditions respond to the source.

Mechanisms that can break this invariance include width-dependent fields or
backgrounds, finite geometry and losses, nonuniform diagnostic weighting, and
local nonlinear chemistry or source-dependent plasma feedback. Their presence
does not automatically make a measurable distinction large enough to identify.
Conversely, a nonlinear functional that depends only on the preserved marginal
is still identical. Width-integrated spectra can discriminate in suitable
regimes; that conclusion must come from the qualified forward calculation and
measurement constraints.

The integrating owners should retain the paired time/state/velocity/current
mixture and recover this equality as a limiting control before attributing a
calculated difference to width-dependent desorption. Then enable a stated
mechanism that breaks the limit, quantify its contribution, and compare the
result with diagnostic calibration, timing, collection and reconstructed-data
uncertainties. A failed invariance control calls for checking the source
marginal, boundaries and response assumptions before interpreting the result.
The joint-source statistics in [joint_source_moments.md](joint_source_moments.md)
provide source diagnostics; they are not measured spectral differences.

The accompanying work-only arithmetic verification uses nonnegative rational
component matrices and unequal patch areas. It checks exact marginal equality
and several generic linear channels, then supplies nonuniform-weight and local
quadratic counterexamples. These are dimensionless manufactured functionals,
with no wavelengths, reaction model, state catalog, collisional-radiative data
or physical spectrum prediction. The verification establishes the stated
algebraic control and its limits, not identifiability in the dissertation data.
