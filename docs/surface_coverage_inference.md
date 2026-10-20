# Inferring surface coverage with an explicit calibration model

Surface coverage can be fitted against spectral observations through externally
supplied forward predictions. The inference layer does not calculate molecular
populations, collisions, radiation or gap transport. Those owners supply a
coverage grid whose predictions include the selected heating/desorption
hypothesis and a stated diagnostic response.

The parameter is coverage in explicitly defined monolayer equivalents:

```math
N_{\rm initial}/A=\theta_{\rm ML}\,n_{\rm ML}.
```

Both the monolayer areal density `n_ML` and its definition/provenance are
required. There is no universal water-on-steel value selected here. A kinetic
law's reference coverage, its initial inventory normalization and a physical
monolayer are distinct. The physical adsorbate inventory must also be allocated
among molecular states without creating a separate full inventory for each
state. Fragment populations do not each receive another parent-water layer.

## Regression contract

`CoveragePredictionGrid` carries nonnegative coverage knots, owner-supplied
nonnegative bin predictions, the monolayer definition, bin edges/support/units
and input status/provenance. `BinnedObservation` carries the same bin contract,
finite observed values, the full covariance and an explicit known background.
Signed dark/background-subtracted values remain usable; discarding negative
noise samples can bias the fit. Bin support and units must agree exactly.

The model is

```math
y=b+s\,m(\theta_{\rm ML})+\epsilon,\qquad
\epsilon\sim\mathcal N(0,C).
```

The known background `b` is not scaled with the source. The supplied covariance
`C` must be symmetric positive definite and is used through Cholesky whitening.
Correlations between reconstructed samples or instrumental bins are retained.
Display readout envelopes are not automatically Gaussian measurement errors;
an actual-data fit needs a justified covariance and calibration model.

The forward grid is interpolated linearly between supplied coverage knots.
Searching every interval and endpoint avoids selecting a minimum solely from
the original knots. This qualifies a finite interpolant, not the accuracy of
the physical forward model. There is no extrapolation outside the grid; a
boundary optimum is reported. Refine the grid independently when its shape or
minimum is not resolved.

The search uses analytic stationary candidates, float64 whitening/coordinates
and long-double scalar contrasts. A large common residual can make different
absolute objectives round to the same float; objective differences therefore
retain centered contrasts instead of subtracting those rounded totals. Numeric
range failures raise explicitly. These safeguards do not supply an error bound
for the prediction grid or a proof of physical identifiability.

`CalibrationScale` requires a stated provenance and one of three modes:

- `fixed`: use the explicitly supplied nonnegative scale, retaining a zero-scale
  boundary as an uninformative response case.
- `gaussian_prior`: profile a nonnegative scale with the supplied linear-scale
  mean and standard deviation. This adds a stated calibration penalty.
- `free`: profile the nonnegative scale without a calibration constraint and
  report the resulting identifiability limits.

For `w=y-b`, a Gaussian scale constraint gives the profiled solution

```math
s_*(\theta)=\max\left(0,
\frac{m^TC^{-1}w+s_0/\sigma_s^2}
     {m^TC^{-1}m+1/\sigma_s^2}\right).
```

The free-scale limit drops the prior terms. Zero response or a scale boundary
requires separate diagnostics. A penalty is additional calibration information;
it is not evidence that the observations alone identify an absolute coverage.
A Gaussian penalty in linear scale is not a lognormal prior. Its nonnegative
support and any zero-scale optimum remain visible.

The result retains the best coverage/areal inventory, scale, data and prior
contributions to the objective, grid profile and objective differences, search
candidates and conditioning/boundary/degeneracy diagnostics. Objective changes
are not automatically confidence intervals. Gaussian-noise assumptions,
finite support, calibration information, interpolation and model discrepancy
must be assessed before attaching physical uncertainty to a coverage estimate.

## Absolute scale and comparison of source hypotheses

If changing coverage only multiplies a fixed template, then a free absolute
calibration scale trades directly against coverage. The product can be fitted
while the separate factors remain unidentified. Increasing the number of bins
does not resolve that exact amplitude degeneracy. Coverage-dependent shape,
independent calibration information or another qualified constraint is needed.

The resolved-width and uniform-width hypotheses should use the same physical
coverage convention, diagnostic treatment and nuisance-parameter assumptions.
Fit and retain both objective profiles before interpreting a best-fit
difference. Use the same calibration record and constraint for both hypotheses.
Their conditional best-fit nuisance values may differ; report those shifts and
their penalties. A joint multi-shot fit must share a calibration parameter
where the instrument conditions require it. Multi-shot shared calibration and
wavelength-dependent response uncertainties need an
explicit extension; a single global scale cannot repair an incorrect response
shape, line-spread function or viewing geometry.

The [observation invariance control](source_identifiability_controls.md) remains
applicable. Matched full retained marginals give identical observations in the
fixed linear launch-width-independent full-width limit. A source covariance
alone does not establish a distinguishable spectrum.

The [absolute-calibration recovery](absolute_spectral_calibration.md) separates
published calibration display data from a corrected instrument response.
Ambiguous radiometric units cannot be converted to an absolute monolayer count
by choosing a convenient multiplicative factor. Keep calibration uncertainty
in an explicit nuisance term or a justified observation covariance, with no
double counting of the same calibration error.

## Explicit radiometric response

`RadiometricResponse` supplies the dimensional detector conversion when its
inputs are known. It accepts photon or energy spectral radiance per nm,
wavelength-bin edges in nm, étendue in m²·sr, effective exposure in seconds and
an explicitly supplied counts-per-collected-photon response. Geometry, exposure
and wavelength-bin integration each enter once. The response may include
camera gain greater than one count per photon and must identify the reference
plane, included optical/detector factors and matching measurement conditions.
The condition identifier is a caller association, not recovered calibration
evidence.

Its bounded assumption is constant mean spectral radiance and photon response
inside each bin. Energy-to-photon conversion uses the exact bin-average
wavelength for that assumption. It is not an arbitrary-spectrum midpoint rule
or a line-spread calculation. A line-integrated prediction must first use a
qualified profile/bin convention, with no duplicate bin-width multiplication.
The dissertation's ambiguous single-`mm` ordinate is not an accepted radiance
kind.

Inverse calibration preserves signed signals and propagates the complete
supplied count covariance by the diagonal conversion factors. It supplies no
missing noise, gain or absolute-response uncertainty. Its numerical PSD check
can retain singular covariance; the coverage likelihood separately requires a
positive-definite covariance and cannot be made valid by silently adding a
noise floor. Correlated calibration uncertainties need their own explicit
treatment.

## Manufactured demonstration

```sh
.venv/bin/python tools/fit_surface_coverage.py \
  --config configs/coverage_regression_manufactured.json \
  --output-dir /tmp/coverage-fit --plot
```

The supplied four-bin grid is an arbitrary numerical fixture. Its monolayer
areal density, count vectors, covariance and calibration constraint are all
manufactured. It assigns no molecular wavelengths or reaction processes.
Recorded truth values validate the regression and are not optimization inputs.
No number of physical water monolayers is inferred from this example.

For publication, replace the fixture with qualified owner predictions and
calibrated or explicitly reconstructed observations. Original detector frames,
per-shot response/gain/gate records and a documented monolayer convention remain
distinct evidence requirements. The module installs no FLASH boundary, changes
no separately owned solver and supplies no Quark implementation.
