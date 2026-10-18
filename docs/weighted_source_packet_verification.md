# Checking weighted launch packets against source moments

The deterministic source report and the paired boundary packets describe the
same emitted population at different levels. A useful check compares their
joint time, width, current-coordinate and energy statistics before any gap
evolution. The source history is reused without integrating release again.
This check does not predict transported populations or spectra.

The work-only verification uses a hash-associated manufactured planar history,
an original source interval and explicit exact-state masses. A dedicated
configuration retains the source, mass and geometry labels. The planar sampler
already uses each physical width/current patch area; a one-dimensional
extrusion length must not multiply those counts again. The hypotheses are
alternative position arrays for one packet population, with every packet's
time, current position, state, mass, temperature, velocity and weight shared.

## Expected local energy and actual kinetic energy

For the implemented stationary, fully accommodated outgoing thermal flux,
write `theta_c = kB*T_c` for source component `c`. Its two tangential velocities
are Gaussian with variance `theta_c/m_c`; its outward normal speed follows the
flux-weighted distribution. Consequently its translational kinetic energy
`K = m_c*|v|²/2` has a gamma distribution of shape two:

```math
E[K\mid c]=2\theta_c,\qquad
\operatorname{Var}(K\mid c)=2\theta_c^2,\qquad
E[K^2\mid c]=6\theta_c^2,\qquad
E[K^4\mid c]=120\theta_c^4.
```

These are outgoing flux moments. Substituting a bulk Maxwellian's `3*kB*T/2`
would change the launch model. The exact-state mass sets the velocity scale;
the mean translational energy at a given temperature does not depend on it.
Masses and state populations are supplied inputs, with no species fallback or
inferred chemical partition.

Let `P_c` be the emitted-particle fraction in component `c`. The complete
population then satisfies

```math
E[K]=\sum_cP_c\,2\theta_c,\qquad
E[K^2]=\sum_cP_c\,6\theta_c^2,
```

```math
\operatorname{Var}(K)=\operatorname{Var}_c(2\theta_c)
                    +\sum_cP_c\,2\theta_c^2.
```

The joint-source interface reports the first variance term: variation of the
local conditional means. Comparing that term alone with the variance of raw
sampled kinetic energies omits the within-component thermal spread. For a
single-temperature population the reported conditional-mean variance vanishes,
while the kinetic variance remains positive.

For independent uniform input columns, time, position and thermal velocity
are independent inside each frozen component. Thus source cross-covariances involving
the local mean energy and time or squared width also give the corresponding
population covariances with actual `K`. Thermal scatter changes energy variance
without adding a within-component time/position–energy covariance. This relies
on the stated independence and frozen local-temperature model.

Expected `2*kB*T`, actual sampled energy, internal-state energy and the bound
removal reference remain separate accounts. A wall transaction uses the actual
weighted launch energy; finite sampling error is not another heat source and
does not justify rescaling velocities. An activation barrier alone does not
determine a fragment's translational energy. Nonthermal release, incomplete
accommodation, breakup recoil and state-dependent energy exchange require
separately qualified channels; they are not inferred by this thermal check.

## Component weights and sampling uncertainty

The existing sampler seeds every positive source component, then refines its
allocation. With `n_c` samples and emitted count `N_c`, each packet in that
component has physical weight `N_c/n_c`. These weights preserve all component
inventories, including low-population components. Counting packets equally
generally changes the mixture and is not an alternative estimator of the same
physical source.

Conditional on the allocation, a weighted raw mean has the form

```math
\widehat\mu_A=\sum_cP_c\,\overline A_c,\qquad
\operatorname{Var}(\widehat\mu_A\mid\{n_c\})
   =\sum_c\frac{P_c^2}{n_c}\operatorname{Var}(A\mid c).
```

These uncertainty formulas assume externally supplied independent uniform
variates across packets and sampling columns, including independence of the
allocation and within-component columns. The sampler validates open variates;
it does not establish their statistical independence. Correlated inputs or
quasi-random designs require a separately derived uncertainty calculation.
Under the stated model, conditional component means are unbiased even though
some allocation counts are random. A sample-based variance estimate needs at least
two samples in each relevant component; one sample does not establish zero
uncertainty. A single effective-sample-count formula for independent identical
draws does not describe this stratified allocation.

For the selected thermal law, conditional-model variances are available even
for singly sampled components:

```math
\operatorname{Var}(\widehat\mu_K\mid\{n_c\})
 =\sum_c\frac{P_c^2}{n_c}\,2\theta_c^2,
```

```math
\operatorname{Var}(\widehat\mu_{K^2}\mid\{n_c\})
 =\sum_c\frac{P_c^2}{n_c}\,84\theta_c^4.
```

Here `84 = 120 - 6²`. These quantify packet sampling under the supplied model,
not physical uncertainty in temperature, desorption parameters or the data.

The usual weighted population-variance plug-in,
`mu_hat_K2 - mu_hat_K²`, is biased downward by the conditional variance of
`mu_hat_K`. Adding that variance gives a conditional-model unbiased estimator
of the population kinetic variance. Preserve both values and the correction;
the corrected estimate can still fluctuate in a finite realization. More
generally a covariance plug-in needs the corresponding covariance of its two
estimated means. Finite-packet variances must not be treated as exact conserved
identities.

## Verification evidence and limits

A bounded study should retain independent fixed seeds, multiple packet
budgets, component allocations, the supplied open seven-column variates and
the sampled arrays. Readback can then recompute physical counts, weighted
statistics, paired differences and output hashes. Exact controls cover paired
retained fields and component inventories; stochastic comparisons cover
normalized moments. Report the actual finite ensemble bias and root mean
square error without requiring every seed's error to decrease monotonically.
Do not select new seeds to obtain a preferred outcome.

For an original interval, compare with that same interval's deterministic
conditional record. Numerical emission-time quadrature changes remain distinct
from packet sampling error. A source emission-time restriction is also distinct
from the recorded photon gate after transport and radiation.

This verification does not determine molecular chemistry, allocate a shared
physical adsorbate inventory among internal states, advance charged particles,
or establish a spectral difference. The
[observation invariance control](source_identifiability_controls.md) still
applies: preserving the full retained marginal gives identical observations
under a fixed linear launch-width-independent full-width response. The
separately owned forward model must quantify a physical mechanism that retains
the width information and compare its signal with measurement uncertainty.
