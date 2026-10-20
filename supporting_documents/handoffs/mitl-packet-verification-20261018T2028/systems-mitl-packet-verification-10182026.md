# MITL weighted launch verification and source-window reference

Date: 18 October 2026

Power-Flow is pinned to `4accc51e1400099378d9ba85ac51a983594048ec` on
`research/flash-publication-foundation`. The public guide and its local commit
are supplied with this packet. Existing source and packet implementations are
unchanged. The preceding frozen heating/enthalpy/joint-source packet remains at
`/workspace/handoffs/mitl-joint-enthalpy-20261018T1947.tar.gz` and passed 262
isolated checks. That unchanged suite was not repeated.

## Weighted source comparison

The study sampled 24 fixed cases: 1,024 / 4,096 / 16,384 packets for each of
eight fixed seeds, totaling 172,032 packets. It reuses the original manufactured
source without reintegrating release. Original interval 6 spans approximately
60–70 ns and has 16 positive source components. Its exact-state mass is an
explicit manufactured input; no water species or molecular catalog is inferred.
Patch area enters once, without a one-dimensional extrusion multiplier.

For outgoing flux-Maxwellian translation, `K|c` has mean `2*kB*T_c`, variance
`2*(kB*T_c)^2`, raw second moment `6*(kB*T_c)^2` and raw fourth moment
`120*(kB*T_c)^4`. The source interface's conditional-mean variance must therefore
be augmented by the within-component spread before comparing with actual
kinetic-energy samples. In this manufactured interval the full variance is
1.79968926e-40 J² and the conditional-mean
mixture contributes 2.15497832e-42
J². Omitting thermal spread would reduce this variance by a factor of roughly
83.51. This is a numerical fixture result, not a measured species
temperature or a calibrated desorption result.

Each component's physical inventory is preserved by unequal packet weights.
The paired hypotheses retain all time, current, state, temperature, velocity,
momentum and weight fields and differ only in width position. Actual sampled
energy is checked separately against the component energy ledgers; no velocity
rescaling or additional wall heat is introduced. Weighted mean uncertainties
condition on the actual component allocation. The reported population-variance
correction removes its conditional-model finite-sample bias.

The supplied pseudorandom inputs assume independent uniform variates across
columns and packets. The open-variate API itself establishes no independence.
Eight seed streams are reused across budgets, with common prefixes verified;
different-budget cases are correlated. Finite ensemble bias/RMSE and
standardized errors are descriptive. The raw-second-moment RMSE is
nonmonotonic across the smaller two budgets, and all fixed outcomes are retained.
There are no certified convergence rates or physical confidence intervals.

## Independent verification

A second reducer reopened every saved NPZ, checked association with the exact
source/reference/configuration and producer report, and required the exact
budget/seed/component sets. It derived counts from source areal release and
geometry, checked equal weights within each component, and reconstructed energy
both from velocities and directly from the saved uniform transforms.
Long-double reductions of means, variances, sampling uncertainties and
time/width/energy covariances agree with the producer in 576 scalar comparisons.
The largest difference normalized by the declared physical scale is
3.823e-17. Sampling was not
repeated. A comparison-output serialization error was corrected; the failed log
and successful retry are retained. The numerical samples were unchanged.

`evidence/packet-figures-final/weighted_source_packets.png` shows one fixed
largest-budget realization, the theoretical energy mixture, its variance
decomposition and the finite ensemble errors. The display choice is the lowest
fixed seed, independent of its errors. A label-placement revision regenerated
only the figure. The earlier plotting log/report remain in the work directory.

## Next task advanced: source-emission windows

The work-only constant-barrier reference under
`evidence/source-window-reference/` supplies 11 manufactured subwindow cases
with 180-digit Decimal references. Ten projected float input domains are
qualified; one exact-decimal window is rejected because its float endpoints
merge. Represented-count underflow, unresolved linear delay quantiles and
lost absolute-time/log-time spread are flagged separately. A zero partition
residual at the stated precision is an arithmetic result, not a certified bound.
This is a bounded reference for the next upstream interface;
it does not add a production source-window implementation or qualify nonlinear
Temkin tails. Read its report for supported cases, representation and budget
limits. The accompanying proposed contract identifies the missing log-survival,
log-mass and gate-relative timing qualification.

An emission window selects launch times. It is distinct from a recorded photon
gate after transport, excitation and radiation. Earlier neutral launches may
contribute photons during an observation gate. The recovered timing-origin and
shot-association limitations remain applicable.

## Scope and reuse

The guide is byte-identical to the pinned public commit. Source inputs, all
sample arrays, configuration, code/input/artifact hashes, exact commands/exits
and logs are included. Work-only scripts retain their recorded workspace paths
and require that pinned public checkout; relocation requires explicit path
adaptation. The included original mass configuration records its provenance;
the planar study's separate configuration deliberately excludes its extrusion
length. Original source papers, scanned pixels and private implementation are
excluded.

This packet develops no Boltzmann or collisional-radiative solver, collision/
radiative dataset or kinetic-gap evolution. No Quark source, runners, builds,
production jobs or external messages were used. Actual FLASH execution still
requires installed licensed source and the configured workstation. No push,
merge or publication occurred.
