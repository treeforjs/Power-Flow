# MITL joint source, enthalpy and conductive-heating handoff

Date: 18 October 2026

Power-Flow is pinned to `c236ef2946484804ecf229117e3d0e912911e283` on
`research/flash-publication-foundation`. Local commits separate source-qualified
enthalpy and its inspection workflow, the periodic conduction reference and its
exports, and joint source statistics. The public patch is
`joint-enthalpy-heating.patch`; the extracted subset contains 39
byte-identical public files plus four named packaging wrappers. The preceding
frozen packet remains at
`/workspace/handoffs/mitl-convergence-heating-20261018T1858.tar.gz`.

## Joint source comparison

`joint_source_moments` integrates law-specific conditional emission times and
actual current/width patch populations. It retains original-width and original-
interval conditional means, variances and covariances, including the analytic
within-patch second/fourth spatial moments. Cov(time,x) can cancel under
symmetric edge heating; Cov(time,(x-face_center)^2) retains that structure.
The matched-width hypothesis removes width correlations while preserving the
same time, current, state and local-temperature mixtures.

In the manufactured example, the order 32 resolved source has
Cov(time,width²) approximately **-1.31395838e-14 s m²**, while the matched-width
value is exactly zero by construction. Cov(local expected translation,width²)
is approximately **2.37014563e-27 J m²**. These are numerical fixture results,
not calibrated water release or experimental discrimination.

The energy variable is each component's expected outgoing `2*kB*T`. Its
variance describes variation of those conditional means, excluding the
within-Maxwellian energy spread. Sampled energy, molecular internal energy,
bound-removal reference costs and the kinetic barrier remain distinct. A
molecular-state owner must allocate the physical adsorbate inventory without
duplicating it among independent pools.

The source was integrated once. Orders 8, 16 and 32 inspect the same hash-associated
HDF5; higher-order inspection does not repeat depletion. The relative change
in time–width² covariance decreases from roughly 5.89e-6 (8→16) to 2.66e-7
(16→32). These are finite quadrature changes, not certified error bounds or
source-grid convergence. `evidence/joint_quadrature_comparison.json` retains
the exact differences and identities.

Split start/local-time terms, dominant-population anchors and centered spatial
polynomials retain resolvable local correlations beside large common origins.
Absolute means still have coordinate quantization. Nonzero statistic products
that round to structural zero reject explicitly; the time kernel's underflow
flags/logs remain attached. Source-time subinterval truncation and arbitrary
diagnostic gates are not implemented. The source statistics do not advance
gap particles, calculate molecular reactions/radiation or predict spectra.

## Bounded solid enthalpy input

The Mills 2004 factory pins the reviewed source-manifest and PDF metadata,
requires explicit austenitic-family/solid-phase qualifications and restricts
the caller-selected interval to 300–1700 K. The caller also selects the reference
temperature and either printed Eq. 11 or an explicitly re-zeroed convention.
Printed rounding/reference offsets remain visible.

Exact rational printed coefficients and exact float input temperatures retain
tiny enthalpy increments and faithfully rounded endpoints. The inverse rejects
outside-domain values without a tolerance or clamp. Independent Decimal,
quadrature and derivative checks cover both conventions, adjacent-temperature
increments, source endpoints and strict exterior rejection. The example exports
source-evaluated Cp/H at illustrative temperatures; it does not manufacture new
material measurements or select actual foil properties.

Mass-specific Cp is not Cv, and enthalpy is not automatically EOS internal
energy. This interface supplies no density, liquid/fusion branch, actual foil
composition/condition/phase validation, transport-property splice or solver
default. The example enthalpy and heating/source fixtures are independent;
their outputs are not silently coupled into one physical foil calculation.

## Resistive heating with thermal conduction

The periodic constant-property reference evolves a decaying magnetic cosine
mode, its signed current and positive Joule heating, together with the exact
forced thermal-conduction mode. It includes zero conduction and the coincident
thermal/resistive decay-rate limit. Point and cell means are distinct; Joule
means integrate J² independently. Magnetic energy decreases into thermal energy
while conduction redistributes heat without changing its spatial mean.

The manufactured export closes its analytic combined energy to about 4.89e-16
relative residual. Twelve independent conservative finite-volume B/temperature
calculations on 32/64/128 cells span zero, subresonant, resonant and superresonant
conduction. Their errors approach second order (observed 1.9940–2.0013); maximum
relative discrete energy residual is 2.12e-10. Numerical artifacts were reopened
and verified. These are bounded independent verification prototypes, not a
production heating solver or a FLASH result.

Material is stationary, coefficients are supplied constants and the domain is
periodic. Full Hall-MHD, flow, foil-boundary geometry, variable properties,
phase change and a material EOS are not qualified. The deliberately large
example thermal conductivity is labelled manufactured; it is not a stainless
steel measurement. The reference uses an explicit volumetric Cv independently
of the mass-specific Cp source contract.

## Owner transaction exercises

The work-only exact-accounting harness tests 15 manufactured one-shot cases:
four accepted outcomes and eleven rejections. Every rejected committed state
remains byte-identical to its pre-call state. Accounts separate actual launch
energy, diagnostic 2kBT expectation, adsorption-reference change, internal
energy, departing advected mass and prescribed-wall support work/impulse.
Cases include stale versions, duplicate IDs, exhausted inventory, a failure
after provisional insertion, deterministic replay and JSON restart.

This is a copy/commit accounting example with symbolic state counts. It tests
neither a live host, molecular nuclei/charge catalog, EOS admissibility,
concurrent ownership nor durable transaction storage. It evolves no
trajectories. The integrating owner receives the script/results as a handoff;
no Quark source, runner, build or production job was used or changed.

## Observational identifiability control

The source's full retained time/state/velocity/current marginal is preserved
by the matched-width construction. A fixed linear response independent of
launch width, with appropriate full-width weighting/boundaries, therefore
produces identical observations. Nonzero source covariance alone does not
establish spectral discrimination. Width-dependent loss, fields/transport or
diagnostic weighting, and source-dependent local nonlinear response can break
the limit. This is a sufficient equality limit, not a claim that integrated
spectra can never discriminate.

The exact-arithmetic work example contracts generic nonnegative component
matrices with fixed channels on unequal patch areas. It verifies the limit
and manufactured counterexamples without supplying a spectral solver,
molecular catalog or physical reaction law. The guide identifies the control
and mechanism checks for the integrating owners. Lower-order mean matching
alone is insufficient to establish the required full-marginal equality.

## Verification and next work

The final independent source extraction passed **262 tests**. Its file
counts, exact commands/exits, log hashes, runtime versions, input/source/output
hashes and resolved findings are in `verification.json` and `evidence/`.
The extraction exercises new import, data and CLI boundaries. Earlier unchanged
whole-project suites were not repeated. Initial failed checks and their
corrections remain recorded; targeted retries do not inflate unique counts.

Use `source/README.md` for independently runnable public workflows. Work-only
finite-volume/accounting scripts preserve their recorded workspace paths and
explicit limits. Original papers, scanned pixels and private implementation
are excluded. Boltzmann, collisional-radiative, state-resolved collision/
radiative data and kinetic-gap development remain separately owned. No remote
push, merge, publication or external message occurred.

Next work compares these deterministic moments with weighted paired packet
samples, adding the correct within-Maxwellian kinetic scatter, then qualifies
partial source-time intervals. A source emission-time gate is not a recorded
spectral gate after transport and radiation. Actual FLASH execution still
requires installed licensed source and the configured workstation.
