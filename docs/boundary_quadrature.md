# Positive joint boundary quadrature

This manufactured interface represents finite emission windows with positive
weighted records of parent state, emission time and surface temperature. It
uses accepted saved polynomial pieces from the
[coupled emission-window interface](coupled_emission_windows.md). It computes no
surface trajectory, source rate, old window export, particle flight or gap
solution. Transport-owner implementation and acceptance remain separate.

The records are deterministic quadrature contributions, not sampled particles.
They retain their full source identities, exact piece support and patch
geometry. No position, molecular mass or velocity is assigned.

## New query support and reuse

The two new queries span exact source-duration fractions `[3/32,11/32]` and
`[13/256,99/256]`. They are source-emission windows, with no experimental
trigger or photon-gate registration. Only accepted complete-support pieces
supply their coefficients; previous window exports are not called again.

Exact de Casteljau restriction defines each new supported piece. A cache keyed
by canonical source identity and exact restricted bounds reuses complete cells
shared by the nested queries. Source-query references and patch references
remain distinct from unique node evaluations. Each query retains its own
exact local-time anchor before moments are combined.

Exact-zero emission pieces receive an explicit catalogue entry. Positive
pieces require finite positive exported number weights and interior nodes.
Unsupported projection or source support rejects; there is no weight floor,
clipping, jitter or renormalization.

## Five-point rule and finite arithmetic

On a restricted coordinate `v` in `[0,1]`, cumulative emission and temperature
are cubic polynomials. The emission measure is their positive quadratic
`dM/dv`. The required integrands are

```text
[1, tau, tau^2, T, T^2, tau*T] * dM/dv
```

Their maximum degree is eight. The ideal five-point Gauss–Legendre rule is
exact through degree nine. Its interior nodes and positive weights therefore
represent all these polynomial moments without replacing the joint mixture
with a mean temperature.

The implementation generates its finite rule at 80 decimal digits. Its rational
decimal approximations are not exact algebraic Gauss–Legendre nodes. Evaluating
polynomials at those approximations and projecting exported fields introduces
two separate discrepancies, each retained against the exact polynomial target.
The restricted derivative already contains the interval scale; node number
weights receive no additional duration factor.

## Actual exported fields and energy

Unprojected areal number contributions are cached once. True geometric area
is applied once to each patch view before its own explicit projection.
Extensive fields are never scaled by area a second time. Reference sums made
from actual exported binary fields remain separate from unprojected target
sums and from projections of the areal display fields.

At each node, expected outgoing translation is `2*kB*T`. Supplied removal
energy `D` and internal energy `I` remain separate. The kinetic activation
barrier is not used as a removal debit. Costs come from explicitly hash-bound
source inputs; no energy channel is inferred as a residual.

If temperature and expected energy are independently projected fields,
`E_export` and the exact value `2*kB*T_export` can differ. Both definitions and
their weighted-sum difference must remain explicit. The finite-field ledger
uses the actual recorded weights and field values. It does not silently
substitute an ideal temperature or energy to close a budget.

Prescribed input and wall thermal change belong to the upstream window
diagnostics. This representation does not assign either to individual parents.
The records supply no internal-state branching,
fragment recoil, molecular-state data or physical velocity distribution.

## Three new discrepancy layers

For each moment or energy quantity, define:

```text
Q_target = exact integral of the saved polynomial measure on the new window
Q_rule   = exact rational sum from finite high-precision rule approximations
Q_binary = sum formed from the actual exported binary fields

finite_rule_defect = Q_rule - Q_target
field_projection_defect = Q_binary - Q_rule
total_defect = Q_binary - Q_target
```

The total equals the sum of the first two defects exactly on the represented
values. No normalization repairs any term. Normalized statistics use each
measure's own positive emitted number, with explicit zero-mass availability.
Variances and time–temperature covariance retain both within-component and
between-component contributions when patch mixtures are combined.

These new discrepancies are separate from the inherited surface-integration,
interpolation and IVP-reference qualifications. Previous accepted windows do
not certify the underlying IVP error on the two new queries. This stage makes
no new-query IVP-accuracy claim.

Independent checks use ordinary power-basis moment targets and a distinct
higher-precision Newton–Legendre rule. They assess gross quantities, resolved
variations, covariance and width contrasts on declared scales. Finite changes
are numerical estimates, not certified enclosures or physical uncertainty.

The independent rule comparison uses Newton–Legendre calculations at 120 and
160 decimal digits. Before calculation, the declared allowances are 10⁻¹⁰⁰
for that reference ladder, 10⁻⁶⁰ for the finite 80-digit rule and 2 × 10⁻¹⁰
for represented-field and total discrepancies on their stated gross,
variation, covariance and width-contrast scales. Structural zeros receive no
numerical floor. These limits qualify a representation of the saved
polynomials, rather than a new surface-kinetic prediction.

This remains a manufactured stationary-patch boundary reference. No physical
H2O coverage, calibrated spectrum, Hall transport, FLASH result, owner solver
implementation or host-budget acceptance is established.

## Interface

`admit_boundary_source` accepts seven named raw-byte inputs and explicit
expected hashes. It binds the accepted window report, coefficients, original
source and window configurations, independent completion/report and review.
The command-line loader owns paths, copying and late byte-identity checks.
It loads this single standard-library module directly so that the repository's
package initializer does not import unrelated scientific solvers.

`BoundaryPiece` retains exact support and source controls; `BoundaryQuery`
declares a finite duration-fraction window. `gauss_legendre5_rule` supplies the
finite rational rule, and `make_piece_quadrature` evaluates each distinct
restricted support once. `summarize_piece_quadrature` applies an explicit area
and time anchor. `combine_boundary_summaries` combines unnormalized quantities
before forming statistics. `build_boundary_quadrature` assembles the source,
patch, aggregate and contrast views without calling any older producer.

The reference currently admits one explicitly pinned manufactured source bundle.
Changing supplied hashes alone does not qualify a different source or physical
parameter set.

## Export

From the complete repository checkout, supply the accepted source-pass directory
and a new output directory:

```sh
.venv/bin/python tools/export_boundary_quadrature.py \
  --config configs/boundary_quadrature_manufactured.json \
  --source-dir /workspace/research_work/mitl-upstream-20261019T0737 \
  --output-dir /path/to/new-boundary-export
```

The loader captures and validates the configuration bytes, checks all seven
source hashes and acceptance associations, copies those input bytes and
rechecks their identity before completing the export. The report binds the
code, configuration, records and flattened CSV views. Reproduction of the
new quadrature uses those saved coefficients; it does not repeat the upstream
source integration or independent IVP calculation.

## Manufactured export results

The two-query export contains 66 distinct restricted pieces and 330 unique
nodes. Those nodes supply 570 source-query references and 1,140 patch references
across six source queries, twelve patch queries, four hypothesis aggregates and
two width contrasts. All source-evolution, RHS, old assembly and old window
export call counts are zero.

The exact polynomial targets give the following nonuniform-minus-uniform
differences. Displayed values are rounded; the export retains exact fractions.

| Source query | Released parents | Expected translation (J) |
| --- | ---: | ---: |
| `shifted_interior` | +2.84778155446 × 10¹¹ | +5.62208575510 × 10⁻⁸ |
| `quarter_cell_interior` | +2.99243736966 × 10¹¹ | +7.47818577654 × 10⁻⁸ |

These are differences between manufactured saved source measures. They do not
establish a measured desorption contrast, physical fragment temperature or
spectral discrimination.

All 52 new analytic and admission tests passed in one targeted run. One
independent verification evaluated 660 new Newton–Legendre nodes across its
120/160-digit levels and recomputed 330 supplied-rule node fields. It checked
2,130 CSV rows and 68,864 cells, including 18,396 floating-field projections.
Neither calculation imported an earlier scientific solver.

The largest scaled width-contrast discrepancy from finite rule construction
was 2.874 × 10⁻⁸¹ against the predeclared 10⁻⁶⁰ allowance. The largest scaled
represented-field or total width-contrast discrepancy was 8.696 × 10⁻¹⁴
against 2 × 10⁻¹⁰. The largest reference-ladder change over source, patch,
aggregate and contrast checks was 8.038 × 10⁻¹²⁰ against 10⁻¹⁰⁰.
These are numerical checks of saved polynomials on declared scales, not
uncertainty bounds for a physical source.
