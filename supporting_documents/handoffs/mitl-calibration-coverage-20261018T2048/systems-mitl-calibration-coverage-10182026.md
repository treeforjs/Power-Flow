# MITL calibration reconstruction and monolayer regression

Date: 18 October 2026

Power-Flow is pinned to `f4bcf3e4d0e402c3ca57b8ef91dd4e09717e2fda` on
`research/flash-publication-foundation`. Five reviewable local commits separate
coverage inference, its CLI/controls, explicit radiometric conversion, published
calibration recovery and primary manufacturer evidence. The patch is
`calibration-coverage.patch`. The extracted subset contains 19 byte-identical
public files and four named packaging wrappers. The previous frozen packet is
`/workspace/handoffs/mitl-packet-verification-20261018T2028.tar.gz`.

## Coverage regression implemented

The fit consumes externally supplied coverage-to-detector-bin predictions. It
requires an explicit monolayer areal density/definition/provenance and matching
observation/prediction bin units and support. Physical parent inventory is not
duplicated among molecular states. Signed observations and backgrounds remain
usable; full correlated covariance is retained through Cholesky whitening.

Coverage is interpolated only inside supplied knots. An analytic search examines
every interval, endpoint, stationary candidate and active zero-scale transition.
Calibration has fixed, Gaussian linear-scale constraint or free modes. Outputs
retain profiled costs, data/prior terms, stable objective contrasts, boundaries,
conditioning and coverage/amplitude degeneracy. There is no automatic confidence
interval or absolute-coverage claim from an unconstrained amplitude product.
Weak-shape cancellation, large common objectives, parameter rounding and numeric
range failures have explicit checks and qualifications.

The CLI example recovers its manufactured 1.4-ML / scale-1.12 fixture. Six
additional manufactured cases show that an incorrectly fixed scale changes the
inferred coverage, while amplitude-only predictions with a free scale leave
coverage unconstrained. The chosen lowest numerical minimum in that degenerate
case is not a measurement. These arbitrary count templates and their monolayer
convention are not physical water spectra or calibrated surface inventories.
The fitting layer calls no separately owned forward solver.

## Radiometric units and evidence recovered

The response interface separates bin width, effective exposure, étendue and
supplied counts per collected photon. It accepts explicitly declared photon or
energy spectral radiance per nm under constant-within-bin assumptions, with the
appropriate energy-to-photon bin integral. Required factor-ownership flags and
measurement-condition association prevent silently adding integration twice.
Signed inverse calibration propagates the full supplied covariance; it adds no
missing calibration uncertainty. The numerical PSD qualification is distinct
from the fitting likelihood's positive-definite covariance requirement.

The dissertation recovery retains 243 original Fig. 4.18 lamp vertices and 340
Fig. 4.19 records: 101 fit vertices, 36 marker centres and 203 closed-band
vertices. Source paths, axes, literal units, exclusions and display envelopes are
pinned. One real lamp wavelength backstep remains unsorted; the curve is not an
interpolation-ready response. The gain fit outside approximately 20–55 is
extrapolated display data; the plotted fit band has no recovered probability
interpretation. Exact rational axis fits reproduce the same CSV bytes in both
checked Python environments. Original source pixels and documents stay local.

The dimensional forward equations distinguish line-integrated energy radiance
from per-nm spectral radiance, normalized line profiles, spectral-bin means and
photon conversion. This repairs the equation/unit distinction without inventing
an experimental normalization.

## Absolute normalization remains unresolved

The user requested rederivation from the dissertation. Its literal plotted
intensity unit is µW/(sr·nm·mm), while the stated lamp étendue is 0.037 mm²·sr.
Original detector counts, the actual lamp certificate, relay acceptance and
throughput, and shot-specific gain/effective gate records are unavailable.
Therefore no justified numerical amplitude correction or physical monolayer
estimate is installed.

Primary follow-up recorded 11 official manufacturer responses and 15 pinpoint
findings, including failed leads. The legacy 632 document itself offers spectral
radiance while labeling its sample graph spectral radiant intensity without an
area denominator. Generic grating efficiency is distinct from coating
reflectance and complete response; its plotted percent/fraction convention is
also unresolved. Photek separates photocathode QE, MCP amplification, phosphor
conversion and gate/readout operation. These documents do not replace a serial
certificate or actual shot records. Current successor-family specifications
remain separate. Read `source/docs/calibration_primary_source_followup.md` and
the hash-associated source manifest.

## Verification and continuation

The final targeted modules passed 174 unique checks. The portable extraction
passed **173 checks with one explicit external-PDF skip**; original-workspace
exact-PDF recovery and changed-source guards passed separately. Earlier project
suites were not repeated. Two documentary files were added afterward; every
previously tested public source byte remained identical, and no unchanged suite
was repeated for that addition. Exact commands/exits, log hashes, runtime/code/
input/output hashes, prior failed attempts and their resolutions are retained.
The earlier standalone optimizer oracle preserves its earlier-code pin; final
tests rechecked the relevant independent oracles on the frozen implementation.

Use `source/README.md` and `source/docs/surface_coverage_inference.md` for the
interfaces. Work-only scripts retain recorded workspace paths. Source papers,
manufacturer response bodies and source-pixel QA are excluded from this archive.
No state-resolved collision/radiative data, CR/Boltzmann solver or kinetic-gap
evolution was edited. No Quark source, runner, build or production job was used.
Actual FLASH remains pending installed licensed source and the workstation.
No push, merge, publication or external message occurred.

Next unblocked work connects explicitly defined coverage to upstream physical
inventory/source-grid inputs without duplicating state inventories, and recovers
published spectral observations suitable for the supplied-grid fit. A numerical
absolute calibration remains conditional on the missing instrument evidence;
shape comparisons must retain that uncertainty. Qualified source-time windows
and the owner response-control handoff remain queued.
