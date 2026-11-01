Pinned H2 dissociative-excitation rate component, host reference v1

This isolated runtime reads exactly the MCCC X(1Sigma_g+), v=0 to b(3Sigma_u+)
dissociative-excitation source. The licensed source stays at the local path
in INPUTS.json and is never included in a repository artifact. It is a
calculated integral cross section, not an experimental differential event law.

load_primary requires exact source, blob, manifest, spec and strict policy
hashes. It validates target/state labels, eV/a0^2 units, row structure,
strictly increasing energy, finite nonnegative data and the threshold anchor.
The comparison registry is prepared v3.2.2; current native dependencies remain
unchanged. This is not a replacement blob or a new accepted native binding.

diagnostic_rate requires an explicit published-source-energy operand and returns
a Decimal area in m^2, using the actual pinned chemistry Bohr-radius constant.
Interpolation is linear in the source energy and cross section, explicitly a
numerical interpolation choice. No extrapolation, endpoint hold, unknown zero,
rate array export or energy-frame substitution is permitted. Number density
output conventions remain cm^-3; cross section units do not change them.

threshold_assessment derives the finite-target incident threshold from the
actual reference-energy and mass operands, using invariant s and a stable
reference-energy difference. It never subtracts GeV masses to obtain Q:
 T_thr = Q [2(m_e+M_H2)+Q] / (2 M_H2).
The source's declared onset is compared separately with Q and T_thr. The
source print resolution is not a physical error bar. No claim that a model's
cross-section uncertainty establishes a thermochemical threshold tolerance is
made. In the present pinned inputs the thresholds differ; their numerical and
frame compatibility remains unresolved. event_rate always refuses.

Resolving the onset alone would not qualify events: a differential electron
scattering/fragment recoil and energy-sharing law is still required. Other H2
channels, energy tails, e+OH/e+O, and unknown recombination remain unresolved.
Any future group import must retain those explicit unknowns, not turn the
minimum encounter-presence check into a complete-network assertion.

Primary data: Scarlett et al., DOI 10.1016/j.adt.2020.101361, and MCCC database
https://www.mccc-db.org/ . The source table hash and public POST selector are
preserved in the earlier water-fragment-primary-data-inventory-20261101 artifact.
No scientific criteria or source estimates changed in this component.

Numerical controls run only through frozen, one-core SimQueue. They exercise
the real pinned source and registry, all table knots, interpolation, unit and
state mistakes, bad tables, source/registry/policy tampering, threshold conflict,
both low/high unsupported ranges and immutable state. Two independent processes
must emit identical assessment reports. No GPU/native event qualification is
claimed. Evidence contains metadata and pass/fail results, not raw rate arrays.

production_qualified=false
physical_model_validated=false
network_complete=false
native_adopted=false

Revision v1.1: explicit Decimal rounding/exponent/trap/flag policy, independent of caller context, and stream-bounded reads. Draft v1 job 433 was canceled while queued after peer review; it never ran. No physical policy changes.
