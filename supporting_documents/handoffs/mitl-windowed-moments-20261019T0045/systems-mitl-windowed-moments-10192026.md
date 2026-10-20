# MITL window-conditioned timing, geometry and thermal statistics

Power-Flow is pinned to `313d4bb128c39a6996c06ce2dc1527f35f7f044e` on `research/flash-publication-foundation`.
Three local commits separate local timing quadrature, joint mixture statistics
and guarded inspection of saved exports. `windowed-moments.patch` contains this
pass's changes. The prior packet is
`/workspace/handoffs/mitl-windowed-boundary-20261018T2348.tar.gz`.

Open probability-space quadrature integrates the qualified source window's log
local delays. Centered scaled variance avoids raw-moment subtraction; logs retain
positive moments outside their SI projection range. Original interval origin,
window start and local mean stay separate in the joint calculation. A fixed
dominant-component anchor supports comparisons between quadrature orders without
first rounding away local timing at a large origin. Global mean time has its own
precision flag. The scalar law's numerical query qualification and finite-order
quadrature resolution remain distinct.

Joint statistics use exact rational arithmetic on represented component counts,
rectangle geometry and local moment estimates. Tiny displayed probabilities do
not delete components. Width redistribution changes only the x law and preserves
the time/current/state/temperature mixture. The full kinetic-energy variance is
`Var(2*kB*T) + E[2*(kB*T)^2]`; reporting only the first term would omit thermal
scatter. Supplied internal and removal costs remain separate from translation
and the kinetic activation barrier. Required unrepresented timing or final SI
statistics reject explicitly; optional fourth-energy diagnostics retain logs.

All 105 new checks passed, and all 105 passed in the portable subset: 48 timing,
27 joint-mixture and 30 inspection checks. One added analytic covariance
assertion reran only its affected existing test and adds no unique test count.
Three new exact-rational finite-allocation cases enumerate 120 outcomes; six
focused aggregate-helper checks retain rare-component effects. The latter's
extreme fixture does not qualify standalone component energy moments. Exact
scope and pins are in `verification.json` and the independent reviews. Earlier
passing suites and saved oracles were not rerun. The original saved source
evolution and packet generation were reused; new tests use fresh small fixtures.

The following bounded verification reused the saved 48-component, 512-packet
manufactured export. New moment estimates use orders 8/16/32. Six weighted means
and six conditionally independent cross-covariances use actual represented
packet-weight targets and exact finite-n allocation formulas. Thirty-five
components have one packet: their conditional variances come from the source
model, rather than a false empirical singleton variance. Signed count/target
rounding differences remain explicit. Maximum standardized deviations are about
1.04 for means and 1.59 for cross-covariances; these are descriptive values,
not a statistical pass/fail or physical validation. Quadrature changes remain
separate from sampling, source-parameter and experimental uncertainty.

`evidence/windowed-packet-qualification/` contains the reviewed JSON/CSVs.
`evidence/windowed-packet-figure-final/` contains the final PNG and vector PDF;
the earlier layout and its original script are preserved separately. All figures
are labeled manufactured. The source subset contains byte-identical public files
plus four named packaging wrappers. Exact commands, exits, input/code/output
hashes and archived payload verification are recorded. Commands after freezing
remain in the local pass's `commands.jsonl`.

The original saved packet export's files are preserved in
`evidence/windowed-moment-example/` with the prefix `packet_export_`. Removing
that prefix into a new directory restores the original export bytes; its report
and `input_export_sha256` bind their original names. Inspection validates the
stored packet schema, units, array dimensions and canonical source association;
it does not independently rerun each saved sample's realization.

Source-emission restrictions are separate from photon recording gates. No gap,
charged trajectory, collision/radiative calculation, molecular state data,
host transaction, wall feedback or actual FLASH run was added. Quark integration
remains an owner handoff. Nothing was pushed or published. Physical water
coverage, raw diagnostics, absolute VUV normalization and configured licensed
FLASH/workstation remain unresolved.

The selected next task is fixed physical SOURCE-window intersection and
resolution qualification across already preserved refined histories. New
partial-release/moment comparisons will separate quadrature changes from source
space/time discretization, while reusing the stored evolution. Read
`source/docs/windowed_source_moments.md` for the exact interfaces.
