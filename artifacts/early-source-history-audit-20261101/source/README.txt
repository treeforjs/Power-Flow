Early source history audit v1, 2026-11-01

This tool extracts existing recorded expectations and sampled births. It does not
evaluate an emission/desorption model, rerun heating, replay a simulation, change
deployed configuration, or interpolate missing per-patch history.

Supported coupled evidence
  Qualified CPU candidate c2 pilot, 2,000 accepted 1e-14 s steps:
  0.749788055492 to 0.769788055491965 ns (approximately 20 ps).
  Recorded expected electrons: 4.541809285965166e-57.
  Recorded expected H2O desorption: 174429.1661760889 molecules.
  Sampled electron/H2O macro births and physical sampled release: all zero.
  These expectations are source-model values, not confidence/accuracy bounds.
  Both nominal macro weights are 1e6. Zero macros do not establish zero source.
  Expected H2O is 0.1744291661760889 nominal macro weights. Recorded source stock
  stays 2.2308e15 molecules; native H2O/H2O+, sampled release and sinks stay zero.
  Every row reports exact water conservation. Arrival v1.1 stochastically rounds
  each step's D/w0 to quanta and only debits sampled n*w0; unsampled material stays
  in stock. Its expected_sum independently retains cumulative expectation.
  The pilot exports binary64 projections and exact-owner conservation flags,
  with no separate per-patch credit/remainder values or final arrival checkpoint.
  No fractional-credit value is invented. Source-quantization convergence remains
  a gate; this audit changes no macro weight and makes no onset-time prediction.

This exact run uses source-v7 benchmark surface and interval settings: clean
constant phi=4.3 eV, theta0=0.33, one flat beta=1 class, wall/photo off, conditional
analytic oxide-water desorption. Eight source owners cover 32 patches of both
foils (four axial slices per foil, four faces each). The source starts from the
declared stock at the first supported drive knot. The scenario string containing
"oxide-water" does not enable the nonbenchmark work-function/enhancement model.
The interval cap is 1e-14 s and error targets are disabled in this benchmark.
No production source-accuracy or particle-statistics convergence is inferred.

Hall is the v3.3 solid-only thermal binding, native is auditfix v8 archive with
v5 compatibility configuration/replay, and chemistry is parent ionization only.
This is not full H2O fragmentation or coupled SS304 phase evolution.

The 20 ns vacuum history is distinct
  Hall v3.4 contains 194026 vacuum/solid-only steps, 1927 aggregate samples and
  20 complete field/heat checkpoints at integer ns from 1 through 20 ns.
  The aggregate samples contain foil mean/max temperature and four absolute
  face-field maxima, not signed per-patch source states at every accepted step.
  Sparse checkpoints contain fields/heat, but not source-stock/RNG/plasma history.
  No interpolation error bound or cumulative omission bound is established.
  restore_as_genesis rebasing does not supply this missing physical history.
  Extending this 20 ps pilot's expected counts to 20 ns is explicitly refused.
  A 20 ns vacuum restart cannot be called continuous plasma evolution.

The earlier 295-time high-resolution temperature cache is retained and inventoried
separately. It belongs to the older 2D heating/supply calculation, including the
earlier symmetry/anode and magnetic-field limitations. It is not used to fill the
corrected Hall/source history gap. The separately qualified SS304 phase owner v1.3
is pinned for context and has not been silently substituted into either history.
Active Hall v3.6 is excluded from this audit and is not promoted.

Provenance
  INPUTS-SHA256.json binds 201 current input files. Pilot outputs are additionally
  checked against the original run receipt's file hashes; compiled dependencies
  and owned sources are compared against its build receipt. Source model settings
  and owner identities are in report.json. Hall checkpoint bytes are SHA256-pinned
  now; their native semantic checkpoint seals are not revalidated by this tool.
  The Hall run receipt does not itself hash all those original output files.
  Read-only metadata extraction uses frozen Python; no numerical workloads run.

Reproduce (output must not exist):
  /home/treefor/projects/runtimes/arch-base-20261028/run-frozen /usr/bin/python -I \
    /home/treefor/projects/runtimes/early-source-history-audit-v1-20261101/audit.py \
    --output /tmp/early-source-audit.json
  /home/treefor/projects/runtimes/arch-base-20261028/run-frozen /usr/bin/python -I \
    /home/treefor/projects/runtimes/early-source-history-audit-v1-20261101/test_audit.py

report.json intentionally has bound=null and refused_unsupported_history for the
omission through 20 ns. The unsupported range starts after the recorded coupled
segment; there is also no supported history before the first C19 knot. A complete
qualified per-patch history with source/stock provenance or a qualified
continuous coupled evolution would be required to remove this gate.
