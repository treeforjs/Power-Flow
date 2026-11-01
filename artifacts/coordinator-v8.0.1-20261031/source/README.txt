Coordinator v8.0.1 audit correction -- 2026-10-31

This version is based on the sealed coordinator-v8-stage1 artifact. Historical
measurements and architecture are preserved in README-v8-stage1.txt; they are
not qualification results for this corrected build. The release artifact's
RESULTS.txt and QUALIFICATION.json record the corrective build and test receipts.

Corrections

* Restore verifies saved owners, receipts and replay state before applying a
  new operational stopping target. --window is relative to the accepted
  checkpoint time; --end-time is absolute; --steps is the absolute accepted-step
  target and clears a saved time limit. With none specified, the saved target
  is retained. Conflicting requests and non-advancing time targets refuse.
  Saved step targets also survive restart. New checkpoint identities prevent
  an older build from silently interpreting these changed semantics.
* A surface with 0 < extraction field < 1 V/m now refuses the step. The former
  path recorded emission bounds but omitted the electron source. This is a
  model-domain refusal, not a timestep refusal; retries cannot bypass it.
  No low-field quadrature or microscopic launch model has been qualified.
* Execution reports say one MPI process. --ranks only changes saved metadata;
  varying it does not demonstrate distributed restart or scaling.
* A successful executable does not qualify a run if the post-run dependency,
  source or executable hashes changed. Integrity failure now sets a failed
  receipt and nonzero launcher exit status; missing inputs also fail closed.
* The corrected native archive is pinned by its source manifest and library
  hash. The legacy parent-ionization path remains a limited chemistry binding;
  a broader native API does not automatically wire full water fragmentation.

Qualification

All numerical builds and tests run inside arch-base-20261028/run-frozen under
SimQueue, with one allocated core and OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=
MKL_NUM_THREADS=1. GPU tests additionally require a scheduled GPU allocation.
The build receipt pins inputs, dependency headers, compiler commands, objects
and the executable; run_v8.py refuses source changes after the build.

  run-frozen /usr/bin/python -I run_v8.py build LABEL
  run-frozen /usr/bin/python -I run_v8.py run LABEL checks NAME
  run-frozen /usr/bin/python -I tools/audit_regressions.py LABEL TAG direct|replay
  run-frozen /usr/bin/python -I tools/restart_tests.py LABEL TAG
  run-frozen /usr/bin/python -I tools/audit_run_integrity.py OUTPUT.json

The audit regressions run coupled short fixtures and fresh-process restores:
direct genesis restore, continuation after a completed time window, a second
restart after extension, saved time and step targets, switching to --steps,
and refusal of past/conflicting targets. The check suite includes hot thermal
emission at nonpositive fields and repeated refusal throughout the unqualified
weak-positive-field interval. These are small validation fixtures, not a
60 ns physical breakdown run. Outputs use the versioned v8.0.1 directory.
Both restore modes require explicit evidence: this coordinator retains the
v5 native configuration even when using the corrected v8 archive. Genesis
restores directly; nonempty checkpoints still use verified replay. The test
asserts both paths, including chained replay, instead of inferring the mode
from the library version. Native v8's own nonempty direct-restore qualification
is separate. Run targets stop after the first accepted step reaching them;
non-grid-aligned targets may overshoot by a step. Exact endpoint landing is not
claimed by this correction.

Limitations

The default Hall binding is still v3.3 and solid-only. The SS304 phase owner,
full reaction network, wall responses and circuit coupling are not integrated
by changing the native library. Starting plasma from a 20 ns vacuum pre-run
still omits earlier charged-particle/source history and requires a physical
qualification. The small field box and inferred drive require further
calibration. Full plasma shielding and breakdown remain unqualified.

Old checkpoints are retained in their original runtime. This corrected owner
identity does not silently admit them; migration needs explicit comparison of
owners, source physics and accepted histories.
