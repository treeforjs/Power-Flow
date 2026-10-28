# Pre-push Python verification, 2026-10-28

The complete suite is **not green**. These runs used the current worktree and
removed only the invoking process's stale `_editable_skbc_mhdlab` finder before
importing pytest. No global Python installation or Git state was changed by
verification. Testing stopped after the publication/pause instruction.

| Run | Passed | Failed | Errors | Skipped | Duration |
| --- | ---: | ---: | ---: | ---: | ---: |
| Full suite, before the approved harness fix | 3111 | 172 | 128 | 7 | 130.15 s |
| Affected-session selection | 631 | 6 | 5 | 2 | 25.15 s |
| Boundary quadrature alone, original fixture | 52 | 0 | 0 | 0 | 0.25 s |
| Fixed boundary fixture plus water-source tests, with `mhdlab.planar_release` deliberately preimported | 61 | 0 | 0 | 0 | 2.50 s |

The affected selection covers MFEM surface/metadata/comparison, planar/Temkin
boundary interfaces, neutral inference/handoff/comparison, optics, EOS,
water-source binding, rendering and conservative transport. It also includes
the existing coverage-inventory and Temkin convergence CLI tests.

## Diagnosed failures

- **Four affected failures:** existing extreme-range coverage-inventory tests in
  unchanged `mhdlab/coverage_inventory.py`. On this Windows/NumPy installation,
  `np.longdouble` has the same 64-bit range as float64. Intermediate positive
  products can underflow before representability checks, or overflow before a
  later factor restores a finite result. Two expected underflow rejections and
  two finite-intermediate-range cases fail. This is a real pre-existing
  portability issue, not a fitted monolayer result.
- **Two affected failures and five setup errors:**
  `tests/test_temkin_planar_convergence_cli.py` invokes the existing convergence
  CLI in subprocesses that retain the stale editable finder. They fail to import
  `mhdlab.emission_time_moments`. Removing the finder only from the parent pytest
  process does not alter those subprocesses. Many additional full-suite CLI
  failures have the same stale-checkout import cause.
- **Full-suite import-boundary fixture:** 52 setup errors arose because the
  original boundary-quadrature fixture required the entire pytest process to
  have imported no `mhdlab` modules. Other test collection legitimately imports
  them. The approved fix now verifies the direct-core import in a pristine
  subprocess, preserving the XML properties and semantic guarantee without
  deleting global modules. The 61-test combined run above verifies this fix.
- **Full-suite reviewed-source hash checks:** Mills source-manifest, Flynn
  transport and calibration data are checked out as CRLF while their reviewed
  digests bind LF bytes. Representative Mills, Flynn and both calibration CSV
  digests match exactly after newline-only CRLF-to-LF conversion. Git reports
  `i/lf w/crlf`. No scientific hash check or source-data bytes were changed by
  this verification task. Other recovered-data hash failures remain in the
  full-suite result.
- **Isolated legacy-retirement subprocesses:** Python `-I` cannot import `h5py`
  in this host's isolated environment. This is separate from the stale editable
  finder and was not repaired here.

The full suite was not rerun after the one approved test-harness fix. No other
production fixes, environment changes, staging or commits were performed by
this verification task. The passing affected tests do not establish that the
full historical suite or the physical spectral model is complete.

## Preserved local evidence

Logs remain in `.cache/prepush_pytest_full_20261028T2255.log`,
`.cache/prepush_pytest_affected_20261028T2258.log`,
`.cache/prepush_pytest_boundary_isolated_20261028T2258.log` and
`.cache/prepush_pytest_boundary_combined_20261028T2301.log`.
The affected selection, per-file failure grouping, JUnit XML and exact
newline-identity comparisons are preserved alongside them as
`prepush_pytest_affected_selection.json`, `prepush_verification_diagnosis.json`,
`prepush_affected_20261028T2258.xml`,
`prepush_boundary_combined_20261028T2301.xml` and
`prepush_source_identity_newlines.json`. These local cache records are not
claimed to be bundled by this concise note.
