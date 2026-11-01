Coordinator CPU performance candidate v8.0.2-cpu -- 2026-10-31

Status: final c2 functional qualification and all five deterministic fixture
comparisons passed (SimQueue316/326). Timing remains unqualified; job327 is
queued. See RESULTS.txt. Do not use study speedups as evidence for this build.

Baseline: ../pic-water-common-commit-v8.0.1-20261031, qualified build c2,
SHA256 3e5d52b2eb93d15fdfbaa79a7a3541911fa987905a9bbb55c14f9411d7e86072.
All source/ configuration/ owner identity strings are unchanged. This candidate
changes host computation only, not physics, precision, particles, schedules,
refusal guards, or the default v5 native configuration compatibility binding.
That binding still uses direct genesis restore plus verified nonempty replay.
It is not the fully integrated native-v8/chemistry-v3.2/thermal/circuit system.

Changes:
- P1: vendored collision_scheduler.hpp uses an exact product-error helper for
  compensated cross products. Dekker's 2^32+1 split is admitted only on x87 with
  64 significand bits, extended-precision nearest rounding, masked exceptions,
  operands in [2^-4000,2^4000] magnitude, and x == a*b. Else the original fma
  is called. The bounds leave every split and product well inside the normal
  extended exponent range; the fixed split matches the 64-bit significand.
- P3: vendored water_parent_bridge.hpp caches immutable parent linear-table and
  kernel metadata per channel construction, avoiding per-query JSON allocation
  and lookup. Expression ordering, threshold/support gates and error messages
  follow the pinned collision_data.h. Unsupported models or malformed cache
  metadata retain the original query path. The dataset remains a const shared
  object bound by the unchanged data/kernel digests; mutating it through a
  separate alias would violate that frozen input contract.
- No P2 string/container rewrite, CPU ISA flags, fast-math, lower precision,
  altered convergence target or particle reduction.

The vendor manifest records baseline, study, and candidate hashes. All vendor
files enter the before/after source pins; sealed_pins verifies their fixed hashes
and that the original source headers remain unchanged. Compiler dependencies
must contain exactly the vendored scheduler/bridge/collision-data headers.

P1 contract: exact finite value and signed-zero parity; NaN classification rather
than NaN payload identity. Splitting can set FE_INEXACT on an exact product, so
floating exception-flag identity is not claimed. There are no floating exception
flag readers in the coupled callers reviewed. Unmasked exceptions and changed
rounding/precision modes use fma. Tests verify the mode is not changed.

Tests (all numerical execution through SimQueue, one core, RAM cap,
frozen arch-base-20261028; functional checks may use an explicitly borrowed
core with --on-preempt freeze; timing requires a normal exclusive owner core;
no shared frozen dependency edits):
  tools/qualify_candidate.py LABEL
    targeted product-error and cross-section boundary/refusal equivalence,
    full generic build, existing 53 checks, 12 continuation cases, and six
    postrun integrity fault cases, including actual .inc source tampering.
  tools/compare_candidate.py LABEL TAG parity
    final-build check evidence and deterministic qual20/dense-adaptive/load100k/
    load1m/pilot2000 comparisons against qualified v8.0.1 c2.
  tools/compare_candidate.py LABEL TAG timing
    three interleaved baseline/candidate pairs, 1e6 manufactured carriers,
    the same 200um mesh, three accepted steps; compare steps 2 and 3 separately
    from genesis. Requires a normal one-core --exclusive-timing job and rejects
    measured SMT-sibling load above the recorded threshold. Existing workers
    are not stopped. This is a CPU fixture timing, not a full-run prediction.

Output root:
/home/treefor/sim-workspace/powerflow/outputs/pic-water-common-commit-v8.0.2-cpu-20261031
