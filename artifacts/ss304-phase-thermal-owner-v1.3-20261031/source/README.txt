SS304 phase thermal owner v1.3 -- transaction corrections, 2026-10-31

Derived only from sealed v1.2. The phase law, phase/melt/ejection kernels and floating-point arithmetic are unchanged.
Namespace: ss304_phase_v1_3. Owner configuration/state/input/quote/checkpoint identities use v1.3; cell and layer
Merkle section definitions stay as in v1.2. v1.2 checkpoint files are refused by this version; no silent conversion.
V1.2-BASELINE-SHA256.json records all copied baseline inputs. Frozen dependencies are pinned before/after each build
and run by run_owner_v1_3.py (including historical material/session/NVRTC and v1/v1.1/v1.2 headers).

Corrections
1. prepare invalidates previous proposals before reusing engine work slots, including when the replacement refuses.
   Only the most recently prepared proposal may commit. Accepted numerical state is unchanged until commit.
2. Back pressure supplied in StepInput stays in that proposal and publishes with commit. A refusal or abandoned
   proposal never changes accepted pressure. quote(dt, pressure) evaluates proposed pressure without publishing it;
   set_back_pressure is an explicit accepted-boundary update and invalidates pending proposals.
3. Every accumulation invalidates prepared proposals, preventing their commit from clearing newer field inputs.
   A pending field window must be consumed before a direct thermal step. A pending window cannot be saved or loaded
   over. State identity binds pressure and every pending per-EM-cell energy word, force-impulse word and invalid-input
   flag. Accumulated-step input identity binds the same full snapshot, including spatial distribution.

GPU accumulator identity currently downloads the small EM-cell accumulator arrays on demand when hashing a nonempty
window (state_sha256 and accumulated prepare). It does not change accumulation arithmetic or introduce a per-particle
loop. There is no performance claim in this correction; future device hashing must preserve this exact encoding.

Qualification
test_baseline_defects compiles against the unmodified v1.2 headers and checks that the three original defects
are reproduced. New tests/transactions.hpp is used by test_transactions (CPU) and test_gpu transactions (CPU/GPU). It covers stale
commit after different prepares and after a failed replacement, pressure refusal/commit/quotes, explicit pressure
invalidation, extra accumulation after prepare, exact combined-window consumption, direct-step refusal, pending save
and restore refusal, pressure and consumed-window checkpoint continuation, and energy-distribution/force/flag digest
identity. Existing law, owner, v1.1/v1.2 and historical v1.1 physical-regression suites remain available unchanged in
scope. Actual results and source manifest are recorded separately after tests; the old README's results do not qualify
this release. See README-v1.2.txt for inherited physics and limitations.

Queue all builds/tests through SimQueue with explicit --mem; run_owner_v1_3.py verifies actual job cgroup, allocation,
affinity and owner/borrow leases. CPU work may borrow cores4/5 through the queue. GPU work needs a queued GPU claim
and a powerflow core6/7. Frozen userspace/dependencies are unchanged. No integration with Hall is performed here.
