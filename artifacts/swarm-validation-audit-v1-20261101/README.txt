Swarm completeness audit -- isolated, not deployed
================================================
Only this directory is owned by this work. Live chemistry v3.2/v3.2.1, the
existing swarm outputs, and gas-target-owner-v1 were not modified. No duplicate
swarm or simulation was launched. Tests use synthetic output records only.

Fixed integrity/completeness defects
-----------------------------------
1. The existing writer prints E/N with six significant digits. The old summary
   matched those labels to unrounded measured coordinates at 1e-9 relative
   tolerance, silently dropping valid measured points. The contract now maps
   each of the 36 expected coordinates to its exact six-digit serialized label;
   it refuses collisions/duplicates/unexpected coordinates. This does not
   loosen any physical measurement tolerance.
2. Missing, malformed, duplicate, wrong-source, wrong-seed, nonfinite, or
   incomplete run records cannot produce a complete summary. Exactly six seeds
   of each set, 36 points per seed, all required quantities and the e/H2O process
   set are required. All 50 metric/data rows must exist; 43 are judged under
   the existing rule. Zero judged rows refuses.
3. Summary source identities, sample counts and reconstructed mean/sample SD
   are validated before comparison. The observed max(2 SD_reference, 0.01)
   rule is unchanged. It is labeled as an observed worker criterion whose
   authorization and independent physical validity this audit does not prove.
4. Complete validation summaries are published through a flushed temporary
   file and atomic replace. Refused summaries do not replace prior output.
5. The staged runner holds an exclusive output claim, writes a unique temporary
   JSON file, and publishes only after exit 0 and exact source/configuration and
   output checks. Resume requires a matching completion receipt; parseable JSON
   alone is insufficient. Cancellation and the enclosing qualification timeout
   terminate and reap this runner's own child process groups.
6. Technical sealing requires all ten stages in order, immutable source/input
   identities and tested binary identities, complete per-job proofs, and the
   canonical 33-job plan. It recomputes the comparison from current completed
   outputs, so a stale passed comparison cannot qualify replacement raw files.
   It explicitly reports production_qualified=false and physics_model_validated=false.

Files
-----
source/swarm_validation.py       strict read-only input/summary validator
source/swarm_multiseed.py        patched existing comparator; scientific rule unchanged
source/swarm-contract.json       pinned source/spec/blob hashes and execution plan
tests/test_validation.py         baseline reproductions and fault regressions
tests/test_sealing.py            stage, build, job and comparison provenance faults
tests/test_runner.py             actual child publication/resume/termination tests
source/qualification_integrity.py complete bound-stage technical receipt validator
source/swarm_job_integrity.py    exact successful job proof/publication helper
source/run_jobs.py               claimed, atomic, cancellation-safe staged runner
source/run_checks.py             stage identity capture and graceful timeout cleanup
source/make_seal.py              fail-closed technical seal publication
adoption/swarm-completeness.patch  exact-base adoption patch (not applied live)
adoption/base-sha256.json        base/patch/output identities
baseline/                       local-only source snapshots and measured extracts
evidence/                       tests and live read-only classification

Adoption
--------
Wait for the live worker's stable source/output boundary. Verify the exact base
hashes of all four existing files and absence of all four new paths against
adoption/base-sha256.json; then review/apply the patch in a new staged version.
git apply --check is part of the queued regression suite. Do not force a changed
base. The contract intentionally pins the existing v3.2 spec/blob and writer;
v3.2.1 regeneration needs a newly reviewed contract and qualification, not a
silent replacement of the hashes. Measured/spec extracts in baseline/ remain
local and must not be committed. The adoption patch includes hashes, not those
source tables.

The staged runner fixes the publication/resume protocol for newly executed jobs.
Existing live writers do not use its output claims, so use a new output tree after
a coordinated handoff; do not mix old and new writers. The live runs were left
untouched. Raw simulation JSON alone still lacks executable/exit provenance.
No retrospective completion receipt is created for those files. Complete content
can be inspected separately, while technical sealing requires fresh valid proofs.
The original leased-core/resource rules remain unchanged. This is not permission
to run an unqueued qualification, take other cores, or adopt scientific changes.

Run the read-only classifier:
  python3 -I source/swarm_validation.py <raw-multiseed-directory> \
    --contract source/swarm-contract.json --source-root <pinned-v3.2-runtime>
It exits 2 for incomplete/invalid input and never repairs or changes live files.

Qualification status is recorded in evidence/test-summary.json when the
scheduled frozen one-core test job finishes. No physical swarm claim is made
by synthetic output-fixture tests.
