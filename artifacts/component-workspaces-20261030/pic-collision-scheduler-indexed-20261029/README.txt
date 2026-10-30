Indexed local collision candidate builder

This isolated header replaces the frozen scheduler candidate enumeration only.
The reference source SHA256 is a4cae6d46fb9b7b99244149b543388439dbf717a62dd69541f7d91b753b2ce6b.
The root scheduler and frozen bridge are unchanged. No native particles, physical
cross sections/product kernels, RNG clocks, fields, GPU, or accepted history are
changed by this package. This is scheduling/data-structure qualification.

Use bash run_checks.sh. It uses testing-core7.lock, CPU7, nice15, one thread,
the frozen Arch userspace, a1GiB address-space cap and60 CPU-second cap. It runs
the unchanged47-check scheduler suite, exact reference comparisons, cap/coverage
checks, and bounded end-to-end builder benchmarks. No dependencies are installed.

Enumeration and semantic contract

All cells, particles and channels are validated in the original relative order,
after the new input-size/byte preflight. Cell volume is found with an ordered
cell index. Per-cell target pointer lists retain the ORIGINAL target order.
Projectiles retain input order, and matching channel lists keyed by reactant
states retain ORIGINAL channel order. Every same-cell pair still consumes the
pair budget and computes the unchanged relative kinematics, even with no eligible
channel. This preserves unmatched-pair rejection and all support-query behavior.
Every matching channel is queried, including zero relative flux. No physical
pair/channel is sampled, discarded or approximated. Final canonical sort is
unchanged. Candidate identities/order and all numeric fields are compared by
binary64 bits; sizeof/padding bytes are deliberately not compared.

The public candidates() retains ordinary existing direct-call forms. It adds
two optional trailing arguments: CandidateInputLimits and CandidateWork*.
The latter is diagnostic only and is reset at entry. The extended function type
is not an ABI/function-pointer compatibility guarantee. Inputs, metadata strings
and callable ownership must remain immutable during the call. The indexes store
string views and particle/channel pointers with call-local lifetimes.

Default resource policy and exact counters

Original limits remain:1,000,000 local pairs,1,000,000 matching-channel queries,
1,000,000 positive-propensity candidate rows. Zero limits permit zero work and
reject the first corresponding item before kinematics/query/row allocation.
The comparison-before-increment form also avoids size_t wrap at SIZE_MAX.

New defaults:1,000,000 projectiles,1,000,000 targets,1,000,000 cells,16,384 channels;
256MiB aggregate input string content;256MiB candidate copied string payload.
Count and aggregate input-byte checks precede internal map/set/vector allocation.
Input bytes include packet ID/state and channel ID/reactant states/data/kernel
digests; callable captures and input vector storage were already caller-owned.
Candidate payload counts the five copied strings plus five NUL terminators per
row, checked before constructing it. Row objects/vector capacity are separately
limited by candidate_limit. All byte sums use remaining-budget comparisons before
addition. These are logical payload/count caps, not allocator/RSS guarantees.
The harness additionally caps process address space and CPU time. A host must
choose explicit limits appropriate to its complete inventory; never truncate
inputs or silently skip a physical channel to fit a cap.

For inputs within new caps, validation/support exception messages, query order
and cap boundaries match the legacy builder. Oversized malformed requests may
now fail at preflight before a legacy validation error. Any cap failure aborts
the whole request; no partial candidate vector is returned. Previously executed
query callbacks may have occurred before a local pair/query/row/output-byte cap,
as in the original builder. Physical cross-section queries must remain read-only.

CandidateWork counts successful input string preflight bytes, cell membership/
index lookups, target/channel index entries, visited projectiles/local pairs,
kinematics calls, channel-group lookups, actual channel queries, and emitted rows/
copied payload. Counters at exceptions describe the permitted work prefix, not
accepted physical state. Input preflight failures report zero index/query work.
Counters do not count map comparator instructions or bound arbitrary callback
cost. Count/byte caps bound validation/index size; process limits remain useful.

Complexity and interpretation

Validation/index setup is bounded by input counts and ordered-map/set operations.
There is no global projectile x target scan or per-pair linear cell-volume scan.
The pair loop visits sum_cells(N_projectile_cell*N_target_cell), with a channel
state-group lookup per local pair and only matching-channel query traversal.
Canonical output sorting remains. A dense cell still requires all local pairs,
and its original pair/query/row caps still apply. Index setup can dominate very
small populations; this package does not promise universal speedup or full-PIC
speedup. Benchmarks use manufactured finite cross sections to isolate enumeration,
not physical water event validation or collision/plasma evolution.

Qualification and integration

The output directory contains results.json, counter-limits.json, benchmark-summary.csv,
resource receipts, input/output manifests, exact unified integration patch and
handoff.json. The patch is checked by applying it to an isolated copy of the exact
frozen legacy source. relative(), Remainder, Clock/Proposal/Scheduler code are
unchanged. The unchanged original contract suite still passes47 checks.
The expanded independent suite checks300 varied cases and their permutations,
embedded-NUL/long IDs, dense/sparse/mixed states, exact candidate bytes/query traces,
zero-flux support errors, unmatched-state kinematic errors, and exact/one-below
input/pair/query/row/payload boundaries. It also exercises a100k-by100k disjoint
population at zero local-pair/query/row budgets, with actual wall/RSS measurements.

Root integration must rebuild/rebind the dependent kernel manifest to this new
scheduler source in a new integration prefix. Do not overwrite the old frozen
bridge input or its receipts. Authoritative stage/state/mesh checks, iterative
compensated reservoir, native mutation/current and atomic publication remain
separate unfinished integration. No density output is generated; any future
number-density display remains cm^-3.
