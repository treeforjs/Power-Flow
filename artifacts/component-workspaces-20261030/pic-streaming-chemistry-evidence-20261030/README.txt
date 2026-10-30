Streaming generic-v6 chemistry evidence, separate primitive v1 — 2026-10-30

This processes evidence from the existing frozen generic-v6 provider. It does
not evaluate or change a physical kernel, add channels, choose a chemistry
policy, change the1keV seam, select collisions, consume RNG or own accepted
particles/fields/history. Generic-v6 and the native validators remain unchanged.

source/streaming_coverage.hpp provides Stream(reference_coverage,bindings,limits),
observe(Query,coverage), commitment(),snapshot(). The reference fixes the exact
40-key inventory and policy. Every observed record must contain that complete
inventory and the unchanged provider policy, including all unapplied/missing
rate rows. A removed key, unknown-to-zero conversion, promoted diagnostic event,
changed policy, invalid operand, ordinal gap or resource limit refuses without
publishing a partial evidence state. Finite source values remain caller/provider
observations: this primitive cannot authenticate a caller or certify a kernel.

Category partition per key (never per distinct energy):
  selected_positive; selected_zero; selected_unsupported;
  unselected_positive; unselected_zero; unselected_unsupported;
  unqualified_no_rate.
Each has an exact uint64 query count and exact input-weight sum, energy extrema,
positive-weight energy extrema, finite known-value extrema, first/last ordinals
and weighted energy/value moments. Zero-weight unsupported queries are counted
and committed, even though their physical evidence weight is zero. Unknown values
stay null; no unsupported weight or positive unmapped rate is renormalized away.
The conditional-model, unapplied-observable and missing-rate lanes are separate,
nonadditive views; in particular the2a1 observable/contract is intentionally not
two physical events. eV*m2 stopping moments stay distinct from m2 cross sections.

Weights are caller-declared evidence measures. Repeated particle/pair queries
count repeatedly, so these totals are not automatically a unique population,
number density or EEDF. In the captured fixture each old electron is queried once
against a stationary ground-state H2O target using its captured pre-kick u. No new
EEDF or physical history is produced. A future owner must choose and bind the
measure and actual query census explicitly; no hidden slab/depth/volume factor
is introduced here.

Exact sums reuse the sealed, unchanged ExactParticleCount header (2098-bit frozen
libcrypto integer representation), including its existing finite range. No new
count backend is introduced. Input-weight sums are mathematically exact for the
admitted binary64 operands, including weights below the displayed total's ULP.
Moment sums mean EXACT SUMS OF PER-TERM BINARY64 ROUNDED PRODUCTS, not exact-real
multiplication: each w*E and w*source_value is formed once in binary64. Positive
product underflow, product overflow and non-finite rounded totals refuse before
publication. Unknown source moments are null. The raw values/weights are still
committed, so a future wider moment representation can be separately qualified;
no floor/clamp or zero is silently substituted. Round-to-nearest and preserved
binary64 subnormals are checked; fast-math builds are unsupported.

Bounded memory:
  No map/set/vector grows with query count or distinct energies.
  At most64 keys,7 fixed category slots each; the actual provider uses40/280.
  Stats are staged in one bounded successor copy and published by no-throw swaps.
  Reference and query text are bounded32KiB; each encoded query <=256KiB.
  Default1M queries and16GiB encoded query bytes per receipt; hard query max100M.
  Caps are explicit refusals, never dropped observations. Header bytes are
  separate from encoded_query_bytes. Process/source-provider allocation overhead
  is measured independently; record_bytes is not a process RSS claim.

Deterministic content commitment, SHA256:
  seed: ASCII CHEM-EVIDENCE-V1 then BE-u64-length + canonical_header_json bytes.
  each query: ASCII Q, BE-u64 ordinal, length-prefixed id bytes and64-character
  context digest, raw IEEE754 energy bits (BE-u64), raw weight bits (BE-u64), then
  length-prefixed full canonical provider JSON bytes.
The canonical JSON is frozen nlohmann dump() with ordered object keys; array order,
source reason text and every supported/unsupported value are included. Exact
component/arithmetic/provider/data/kernel/mechanism/missing-ledger/query-input
hashes and limits are in the header. The input sequence is ordinal0..N-1.
Chunk boundaries do not affect results; stream order does. Snapshots include the
canonical header, so independent framing can be reproduced. Digests provide
integrity/content identity, not external provider authentication.

A receipt stores bounded summaries plus the full-stream commitment, not a growing
list of per-energy JSON. Exact unsupported-query retrieval requires replayable
input operands and the pinned provider/generator. This fixture retains its full
original CSV and generation manifest. A future dynamic owner must retain such
operands externally or provide equivalent checkpoint/event replay; a digest
alone is not a replacement for missing history. This is not a rolling native
checkpoint or accepted chemistry completion and cannot weaken those validators.

Tests compare the seven existing frozen pointwise fixtures and a16-projectile
legacy-admitted candidate batch; candidate scalar bits and scheduler RNG remain
unchanged. An independent Fraction oracle checks exact weighted categories,
moments and extrema; an independent hashlib oracle checks wire commitments.
The46,726-query captured pass demonstrates constant retained state and records
actual peak RSS/time. A second process replays the same list and compares the
complete receipt. Full-provider evaluation is included in timing; this is not a
production-throughput qualification or a new coupled test.

Run run_checks.sh under the exclusiveCPU7/nice15/frozen/one-thread wrapper it
supplies. It uses existinglibcrypto and no new dependencies. The source/receipt
integration proposal and exact remaining limits are in the output report.
