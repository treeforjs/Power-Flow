Generic v6 water bridge with streaming evidence, separate wrapper v1

Scope
The new entrypoint is raikou::streaming_water::Bridge::candidates in
source/streaming_water_bridge.hpp. Generic v6 and stream-evidence v1 are copied
unchanged under provider/ and inputs/. No accepted particle owner, native event
admission, transport, field, trajectory, source, heat or common-time evolution
is implemented or exercised here. Complete-mechanism readiness remains false.

The returned rows are the original physical Candidate objects. Evidence uses a
new schema; it is not the old per-distinct-energy coverage array and cannot be
silently substituted into a native completion/history schema. The old generic
entrypoint, including its default 256/hard 1024 evidence ceiling, remains intact.

Preserved physical operations
Physical channel creation and both product-sampling overloads delegate directly
to the immutable generic v6 Bridge. The indexed builder in observed_candidates.hpp
has the original arithmetic, validations, ordering, work/allocation caps and
sorting body byte for byte, with only begin-pair/channel-return/end-pair evidence
hooks. The runner verifies that equality and the exact coverage projection and
support predicate before compiling. The physical input checks remain unchanged.
No extra RNG draws, rate changes, channel splitting, seam smoothing, invented
zeros, attachment relocation or weakened unsupported-kernel guard is added.

Evidence processing
All eight physical channel callbacks are evaluated in original order. Actual
active sigma returns seed a bounded lookup memo; an inactive comparison zero is
never substituted for a raw positive/unsupported source observable. Remaining
raw diagnostics are queried once per raw ID in a pair. One last exact-bit-energy
coverage record can be reused, after its active values are checked against the
current physical returns. No container grows with distinct energy. Stream v1
retains the same 40 keys / 280 category slots and hashes every complete provider
record, including reasons, nulls, units, missing contracts and unapplied positives.

The per-query weight is projectile-carrier weight for each evaluated e-H2O pair.
Repeated pairs are not a unique population census, density, EEDF, event count or
rate-weighted integral. Raw source moments retain their units (including eV m2
where appropriate). Unapplied diagnostic views are nonadditive and cannot be
summed into physical rates. Evidence storage is bounded independently of query
count; indexed input/output candidate storage still scales within original caps.

Publication and replay
The recorder, memo, last record and candidate rows are private to one call. A
failed physical or evidence preparation returns no receipt or rows and changes
no Bridge state. The caller must keep a successful batch private until its own
accepted transaction, if integrated later. accepted_native_completion is false.
The input descriptor binds the host query-input SHA, scope-context SHA, actual
wrapper/observed-builder SHAs and explicit weight meaning. Each pair query ID
hashes packet ID/state/cell/weight/mass/proper velocity, cell volume and actual
relative energy/flux. Physical positions and full authoritative census/stage/
geometry must also be bound by the host input/context; the pair ID is not their
substitute. Retain the exact input and generator/provider for replay; a digest
does not reconstruct individual queries. The frozen test generator and complete
copied provider inputs are included in this release.

Limits
Default pair/query/row caps remain 1,000,000 / 1,000,000 / 100,000. Each matched
pair performs all eight physical queries; all work remains counted. Generic
populations are capped at 100,000 each, packet IDs at 4096 bytes and their total
ID/state payload at 64 MiB. Original indexed caps include 1,000,000 cells,
16,384 channels and 256 MiB each for input/candidate string payload. The fixed
physical channel set has eight entries. These are logical caps, not an allocator
RSS bound. Stream v1 adds default 1,000,000 pair records, 16 GiB encoded bytes,
256 KiB per record, exact finite count/moment guards, FE_TONEAREST and preserved
subnormals. These evidence limits are explicit additional refusal domains.

Two tested examples physically admitted by v6 are refused by stream evidence:
20 eV with projectile weight denorm_min (positive source moment underflow), and
20 eV with DBL_MAX (moment overflow). They are not physical-builder parity.
Other stream finite-total/environment/binding/budget refusals remain applicable.
No lossy moment, unsupported-positive omission or normalization is used to admit
those inputs. A later policy/backend change would need a separate version.

Qualification and timing
3,781 assertions pass, including complete candidate field bytes/work counters,
14 exact physical-error messages, support edges, moving target, 15 candidate
product/RNG plus explicit-variate comparisons, private-state failure controls,
and 32 competing scheduler seeds for each of 42 admitted fixtures. Complete
evidence matches original pointwise coverage for 4,139 pair projections.
The 2048-distinct-energy case passes with 7,687 rows / 16,384 channel queries;
old v6 refuses it at its hard 1024 evidence cap. A second process gives identical
candidate/evidence receipts for both broad cases. The evidence has 280 category
slots for both distinct and repeated energies.

Measured wrapper call time on CPU7 was about 0.500 s for 2048 distinct energies
and 0.150 s for 2048 equal energies. Full broad qualification was 1.97 s and
56,836 KiB peak RSS; its independent references, stored test outputs and RNG
comparisons are included. These timings are bounded fixture observations, not
production throughput. The earlier 46,726-query full diagnostic pass took
8.9 s and is NOT suitable for every production step. Hashing complete records
still costs work even when source values are reused. No production timing or
chemistry/native/joint qualification follows from this storage improvement.

Reproduction
Use ./run_checks.sh from this runtime (or its extracted archive). It uses the
existing CPU7 lock, affinity 7, nice 15, frozen userspace, one thread, no visible
CUDA devices, 1.5 GiB address and 180 CPU-s limits, and subprocess timeouts.
The copied provider manifest is verified before compilation; no live source is
needed. Optional PIC_STREAMING_WRAPPER_OUTPUT selects a separate output path.
Run replay_and_summarize.py under the same lock/frozen/CPU contract to reproduce
the independent-process receipt comparison. No dependencies are installed.

Remaining physical limits
The conditional eight-channel comparison, missing positive 2a1/core/multiple/
fragment-continuation/low-energy angular/inelastic coverage, conditional OHv
family interpretation and attachment exact-site requirement are unchanged.
The 1 keV model seam is unchanged. There is one competing candidate set; no
sequential full-duration channel loops are introduced. Native multi-channel
family/completion admission and the full chemistry/source/field model remain
separate, unimplemented or unqualified integrations under their existing audit.
