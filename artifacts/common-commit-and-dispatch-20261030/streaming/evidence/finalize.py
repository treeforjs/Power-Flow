from pathlib import Path
import json,hashlib,tarfile,os
P=Path('/home/treefor/projects/runtimes/pic-streaming-chemistry-evidence-20261030');O=Path(__file__).parent
if os.sched_getaffinity(0)!={7} or os.getpriority(os.PRIO_PROCESS,0)!=15:raise RuntimeError('CPU7/nice15 required')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for n,want in json.loads((P/'inputs/frozen-source-manifest.json').read_text()).items():
 if sha(P/n)!=want:raise RuntimeError('sealed input/source changed: '+n)
small=json.loads((O/'small-results.json').read_text());large=json.loads((O/'capture-results.json').read_text());oracle=json.loads((O/'independent-oracle.json').read_text());replay=json.loads((O/'replay-results.json').read_text());receipt=json.loads((O/'capture-receipt.json').read_text());memory=json.loads((O/'memory-checkpoints.json').read_text())
if not all(r['status']=='passed' for r in [small,large,oracle,replay]):raise RuntimeError('checks incomplete')
if receipt['queries']!=46726 or receipt['observations']!=1869040 or receipt['retained_category_slots']!=280:raise RuntimeError('capture inventory differs')
if not replay['whole_receipt_byte_exact']:raise RuntimeError('replay differs')
proposal={'schema':'streaming-chemistry-evidence-integration-proposal-v1','role':'diagnostic rate/evidence primitive only; no accepted chemistry owner','implemented_source':str(P/'source/streaming_coverage.hpp'),'proposed_fork_only':True,'generic_v6_unchanged':True,'current_generic_entrypoint_still_has_old_map_cap':True,'integration_steps':['Create a new versioned generic wrapper; preserve the sealed v6 source and native validators.','Replace only the unique-energy JSON map/evidence cap. Invoke one evidence hook per declared source query or physical pair, including zero-flux/unsupported queries. Keep all existing sigma/kinematics/rate/weight/RNG/selection operations unchanged.','Bind actual immutable input/provider/component/arithmetic/data/mechanism/missing-ledger/phase contexts and declare the query-weight measure. Repeated pair queries are not a unique electron census.','Keep the evidence stream private to the candidate/exposure. A failing physical builder discards its evidence too. Freeze the finished evidence receipt/hash before a common preflight; no independent scientific-state publication is introduced.','Retain replayable query operands externally, or qualify exact checkpoint/event replay. Full-stream hashes alone cannot reconstruct missing dynamic history.','Prove complete candidate/rate/RNG equivalence and error rollback on the future fork. Add production-scale cost admission before using full-provider diagnostics inside every step.'],'unchanged_physics':['eight conditional v6 kernels','no new channels or missing-rate zeros','exact1keV subset policy and its known seam','single-channel native completion and other validators','collision-site/current restrictions','full-mechanism readiness false'],'remaining_limits':['The old generic-v6 candidates() remains capped until a separately reviewed fork integrates this primitive. The46726-query test evaluates scalar evidence, not a full large candidate population.','Default1M queries/16GiB encoded query bytes per receipt; hard100M queries; <=64 keys; <=256KiB encoded record and32KiB text. All caps refuse without partial publication.','Nonnegative binary64 evidence weights/source values; exact weight sums reuse2098-bit ExactParticleCount. Moment sums are sums of individually binary64-rounded products. Positive product underflow, overflow or nonfinite rounded totals refuse.','The receipt commits all provider JSON, including unsupported reasons/values. It keeps bounded summaries rather than all raw pointwise JSON; replay requires bound operands/provider/generator.','Evidence validates record shape and content commitment, not an external provider or native epoch. Caller digests must be verified by the host.','Only the frozen v6 provider/layout and recorded floating environment are qualified here. No full chemistry, EEDF, field, trajectory, native mutation, restart-owner or host/GPU deployment qualification.','Full-provider evaluation included in measured8.81CPU-s pass is not production-throughput qualification; typed hooks should reuse already computed query results where possible without changing any value.']}
(O/'integration-proposal.json').write_text(json.dumps(proposal,indent=2)+'\n')
archive=O/'streaming-chemistry-evidence-source-v1.tar.gz'
with tarfile.open(archive,'w:gz') as tar:
 for f in sorted(P.rglob('*')):
  if f.is_file() and 'build' not in f.parts:tar.add(f,arcname=str(f.relative_to(P)),recursive=False)
(O/'source-manifest.json').write_text(json.dumps({str(f.relative_to(P)):{'sha256':sha(f),'bytes':f.stat().st_size} for f in sorted(P.rglob('*')) if f.is_file() and 'build' not in f.parts},indent=2)+'\n')
(O/'report.txt').write_text('''Streaming chemistry coverage evidence primitive v1 — 2026-10-30

Result: bounded evidence processing passes all existing small-fixture controls,
independent exact-sum/framing checks and the full captured46,726-query replay.
This is a separate component/proposal. Generic-v6 remains immutable, including
its old256/1024-energy map cap until a new fork integrates this primitive.
No physical kernel, channel, chemistry policy,1keV seam, RNG, native validator,
field, source or accepted owner was changed. No new EEDF or physical evolution.

Source/runtime: /home/treefor/projects/runtimes/pic-streaming-chemistry-evidence-20261030.
source/streaming_coverage.hpp and README.txt specify the API and wire format.
integration-proposal.json lists the exact future fork boundary and remaining
limits. The new receipt schema is streaming-generic-v6-evidence-receipt-v1; it is
not a native event/completion/restart schema and cannot weaken existing admission.

Semantics and support partition
Each declared query has a contiguous ordinal, opaque bounded ID, context hash,
binary64 energy and nonnegative evidence weight. The frozen provider's complete
40-key inventory is required for every query. Seven mutually exclusive bins per
key distinguish selected positive/zero/unsupported, unselected positive/zero/
unsupported and unqualified-without-rate. Counts are exact uint64; weight sums
are exact admitted binary64 sums using the unchanged sealed ExactParticleCount
header, not a new count backend. Zero-weight unknown queries remain counted.
There is no support-weight/tail renormalization and no unknown-to-zero conversion.

The40 keys occupy conditional-model, unapplied-observable and missing-contract
lanes; these are nonadditive evidence views. The2a1 observable and its explicit
missing-rate contract are not two events. Stopping moments retain eV*m2 units,
separate from m2 cross sections. Each category records finite energy/value extrema,
positive-weight energy extrema and first/last query ordinals. Unknown values and
their source moments remain null. Query weights are not automatically a unique
population or density: repeated queries count repeatedly under the declared
measure. No observation slab/depth/volume convention is silently introduced.

Weighted moments are exact sums of PER-TERM BINARY64 ROUNDED PRODUCTS, not exact
real multiplication. Positive product underflow, product overflow, count-range
or finite-total overflow refuses atomically before evidence publication. The
rounding/subnormal environment is checked; fast-math builds are unsupported.
Raw values/weights are still committed, and any wider moment representation
would require a separately qualified change. These explicit numerical admission
limits do not silently floor a physical rate or drop an unsupported query.

Commitment and reproducibility
The SHA256 commitment includes a versioned, length-framed canonical header with
component/arithmetic/provider/data/kernel/mechanism/missing-ledger/query-input
hashes, policy, layout, weight meaning and resource limits. Each ordered query
then commits its exact ID/context bytes, ordinal, raw IEEE754 energy/weight bits
and the full canonical provider JSON. All source reason strings, nulls, positive
unmapped values, inactive contributions and array order are included. The wire
format is documented in README.txt. Chunk boundaries leave it unchanged; input
order and contents are bound. Hashes give integrity, not external authentication.

Only fixed summaries and a commitment remain in memory. Exact per-query retrieval
requires the retained immutable operands/provider/generator or a separately
qualified dynamic replay source. This fixture retains its complete original CSV
and input manifests. A digest cannot replace missing accepted native history;
no rolling owner/checkpoint implementation is claimed here.

Evidence
129 C++ checks compare all7 frozen pointwise fixtures, a16-electron legacy-admitted
candidate batch, repeated-energy queries and malformed/omitted/zero-substituted
observations, resource caps, numeric under/overflow and changed rounding mode.
Source coverage and candidate scalar bits remain unchanged; scheduler RNG parity
passes. Legacy unique-energy metadata is equivalent pointwise; the new stream
additionally counts all query multiplicities/weights rather than deduplicating
repeated energies. A separate Python Fraction/hashlib oracle passes6,931 checks
of exact category sums, per-term moments, extrema, support partitions and SHA
framing. No implementation-mirroring floating tolerance is substituted for exact
weight checks. The captured exact total query weight is
210953412474120907543/144115188075855872 represented electrons.

The captured pre-kick list is evaluated against a declared stationary H2O target:
46,726 queries,1,869,040 evidence observations,40 retained keys/280 category slots.
It hashes501,032,796 query bytes, maximum11,194 bytes for one encoded query, while
retaining no per-energy container. Processing peak RSS is9,160KiB at256/1024/4096/
16,384/46,726 checkpoints; final receipt serialization reaches10,520KiB process
peak. The first full pass used8.81CPU-s/8.89wall-s including all frozen provider
value/coverage evaluations. A second process produced a BYTE-IDENTICAL full
receipt. Exact timestamps/resources are in logs. This proves bounded memory on
this list, not production throughput or a new kinetic/EEDF result.

Integration scope
A new wrapper fork may replace only the evidence map with a private stream hook,
once per declared query/pair while preserving all numeric callbacks and caps on
physical work. Bind actual phase/census and weight semantics; discard the private
evidence if candidate preparation fails. Freeze the finished receipt before a
common owner preflight. Reuse already computed provider results where possible;
this prototype's full diagnostic evaluation cost is not approved for every
production substep. Preserve raw/replayable query input history externally.
The new component does not fix missing physics, native channel-set/product-family
admission, source counts, field/current/heat ownership or the native history
ceiling. All those review gates remain separate.

Resource limits are finite: <=64 keys,7 bins/key, <=32KiB validated text and256KiB
encoded query; default1M queries and16GiB encoded query bytes per receipt, hard
query maximum100M. Shapes/caps are checked before copying admitted policy/identity
strings; each update stages a bounded successor and commits by no-throw swaps.
The caller's already allocated input JSON/provider memory is a separate resource.
Overflow or limit refusal leaves counts, sums, digest and extrema unchanged.

All execution used the exclusiveCPU7 lease, nice15, frozen userspace, one thread,
1.5GiB address and180s CPU limits, existinglibcrypto, no new dependencies or GPU.
Input/source manifests and the source archive preserve the exact implementation.
No accepted scientific state, production process or full integration was advanced.
''')
summary={'status':'streaming_primitive_passed_fork_integration_pending','cpp_assertions':small['assertions'],'independent_assertions':oracle['assertions'],'queries':receipt['queries'],'observations':receipt['observations'],'retained_keys':40,'category_slots':280,'processing_peak_RSS_KiB':max(x['peak_RSS_KiB'] for x in memory),'whole_receipt_replay_byte_exact':True,'exact_input_weights_and_partitions':True,'content_sha256':receipt['content_sha256'],'generic_v6_modified':False,'current_generic_map_cap_removed_in_place':False,'physical_kernels_or_policy_changed':False,'native_or_field_or_joint_evolution':False,'complete_mechanism':False}
(O/'validation.json').write_text(json.dumps(summary,indent=2)+'\n')
h={'status':summary['status'],'runtime':str(P),'report':str(O/'report.txt'),'source_header':str(P/'source/streaming_coverage.hpp'),'source_header_sha256':sha(P/'source/streaming_coverage.hpp'),'source_archive':str(archive),'source_archive_sha256':sha(archive),'capture_receipt':str(O/'capture-receipt.json'),'capture_receipt_sha256':sha(O/'capture-receipt.json'),'integration_proposal':str(O/'integration-proposal.json'),'validation':str(O/'validation.json'),'full_integration_authorized':False,'cpu7_lease_released_on_completion':True}
(O/'handoff.json').write_text(json.dumps(h,indent=2)+'\n')
(O/'output-manifest.json').write_text(json.dumps({str(f.relative_to(O)):{'sha256':sha(f),'bytes':f.stat().st_size} for f in sorted(O.rglob('*')) if f.is_file() and f.name!='output-manifest.json'},indent=2)+'\n')
print(json.dumps(h))
