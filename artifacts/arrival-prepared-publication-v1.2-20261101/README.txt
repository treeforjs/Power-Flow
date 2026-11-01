Arrival prepared publication v1.2
================================

The v1.1 namespace, sampling policy, configuration identity, ledger bytes and
checkpoint format are retained. Only the transaction API changes. The original
v1.1 runtime remains untouched.

Owner::prepare_commit(Staged&&, const NativeBirthReceipt&) verifies the actual
mapping and prepares all ledger/count/hash/state storage while accepted state
is unchanged. It consumes the staged input only after fallible work succeeds;
failure leaves that input reusable. The returned move-only PreparedCommit
contains an immutable staged view, ledger record and candidate state digest.
It is bound to this exact owner lifetime and accepted generation. It neither
publishes nor modifies the owner's accepted state.

can_publish and try_publish are noexcept and perform no hashing/allocation.
Publication consists of vector/string swaps, scalar assignments and token
invalidation. Multiple private proposals may coexist for the same accepted
state; publishing one invalidates all competing old-generation tokens. A
restored or equal-content foreign owner cannot publish another owner's token.
Discarding a prepared token leaves accepted state unchanged. Use serial owner
access; these APIs do not provide concurrency synchronization.

The compatibility commit method delegates to prepare and publish, returning
the same ledger record. A common coordinator must prepare arrival before its
first accepted owner swap and include can_publish in its joint preflight.

Qualification: frozen SimQueue CPU 447 and ASan/UBSan 449 each pass 295 checks,
including 73 C++ allocation-failure positions. The real original v1.1 header
is separately compiled: populated and empty histories, record hashes and
checkpoint bytes match; nonempty restored continuations agree. Prepared hashes
equal actual accepted hashes, and publication passes under immediate allocation
failure. Exact-count OpenSSL internal allocations are not independently fault
injected; the publication path performs no exact-count or crypto operations.
The first job446 failed because the zero-exposure fixture claimed a finite log
while supplying -infinity. The fixture was corrected; the admission guard was
preserved. The corrected runner records per-translation-unit dependencies at
their actual frozen-root paths. All source, dependency and binary hashes are
verified unchanged.

This is a component qualification. Actual all-owner publication, source state
projection, native/gas coupling and physical breakdown remain separate work.
The development coordinator successor is
pic-water-common-commit-v8.0.3-publication-20261101. It must not be promoted until
its full-window failure, restart and baseline comparisons pass.
