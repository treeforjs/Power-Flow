Common coupled PIC integration support
=====================================

gas_step_scope.hpp owns rollback for one gas candidate during a coordinator
trial. Construction retains a shared gas owner and begins a new attempt.
Destruction/reset abandons only that attempt, even after gas mutation, failure
or publication preparation. It cannot cancel a newer retry, a foreign owner,
or a step that has already published. Move assignment first cleans up its old
attempt. The accepted state remains unchanged when rollback occurs.

Dependency: gas-target-owner v1.5 rollback API (AttemptBinding, begin_attempt,
matches_attempt, abandon_attempt). CandidateBinding remains the separate
mutation-revision/content preflight mechanism; PreparedCommit remains the
separate publication mechanism. Never substitute the rollback token for either.
Use one owner serially. This helper supplies no thread synchronization.

Intended use in a common trial:

  auto gas_trial = pf_coupled::GasStepScope::begin(owners.gas, step);
  // Prepare all source, field, thermal, native and gas candidate operations.
  // Seal and jointly preflight every owner before publishing any of them.
  // Publish through the existing no-fail serial commit window.

If preparation throws, normal stack unwinding rolls back this gas attempt.
After successful gas publication the guard is inert. The other owners still
need their own equivalent rollback and a common preflight/publication protocol.
This helper alone does not implement or qualify the all-owner integration.

The qualification driver includes tests/gas_step_scope.inc against the actual
gas owner. Controls exercise populated mutation, failure after preparation,
same-digest retry, foreign/dead owner, move construction/assignment, allocation
failure during begin, poisoned candidate cleanup, nonallocating rollback,
accepted-byte equality and nonempty checkpoint continuation. Device tests also
inject failures before/after the actual device begin copy; this is an explicit
fault-injection test and does not establish recovery from driver failure.

Reproducible runtime:
  /home/treefor/projects/runtimes/gas-target-owner-v1.5-rollback-20261101
Run its tools/qualify.py normal|sanitized|device LABEL using the pinned
arch-base-20261028/run-frozen launcher inside SimQueue with one assigned core,
explicit memory, single-thread library settings and a GPU claim for device.
The runtime snapshots this maintained header/test and records source hashes.
Refer to its receipts for passed modes; no speedup or full physical simulation
is claimed by these lifecycle controls.

native_gas_motion.hpp is a pure checked records converter for the source and
motion phases. It checks before/post-birth/after carrier partitions, native
mass/charge/weight, stored internal energy and actual sink incident velocity.
It produces gas flow/census inputs without granting publication authority.
The exact native prepared-motion snapshot and its live preflight are still
required at integration. Current v2 uses gas-v1.7.1 census-conservative mode;
CPU565 and ASan/UBSan569 each pass71 controls. Source and receipts are in
artifacts/native-gas-motion-census-v2-20261101.

Preserved v1 tests reproduce a gas-v1.6 carrier-share limitation: nonzero exact
correction can remain in an emptied cell. Its wall/* aggregate also combines
conductor and exterior correction losses. V2 checks that the corrected gas
owner leaves zero stock in emptied cells and keeps exact exit destinations.
The motion integration gate remains closed pending native capability and
common field/thermal accounting; conductor correction is not automatically
electrode heat. See the v2 README for the candidate sequence and limits.

The native_genesis.hpp builder admits the actual native/Hall/source/gas startup.
Normal and UBSan job606 each pass20 controls, including all8 foil source ports,
actual blob process operands, grid/mass identity, seeded RNG consistency and
populated census parity. See artifacts/joint-native-genesis-v9-20261101 for the
pinned development SourceStage snapshot and launcher. This startup qualification
does not open connected-evolution or checkpoint/continuation gates.

final_charge_guard.hpp checks actual post-stage2 native charge against Hall
motion endpoints with the unchanged per-node Hall bound. ASan/UBSan job616
rejects an actual8-event recombination counterexample that conserved total
charge but changed8 CIC nodes. Candidate authentication and this check must
both precede common publication; this helper alone grants no authority.
