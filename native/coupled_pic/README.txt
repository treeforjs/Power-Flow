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
