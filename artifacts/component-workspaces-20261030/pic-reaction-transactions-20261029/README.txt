Conservative internal reaction transaction reference

This isolated candidate is manufactured engineering evidence. It contains no
physical H2O rates, thresholds, branching probabilities or collision sampler.
It does not advance GapPic, Raikou electron F0, fields, material or production.
The shared repositories and sealed particle/checkpoint ABI are unchanged.

Source
  source/include/reaction_transactions.hpp: immutable catalog, packet/parent
    types, forced event, prepared transaction and serialized owner publication.
  source/src/reaction_transactions.cpp: validation, conservative staging,
    replay protection and standalone checksummed/atomic restart.
  source/include/integration_ports.hpp: typed, prefix-bound proposals for a
    common accepted interval: charged removal/insertion terms, one gas-parent
    debit reservation, neutral products, current moment and energy transfers.
  source/include/native_particle_view.hpp: value conversion to the actual
    frozen GapParticle type. It does not call authored add_particle or commit
    into native fields/particles. Native IDs must be reserved by the host.
  source/tests/test_reaction_transactions.cpp: manufactured forced-event and
    restart/rollback tests, including the actual frozen GapParticle value type.
  source/vendor: byte-identical existing Raikou nlohmann JSON header/license.
    OpenSSL SHA256 is the already installed frozen libcrypto, not new software.

Ownership and numerical contract
  The catalog has immutable species/state/channel IDs and a SHA256 fingerprint.
  Accepted event IDs bind copied payloads, a sequential journal and hash chain.
  Previously consumed particle IDs cannot be reused. Catalog/reference-energy
  changes invalidate restart. PIC owns free electron translation/number; an
  alternate Raikou F0 owner is rejected. There is no second electron integrator.

  Each forced reaction acts on one co-located gas-parent lot and an entire
  primary macro with matched represented event weight. Ionization and neutral
  dissociation replace the primary and insert products; attachment consumes it.
  Every product has the same physical event weight, immutable state and lineage.
  One gas-parent debit is staged. This is not a second adsorbed-surface debit;
  surface release and its energy/support exchange belong to their own owner.

  Local ledgers cover H/O nuclei, signed charge, vector momentum, electron/ion/
  neutral relativistic translation, explicit internal energy and species
  reference energy. Closed nonradiative internal events have zero external work
  and zero photon port. Authored external source counters are unchanged.
  The fixed validation budget is256 binary64 eps times the local sum of absolute
  terms, with a1e-300 numerical floor. Numeric residuals are recorded, never heat.
  The manufactured fixtures conserve nuclei/charge exactly and energy/momentum
  within that budget. No global field/driver/material energy claim follows.

  All potentially throwing work is private. prepare() returns an immutable
  proposal without publication; can_publish()/try_publish() reject a stale or
  foreign input snapshot and perform the local no-throw publication. apply() is
  a convenience for this standalone owner. Owner objects are noncopyable and
  single-threaded; this is an all-or-nothing state contract, not a lock-free
  multithreaded application commit.

  Standalone restart stores genesis, immutable catalog binding, full state and
  accepted event journal in one file. Restore replays every event and compares
  all receipts and final state, preserving duplicate-event protection. Save uses
  a same-directory temporary, file fsync, rename and directory fsync. An I/O
  error after rename can have an uncertain durability acknowledgement; inspect
  the valid checkpoint before retrying. SHA256 detects corruption/provenance
  mismatches; it is not an authentication signature for hostile writers.

Manufactured-only content
  Three stoichiometries are tested: e+H2O ->2e+H2O+, e+H2O ->e+OH+H, and
  e+H2O ->H-+OH. Test states deliberately carry arbitrary reference/internal
  energies, including a3eV ion reference and0.5eV internal energy. These are
  intentionally not physical ionization thresholds. The test's transverse
  equal/opposite recoil construction is an arbitrary energy partition, not a
  measured differential cross section or product distribution.

  Event clocks, positions, counts and source inventories are declared fixtures.
  No actual population or neutral scenario is sampled. Whole-primary reaction
  weights are required; partial macro splitting is rejected. Gas products are
  retained as neutral packets but arbitrary heavy-reaction targeting is absent.
  Radiation and wall/support exchanges require separate future transaction kinds.

  The reference stages a full bounded snapshot, at most4096 records/events. It
  is a correctness implementation, not a qualified high-throughput MCC operator.
  Production batching, memory reservation and weighted splitting remain work.

Build/checks
  Run build_and_test.py under taskset CPU7, nice15, one thread/job and
  flock <workspace>/outputs/pic-profile/testing-core7.lock.
  It uses the frozen userspace, -std=c++17 -O2 -Wall -Wextra -Wpedantic -Werror,
  pinned MFEM headers for a value-type check, existing JSON and frozen libcrypto.
  No MFEM/PIC solver library is required at runtime. run_sanitizers.py repeats
  with address/undefined behavior/leak checks; verify_ledgers.py independently
  reconstructs the three closed event ledgers with80-digit decimal arithmetic.
  Full command arrays, inputs, outputs and results are in
  <workspace>/outputs/pic-profile/reaction-transactions-20261029.

The exact native integration boundary is native-commit-contract.txt in that
report directory. Typed ports are staged proposals. They do not implement native
charge/current deposition, joint EM/thermal commit, circuit/field work, collision
selection, measured water chemistry or physical convergence.
