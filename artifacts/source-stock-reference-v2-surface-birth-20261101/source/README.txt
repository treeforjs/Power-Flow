Immutable source-stock chemical-reference transfer, version 1

Status: implementation candidate, awaiting component qualification; no deployed
simulation or coupled physical result is claimed.

The source remains the sole finite surface inventory owner. For every enabled
neutral patch this adapter takes the actual pinned native registry reference R,
computes exact N_before*R and N_after*R from the source's binary count carriers,
and verifies (N_after-N_before)*R + sum(actual birth weights*R) == 0 exactly.
The heat/source Trial loads, kinetic energy, desorption binding, work function and
velocities are unchanged. No chemical reference offset becomes electrode heat.
Zero, positive and negative reference origins are supported; water is negative.

Adapter configuration uses the real full_time_chemistry_registry/states from the
pinned v3.1.1 blob and actual native options. It validates every source patch,
kinetic mass and declared electron mass policy; all enabled wall, phase/metal,
writeoff and nonzero internal-energy source paths are refused explicitly. A
reference stock checkpoint is a derived immutable value of an actual source
checkpoint, never a mutable second inventory. Restore recomputes and compares it.

Adapter::prepare_group accepts actual source ports (owner, live Trial, immutable
field/final view, stream/patch registration, neutral/electron links). Every full
native credit stream and birth belongs to exactly one port; duplicate owners,
physical patch identities, streams, source packets and native births are refused.
Each neutral is matched to the actual sealed Trial packet's event ID, weight,
state, mass, time and phase space. Electron entries match genuine Candidate and
Handoff objects from that same Trial, registry electron mass/zero reference and
actual birth/handoff-start phase space. An already-advanced Handoff particle is
not used as the native birth state. Normal-coordinate face projection follows
the existing explicit field-to-native translation with a 1e-12 m mapping check.

Only after all validation does the adapter call the actual BirthScope::prepare
once, on the full owned immutable credit. The returned move-only
PreparedReferenceTransfer owns that genuine native PreparedBirths and credit,
plus exact per-port stock reference ledgers. Only const views are exposed.
Native original outer/base/scope lifetime guards remain authoritative. The
coordinator must still complete source/Hall/arrival/gas/native preflight and
common publication; this receipt independently authorizes no accepted mutation.
Discarding it performs no stock or thermal commit.

Initial electron support is an actual zero-microscopic-displacement thermal
surface handoff. Positive-field microscopic bridge, persistent-layer and
returned-source mappings are explicitly refused pending joint native/Hall
current/work ownership. This is an intermediate gate, not a reduction of the
fully coupled goal. Existing Handoff current is a preview of the remainder:
its current/work must not be deposited a second time if native motion owns that
same flight. Strong-field integration must settle this ownership explicitly.

Qualification uses actual source-neutral and thermal-electron packets and actual
native chemistry/birth scopes against v8.7, but fields/stock/arrival credit are
manufactured component inputs. Empty chemistry and birth-transport forks are
not coupled plasma evolution. Component source/thermal publication is used only
for exact stock depletion and checkpoint restoration. No accepted common native
step or physical reachability is claimed. Normal and UBSan boundary checks use
frozen numerical libraries, explicit memory and one SimQueue CPU; no GPU work or
performance claim. Source v7.2 must finish its gate before these tests start.
