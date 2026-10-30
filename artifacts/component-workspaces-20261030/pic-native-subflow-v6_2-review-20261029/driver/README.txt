Native finite-rate water parent-ionization subflow

The 442-check normal suite passes against the corrected v6.1 native arena.
Independent review status and exact source identities are in handoff.json.
No physical trajectory time has advanced and no full coupled run has launched.
Earlier review failures, fixes and intermediate receipts remain preserved.

This driver joins the actual parent-yield table, the conditional equal-pair CM
product model, indexed local weighted SSA, and the canonical native particle
arena. Rates are rebuilt after every event from positive exact-weight carriers.
The collision grid covers the actual native kinetic domain in three dimensions;
its cell volume includes the finite axial extent once. Foil boundaries must align
with the collision bins: fully solid cells are excluded, and partially cut cells
are refused until an exact vacuum-volume treatment is implemented. A particle
on a material surface needs a side-resolved transport stage before collision
admission. This grid is distinct from the field mesh and an observation slab.

Each prepare stages all reaction subevents privately, then restores the accepted
particle/field state. The single resulting token publishes the native inventory,
scheduler RNG and product RNG together. Every accepted exposure, including one
with zero events, records its common clock and before/after random-stream hashes
in the native history. Mixed-age restart combinations are rejected. Rejection
after an earlier supported event restores all owners. Waiting time is subflow time at a
frozen physical stage; it must never be used as particle flight time. The current
native placement restrictions remain initial-state proof or provisional equal-dt
mid-x splitting. A common variable-step EM/heat/source driver is not supplied.

The 100 eV fixture generates seven actual finite-rate parent events, with water
depletion and secondary production. A second fixture retains sub-ULP depletion
from a 1e20-water reservoir in its two-limb count. The one-target ensemble tests
the known first-event Bernoulli probability. Long artificial collision durations
hold fixture coordinates fixed solely to test this operator. They are not foil
simulation histories or physical times attained by the requested calculation.

Coverage remains the original conditional parent channel: 13.5..1000 eV, zero
below the physical threshold, and explicit refusal in the unresolved onset gap
and above support. Other mechanisms are not silently represented by this parent
channel. The 38.3 eV fixture has a valid first event whose products enter the
unsupported onset gap; the entire trial rejects rather than zero-filling it.
No electrode heat, surface emission, neutral source, Maxwell transport, or
macroscopic shielding is qualified by these tests.

The gas-only fixture uses H2O reference zero and ion reference 12.621 eV. Adsorbed
binding and electron/image reference reservoirs must be reconciled explicitly
when adding source and thermal transactions. All species use one particle arena.

The bulk native census validates all identities and exact counts once; it does
not repeat whole-population scans per particle. Resource caps reject a trial
without truncating its population, rates or history. This prototype admits at
most 100,000 native carriers, 4,096 native journal entries and 16 MiB of journal
receipts. Default subflow limits are 128 reactions and one million pair/query
visits. Production needs a qualified bounded history/checkpoint implementation;
the current full journal must not be grown through millions of steps. Two-limb
inventories also reject an exact debit that needs a third residual; usable source
weight ranges must be checked before a physical pilot.

Reproduction uses build_and_check.py under the CPU7 exclusive lease,
nice15 and one thread. It compiles and runs in the pinned Arch numerical runtime;
it does not use GPU or signal production. The kernel-binding.json includes the
new indexed scheduler bytes, leaving the older frozen bridge release unchanged.
Extract native-source.tar.gz into native/ first. The build checks the archive and
every extracted source file before compiling; no system packages are installed.
See outputs/pic-water-native-subflow-20261029 for commands, failure attempt,
normal results and the exact native dependency archive identity.
