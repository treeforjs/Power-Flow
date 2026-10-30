Native finite-rate water parent-ionization subflow

Status: provisional integration. The 428-check normal suite passes; independent
arena review found restart/ownership/stoichiometry deficiencies in its sealed
dependency. Rebase onto the audited fix before qualification or use in a pilot.
No physical trajectory time has advanced and no full coupled run has launched.

This driver joins the actual parent-yield table, the conditional equal-pair CM
product model, indexed local weighted SSA, and the canonical native particle
arena. Rates are rebuilt after every event from positive exact-weight carriers.
The collision grid covers the actual native kinetic domain in three dimensions;
its cell volume includes the finite axial extent once. It is distinct from both
the field mesh and a visualization observation slab.

Each prepare stages all reaction subevents privately, then restores the accepted
particle/field state. The single resulting token publishes the native inventory,
scheduler RNG and product RNG together. Rejection after an earlier supported
event restores all of them. Random waiting time is collision-subflow time at a
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

Reproduction currently uses build_and_check.py under the CPU7 exclusive lease,
nice15 and one thread. It compiles and runs in the pinned Arch numerical runtime;
it does not use GPU or signal production. The kernel-binding.json includes the
new indexed scheduler bytes, leaving the older frozen bridge release unchanged.
See outputs/pic-water-native-subflow-20261029 for commands, failure attempt,
normal results and the exact native dependency archive identity.
