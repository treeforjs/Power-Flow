Actual native recombination can change local charge

SimQueue609 reproduces this in the frozen nativev8.7 implementation, using16
manufactured carriers (8 electrons and8 H2O+ ions, eachweight1e15) at different
positions in the same native collision cell. This deliberately large weight
makes the admission defect easy to measure; it is not a density prediction.

The actual first-stage candidate admits8 recombination events and8 charged
tracked-target debits. All charged carriers are consumed, leaving neutral
products in the actual pending reservoir. Global charge residual is exactly0.
The canonical Hall CIC deposition changes at8 nodes: maximum
3.1530836157120444e-4C, L1change1.281741307200009e-3C.
Hall same_net_charge returnsfalse with its unchanged gross-charge rounding bound.
The candidate is discarded; accepted native state and RNG are unchanged.

Mechanism: actual native BIRTH products use projectile position; tracked target
debits remove particles at their own positions. Global event charge conservation
does not imply local Maxwell charge continuity for charged-target consumption.
Neutral-target ionization at the projectile location is a different case.

This diagnostic explicitly enables the existing diagnostic target-energy mode
only to inspect candidate spatial charge. It is neither an admitted coupled step
nor a claim of energy or physical model closure. No production guard is changed.

The first Hall admission catches this before motion. Its preflight_final checks
motion endpoints and the transport-start census, not actual post-stage2 native
carriers. The new coordinator must therefore compare the authenticated final
native census with Hall proposal endpoint charge before any owner publication.
Reaction relocation current and its field-work/energy ownership require a
separate physical/numerical treatment; no artificial field repair is supplied.

The reproduction runner requires one SimQueue powerflow CPU with6GiB explicit
memory and the frozen userspace. Its compiler headers, source/archive identity
and result are recorded in receipt.json. Job605 was a missing Hall-motion object
at link time;609 links the actual pinned implementation and passes.
