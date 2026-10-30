Weighted finite-rate collision scheduler

This component supplies the event clock missing between supported Raikou
cross-section queries and conservative native reaction transactions. It does
not define a water mechanism or any missing differential reaction data.

For an explicit projectile packet of represented count wp and target packet
of count wt in physical volume V, a local binary event affects w=min(wp,wt).
Its propensity is max(wp,wt)*sigma(E_target)*g_Moller/V. The represented
reaction rate is therefore w*rate=wp*wt*sigma*g/V. The larger packet's
unreacted remainder must survive at identical state, position and velocity.
Both reacting pieces and all products carry w. This is a directly derived
weighted jump-process construction, not a claim to reproduce a named
published unequal-weight algorithm. It has exact stoichiometric multiplicity
when joined to a conservative transaction. Finite-weight fluctuations and
spatial/operator splitting still need convergence. No merging is performed.

Incident energy is evaluated in the target rest frame by Lorentz transforming
the projectile proper velocity. The lab collision flux is the Moller value
sqrt(|v1-v2|^2-|v1 cross v2|^2/c^2); it need not be a particle velocity below c.
The neutral density/packet count and cell volume are in the simulation frame.
All input mass, velocity, volume and weight quantities use SI and physical
represented counts. An extrusion depth enters V once. A plotting slab never
defines a second reacting population.

Channel functions must reject gaps and unsupported energies themselves. The
scheduler calls them even at zero relative flux, and rejects nonexclusive
product yields, aggregate inelastic observations, absent immutable data or
kernel IDs, invalid values, duplicated IDs and excess pair work. The clock
requires the actual state, mesh and mechanism hashes. Hash presence does not
qualify product data: the parent must bind the actual differential kernel
and all stoichiometric/energy conventions before preparing a physical event.

The direct exponential clock holds particle positions/velocities fixed during
one collision stage; after each accepted event the parent must rebuild rates
using the changed populations and remaining interval. Proposal creation copies
the RNG and is side-effect-free. Global publication must preflight the native
particles, neutral inventory, reaction, EM, heat and this RNG together, then
publish a single accepted state. Publishing the RNG alone is not an atomic
physical step. Rejected proposals are simply discarded. A no-event proposal
can be accepted through the requested endpoint; across intervals the process
is memoryless only for the frozen local-rate approximation. Step refinement
must bound field/time variation and the splitting error.

The scheduler checkpoint stores only its versioned deterministic RNG and
generation, not a complete physical restart. The parent must put it in the
same hashed accepted epoch as every other owner. Direct serialization here
uses the frozen C++ standard library mt19937_64 implementation. Input ordering
is canonicalized, and uniforms are strict open-interval 52-bit midpoints.

The native D-K-D field solver and compatible EM adapter must supply an admitted
collision insertion stage and compatible pre/post segment current deposition.
Do not simply append MCC to a velocity push. Vay et al., PRE111025306(2025),
https://doi.org/10.1103/PhysRevE.111.025306, is the primary integration reference;
its analysis is not automatic validation of ionization or newly born particles.
Morris et al., PoP29,123907(2022), https://doi.org/10.1063/5.0126336, identifies
the importance of preserving real-event counts for unequal macroparticle
weights. This implementation does not claim its differential atomic model
provides gas-phase water branching.

Initial verification uses manufactured constant cross sections: analytical
rate/weight identity, target-frame boost invariance, unsupported channels,
finite-volume locality, rollback/replay/restart, deterministic ordering, and
50,000 two-channel exponential clocks with specified six-sigma ensemble checks.
This establishes scheduling machinery only. It does not evolve native PIC,
create physical water ions, prove screening, or constitute the requested
continuous coupled foil simulation.

Run: bash run_checks.sh
Frozen userspace; CPU7 nice15 and shared CPU7 lock; one thread, 1GiB virtual
memory and 60CPU-second limits. No GPU or production action.
