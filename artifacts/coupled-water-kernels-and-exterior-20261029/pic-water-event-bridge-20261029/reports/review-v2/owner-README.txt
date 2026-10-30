Finite-rate conditional water parent-ionization bridge

ParentBridge connects the scheduler's local binary event candidate to the
frozen gas-phase H2O parent-ionization table and explicitly selected
parent_ground_equal_pair_CM_v1 product model. The initial table supports
13.5..1000 eV, with a physical zero below 12.621 eV and an unresolved
12.621..13.5 eV onset interval. Above 1000 eV is refused; total BEB ionization
is never substituted for the parent-ion yield. This channel alone is not a
complete collision mechanism or full water breakdown.

The conditional model retains one intact ground parent ion, puts an equal
back-to-back electron pair in the center-of-momentum frame and boosts all
products to the simulation frame. It is a conservation comparison model,
not measured differential branching, secondary energy or angular data.
The host supplies an isotropically sampled unit axis and owns its staged RNG.
Thermal and field electrode emission remain separate physical sources.

The scheduler uses the target-rest incident energy and lab Moller flux, one
actual cell volume in m3 (axial depth once), and represented particle weights.
Each event reacts min(projectile weight,target weight); its propensity uses
max(weight)*sigma*flux/volume. Neither the bridge nor the rate scheduler
debits inventory, inserts particles, deposits current, advances trajectories,
changes fields or publishes thermal state. SSA elapsed collision-subflow
time is distinct from the fixed physical stage of the operator split.

Products contain masses, signed charge numbers, proper velocities, represented
weights and state-reference energies. The neutral mass is bound explicitly
from the species catalog at bridge construction, then checked on each event.
Mion=Mneutral-me+12.621 eV/c2. Kinetic
energy plus the single reference increase is conserved, including target
translation. Do not additionally deposit the ionization threshold as heat.
The kernel digest binds the exact scheduler, loader, rest-frame kernel and
bridge sources through the canonical manifest in kernel-binding.json.
Inputs are copied and SHA256-checked before every frozen-runtime check.

Exact high/low weight remainders are retained. A downstream single-double
inventory conversion is refused when it would discard the low word. This
bridge does not implement an iterative compensated neutral reservoir or a
partial-primary arena. The native accepted-state owner must provide that
storage or refuse the entire event without mutating any owner. It must bind
actual positions, identities, revisions, states, cell geometry, source lineage,
accepted clocks and fresh product IDs, and preflight all owners before one
no-fail publication. A candidate hash alone does not perform this transaction.
In particular a coherently changed volume/rate or same-energy rotated input
requires authoritative host epoch/state preflight; pair consistency alone
does not certify that it is the currently accepted physical state.

The stable inverse boost uses u_rest=u_e-(gamma_e+gamma_relative)/(Gamma+1)*U
with long-double FMA. Numerical coverage requires a conservative normalized
direction roundoff bound <=1e-12, and the unchanged lab closure bound must
resolve target-rest reaction energy to1e-10. These explicit precision limits
refuse unresolved large common boosts; they are not water coefficients.
Unweighted long-double conservation is checked before scaling by macro weight.
Positive weighted receipt bounds must remain representable; underflow and
overflow are refused. Individual reported signed defects remain binary64
rounded diagnostics and can quantize to zero; acceptance never relies on that
quantization. Initial extreme-boost/tiny-weight failures are retained in review/.

Tests cover actual-table rates and conditional products at 14,20,100,500,999 eV,
stationary/thermal/moving targets and admitted common boosts up to0.001c;
500 product events;
20000 constant-population optical-depth samples; support gaps, unsupported
8 keV tail, stale candidates, mass and axis validation, finite parent debit,
and a sub-ULP debit that must remain compensated. The constant-population
clock ensemble is not trajectory evolution or a plasma-growth prediction.
No native field solve or electrode feedback is claimed by this component.

Execute bash run_checks.sh. CPU7, nice15, one thread, existing CPU7 lease,
frozen Arch runtime, 1GiB address cap and 60 CPU seconds. No GPU or production
actions. Number density is not generated; future displayed densities cm^-3.
