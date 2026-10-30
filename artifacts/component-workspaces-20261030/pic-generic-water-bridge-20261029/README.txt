Generic conditional water bridge, isolated component — 2026-10-29

This implements eight sealed Raikou v6 conditional event models in ONE competing
SSA candidate set. It is a supported comparison subset, not a complete water
mechanism or a physical evolution driver. RequireCompleteMechanism refuses.
No production/native sealed source, source law, field/current path or GPU is
changed. Exact source/data bindings are in outputs/pic-generic-water-bridge-20261029.

source/generic_water_bridge.hpp:
  Bridge(data,dataset_sha,kernel_sha,Mwater,MH,PreserveReactantSites)
  candidates(electron PhysicalPackets,H2O PhysicalPackets,physical 3D cells,
             CompareSupportedSubset, pair/query/row/coverage caps)
  sample(candidate,physical e,physical H2O,private mt19937_64 OR replay variates)
  catalog(), specs(), coverage(E), OH_family_metadata()
Every batch has an explicit coverage ledger with source values for unapplied
observables, positive unmapped2a1 and null unknown rates. No tail or branch rate
renormalization occurs. At exactly1keV outer3 replaces the overlapping empirical
parent model; IAM elastic competes in the same set. Outside conditional support,
the contribution is excluded only from this declared subset and its missing
coverage remains visible. Unknown rates are not physical zeros.

Products have variable count2/3/4, canonical state/mass/internal/reference data,
proper velocity, represented event weight, and explicit primary/target placement.
A surviving scattered electron stays at the primary site. Other products stay
at the target. Hminus/OHminus capture requires EXACT full3D co-location, including
z; separated capture refuses without committing the private RNG. The bridge
never moves an attachment ion to evade this condition or invents a reaction
current. Hminus has ONE base event rate followed by ONE source-bound conditional
OH(v=0..4) draw, not five competing full-rate channels.

The common source-rounded mass basis is Mwater=18.01056468403u and H=1.00782503223u,
u=1.66053906660e-27kg, with the frozen conditional binding/internal energies.
These are explicitly supplied constructor operands and must match the registered
native catalog. They are not precision measured species masses. Equal-pair CM,
ground/internal endpoint, factorized DEA and IAM free-atom angular assumptions
remain the source model's conditional approximations. IAM is not molecular
angular validation. Low-energy angular channels, core/multiple/Auger processes,
positive2a1 completion, other inelastic channels and fragment continuation remain
unqualified. Non-H2O target continuation is refused, not silently ignored.

source/generic_water_native_adapter.hpp:
  native_catalog(bridge), native_event(outcome,...)
  fixed_native_comparison_channel(outcome,binding_file)
  mechanism_contract(bridge,binding_file)
The fixed native channel helper is ONLY for independently forced admission
comparisons with a product multiset known for that fixture. It must not be used
to create duplicated OHv rates or sequential full-duration channel loops. Actual
multi-channel native SSA completion still needs an immutable mechanism/channel
set and conditional product-family contract. Existing single-channel completion
checks have not been weakened. Native phase/state/field/geometry authority and
accepted RNG publication belong to the common owner; this sampler authenticates
neither a caller's external positions nor physical3D cell volume.

Validation exercises actual sealed kernels, all five OHv states, moving targets,
unequal weights, deterministic/restarted private RNG, actual-rate competing SSA
clocks, exact Fraction nuclei/charge/weight and140-digit independent energy/
momentum. Native v6.2 forced comparisons test variable products, full/partial
capture, exact restart/rollback and nonlocal refusal. Those tests have no physical
time advance and do not establish coupled chemistry, EEDF or a full mechanism.
The clock fixtures use V=1e-12m3, ne=2e6cm^-3 and nH2O=7e6cm^-3. Those are declared
manufactured populations, not inferred foil densities.

Bounds:100000 projectiles/targets,4KiB each ID,64MiB copied input identity payload;
underlying indexed builder retains separate pair/query/candidate and string caps.
Coverage evidence defaults256 distinct energies (hard maximum1024). Any cap is a
refusal, never silent pair/channel omission. The runner uses exclusiveCPU7,
nice15, frozen userspace, one thread,1.5GiB address limit and180s CPU limits.
All supplied source/data dependencies are byte-verified. No new dependencies.

Run: bash /home/treefor/projects/runtimes/pic-generic-water-bridge-20261029/run_checks.sh
Frozen native objects are reused only at their pinned hashes. A clean checkout
can rebuild the copied v6.2 source with PIC_GENERIC_REBUILD=1 before the command
above; run_checks.sh keeps the SAME CPU7/nice15/frozen lease for that build. Its derived object receipt is written
under build/, never over the frozen input reference. Source/debug object bytes
need not match across paths; archive, command, source identity and tests bind it.
Source and result manifests in the output directory define the sealed receipt.
