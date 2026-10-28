# Research handoff for MITL kinetic chemistry and energy coupling

Date: 18 October 2026

Status: source port map and proposed verification work. No Quark source changes, builds, runners, tests or simulation jobs are included. The Quark owner integrates contributions through handoffs. Quark runners remain excluded until a later setup decision.

## Pinned repositories and scope

| Repository | Reviewed commit | Local reference |
| --- | --- | --- |
| electrothermal | `ff21224f975bb382fb7a09f4da72409b50f43ebb` | `/workspace/electrothermal-review` |
| Quark | `2546ecf54daa4284a96509da6dcda8800e65788b` | `/workspace/quark-review` |

Source paths below refer to the pinned electrothermal commit unless explicitly labeled Quark. Line numbers identify starting locations in that source snapshot. No `AGENTS.md` appears in either pinned recursive tree. Quark handoff conventions were checked in `docs/handoff/2026-10-15/README.md`; contribution and ownership conventions are in `docs/craft/contributing.md`.

The useful inheritance is explicit population ownership, same-event species and energy transfers, conservative physical histories, signed flux/work accounting, and provisional updates with rollback. Existing nitrogen and oxygen data do not supply the complete water-contaminant model needed for a magnetically insulated transmission line. Public Power-Flow reference work remains separate from this private source map. No private source has been copied into that repository.

## Source components and proposed native responsibilities

### Reaction events and reservoir accounting

- `zapdos-multibolt/apps/pore_plasma/src/materials/ADMassActionEvent.C:9` defines one shared molar event rate. Rate orders are separate from signed species stoichiometry, and third-body participation must be explicit.
- `zapdos-multibolt/apps/pore_plasma/src/materials/ADReactionEnergyLedger.C:12` requires named reservoirs, signed energy allocations and a provenance record. Lines 41 and 53 check the balance and use the same event rate for all powers. Formation and excitation storage already represented by evolving populations must not become duplicate energy sources.
- `gas-kinetics/gas_kinetics/pathways.py:573` supplies source proposals, commit/replay guards and checkpoint ownership. Its frozen-rate increment is an accounting operation, not a general stiff kinetic integrator.
- `gas-kinetics/gas_kinetics/coupling/nitrogen_excitation_relaxation.py:108` provides an adaptive population-relaxation reference with changing hazards for one-to-one transitions.

Proposed native responsibility: one evaluated event produces the species extents and all supported electron, heavy, internal, chemical, excitation, photon and surface transfers. Missing product or energy branches remain explicit missing data. A balanced numerical remainder does not determine a physical heat partition.

### Physical population and energy histories

- `zapdos-multibolt/apps/pore_plasma/src/kernels/ADConservativeLogTimeDerivative.C:35` advances differences of physical concentrations when the nonlinear unknown is logarithmic. It supports implicit Euler and variable-step BDF2 with accepted history.
- `zapdos-multibolt/apps/pore_plasma/include/utils/ConservativeLogMath.h:15` preserves small concentration differences without adding a density floor.
- `zapdos-multibolt/apps/pore_plasma/src/kernels/ADMaterialEnergyTimeDerivative.C:33` differences the full stateful energy density, including changing composition. Its implemented scope is fixed geometry with implicit Euler or BDF2.
- `zapdos-multibolt/native/n2_driven/RawPopulationWords.h:6`, `TranslationalBathState.h:10` and `RotationalLevelState.h:117` retain small increments against large inventories through numerical state representations. Those numerical words do not introduce extra physical reservoirs.

Proposed native responsibility: the accepted state owns physical species inventories, distribution populations, reservoir energies and time history together. A host using logarithmic coordinates must retain a conservative physical-inventory residual. Restart and rejected-step recovery must preserve the complete history needed by the selected integration method.

### Electron population and internal-state coupling

- `zapdos-multibolt/native/ionization_population/ionization_population.hpp:13` declares an origin-inclusive electron energy-bin ionization operator. Supplied thresholds, event frequencies and selected daughter sharing determine population gain/loss vectors. It does not determine heavy-product energy or angular scattering.
- `zapdos-multibolt/native/n2_driven/OriginPopulationStep.h:91` provides a sparse backward-Euler population step with supplied cross sections, scalar electric field, ionization and optional forward excitation. It rejects negative proposals rather than clipping or renormalizing them. Its electron-birth restriction must remain part of the timestep contract.
- `zapdos-multibolt/native/n2_driven/CoupledKineticBathStep.h:110` and `CoupledRotationalStep.h:46` couple the electron distribution, translational bath and explicitly resolved rotational populations at the same proposed time level. Only signed elastic exchange is automatically placed in the translational bath. Other inelastic losses and stored internal energy remain distinct.

Proposed native responsibility: retain distribution normalization and its relation to electron density, explicit level populations, collision event extents and energy moments in one local kinetic proposal. The host decides acceptance; failure leaves the prior state intact. Use the existing distribution and bath calculations as bounded references while implementing the required magnetized and spatial operators separately.

### Heavy kinetic and internal energy ownership

- `gas-kinetics/gas_kinetics/coupling/common_flow_random_energy.py:83` separates common bulk flow from neutral and ion random translational energies and passive internal sensible energies. Its common-velocity assumption does not represent neutral-ion drift. It does not supply electron kinetics, transport or global acceptance.
- `gas-kinetics/gas_kinetics/caloric.py:30` contains explicit caloric conversion interfaces. Three translational degrees remain applicable when the host geometry is planar.
- `zapdos-multibolt/apps/pore_plasma/src/materials/ADIonTemperatureSources.C:112` evaluates signed ion-neutral energy relaxation, explicitly supplied physical-work allocations and same-event carried creation/removal energy. Lines 132 and 133 pair the ion and gas gains.
- `zapdos-multibolt/plasma_models/pore_ion_temperature.py` requires coordinated replacements of gas/ion storage, sensible transport, pressure and the preexisting ion-work contribution before a separate ion temperature can be activated.
- `zapdos-multibolt/apps/pore_plasma/src/materials/ADResolvedHeavyEnergyState.C:122` performs sensible-energy inversion with composition-dependent capacity. Lines 156 and 169 retain implicit derivatives and composition gradients. Its present sourced registry is a bounded N/O subset with a hard 300–1000 K domain, not a hot H/O/OH/H2O/H2 equation of state.

Proposed native responsibility: declare which modes belong to each energy variable and which are represented by separately evolved species. A new ion or molecular internal-energy variable replaces its previous ownership everywhere; it must not be added alongside an unchanged common-heavy inventory.

### Transport, electrical work and boundary fluxes

- `zapdos-multibolt/apps/pore_plasma/src/materials/ADSpeciesFlux.C:69` supplies scalar electrostatic drift/diffusion. Lines 141 and 145 preserve physical and numerical flux contributions. This is not a magnetized MITL transport law.
- `zapdos-multibolt/apps/pore_plasma/src/materials/ADFluxElectricalWork.C:35` calculates current and signed bulk electrical work from the actual continuity flux properties. Bulk work and terminal circuit power have different owners.
- `zapdos-multibolt/apps/pore_plasma/src/materials/ADSpeciesSensibleHeatFlux.C:8` shares species fluxes with sensible enthalpy transport. Its constant-capacity closure excludes nonthermal ion energy and independently carried chemical/excitation energy.
- `zapdos-multibolt/apps/pore_plasma/src/materials/ADSpeciesStoredBoundaryPower.C:49` diagnoses signed formation/excitation energy carried by the normal flux actually applied in the species boundary condition. It does not apply a second heat source.
- `zapdos-multibolt/native/n2_driven/KineticGridAlignedTransport.h:12` and `KineticPotentialTransportReference.h:13` provide bounded one-dimensional fixed-potential transport references with reflecting boundaries and scalar face frequencies. They do not provide a magnetic-field or Lorentz-force interface.

Proposed native responsibility: expose actual face particle, momentum and energy moments, including chemical/excitation transport and separate physical/numerical contributions where applicable. Electrode Joule heating, gap work, directed kinetic work and collisional heat exchange must remain distinguishable. Face constitutive traces are not substitutes for the applied boundary law.

### Desorption and retained surface inventories

- `zapdos-multibolt/apps/pore_plasma/src/materials/ADSurfaceReactionNetwork.C:108` checks elemental and finite-site consistency. Lines 163–189 use one current-state surface event for gas flux, retained inventories, chemical release and wall energy.
- `zapdos-multibolt/apps/pore_plasma/src/bcs/ADWallPairedFluxBC.C:44` assembles the gas flux and retained patch gains from the same face evaluation with their coupled derivatives.
- `zapdos-multibolt/docs/wall_coupling.md` records finite-site, competing-species, charged-collection and integrated wall verification and its constitutive limits.

The existing surface is a uniform scalar patch. Its carried `Cv*T` sensible-energy prescription is explicitly a synthetic reference, not a validated effusive kinetic distribution or sheath accommodation law. It should not be used unchanged for foil-width heterogeneity. The MITL interface needs independently evolving per-face inventory, local foil surface temperature and an emitted distribution whose particle, momentum and energy moments agree with the surface/binding-energy debit. No SS304 water coverage, desorption prefactor, barrier distribution or accommodation coefficient is supplied by this source map.

## Native port sequence

1. Agree the owning Quark multiphysics paths and handoff interfaces before editing. The current contribution rules assign `physics/`, `solver/`, `discretization/`, `operator/`, `model/`, `runtime/`, `charon_native/` and `gpu/` to the multiphysics owner. Concrete target classes and insertion points remain an owner integration decision.
2. Define a source-independent native species/state registry, units, energy reference and reservoir ownership. Preserve explicit identities for H angular-momentum states and H2 rovibrational levels where the spectral model requires them; do not substitute total principal-shell or vibrational population for every individual radiating state.
3. Implement local event and population proposals with consistent source/energy moments and host-controlled acceptance. Use bounded synthetic networks and independent analytic references first. Preserve source identity, unsupported branches and rollback semantics.
4. Connect separate electron, heavy translational, ion translational and resolved internal modes only after their prior storage/transport/work ownership is replaced consistently. Retain the Python source and caloric interfaces as comparison references where they cover the selected equations.
5. Add per-face desorption and surface inventories driven by the externally supplied electrode thermal history. Preserve finite inventory and particle/energy moments. A uniform-surface reference must precede the foil-width-dependent case.
6. Connect actual gap transport and electromagnetic work. The new charged transport must implement the selected magnetization approximation and field coupling; scalar zero-magnetic-field references remain separate verification cases.
7. Add the sourced contaminant mechanism, state-resolved excitation/decay and radiation/instrument projection. Validate local chemistry and spectra before using spatially projected spectra to infer transverse desorption.
8. Integrate into Quark only through a reviewed owner handoff. No Quark runner or production-queue execution is part of this sequence at present.

## Required native interfaces

- **State:** stable species/level identities, physical concentrations or number densities, electron distribution and density, declared translational/internal reservoirs, chemical/excitation references, surface inventories, fields and accepted history. Convert molar versus particle SI at one declared interface.
- **Event:** current-state event extent, signed stoichiometry, electron redistribution and energy moment, product momentum/energy branches, photon channels and explicit reservoir transfers. Rates and all transfers use the same state and event.
- **Proposal:** candidate states and moments, domain/positivity failures and separate conservation defects. An unsuccessful proposal cannot modify accepted state or cumulative exchange histories.
- **Face exchange:** applied species flux and consistent momentum/energy moments, signed boundary storage transfer, wall inventory/energy counterpart and local surface state. Declare geometric area/depth units.
- **Commit and recovery:** one host acceptance decision covers coupled populations, reservoir energies, event history and field/source state. Checkpoints include replay identifiers and integration history; rejection restores caches and all provisional changes.
- **Diagnostics:** retain species/nuclei/charge inventory, every physical reservoir, boundary escape, imposed work, radiative transfer and numerical transport separately. Diagnostic cancellation alone cannot establish that the corresponding PDE source was applied.

## Existing verification evidence to preserve

These are archived numerical component results described by the pinned source. They were not rerun for this handoff and do not certify the pinned whole repository or a Quark implementation.

| Reference | Relevant verification to reproduce in the native host |
| --- | --- |
| `zapdos-multibolt/docs/material_energy_verification.md` | Changing composition with conserved physical energy; variable-step BDF2; source-driven energy; composition derivatives and sensible transport |
| `zapdos-multibolt/docs/component_qualification_status.md` | Conditional heavy channels and an evolved electron/excited/heavy/photon cascade; missing applied heat/photon counterparts must produce failures |
| `zapdos-multibolt/docs/species_stored_boundary.md` | Inward/outward species transport with signed formation/excitation boundary power and omitted-species negative control |
| `zapdos-multibolt/docs/open_thermal_budget_verification.md` | Open conduction budget with actual storage evolution; omission or sign reversal of outward heat must fail |
| `zapdos-multibolt/docs/ion_temperature_integration_20260927.md` | Changing ion capacity, paired relaxation, carried birth/removal energy, coordinated storage/flux/work replacement and omitted gas counterpart |
| `zapdos-multibolt/docs/wall_coupling.md` | Finite coverage, competing events, nuclei and energy balances, current-state gas/patch coupling and charge accounting |
| `zapdos-multibolt/native/n2_driven/accepted_collision_events_probe.cpp:19` | Same-state electron number and threshold-energy event increments; this probe does not advance or commit a trajectory |
| `zapdos-multibolt/native/n2_driven/accepted_local_energy_probe.cpp:72` | Actual local-generator moments with elastic, forward, rotational and excluded reverse contributions distinguished |
| `zapdos-multibolt/native/n2_driven/kinetic_bath_probe.cpp:41` and `coupled_rotation_probe.cpp:85` | Coupled bath/level states, signed exchange, work, births and residuals |
| `gas-kinetics/tests/test_common_flow_random_energy.py:191` | Force-only bulk work, frame consistency, exchange signs, trace inventories and small-increment loss detection |

Additional MITL verification is required: zero-field ballistic neutral release, Maxwellian boundary moments, finite inventory exhaustion, zero-magnetic-field and magnetized charged transport limits, level-resolved detailed balance where valid, quenching/radiative competition, resolved energy reference consistency, rejected-step/restart identity, and the actual Quark spatial/time/distributed discretization. Analytic conservation and source-removal controls should precede comparison against another solver using the same tables. Numerical verification does not establish the physical validity of those tables.

## Remaining physical and source limits

- The reviewed generic source machinery does not supply a complete H2O/OH/H/O/H2 contaminant mechanism, fragment energy distributions or complete recombination/superelastic partitions.
- Nitrogen finite-J metadata contains `N2(X,v=0,J=0..32)`, effective Biagi products and electrons. The bound source provides a limited reversible rotational set. This is not a full H2 rovibrational population model and does not qualify state-dependent excitation/ionization, scattering or product branching.
- `zapdos-multibolt/docs/remote_water_decay_chains_review.md` and `plasma_models/remote_water_decay_chains.py` contain band/decay structure work with explicit production limits. Band origins and nearby hydrogen wavelengths do not establish rotational line coincidences, oscillator strengths, opacity, continuum branching or a complete CR mechanism. The documentation includes acquisition-stage history; current data fields and qualification records take precedence over an isolated early access statement.
- No reviewed native transport interface establishes magnetic insulation, magnetized electron kinetics, sheath losses or MITL field closure.
- The bounded N/O thermal inversion cannot be extrapolated to the approximately 5000–9000 K spectral-inference range reported in the dissertation.
- `gas-kinetics/gas_kinetics/adapters/flash.py` is offline cgs conversion for an explicitly defined Cartesian state. It does not register native FLASH physics or define Quark integration.

## Source reuse and distribution boundary

The reviewed native source includes GPL-3.0-or-later notices. Quark `docs/adr/ADR-040-supply-chain-and-licensing.md`, section 2.2, excludes strong-copyleft code from compiled Quark targets under its current dependency policy. Before any direct reuse, the owner must establish a compatible license for the exact files and their dependencies, or select an independent implementation from appropriately licensed public specifications and scientific references. Private ownership of a repository alone does not determine the license of every included or derived file.

This handoff contains behavior and interface observations, proposed work and source references only. It carries no source patch, no copied implementation, no integration approval and no claim of native Quark execution. Keep this source map within the private handoff packet.
