# Starting the full neutral collisional-radiative model

The initial implementation now exports a state-resolved **H2 radiative
subgraph** to the existing Raikou mechanism format and exercises its native
transient chemistry kernel. It is not yet a complete H/H2/H2O/OH CR mechanism
or a fit to the dissertation. The explicit
[requirements manifest](../configs/neutral_cr_requirements.json) records the
remaining state, collision, source, electron and observation inputs.

## Current runtime audit

The locally inspected Raikou source is clean at
`d58dca88c5e9f93043a7c48852f6f3c1b55b37c2`, newer than the reference in the
earlier handoff. The external repository was read and executed without source
changes. The installed `build-v1-install/bin/raikou-native.exe` advertises
`chemistry.advance`, per-cell batch chemistry and EEDF operations. Its exact
binary hash is retained separately from the current source revision; the probe
does not attest that this older installed binary was built from that revision.

The generic chemistry API supplies the useful production path: arbitrary
state-as-species mechanisms, mass-action rates and analytic concentration
Jacobians, external EEDF coefficients, temperature tables with bounded support,
adaptive backward Euler/Newton/step doubling, positivity rejection and
transactional failure. Electron and heavy temperatures/rates remain frozen
during each call. A transport/field host must update them and check coupled
energy convergence. Internal-state labels can be encoded in species IDs; their
physical mapping, energy and degeneracy data remain caller responsibilities.

The separate legacy `CollisionalRadiativeModel` is an approximate five-level
principal-quantum-number hydrogen QSS model. It aggregates states needed to
distinguish metastable H(2s) from emitting H(2p) and uses a molecular-band proxy.
It therefore cannot supply the requested full neutral model. Its legacy slab
escape approximation and bounded inverse fit are not adopted here.

Raikou's `h2o_research` scaffold has 19 species and 13 reactions. Native
validation accepts its elemental/charge structure, but reports eight missing
EEDF coefficients, five missing heavy-particle parameters and all 13 reaction
energy partitions unchecked. It has no resolved H2 B/C ladder, no atomic H
2s/2p populations and no complete OH/H2O electronic/rovibrational model.
Adding species labels alone does not fill those rates.

## Executed molecular radiative foundation

The existing curated MOLAT H2 files contain 28,466 B/C-to-X transitions. The
new [probe/export tool](../tools/neutral_cr_probe.py) preserves all branches,
including wavelengths outside the recorded band, and creates canonical
`H2:electronic:v#:J#` state IDs. C-plus and C-minus retain their source labels.
The exported graph has:

- 1,378 upper B/C levels and 295 lower X levels, or 1,673 states in total.
- 28,466 spontaneous radiative reactions with sourced Einstein A values.
- 16,609 transitions within 115-165 nm.
- Matched per-event internal-energy removal and photon-energy gain, derived
  from each retained wavelength. This does not establish a complete state
  energy table, excitation library or thermodynamic closure.

The complete graph passed the installed runtime's `mechanism.validate`.
Its `rates_ready` means that this restricted set of A values is supplied; it
does **not** mean the physical CR mechanism is complete. Absolute molecular
level energies, statistical/nuclear-spin conventions and mapping to a compatible
collisional basis remain required.

A small native control initializes an **authored** 1e12/m3 in B(v=0,J=0), with
all receiving levels empty, and evolves its 15 sourced radiative branches.
The total A is 1.865137926e9/s, giving a 0.536153 ns radiative lifetime. Over a
saved history through 60 ns:

- Maximum population error against independent exponential decay/branching is
  3.98e-5 of initial density.
- Molecular-number balance closes to 5.0e-15 relative.
- Integrated native photon-event energy agrees with the independent branch
  energy to 2.3e-16 relative.

The existing observation projector then produces 0-3 ns and 3-6 ns spectra
over 115-165 nm using an authored 1 mm path and 0.1 nm Gaussian FWHM. Missing
wavelengths are not renormalized into the recorded band. This verifies native
state readback, branching and gated projection. It is not an excitation model,
physical source admission, calibrated line width or synthetic shot prediction.
The projector integrates its declared piecewise-linear saved-population
interpolant; accuracy of that interpolant is separate from the native kinetic
and analytic checks above.

Outputs and the visually inspected control figure are in
`runs/neutral_cr_radiative_start_v3`. Reproduce with a fresh directory:

```powershell
python tools/neutral_cr_probe.py --raikou-root C:/Users/tsmit/Documents/Raikou --executable C:/Users/tsmit/Documents/Raikou/build-v1-install/bin/raikou-native.exe --output-dir runs/new_neutral_cr_start
```

The output includes the full mechanism, state catalog, source/binary hashes,
native validation, every control request/response, population arrays, gate
spectra, figure and verification report. No replacement CR solver is added.

## Physical mechanism to complete next

Atomic H must distinguish 1s, 2s and 2p, with an adequate nl-resolved cascade
ladder and ionization/recombination channels. Atomic O from water breakup also
needs explicit VUV emitting states. NIST ASD supplies evaluated atomic levels,
wavelengths and transition probabilities; it does not by itself supply a full
electron-collision model. [NIST ASD](https://www.nist.gov/pml/atomic-spectra-database)

H2 needs a source-resolved X(v,J) distribution feeding B/C levels through
electron excitation, followed by radiative branches, excited-state transfer,
quenching, predissociation, ionization and omitted-state cascades. A single
band population must not be copied into every line. The source's parity and
nuclear-spin conventions must match those of collision data. The molecular
CR study by Fujii et al. identifies density-dependent sensitivity to both
ground-state and excited-state excitation data, with sparse excited-state
rates remaining a limitation. [Fujii et al.](https://arxiv.org/abs/2405.10227)

The MCCC database is a candidate for compatible electron-molecule state data.
Its provider specifies citation and redistribution conditions; the tables
have not been imported or redistributed here.
[MCCC](https://mccc-db.org/)

OH and water require their own electronic and rotational/vibrational bases,
formation pumping, fragmentation/recoil and state-specific losses. OH A-X near
308 nm is a useful separate optical diagnostic, outside the present VUV band.
Water is a polyatomic rotor, so a diatomic `(v,J)` label is insufficient. The
existing [all-molecule state specification](rovibrational_data_requirements.md)
also applies to every charged or intermediate molecular product retained by
the accepted chemistry.

The retained Song2021/NIST total-ionization and attachment/fragment tables are
useful partial inputs. Total ionization and its partial fragment channels must
not be summed as independent losses. Total ionization is not a replacement
for resolved excitation into the emitting states.
[NIST electron-impact database](https://www.nist.gov/pml/electron-impact-cross-sections-ionization-and-excitation-database)

## Source, electron and optical connection

The current four-face maps supply a provisional neutral spatial/velocity
distribution, not excited-state populations. The `assumed_source_state` labels
must receive justified electronic/rovibrational distributions. Independent
assumed H/H2/H2O/OH reservoirs do not establish water-fragment branching or a
shared parent inventory. Native source acceptance must conserve the chosen
parent nuclei, number, internal energy and actual sampled kinetic energy.

The physical source basis now has a separate
[common-water contract](../configs/neutral_cr_water_inventory_source.json).
For each actual material-face patch,

```math
\Sigma_{\mathrm{H2O},0}=N_{\mathrm{ML}}\,\sigma_{\mathrm{ML}}\,f_{\mathrm{patch}},
\qquad N_{\mathrm{H2O},0,\mathrm{patch}}=\Sigma_{\mathrm{H2O},0}\,A_{\mathrm{patch}}.
```

`N_ML`, named **`initial_h2o_monolayers`** in the fitting interface, is the
primary free nonnegative parameter for spectral optimization/regression. It
is bound anew for every forward-model candidate, rather than fixed at an
assumed experimental value. `sigma_ML` is the separately defined number of H2O
molecules per square metre in one monolayer. No numerical best-fit monolayer
count has been selected or inferred, and there is no universal water-on-steel
areal density installed here. The
definition must state its area convention; microscopic roughness must not be
silently multiplied into a geometric material-face area. A common scalar
coverage may apply to all four faces while their actual different areas remain
in the inventory calculation.

Reuse the existing `MonolayerDefinition`, `ParentSurfaceInventory` and
`prepare_coverage_inventory` in `mhdlab.coverage_inventory` when those inputs
are supplied. One H2O pool, or normalized H2O internal-state fractions, owns
this parent inventory. Each desorption event removes one adsorbed H2O and
inserts one gas H2O; it does not also insert another copy of H, H2 and OH.
Water-derived nonparent fragments initially have zero independent source
inventory. Any genuinely separate contaminant requires its own declared amount
and provenance.

H, H2, OH, O and supporting ions then arise from accepted dissociation,
dissociative excitation and subsequent chemistry. For example,
`e + H2O -> e + OH + H` conserves two H nuclei, one O nucleus and charge.
That stoichiometry specifies no rate or fixed branch probability. Alternative
fragment channels require their own supported state-resolved rates and product
energy/recoil distributions. Closed-control nuclei budgets include the
remaining adsorbed water, every gas state/product and accounted boundary
exits. Native chemistry's atom/charge checks remain mandatory.

The historical four-pool GIF and the authored initial H2 radiative control
above remain separate controls; neither is reinterpreted as this water-derived
physical source. The common-water source cannot run physically until each
candidate is bound, its monolayer definition and physical release/state inputs
are supplied, and the existing heating qualifications are addressed. A globally
fixed monolayer value is not required: the optimizer supplies it at each
evaluation. The requirements manifest lists these numeric source gaps beside
the missing electron-density, EEDF, rate and clock inputs.

The new [parameter binding](../mhdlab/neutral_cr_water_source.py) accepts a
scalar candidate or a strictly increasing nonnegative scan and returns the
existing validated `CoverageInventoryPlan`:

```python
from mhdlab.neutral_cr_water_source import (
    bind_initial_h2o_monolayers, build_common_water_source_sweep,
)

plan = bind_initial_h2o_monolayers(
    surface, candidate_monolayers,
    monolayer=explicit_monolayer_definition,
    allocations=water_only_pool_or_state_allocations,
    state_species={water_source_state_id: "H2O"},
    provenance="Declared physical source inputs", status="unvalidated",
)
```

Each allocation's source state must explicitly map to H2O, and all fractions
must partition one parent inventory. Actual patch area enters once.
`build_common_water_source_sweep` repeats the qualified depletion solve at
every candidate with fixed supplied kinetic reference coverages, rather than
scaling an existing fragment history. The complete transport/CR/optical owner
then produces one spectrum per candidate for the existing finite-grid or a
qualified regression interface. No arbitrary physical upper bound is selected;
the supplied reservoir and rate-validity domains bound admissible candidates.
Nine focused binding tests pass, including zero coverage, real area accounting,
fragment rejection and nonlinear source recomputation.

Each cell and time needs an electron-density/EEDF history or an explicitly
qualified electron closure. MFEM wall temperature, neutral translational
temperature, rotational/vibrational populations and electron energy are
different inputs. A neutral temperature fit cannot silently prescribe Te, and
unknown electron excitation can be degenerate with neutral density.

Changing coverage is not universally a post-processing multiplier. A source
can scale linearly with `N_ML` under a prescribed wall temperature and a
coverage-independent first-order law. Light can also scale linearly only when
the transport and state operator remain linear, the electron density/EEDF is
prescribed, neutral-density-dependent losses are negligible, transfer is
optically thin and detector response remains linear. Coverage-dependent
binding, neutral-neutral chemistry/quenching, electron or wall-energy feedback,
opacity and saturation require new coupled cases. Even in an admitted linear
regime, an undefined `sigma_ML` leaves only the product `N_ML*sigma_ML`
identifiable; it cannot establish a unique physical monolayer count.
Unknown absolute optical/detector throughput supplies a further multiplicative
degeneracy with coverage in that linear regime. A regression must constrain
that scale independently, include justified calibration information as a
nuisance parameter, or report the unresolved coverage-throughput combination.
It must not return a unique absolute monolayer count solely by fixing an
unmeasured optical scale.

The state populations then enter the existing
[receiver/observation contract](raikou_receiver_contract.md). Simulation-derived
ray-cell paths and collection weights apply once. Detector etendue, throughput,
pixel allocation and exposure belong to the downstream radiometric response,
not the chemistry. Ground-state optical depth and ordered transfer must be
assessed before treating Ly-alpha or molecular bands as optically thin. Actual
emitting-state velocity and recoil distributions must support the line profile;
0.03 nm/pixel is sampling, not a measured line-spread function.

Finally, the dissertation's recovered derived transition-density displays
still need their explicit reduction/radiance mapping and registered gate before
an absolute spectral comparison. The radiative control above supplies none of
those missing experimental factors.
