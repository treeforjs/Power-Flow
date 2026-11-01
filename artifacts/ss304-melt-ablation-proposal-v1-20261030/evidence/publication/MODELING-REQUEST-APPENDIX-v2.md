This versioned publication successor corrects the surface-model wording only; its earlier audit and appendix remain sealed. This separate appendix records the requested **melt-triggered, energy-limited
ablation approximation** and a finite surface-layer inventory. It does not alter
the [completed readiness audit](/home/treefor/Documents/Codex/2026-10-28/we/outputs/ss304-material-geometry-readiness-audit-20261029/READINESS-AUDIT.md),
its findings, any sealed material/source code, or the existing guards. No model
parameters are selected here and no integration tests are run. The new request
does not resolve the current correctness or integration defects.

The proposed metal rule is a reduced ablation closure, not a statement that alloy
melting and thermodynamic vaporization are the same transition. Melt is an
eligibility condition; finite available energy and finite removable inventory
must determine the amount released. The existing water/contaminant desorption
processes remain active below metal melt. Adsorbed water is a surface inventory;
the emitted water becomes gas after release.

**A consistent local removal budget**

The trigger needs an explicit definition: a temperature, liquid fraction or
enthalpy threshold, and its material/reference identity. The sealed metal law
uses a declared1672–1727K phase band and includes276320J/kg fusion energy once.
That band is available model context, not a new choice of trigger for this foil.
The alloy solidus, liquidus and an arbitrary single melt temperature are different
choices. Neither the trigger nor the requested removal enthalpy is assigned here.

Use the accepted local mass and energy after other admitted heat/work transfers.
A lifetime counter of deposited Joule energy is insufficient: energy may already
have left through conduction, desorption, emission, mechanical exchange or escape.
The energy counted as available must still be present in the same local accepted
stock. Do not use global foil headroom to remove a locally unsupported cell.

Let `E_tilde=E_accepted+Q_net` include the signed, once-only non-ablation thermal
transfers for the interval. A metal-removal candidate must obey

```
m_after = m_before - dm
E_after = E_tilde - E_actual_emitted_packets - X_other_out
E_actual_emitted_packets = sum_p weight_p*(R_p + I_p + K_p)
```

Here R and I are the emitted state's reference/internal energies, K is actual
native per-particle kinetic energy, and X_other_out is the separately identified
signed energy credited to other owners during removal. It must not duplicate a
term already inside Q_net or the packets. An already weighted packet/LayerExchange
operand is not multiplied by weight again. The remaining temperature derives from
`E_after/m_after`; no temperature or energy clamp completes a failed transaction.

For illustration only, an **isothermal-front variant**, if separately chosen,
could retain the remaining metal at T_star with condensed specific energy u_star.
Its excess energy would be `E_excess=E_tilde-m_before*u_star`, and conservation
would require

```
E_excess = E_actual_emitted_packets - dm*u_star + X_other_out.
```

If X_other_out is zero and the declared emitted state costs
`u_out-u_star=Lambda_rem>0` per kilogram, the ideal constant-cost relation is
`dm=E_excess/Lambda_rem`, capped by the actual exposed inventory. This equation
defines the required meaning of Lambda; it is not an assigned ablation law or
parameter. Actual samples must satisfy the energy equation with their actual K.
An implementation cannot substitute mean2kT, resample until energy fits, inflate
weights, or change a registry reference per event to hide a deficit. Packet
allocation and retry must retain the existing sampling/RNG ownership policy.

The removal enthalpy specification must state whether it means a total energy
from the300K reference, an incremental cost from the melt state, or a gas-minus-
condensed enthalpy difference. The thermal law has already paid sensible/fusion
energy on reaching melt. Charging that threshold again and also charging the
full packet reference would double count. A thermodynamic enthalpy input also
needs an explicit internal-energy/pressure-work convention: native packets carry
reference/internal energy plus actual kinetic energy, not an automatically added
h, p/rho or latent term. Mechanical, interface/fracture, escape and chemical
energy must have named receiving owners if included.

The existing conditional vapor model is supported only at1800–3000K; its latent
estimate is anchored at3000K by an explicit model assumption. A melt trigger in
the declared1672–1727K band is below that vapor admission. A future reduced branch
therefore needs its own scenario identity, product reference and applicability
statement. It must not silently relax the sealed vapor sampler, saturation law or
3000K guard. Kim's estimated vaporization energy and the model fusion energy are
different quantities; neither is automatically the requested removal cost.

Melting alone cannot vaporize a complete thermal cell. Complete removal requires
enough energy for all actual products and a supported empty-cell/geometry handoff.
Partial removal must remain possible, with positive remaining mass and local
recession limits. The energy-limited approximation supplies neither a validated
physical rate nor a bulk plasma/droplet distribution. If products are designated
neutral vapor, droplets, fragments or charged plasma, their state/phase identities
and energy/current obligations must be explicit.

**Finite layer inventories and their interfaces**

The requested oxide/intermediate-layer model is recorded as a requested surface
model. A useful inventory taxonomy is bulk alloy, oxide/interphase, other bound
contaminants, hydroxyl groups, removable adsorbed water and optional bulk-bound
hydrogen. This is not an assertion that all are uniform, separately ordered
films on this foil. Different patches can have different coverage and exposure.
The source overlayer illustration includes carbon outside its chemisorbed-water
region; it does not establish one universal layer order.

Each represented region needs a finite mass or elemental/species areal inventory,
chemical/reference-state identity and a mapping to physical exposed area. Keep
adsorbed molecular water, structural OH, oxide oxygen, carbon/nitrogen and bulk
hydrogen separate so the same H/O/C/N atom is not counted in two pools. An oxide
reduction or recombination source cannot create water without the corresponding
elemental debit and a defined reaction-reference energy. An oxide inventory is
not an additional SS304 inventory merely because it contains a metal element.

Inventories can be separate while heat uses one canonical control-volume energy.
If separate layer temperatures are introduced, each layer also needs its heat
capacity, geometry and interface thermal conductance; the same interlayer heat
must be added and subtracted once. Do not create independent copies of deposited
energy for each layer. Instantaneous common temperature is itself a declared
reduction, not a consequence of the word “thin.” Metal k/sigma are not oxide,
hydroxyl, carbon or water properties.

Removing a covering layer changes the exposed patch material, adsorption sites,
normal, physical area and source reference. Do not reset a fresh monolayer or
delete retained OH/water when the underlying support disappears. Any entrainment,
delamination, permeability or new exposure rule needs an explicit transfer and
energy/geometry treatment. Water desorption need not wait for oxide or metal
melting; conversely, bulk metal melting does not by itself prove that a coating
has disappeared. A DFT simulation slab dimension is not a measured coating
thickness, and an XPS monolayer count is not a film thickness.

Electron emission must use the actually exposed material and the same trial
one-sided field. The clean-metal4.3eV component hypothesis does not become an
oxide-covered work function or dielectric-barrier law automatically. Oxide charge,
electrical conductivity/permittivity, interfacial current and field screening
remain inputs/requirements if the chosen reduction intends to model their effects.

**What the existing primary sources recover**

| Source and location | Recoverable information | Limit for the requested foil model |
|---|---|---|
| Lamppa et al., SAND2021-12691, p.11 Fig.2; pp.48–49 | Qualitative304L contaminant overlayer: oxide, hydroxylated, chemisorbed-water and adventitious-carbon regions; discussion of H2O, H2, hydrocarbons, hydroxyls, nitrides and oxides | Schematic surface classes, not measured order, thickness, stoichiometry or release energy for this foil |
| Same report, p.49 Table4, visually verified from the cached PDF | XPS coverage estimates before/after250°C bakeout, in the report's convention1ML≈10^15surface atoms/cm²; values reproduced below | Reference304L/UHV coupon context; not this foil, not ns-rate release kinetics, not separate volumetric layers |
| Same report, pp.48–49 Table3 | TPD release totals for a304L coupon: H2O6.8, CO2 2.3, CO50.7, H2 4.5, total64.3 reported ML | TPD cannot locate the original reservoir; release products can involve chemistry/bulk supply. Do not add these totals to the XPS pools |
| Lane et al.2018, accepted manuscript pp.3–5,22,25–26 | Ideal Fe2O3/Cr2O3 facet context; physical water monolayer convention about13molecules/nm²; molecular desorption/coverage model and scoped MD comparisons | Not actual passive-oxide composition/thickness. The0.87eV low-coverage ClayFF value and0.695eV effective kinetic fit have different roles; neither is a new oxide-removal or foil binding-energy default |
| Mills et al.2004 and existing admitted metal record | Condensed caloric/reference treatment, fusion estimate, liquid cp; declared alloy melting interval/context | Does not specify oxide/hydroxyl/contaminant thermodynamics or this foil's phase calibration |
| Kim1975 and sealed material-v2 evidence | Estimated304L liquid/vapor inputs and the existing conditional vapor reference alternative | Not a measured melt-trigger removal enthalpy, layer energy, accommodation or full EOS |
| Bennett et al.2023, pp.040401-2–3 and melt discussion | Contaminant plasma before bulk melt, H2-based release model and a distinct dense conducting melt/plasma treatment | Does not establish instant dilute metal vapor at melt, nor an H2O/oxide reaction-energy set |
| Smith dissertation, printed pp.17–18 | Initial monolayers remain an important unknown; substrate/batch dependence and the cited thermal-release context | No measured layer thickness or shot-specific water inventory follows from those equations; the≈1eV H2O attribution remains unresolved in the existing primary-source audit |

The visually checked Table4 values are:

| Reported adsorbed species/group | UHV, room temperature (ML) | UHV after250°C (ML) | Reported desorbed20–250°C (ML) |
|---|---:|---:|---:|
| CO |1|0.5|0.5|
| CO2 |2|1.7|0.3|
| H2O |0.3|0|0.3|
| OH |2.2|0|2.2|
| C, C-C or C-H grouping |8|6.4|1.6|

The report totals desorbed coverage as4.9ML and describes8.6ML carbon compounds
remaining. These are source observations for their stated samples and definitions;
**none is adopted as an initial condition here**. The report's approximate
surface-atom normalization differs from Lane's water-molecule/facet monolayer.
Coverage, molecular inventory, multilayer thickness and bulk composition cannot
be interchanged. Slow bakeout/TPD results also do not establish the ns response.

The existing sealed water scenario selects0.33 coverage, an ideal-oxide surrogate,
an Eyring prefactor and a particular adsorbed-reference closure. They remain
explicit component assumptions, not measurements of the new layer stack. The
new appendix neither replaces them with Table4 nor adds all literature pools
together. Cited quantitative-XPS and earlier conditioning reports are documented
follow-up sources; their full measurement/attenuation records and this foil's
preparation data have not been recovered by this appendix.

**Inputs that remain unknown or require an explicit selection**

- Metal trigger definition and value, removable-front mass/area, and the positive
  removal-energy value with its reference/enthalpy convention.
- Actual oxide/interphase composition, stoichiometry, thickness or areal mass,
  spatial coverage, roughness and chemical state; intermediate carbon/nitrogen/
  hydroxyl inventories and whether any are separate films or mixed regions.
- Initial removable-water amount and monolayer definition for each physical
  patch, retained structural OH, multilayer water if present, and any bulk H
  inventory/transport separate from the surface pools.
- Layer caloric/transport properties and thermal contact if separate temperatures
  or thermal barriers are resolved; electrical/charging properties and actual
  exposed work function if their field/source effects are represented.
- Layer desorption, reaction, dissolution, fracture or entrainment rules and
  their reference energies; kinetic barriers and thermodynamic costs must remain
  distinct unless their equality is explicitly justified.
- Product composition/state/charge, internal/reference energy, gas temperature
  and actual velocity/angle sampling, accommodation/return or redeposition, and
  whether the reduced product is gas, droplets or fragments.
- Recession/remap and exposure rules, current path when a conductive layer is
  removed, remaining mechanical applicability and the common accepted-boundary
  checkpoint schema for all inventories and geometry.

The machine-readable [input register](/home/treefor/Documents/Codex/2026-10-28/we/outputs/ss304-melt-ablation-proposal-v1-20261030/publication/model-input-register-v2.json)
keeps unresolved numerical values null; zero is not used to conceal missing data.

**Mass, current, geometry and common-transaction requirements**

Use the existing exact native count/positive-carrier representation for actual
births and actual source-pool debits. Retain finite inventories and the current
physical packet-weight policy. Neutral release changes no net charge; charged
products require charge-balanced particle/electrode currents and their actual
ionization/electron-reference energies. Do not create charge or an external
authored-energy source simply because a material cell changes label.

The removed physical volume follows the selected material density/state, not an
invented oxide thickness or a nominal monolayer conversion. Geometry, physical
areas/normals, electric/magnetic fields, conducting current paths, first-hit
collection and source identities must use the same trial. Field energy changed
by loss of conductor is reconciled with circuit/boundary work, not credited as
arbitrary wall heat. Metal/gas transport and moving-wall support/work remain
separate from the already included caloric fusion and source references.

Every proposal is prepared from one accepted state and the same owned draws.
Joint preflight covers layer/metal mass and reference energy, thermal/elastic/
mechanical state, actual native IDs/counts, fields/current, geometry and RNG/
event watermarks. Depletion, an unsupported product/layer or failed energy budget
rejects the complete candidate. Restart retains the remaining inventories,
exposure state, cumulative removal and actual receiving gas states under the
same references, rather than reinitializing a fresh coating or monolayer.

This refinement therefore supplies a modeling direction and a precise input/
ownership contract. It does not close MG-01–MG-08 or any later source/native/EM
correctness findings, and it does not authorize extrapolation or new integration
tests in the present audit step.

Primary references reused, with no new literature acquisition:

- [Lamppa et al., Assessment of Electrode Contamination Mitigation at0.5MA Scale,
  SAND2021-12691,2021](/home/treefor/Documents/Codex/2026-10-28/we/work/field-comparison/sand2021-12691.pdf),
  pp.11,48–49; Fig.2 and Tables3–4. The cached PDF pages11 and49 were rendered
  and visually inspected for this appendix.
- [Lane et al.,2018, doi:10.1088/1361-648X/aae4af](https://doi.org/10.1088/1361-648X/aae4af),
  [accepted manuscript](https://www.osti.gov/servlets/purl/1595030), cached text
  and the prior sealed page/parameter audit.
- [Mills et al.,2004, ISIJ International44,1661](https://www.jstage.jst.go.jp/article/isijinternational1989/44/10/44_10_1661/_pdf),
  [Kim1975](https://www.osti.gov/servlets/purl/4152287), and the sealed material
  source records. No new property values were fitted or selected.
- [Bennett et al.,2023, PRAB26,040401](https://doi.org/10.1103/PhysRevAccelBeams.26.040401),
  supplied/cached full text; [Smith dissertation](/home/treefor/Documents/Codex/2026-10-28/we/work/dissertation-layout.txt),
  printed pp.17–18 in the cached layout text.
