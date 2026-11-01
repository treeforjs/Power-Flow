# Water fragment reference and source-window audit

The existing H2/OH/O/O−/OH− registry uses the retained ATcT v1.202 **0 K** formation energies correctly. The same species masses, charges, atoms and references are identical in native-bound chemistry v3.1.1 and prepared v3.2.2. This audit changes neither registry nor production admission.

SimQueue 451 passed eight primary-source/invariant controls in frozen userspace on one allocated core. Two fresh processes produced identical reports. Exact dependency and source hashes are in INPUTS.json and the receipt. The independent invariant calculation uses explicit 70-digit arithmetic, not a threshold guessed from a source grid.

| Channel | Reference/kinematic finding | Consequence |
|---|---|---|
| H2 → H + H | Registry cost 4.47806926184 eV; target-rest electron threshold 4.47928801982 eV; MCCC source zero anchor 4.47713 eV | The first positive source knot, 4.5 eV, is above both. A prospective unshifted 4.5–150 eV model window avoids interpolating through the inconsistent onset, but fragment energy/angle sharing still needs an explicit model. |
| OH → H + O(3P) | Registry cost 4.41103316571 eV; the 2019 calculation assumes D0 = 3.88 eV | A 0 K versus 298 K change does not explain the difference: the latter gives 4.453630314 eV. Do not replace the registry reference with the source assumption. |
| e + OH → H + O− | Registry cost 2.94991988664 eV; finite-target threshold 2.95001506374 eV | Five positive source knots fall below this threshold. A prospective 3.01–9.91 eV window must explicitly refuse the excluded interval; no shift, clamp or invented zero. |
| OH → H + O(1D) | Product state is absent from the registry, but primary ATcT singlet O and independent NIST level data exist | O(1D) lies about 1.96736 eV above O(3P). Adding it needs an explicit species/state schema and branch assignment; the unlabelled workbook dissociation column does not supply that assignment. |

The OH paper's individual a4Sigma− and 1²Sigma− excitation channels have a clearer product interpretation: Section 4.2.1 identifies their repulsive curves with ground O(3P)+H(2S). An explicit vertical-transition breakup model can use the published calculated integral cross sections, fixed vertical loss, and exact recoil. The a4Sigma− proposal uses 7.36 eV loss and gives 2.94896683429 eV relative fragment kinetic energy; its first source knot above the finite-mass threshold is 7.3605 eV. This has now been implemented separately in `../oh-a4sigma-diagnostic-v1.1-20261101` as a diagnostic model, with no production adoption.

Source quality matters. Excitation!B546 is missing in the A²Sigma+ column. Excitation!D1350 contains the literal value **1974 bohr²** at 9.2993000000000006 eV in the 1²Sigma− column, inconsistent with the plotted scale and neighboring values. No value is repaired or imputed here. A future D-column importer must refuse the interval spanning its neighboring knots, 9.2925000000000004–9.3061000000000007 eV, unless corrected primary data are obtained. The a4Sigma− C column is complete and has maximum 0.47727999999999998 bohr², consistent with the published Figure 3 scale. These checks do not establish experimental accuracy.

Primary sources and provenance:

- [ATcT v1.202](https://doi.org/10.17038/CSE/2440256): retained complete HTML and selected original 0 K entries, both pinned by INPUTS.json.
- [NIST O levels](https://physics.nist.gov/PhysRefData/Handbook/Tables/oxygentable5.htm): O(1D) 15867.862 cm⁻¹, checked against the ATcT singlet-triplet difference.
- [MCCC H2 dissociative excitation](https://doi.org/10.1016/j.adt.2020.101361): calculated BO/adiabatic-nuclei data; raw source stays local under its redistribution restrictions.
- [Chakrabarti et al., PSST 28 085013 (2019)](https://doi.org/10.1088/1361-6595/ab364c), [author manuscript](https://discovery.ucl.ac.uk/id/eprint/10081520/), [publisher supplement](https://stacks.iop.org/PSST/28/085013/mmedia): exact source workbook SHA f3e96e924b3b6fde7583416f8ec483d708a52a982bfa80a07e61f0532c20b38e. The source calculations and DEA approximation are not measured event kernels. The paper's OH molecular electron affinity and the O atomic electron affinity are different quantities.

All windows in this artifact are proposals or kinematic compatibility checks. `production_qualified`, `physical_model_validated`, `network_complete` and `native_adopted` remain false. Missing encounter groups, unsupported tails and existing coverage refusals are unchanged. Raw tables/PDFs remain local; this artifact contains code, hashes and findings.
