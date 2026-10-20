# MITL dissertation publication handoff

Power-Flow contains the heating, desorption, boundary-source and experimental
evidence work for testing whether foil-width heating changes H₂O release and
fragment transport enough to distinguish VUV spectra from uniform desorption.
Use the [README](README.md) for the current FLASH workflow.

## Current architecture

FLASH is the intended production MHD/heating backend. The homegrown reduced MHD
solver, coupled run entrypoint, surrogate material defaults and their exclusive
configuration/test scaffolding have been removed. Shared geometry and two
compatibility helpers remain; there is no coupled prototype run API.

The [retirement verification record](supporting_documents/publication/mhd_retirement/completion.json)
records seven focused checks passed once. Three FLASH/source command-line import
checks also passed. Owned plasma/kinetic implementations and the qualified
upstream numerical modules retain their recorded hashes.

This repository supplies explicitly configured FLASH-output extraction, SI
surface histories, finite-inventory desorption, energy ledgers and conservative
boundary-source interfaces. Independent diffusion/conduction/Hall references
remain for verification. No licensed FLASH implementation or actual FLASH run
is supplied here.

The plasma CR model will be maintained in a separate repository. Its URL and
versioned receiver interface are pending. Boltzmann, CR, kinetic gap evolution
and state-resolved collision/radiative data remain with their owners. Quark
integration uses owner handoffs; this work has not used its runners or modified
its implementation.

## Completed work and evidence

- [Heating verification](docs/flash_heating_verification.md),
  [conductive heating](docs/conductive_heating_verification.md) and
  [Hall reference](docs/hall_induction_verification.md) define independent
  mathematical benchmarks and energy checks.
- [Conservative surface release](docs/surface_release_pipeline.md),
  [coverage-dependent release](docs/coverage_dependent_release.md) and
  [width-source comparison](docs/thermal_width_power_reference.md) preserve
  finite inventory, geometry and separate energy channels.
- [Coupled emission windows](docs/coupled_emission_windows.md) and
  [positive boundary quadrature](docs/boundary_quadrature.md) retain joint time,
  temperature, state and number mixtures, including covariance. Polynomial,
  finite-rule and actual binary-field discrepancies remain distinct.
- [Dissertation recovery](docs/dissertation_recovery.md),
  [current drive](docs/shot_17133_current_recovery.md),
  [depth/current recovery](docs/dissertation_current_depth_recovery.md) and
  [gap magnetic recovery](docs/dissertation_gap_magnetic_recovery.md) keep
  measured, reconstructed and simulated quantities separate.
- [Calibration rederivation](docs/absolute_spectral_calibration.md),
  [manufacturer follow-up](docs/calibration_primary_source_followup.md),
  [observed-spectrum recovery](docs/observed_spectrum_recovery.md) and
  [coverage inference](docs/surface_coverage_inference.md) record what can be
  compared and which absolute factors remain unidentified.

The [supporting-document catalog](supporting_documents/README.md) includes all
25 completed historical handoff packages, readable handoff reports, source/data
version records, verification commands, results and checksums. Twenty-four
archives are exact copies; the 0150 publication derivative removes generated
Python bytecode and records its original hash and every omission. Earlier
baseline requirements and failed-command/repair lineage remain visible.

Historical snapshots may contain the retired prototype and its old instructions.
They preserve provenance; current work uses the source and workflow at branch
HEAD. Third-party source papers remain cited with recorded hashes and access
details; the packages contain project reports and generated verification figures.

## What the results establish

Manufactured numerical examples qualify upstream source/energy representations
and width contrasts on stated scales. They are not physical H₂O parameters or
experimental uncertainty estimates. Published-figure recovery retains its
readout, normalization, timing and quantity limitations. No calibrated physical
full-water rovibrational spectrum or experimental discrimination result has been
demonstrated.

Absolute VUV normalization and physical monolayer inference still need the
original lamp certificate/reduction, shot acquisition/gain/gate/binning and
response evidence, together with qualified surface inventory/rate/energy inputs.
No single correction factor or physical coverage has been silently adopted.

## Next tasks

1. Agree the external plasma receiver's schema, state identities, units,
   time/precision and energy conventions. Obtain one hash-bound readback of
   the actual received source fields, including covariance and area applied once.
2. Follow the [Yacora and rotational-band comparison plan](docs/external_plasma_yacora_comparison.md).
   Use the identified 2023/2025/2026 references to fix model/data variants and
   overlapping states before freezing
   common-data numerical and native-model prediction cases. Use one observation
   model for both predictions.
3. Keep the spatial-localization null separate from a heating-model comparison:
   redistribute the same joint source for a fixed linear full-width equality
   check. Equal prescribed heating power does not make emitted mixtures equal.
4. Recover new experiment-linked geometry/reduction evidence. An unassessed
   official source lead is Lamppa et al., SAND-2021-12691,
   [DOI 10.2172/1825219](https://doi.org/10.2172/1825219),
   [official full text](https://www.osti.gov/servlets/purl/1825219).
   It is linked in the dissertation to the Mykonos MP3 platform; actual sample
   identity, geometry and diagnostic implications still need inspection.
5. Run the FLASH benchmark ladder only after the licensed source and configured
   workstation are available. Verify actual fields, units, material and Hall
   interfaces before inferring a physical width-dependent source.

The [preserved unfinished NPL acquisition](supporting_documents/in_progress/mitl-upstream-20261019T0953/README.md)
has request logs, source hashes, citation review and draft code. Its new-source
page review and packet validation are incomplete. No result or property from it
has been promoted to a completed source assessment.

## Continuing work

The [scientific checkpoint](supporting_documents/research_checkpoint.json) and
each package's verification record retain prior milestones and the saved queue.
Local operational status is maintained separately in
`/workspace/research_work/mitl-upstream-status.json`; check its lease and reserved
paths before an editing pass. Preserve unchanged passing evidence, record exact
commands/exits/hashes and use bounded local commits. Historical package replay
requires its declared complete baseline and source dependencies.
