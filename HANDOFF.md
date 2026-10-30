# MITL dissertation publication handoff

Power-Flow contains the heating, desorption, boundary-source and experimental
evidence work for testing whether foil-width heating changes H₂O release and
fragment transport enough to distinguish VUV spectra from uniform desorption.
Use the [README](README.md) for the current MFEM Joule workflow.

## Current architecture

MFEM `joule.cpp` is the selected electrode-heating backend, replacing FLASH.
It covers stationary electromagnetic diffusion and heat conduction; the stock
miniapp has one-way coupling and requires an experiment-specific driver for the
physical foil. See the [current plan](docs/mfem_h2o_publication_plan.md).
The current finite-width driver is a custom MFEM 2D A_z current-constrained
extension with a return conductor and true nonconducting exterior. It uses the
recovered shot 17133 current, remains capped at 60 ns for early VUV work, and
does not reuse the stock Joule kernel. Its constant-property solid temperatures
and independent assumed neutral reservoirs remain provisional; edge heating
comes from the solved current distribution. See the
[published ALEGRA comparison](docs/mfem_alegra_comparison.md).
The [current four-face GIF package](artifacts/neutral_desorption_perimeter_20261028/README.md)
contains 65 frames over 0.749788-60 ns, the actual ALEGRA comparison, current and
surface inputs, number ledgers and replay code. Its native heating source is
`runs/mfem_edge_current_source_four_faces`; its neutral maps are
`runs/neutral_desorption_current_17133_perimeter_v2`. All four material faces
emit with actual areas and outward normals. Both the overview and edge close-up
use equal x/y scales, replacing the earlier broad-face-only, y-enlarged view.
The actual first/middle/last GIF frames were
decoded and inspected. Field-derived edge concentration persists under mesh
and domain sensitivity checks, but corner temperatures are unresolved and
overshoot the published ALEGRA profiles. Absolute spectra fitting and physical
desorption/fragment yields remain outside this provisional demonstration.
The [perimeter-source record](docs/neutral_perimeter_release.md) retains the
original broad-face law histories and separately bounds negligible packet-tail
omissions. Sides supply 9.11% of released particles with the case's equal
per-area coverage/rate assumptions. Their H RMS x displacement is about150um
by60ns, without velocity scaling. Number and histogram residuals remain below
6e-15 of initial inventory; source/frame, packet-tail, flight and renderer
checks pass.
The newer [collection/CR record](docs/target_collection_and_cr.md) makes
`initial_h2o_monolayers` the free primary spectral-regression parameter through
the existing parent-water inventory. The historical GIF's four independent
pools are not the physical water source. Finite-pupil target collection,
geometric-edge and source-face overlap, exact supplied-population photon
projection, a sourced Raikou H2 radiative graph/control and bounded SS304 EOS
table admission are implemented. New geometry uses the stated approximately
19-inch target-to-first distance,8-inch mirror separation and50um by1inch slit;
the first-mirror-to-slit few-inch distance is straight-line geometry. Actual
relay prescription/throughput, full collision data and numeric alloy EOS
remain explicit inputs. No calibrated detector or monolayer fit is claimed.
The homegrown reduced MHD
solver, coupled run entrypoint, surrogate material defaults and their exclusive
configuration/test scaffolding have been removed. Shared geometry and two
compatibility helpers remain; there is no coupled prototype run API.

The [retirement verification record](supporting_documents/publication/mhd_retirement/completion.json)
records seven focused checks passed once. Three FLASH/source command-line import
checks also passed. Owned plasma/kinetic implementations and the qualified
upstream numerical modules retain their recorded hashes.

This repository supplies an MFEM surface-sample adapter, optional earlier FLASH extraction, SI
surface histories, finite-inventory desorption, energy ledgers and conservative
boundary-source interfaces. Independent diffusion/conduction/Hall references
remain for verification. The MFEM installation is available locally; FLASH
acquisition is no longer required. Earlier records retain their original scope.

The plasma CR model is maintained in [Raikou](https://github.com/treeforjs/raikou).
The [neutral-spectrum observation handoff](docs/raikou_observation_handoff.md)
records its source reference and population-output interface, with a manufactured
native output/projection control. Physical source admission is still pending.
Boltzmann, CR, kinetic gap evolution
and state-resolved collision/radiative data remain with their owners. Quark
integration uses owner handoffs; this work has not used its runners or modified
its implementation.

The current [neutral-spectroscopy implementation](docs/neutral_spectral_workflow.md)
adds hash-bound shot 15750/15765 readouts, categorical reconstruction scenarios,
a conservative receiver/population/geometry validator and finite neutral-ensemble
comparison. [Implementation evidence](supporting_documents/neutral_spectral_comparison/README.md)
records the tests and local outputs. All physical admission gaps remain explicit.
The newer [MFEM backend evidence](supporting_documents/mfem_joule/README.md)
records native MPI execution, saved-field readback and surface/release compatibility.

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

Absolute VUV normalization and physical monolayer inference still need a
defensible lamp/shot reduction and response reconstruction, registered acquisition
and qualified surface inventory/rate/energy inputs. Original detector files are
unavailable; bounded scenarios must preserve unresolved factors.
No single correction factor or physical coverage has been silently adopted.

## Next tasks

1. Follow the [MFEM qualification ladder](docs/mfem_joule_heating.md): native
   Joule readback, analytic/conservation checks, then a physical finite-width
   SS304/drive case. Use the [surface adapter](docs/mfem_surface_contract.md) to
   retain explicit absolute Kelvin, face support and source provenance.
2. Use the [receiver bundle contract](docs/raikou_receiver_contract.md) to obtain
   a physical hash-bound readback of the received source fields, joint moments,
   covariance, geometry and energy. Schema controls do not establish acceptance.
3. Follow the [Yacora and rotational-band comparison plan](docs/external_plasma_yacora_comparison.md).
   Use the identified 2023/2025/2026 references to fix model/data variants and
   overlapping states before freezing
   common-data numerical and native-model prediction cases. Use one observation
   model for both predictions.
4. Keep the spatial-localization null separate from a heating-model comparison:
   redistribute the same joint source for a fixed linear full-width equality
   check. Equal prescribed heating power does not make emitted mixtures equal.
5. Recover new experiment-linked geometry/reduction evidence. An unassessed
   official source lead is Lamppa et al., SAND-2021-12691,
   [DOI 10.2172/1825219](https://doi.org/10.2172/1825219),
   [official full text](https://www.osti.gov/servlets/purl/1825219).
   It is linked in the dissertation to the Mykonos MP3 platform; actual sample
   identity, geometry and diagnostic implications still need inspection.
6. Qualify actual MFEM fields, drive, units and material response before inferring
   a physical width-dependent source. Hall/plasma/motion physics is outside the
   stock Joule miniapp and remains a separate model requirement where relevant.

The [preserved unfinished NPL acquisition](supporting_documents/in_progress/mitl-upstream-20261019T0953/README.md)
has request logs, source hashes, citation review and draft code. Its new-source
page review and packet validation are incomplete. No result or property from it
has been promoted to a completed source assessment.

## Continuing work

The [scientific checkpoint](supporting_documents/research_checkpoint.json) and
each package's verification record retain prior milestones and the saved queue.
Historical cloud operational status was maintained separately in
`/workspace/research_work/mitl-upstream-status.json`; that path is not the local
Windows worktree's status file. The new implementation evidence records the
local baseline and pending work. Preserve unchanged passing evidence, record exact
commands/exits/hashes and use bounded local commits. Historical package replay
requires its declared complete baseline and source dependencies.
