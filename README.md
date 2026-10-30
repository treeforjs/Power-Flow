# Power-Flow

**Restart handoff:** [Arch setup, selected 2048x256 timed run, temperature
views and deferred physics questions](docs/ARCH_RESTART.md). Work is stopped;
the new native run has not been launched.

MFEM Joule based preparation for studying electrode heating, water desorption and VUV
spectroscopy in magnetically insulated transmission lines (MITLs).

The publication question is whether heating across a foil's width changes H₂O
release and fragment transport enough to distinguish the recorded spectra from
a uniform-desorption model. The [research specification](docs/mfem_h2o_publication_plan.md)
defines the comparisons and the required energy and rovibrational accounting.

## Current workflow

1. Export explicitly evaluated MFEM Joule surface-temperature histories in SI units.
2. Evaluate conservative finite-inventory desorption with separate removal,
   internal and translational energy channels.
3. Hand the complete time/temperature/state/number/geometry mixture to the external
   plasma and kinetic-transport models.
4. Compare predictions through a common observation model using recovered
   experimental evidence and justified calibration and uncertainty inputs.

The homegrown production MHD solver, coupled prototype runner and surrogate
material defaults have been removed. MFEM's `joule.cpp` is the selected heating backend.
Independent magnetic, conductive-heating and Hall benchmarks remain as
verification references. The [MFEM backend](docs/mfem_joule_heating.md) targets
stationary solid-electrode heating; physical foil/drive/material qualification
remains required. FLASH acquisition is no longer an active dependency.

The current foil calculation uses a custom MFEM 2D magnetic-vector-potential
extension with a finite-width return conductor and nonconducting vacuum. It
enforces the dissertation shot 17133 current from Fig. C.19 and stops at 60 ns
for the early VUV window. This is a provisional constant-property calculation;
it does not run the unchanged stock Joule kernel. The
[ALEGRA comparison](docs/mfem_alegra_comparison.md) preserves the published
40/60 ns profiles, unresolved edge probe placement and clock assumptions.
The [four-face delivery](artifacts/neutral_desorption_perimeter_20261028/README.md)
contains equal-scale overview and edge-close-up GIFs, the comparison, all four
material-face inputs and replay code. The
[perimeter source](docs/neutral_perimeter_release.md) explains side areas,
physical launch directions and why the earlier y-enlarged view looked different.

The [collection/CR work](docs/target_collection_and_cr.md) adds finite-pupil
photon collection, foil shadowing and edge-overlap calculations using the
stated 19-inch mirror distance, approximately8-inch separation and
50-micrometre by1-inch slit. Initial H2O monolayers is the primary free spectral
regression parameter in the new common-water source. Fragments come from
chemistry; the earlier independent-pool movie remains a separate control.
Raikou's sourced H2 radiative graph and native transient control provide the
CR starting point. A bounded [SS304 EOS interface](docs/ss304_ethos_eos.md)
is prepared; physical collision rates, optical response and numeric alloy EOS
inputs still require qualification.
The [optics/CR delivery](artifacts/neutral_optics_cr_20261028/README.md)
contains the computed overlap maps, state-resolved radiative control, water
parameter contract, EOS assessment and replay inputs.

The plasma collisional-radiative model lives in the separate
[Raikou repository](https://github.com/treeforjs/raikou).
Boltzmann, CR and kinetic-gap development remain with their owners. The
[external-model and Yacora comparison plan](docs/external_plasma_yacora_comparison.md)
includes staged common-data comparisons and published rotational-band extensions.
The plan identifies the rotational-band references. The new
[Raikou observation handoff](docs/raikou_observation_handoff.md) registers its
population-output interface and a manufactured native readback control. Physical
state mapping, source admission and the qualified excitation/transfer model
remain separate owner inputs.

The [neutral comparison implementation](docs/neutral_spectral_workflow.md)
connects versioned dissertation readouts, separate reconstruction/calibration
scenarios, a resolved receiver contract and finite-ensemble comparison. It
includes the shot 15750/15765 wavelength recoveries and synthetic recoverability
checks. MFEM 4.10 with MPI/Hypre/Metis is installed in Ubuntu 24.04 WSL2.
The older [WSL setup record](docs/flash_wsl_setup.md) preserves the private HDF5
build and MPI evidence as history.

## Getting started

Use Python 3.11 or later and a virtual environment:

```sh
python -m pip install -r requirements.txt
```

The native stock Joule demonstration can be built against the installed WSL MFEM:

```powershell
python tools/probe_mfem_joule.py --wsl-distro Ubuntu-24.04 --source-root /home/tsmit/src/mfem-4.10 --build-dir /home/tsmit/.local/share/power-flow/mfem-joule-build --output-dir runs/new_mfem_joule_demo
```

This verifies native execution/readback in demonstration units. Physical foil
temperature requires the [qualified setup](docs/mfem_joule_heating.md).

These commands exercise explicitly manufactured inputs:

```sh
python tools/export_mfem_surface.py --demo --output-dir runs/mfem_surface_demo
python tools/run_surface_release.py --demo --plot
python tools/validate_pilot_comparison.py --report runs/pilot_readiness.json
python tools/compare_neutral_ensemble.py --config configs/neutral_ensemble_manufactured.json --output-dir runs/neutral_comparison --plot
```

For explicitly evaluated MFEM face samples, configure field semantics, SI units,
temperature reference, geometry and source hashes, then use:

```sh
python tools/export_mfem_surface.py --config configs/mfem_surface_example.json --output-dir runs/mfem_surface
```

Read the [MFEM surface contract](docs/mfem_surface_contract.md) and
[heating-to-release workflow](docs/surface_release_pipeline.md) before configuring
an experimental case. The example configuration is a template, not an actual shot.
The earlier FLASH readers remain available for explicitly labeled saved inputs.

## Results, supporting documents and handoff

- [Start with the project handoff](HANDOFF.md): completed work, ownership,
  verification evidence, limitations and next steps.
- [Supporting-document catalog](supporting_documents/README.md): all 25 completed
  handoff packages, readable reports, checksums and preserved unfinished work.
- [Heating verification](docs/flash_heating_verification.md),
  [Hall reference](docs/hall_induction_verification.md),
  [width-source comparison](docs/thermal_width_power_reference.md) and
  [joint boundary quadrature](docs/boundary_quadrature.md).
- [Recovered dissertation evidence](docs/dissertation_recovery.md),
  [observed-spectrum recovery](docs/observed_spectrum_recovery.md),
  [absolute calibration](docs/absolute_spectral_calibration.md) and
  [coverage inference](docs/surface_coverage_inference.md).

Manufactured examples, reconstructed figure values and sourced measurements keep
distinct provenance labels. Physical monolayer coverage and an absolute spectral
correction remain unresolved without the original calibration/reduction records
or independently justified reconstruction scenarios and qualified surface inputs.
The original detector data are unavailable. No physical full-water rovibrational spectrum or
experimental discrimination result is claimed.

Historical packages preserve the exact source versions and baseline requirements
of earlier work. Their prototype instructions are historical records; use the
current workflow above for new work.
