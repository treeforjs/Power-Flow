# Power-Flow

FLASH-based preparation for studying electrode heating, water desorption and VUV
spectroscopy in magnetically insulated transmission lines (MITLs).

The publication question is whether heating across a foil's width changes H₂O
release and fragment transport enough to distinguish the recorded spectra from
a uniform-desorption model. The [research specification](docs/flash_h2o_publication_plan.md)
defines the comparisons and the required energy and rovibrational accounting.

## Current workflow

1. Export explicitly configured FLASH surface-temperature histories in SI units.
2. Evaluate conservative finite-inventory desorption with separate removal,
   internal and translational energy channels.
3. Hand the complete time/temperature/state/number/geometry mixture to the external
   plasma and kinetic-transport models.
4. Compare predictions through a common observation model using recovered
   experimental evidence and justified calibration and uncertainty inputs.

The homegrown production MHD solver, coupled prototype runner and surrogate
material defaults have been removed. FLASH is the intended heating backend.
Independent magnetic, conductive-heating and Hall benchmarks remain as
verification references. Licensed FLASH and its target workstation must be
configured separately; no FLASH calculation has been run in this workspace.

The plasma collisional-radiative model will live in a separate repository.
Boltzmann, CR and kinetic-gap development remain with their owners. The
[external-model and Yacora comparison plan](docs/external_plasma_yacora_comparison.md)
includes staged common-data comparisons and published rotational-band extensions.
The plan identifies the rotational-band references; the external repository and
its executable interface are pending.

## Getting started

Use Python 3.11 or later and a virtual environment:

```sh
python -m pip install -r requirements.txt
```

These commands exercise explicitly manufactured inputs:

```sh
python tools/export_flash_surface.py --demo
python tools/run_surface_release.py --demo --plot
```

For real saved FLASH outputs, configure field names, units, geometry and snapshot
paths explicitly, then use:

```sh
python tools/export_flash_surface.py --config configs/flash_extract_example.json --output-dir runs/flash_surface
```

Read the [FLASH export contract](docs/flash_surface_contract.md) and
[heating-to-release workflow](docs/surface_release_pipeline.md) before configuring
an experimental case. The example configuration is a template, not an actual shot.

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
and qualified surface inputs. No physical full-water rovibrational spectrum or
experimental discrimination result is claimed.

Historical packages preserve the exact source versions and baseline requirements
of earlier work. Their prototype instructions are historical records; use the
current workflow above for new work.
