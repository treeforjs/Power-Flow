# Legacy prototype audit

Reviewed baseline commit `d69a232cb979e7e27b9a1af8422624b0c3cc4c4a` on 18 October 2026.
This identifies blockers to a publication calculation. Existing outputs must not be
relabelled as validated FLASH, full water chemistry, or resolved CR results.

| Area | Evidence in baseline | Consequence |
| --- | --- | --- |
| Inputs | `configs/mykonos.yaml`: parametric current, zero voltage, generic gate/calibration, provisional rates enabled | Not a reconstructed measured shot |
| Electrode | `mhdlab/runner.py` instantiates `ReducedMHDSolver`; `mhdlab/mhd.py` provides a reduced induction/thermal model | FLASH is not the active heating backend |
| Desorption | `mhdlab/desorption.py`: emits rate times dt, then clips remaining coverage | Can emit more particles than the surface inventory |
| Source geometry | `mhdlab/neutrals.py`: injection into all angular bins; flux conversion uses min(dx,dy) | Does not implement correctly normalized outward thermal face emission |
| Sticking | `mhdlab/neutrals.py`: fixed per-step loss in wall cells | Wall behavior depends on timestep rather than impact statistics |
| Chemistry | `configs/water_fragment_mechanism.yaml` and `mhdlab/chemistry.py`: provisional channels, missing numeric rates | Species names do not establish a complete water reaction network |
| Neutral density | `total_neutral_density_m3` sums stored species including ions | Can give the wrong neutral density and E/N |
| CR | `mhdlab/cr_data.py` loads reference/manifests; `line_emissivity` uses an arbitrary scale times A times ne times donor density | No excited-state or rovibrational population solution |
| Spectral synthesis | `mhdlab/spectra.py`: Gaussian peak heights, single scalar radiometric response | No calibrated, unit-area transition-to-pixel forward model |
| Gate | `mhdlab/runner.py`: spectra use the final diagnostic state | Multiplying by a gate duration does not integrate an evolving gate |
| LOS | Default chord exceeds the domain; `integrate_los` clips samples to edge cells and weights the full length | Out-of-domain emission can be spuriously counted |
| Timing | `mhdlab/runner.py` advances a step before saving a sample labelled t | Potential dt offset relative to drive/gate |

Useful pieces include the repository layout, measured-trace import concepts, geometry
organization, reference manifests and diagnostic storage patterns. The new publication
path should use independently verified FLASH exports and conservative kinetic/CR
components. Passing the old smoke tests establishes software execution, not the
physical validity of the items above. The homegrown production MHD solver and its
coupled orchestration were removed on 20 October 2026. This audit records the
reviewed baseline; earlier source snapshots remain in Git history and the
versioned handoff packages. New heating inputs use the FLASH export contract.

See [the research specification](flash_h2o_publication_plan.md) for the replacement
architecture and [the recovery record](dissertation_recovery.md) for experimental
constraints.
