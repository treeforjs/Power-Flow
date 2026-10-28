# Neutral spectroscopy implementation evidence

This record preserves the earlier FLASH-oriented implementation checkpoint.
The heating backend has since changed to MFEM Joule; see the
[current backend record](../mfem_joule/README.md). FLASH download is no longer
an active requirement. Earlier JSON reports keep their original source hashes
and scope; they are not rewritten as evidence for the new backend.

This local change continues the FLASH research baseline
`18d155b919e70b0b20165e22ce1eb4b7f48ac5a9`. The
[workflow](../../docs/neutral_spectral_workflow.md) describes the interfaces and
remaining physical inputs. [completion.json](completion.json) records commands,
outcomes, artifact hashes and explicit limits. No new physical spectrum, absolute
calibration correction or neutral fit is claimed.

Implemented components:

- Hash-bound pilot manifest, 157 new wavelength samples for shot 15750 and 183
  for shot 15765, readout envelopes and explicit omitted segments.
- Separate dimensional/reduction/response scenarios and the 41-point radiance
  sensitivity scan; derived transition-density displays remain unscaled.
- Literal readout request packets and a manifest-bound comparison guard.
- Receiver receipt, source placement, number/energy/momentum/moment accounting,
  state/rate/clock/velocity/opacity checks and ray-cell projection.
- Finite owner-ensemble comparison with asymmetric maximum-envelope scores,
  alternative inventory/thermal solutions and synthetic degeneracy checks.
- Independent projected-feature resolution/uncertainty qualification.
- Private WSL parallel HDF5 prerequisite build and passing two-rank MPI test.

Accepted focused verification was preserved and only affected checks rerun:

| Check | Observed outcome |
| --- | --- |
| Pilot evidence suite, including exact external PDF recovery | 20 passed before the coverage-only admission addition; its affected 15 admission cases subsequently passed |
| Neutral ensemble module and CLI | 47 passed before sentinel hardening; 32 added sentinel cases subsequently passed |
| Receiver/geometry/projection contract | 52 passed |
| Final convergence, pilot bridge and affected pilot admission guards | 37 passed |
| Pilot evidence validation | 14 referenced hashes valid; physical-ready mode exits 2 as expected |
| Private MPI/HDF5 Fortran compile and collective two-rank file creation | Passed with scoped library runtime path |
| Whitespace and Python parsing | Passed |

The current local runtime outputs remain under ignored `runs/`:

- `runs/dissertation_pilot_readiness_final.json`
- `runs/pilot_requests_final/15750_wavelength.json`
- `runs/pilot_requests_final/15765_wavelength.json`
- `runs/neutral_ensemble_manufactured_20261027_guarded/`
- `runs/raikou_receiver_contract_control/projection.json`
- `runs/projected_convergence_manufactured_final/report.json`
- `runs/flash_preflight_hdf5_runtime/report.json`

Compact accepted JSON reports are copied alongside this record. Exact command
stdout/stderr for the final readiness/request/qualification checks is retained in
`final_commands.json`. Third-party source pixels, the dissertation PDF, installed
libraries and build trees remain external. Their identities and access paths are
recorded without redistributing those files.

FLASH source has not been downloaded. Its version-specific build and physical
benchmarks await acquisition. Physical Raikou source acceptance, complete rates
and populations, reduction/collection/clock evidence, and a sourced monolayer
convention remain pending. The original detector data are unavailable. No remote
push, merge, publication, external message or production model run occurred.
