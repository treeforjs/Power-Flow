# SS304 thermal property handoff

Ready as an isolated conditional nominal model; parent owns integration and launches.

- Source: `/home/treefor/projects/runtimes/ss304-thermal-properties-20261029/source`.
- Native header `include/ss304_thermal.hpp`; compiled scalar C ABI `build/libss304_thermal.so`.
- Contract and exact source/approximation statements: `source/data/scenario.json` and `source/README.md`.
- MFEM callback mapping and restart requirements: `thermal-material-port-v1.json`.
- Frozen hash bindings: `source-manifest.json`; source acquisition/recovery: `source-evidence.json`.

## What is implemented

Solid-only mode supports300-1500K. Explicit nominal-synthesis mode supports300-1800K, including a reversible estimated1672-1727K mushy interval and a liquid branch. Inputs outside the selected temperature/energy range fail without clamping. Initial/cooling states below300K also stop.

Thermal Cp/H use Mills evaluated austenitic correlations. Ho-Chu Table2 provides an18-row nominal304 resistivity/thermal-conductivity table, with sigma=1/interpolated_resistivity. Its source pools304/304L in a quench-annealed recommendation and explicitly extrapolates above1500K. Those facts remain visible in the mode, scenario and each hot evaluation.

The mushy fraction/sensible-Cp bridge are declared model assumptions. The chosen latent fusion contribution276320J/kg is added once and differentiated in apparent heat capacity. Phase endpoints are estimated compilation markers, not actual-foil calibration. A liquid-model parameter1124.3K from the EOS manifest is never treated as a melting point.

The energy model is explicitly incompressible/fixed-reference-volume, using publishedCp as a caloric Cv proxy. Relative u_model=h, both zero at300K; rho_ref is caller supplied and constant. The example8000kg/m3 is a nominal304L EOS reference-density parameter adopted only as a scenario approximation. No real Cv or p(rho,T) is claimed. Moving mechanics requires separate admitted pressure-work/remapping. The mechanics component's773.15K limit is unchanged.

## Evidence and limits

29,912 assertions passed across source anchors,60-digit Decimal comparisons,2012 state samples, bounded inverses, phase continuity, dH/dT, integratedCp/latent-energy checks, signed heat/cooling transactions, subresolution residual bookkeeping and invalid-input/output-atomicity checks. Numeric limits are recorded inverification.json; they are arithmetic checks, not physical property uncertainties or changes to scientific guards. No GPU/MFEM/PIC library is linked. Source and nine shared input files stayed hash-identical during verification.

Both acquired primary PDFs exactly match the prior source manifests. Raw PDFs, text and rendered pages remain outside the source/Git packet. The source README provides primary citations and uncertainty qualifications.

No thermal solve or coupled run was launched. Vaporization, liquid motion, density/geometry changes, surface release/internal energies and wide-range EOS/transport remain absent. The1800K cap and independent mechanical guard must stop an inadmissible coupled trial. No physical-breakdown or foil-calibration claim follows from this verification.

Production was T with original identity/start_ticks and cgroup frozen1 before/after. No production signals, guard edits, dependency installs or shared-source changes were issued here.
