# Conductivity Table Direction

Primary design reference:

- L. J. Stanek et al., "ETHOS: An automated framework to generate
  multi-fidelity constitutive data tables and propagate uncertainties to
  hydrodynamic simulations," Phys. Plasmas 31, 102707 (2024),
  `https://doi.org/10.1063/5.0237197`.

ETHOS motivates three pieces that are useful here:

- Treat electrical conductivity as a density-temperature table, not only a
  one-line analytic closure.
- Keep ensembles of valid tables so uncertainty in conductivity can be
  propagated through MHD outputs.
- Make the table-generation path able to incorporate sparse high-fidelity
  measurements/calculations plus lower-fidelity broad-coverage models.

Current implementation:

- `mhdlab.conductivity_tables.ConductivityTable` reads CSV, NPZ, and HDF5
  tables with axes `density_kg_m3`, `temperature_k`, and
  `conductivity_s_m`.
- Tables may be 2D or 3D. A 3D table is interpreted as
  `ensemble_index, density, temperature`.
- Interpolation is log-density/log-temperature/log-conductivity and clamps to
  table bounds.
- `mhd.conductivity.model: table` selects this path.
- `tools/build_conductivity_table.py` generates starter tables from the current
  Knoepfel closure with optional lognormal ensemble perturbations.

Not implemented yet:

- The Lee-More-Desjarlais model used by ETHOS.
- Bayesian posterior sampling from multi-fidelity data.
- Consistent EOS/transport melt-transition enforcement.
- Optimal ensemble selection for uncertainty propagation.

Near-term path:

1. Use direct EC Knoepfel for fast checks below melt.
2. Use generated table ensembles to exercise the table interface and run
   sensitivity studies.
3. Replace the generator internals with curated SS304/LMD/DFT/AA or literature
   conductivity data as those tables are added.
