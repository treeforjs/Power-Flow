# Conductivity Tables

This folder is for electrical-conductivity tables used by the MHD solver.

The table interface is ETHOS-inspired: tables are functions of material density
and temperature, may contain ensemble members for uncertainty propagation, and
are interpolated in log-density/log-temperature/log-conductivity space. It does
not yet implement the full ETHOS Bayesian/LMD workflow.

Accepted formats:

- CSV long form with columns `density_kg_m3`, `temperature_k`,
  `conductivity_s_m`, and optional `ensemble_index`.
- NPZ with arrays `density_kg_m3`, `temperature_k`, `conductivity_s_m`.
- HDF5 with datasets of the same names.

Generate a starter SS304 table from the current Knoepfel closure:

```powershell
python tools/build_conductivity_table.py --output data/conductivity/ss304_knoepfel_table.h5 --ensemble-count 16 --uncertainty-fraction 0.2
```

Use a table in a config:

```yaml
mhd:
  conductivity:
    model: table
    file: ../data/conductivity/ss304_knoepfel_seed.csv
    ensemble_index: 0
    minimum_s_m: 1.0e5
```
