# SS304 thermodynamics and an ETHOS-informed extension

The material work now has a bounded EOS-table handoff, a recovery of published
304L solid/liquid model parameters, and an explicit generation request. It has
**not generated or installed a physical SS304 EOS**, changed the MFEM heating
fields, or extrapolated the existing solid heat-capacity correlation into the
hot corners.

## What ETHOS supplies

[Stanek et al., ETHOS (2024)](https://doi.org/10.1063/5.0237197) develops a
multi-fidelity framework for transport-table ensembles. Its demonstrated
application fits a parameterized Lee-More-Desjarlais electrical-conductivity
model to beryllium data and propagates parameter uncertainty into MHD results.
This is a useful method for combining sparse measurements, QMD and average-atom
calculations. It is not itself an acquired SS304 pressure/internal-energy table.
No matching public ETHOS source repository or SS304 table was located in the
sources checked for this assessment. No unrelated software named Ethos is used.

EOS and transport remain separate inputs. An electrical-conductivity model
does not provide pressure, internal energy, latent heat or phase coexistence.
The existing Flynn transport data and Mills solid enthalpy branch therefore
retain their original bounds. Mills supplies mass-specific Cp; it is not
automatically Cv, a mass density, or a multiphase EOS.

## A directly relevant alloy EOS

[Sjostrom (2023)](https://www.osti.gov/servlets/purl/2246835), DOI
[10.1063/12.0020364](https://doi.org/10.1063/12.0020364), describes SESAME
**4273** for 304L. It combines solid/liquid Helmholtz models with an electronic
thermal contribution computed using Tartarus. Its calculations use nominal
mass fractions Fe=0.71, Cr=0.19 and Ni=0.10. The author's summary recommends
the EOS for 304 as well as 304L within its stated uncertainty; this does not
establish the actual dissertation foil composition or condition.

The 17 rows of Table 1 are recovered in
[the parameter CSV](../data/materials/ss304_eos_4273_published_parameters.csv)
with a [source manifest](../data/materials/ss304_eos_4273_published_parameters.json).
The table was visually checked on PDF page 7, printed page 310009-5, of the
author/LANL accepted manuscript. Original printed values and units are retained.
Its reference density is 8.0 g/cm3, not a measured density for our foil. The
solid energy convention is u=0 at 298.15 K and one atmosphere. The printed
liquid T_melt=1124.3 K is a model parameter, not an inferred ambient melting
temperature. These parameters do not include the full electronic free-energy
grid or all equations needed to reconstruct the EOS.

The paper directs requests for numeric EOS data to the LANL SESAME contact.
No request has been sent. The numeric 4273 table has not been recovered into
this workspace. Acquiring that table is preferable to recreating an alloy EOS
from pressure plots or assigning independent analytic tails to p and u.

[Mattsson and Desjarlais (2004)](https://www.osti.gov/servlets/purl/920130),
SAND2004-5253, is also directly relevant to ALEGRA. It fits an alloy QMD
soft-sphere model and blends it with SESAME 4272. Printed pages 16-17 explicitly
enforce the identity

```
rho^2 (du/drho)_T = p - T (dp/dT)_rho
```

by thermodynamic integration after blending pressure. Thus smoothly joining
independent pressure and energy curves is insufficient. The report's plotted
model and the recovered 4273 parameters are evidence, not executable complete
EOS tables.

## Explicit generation request

[ss304_ethos_eos_request.json](../configs/ss304_ethos_eos_request.json) records
the nominal Fe-Cr-Ni candidate, source/model references, required fields,
energy conventions, source gaps and acceptance checks. Its initial requested
grid has 15 density points over 100-16000 kg/m3 and 22 temperatures over
300-30000 K. These are requested coverage/refinement targets, not measured
or admitted state points.

The authorized extrapolation remains a generation step: extend a
source-conditioned physics model or use an existing wide-range EOS, retain
the data-support versus model-extension regions, and propagate uncertainty.
No numerical extrapolation rule is silently installed. A free-energy model or
thermodynamic integration must keep the fields mutually consistent. Phase
coexistence and latent heat must come from the selected model. Conductivity
and thermal-conductivity feedback require their own admitted inputs.

For the current fixed-geometry heating path, the first useful coupling would
be u(rho,T) and its bounded T(rho,u) inverse after source qualification. Adding
an EOS alone does not turn the stationary MFEM driver into moving MHD or
resolve gas/plasma transport. Its current constant-property run is unchanged.

## Bounded table interface

`mhdlab.eos_tables.load_eos_table` accepts a complete rectangular SI CSV:

```text
density_kg_m3,temperature_k,pressure_pa,internal_energy_j_kg
```

The adjacent JSON must use schema `mhdlab-bounded-eos-table-v1`, exact CSV
SHA256, explicit SI units, material grade, table version, source reference,
energy zero, phase policy/reference and a thermodynamic-consistency reference.
The caller explicitly selects `expected_material_grade`; no table is selected
by default. `parameter_status` distinguishes `source_table` from
`manufactured`. A phase policy must name a single branch or an equilibrium
multiphase table constructed by its source.

```python
from mhdlab.eos_tables import load_eos_table

# csv_path and metadata_path must name an actually acquired, assessed table.
eos = load_eos_table(csv_path, metadata_path, expected_material_grade=selected_grade)
values = eos.evaluate(density_kg_m3, temperature_k)
temperature = eos.temperature_from_internal_energy_k(density_kg_m3, specific_energy_j_kg)
```

The offline evaluator uses declared bilinear interpolation in density and
temperature. It preserves signed pressure and reference-dependent internal
energy. It rejects missing/duplicate grid cells, nonfinite inputs, inconsistent
units, nonmonotone isochoric energy, and every query outside the supplied
support. It never clamps an out-of-range query. The inverse includes the
forward evaluator's rounded endpoint energies. Arrays and nested metadata are
immutable. Exact rational interpolation weights and sums are rounded once to
binary64.

Interpolation is a numerical observation of supplied p/u fields; it does not
establish thermodynamic consistency between them. A source-consistent table
can also acquire discretization errors under bilinear interpolation. This
adapter does not supply entropy, sound speed, Cv, a Helmholtz derivative
interpolant, mechanical stability or an MFEM constitutive callback. Those
remain explicit coupling qualifications, not claims implied by loading a CSV.

Eight focused tests pass for an explicitly manufactured ideal-gas table,
including the p/u inverse, negative reference energy, adjacent-float domain
rejection, phase/unit/provenance checks, rectangular coverage and immutability.
These are interface checks, not SS304 physical validation:

```powershell
python -c "import sys; sys.meta_path[:] = [f for f in sys.meta_path if type(f).__module__ != '_editable_skbc_mhdlab']; import pytest; raise SystemExit(pytest.main(['tests/test_eos_tables.py','-q']))"
```

Downloaded primary reports remain in `.cache/ss304-eos-primary`. Their byte
hashes are recorded in the published-parameter manifest; no copyrighted full
report or numeric SESAME table is added to the public data tree.
