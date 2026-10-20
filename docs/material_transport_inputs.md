# Bounded material transport inputs

`mhdlab.material_transport` provides an explicit input table and evaluator for
one solid specimen context. It does not select material parameters for a
runner, supply a FLASH material model, or fill missing thermodynamic inputs.

`BoundedTransportTable` requires temperature, electrical resistivity, thermal
conductivity, material grade, sample condition, phase, parameter status, source
reference and source-data SHA256. Arrays and nested source metadata are copied
into immutable storage. The current contract accepts a single `solid` phase;
it cannot interpolate across melting or a change of phase. Manufactured tables
must carry `parameter_status="manufactured"`.

`evaluate(temperature_k)` accepts a scalar or a nonempty one-dimensional array
inside the table's closed temperature interval. It raises for extrapolation,
invalid values and an unrepresentable reciprocal. Its declared representation
is `piecewise_linear_in_temperature`: electrical resistivity rho and thermal
conductivity kappa are interpolated separately. Electrical conductivity sigma
is then derived as `1/rho`, in S/m. Interpolating sigma directly would produce
a different model. The interpolation is a modeling choice, not an assertion
that the primary measurements lie exactly on straight line segments.
The evaluator forms exact ratios/products of its supplied float inputs and
rounds the final properties to float64. This retains a representable property
contribution when its interpolation weight alone would underflow. It is a
bounded input evaluator, not an optimized solver material callback.

```python
from mhdlab.material_transport import load_flynn_304_transport

table = load_flynn_304_transport(
    "data/materials/flynn_1968_304_transport.csv",
    "data/materials/flynn_1968_304_transport.json",
    expected_material_grade="AISI Type 304",
)
values = table.evaluate([253.15, 303.15, 423.15])
rho = values.resistivity_ohm_m
sigma = values.electrical_conductivity_s_m
kappa = values.thermal_conductivity_w_m_k
```

The named loader accepts the audited ten-row recovery from D. R. Flynn,
*Thermal conductivity and electrical resistivity of a sample of AISI type 304
stainless steel*, NBS Report 9869 (1968),
[DOI 10.6028/NBS.RPT.9869](https://doi.org/10.6028/NBS.RPT.9869), printed p.12
(PDF p.16). These are published **smoothed results from measurements** of one
supplied solid plate sample. They are not raw acquisition points, dissertation
foil measurements, or grade-wide defaults. Chemical composition and thermal
history are unreported; the loader retains those unknowns and does not infer
annealing, cold work, oxide composition or agreement with a particular foil.

The supported interval is **253.15–423.15 K**. The CSV retains original Celsius
temperatures and resistivity mantissas. The loader checks decimal
`T_K = T_degC + 273.15` and the printed resistivity scale `1e-7 ohm m` before
converting to float64. Thermal conductivity has the same numerical value for
temperature gradients expressed per degree Celsius or per kelvin. Source units,
row count, grade/context and recovery status are required. Source-PDF and
canonical recovered-CSV digests are pinned, and the full recovery metadata is
retained. Changing table values and rehashing a sidecar does not inherit this
source identity; a revised recovery requires review and a versioned adapter.

The source's stated error limits remain source-specific prose/metadata. This
interface creates no error bars, Gaussian distribution, independent uncertainty,
confidence level or between-foil variability from them. It supplies no mass
density, heat capacity, enthalpy, EOS or phase-transition data. Heating beyond
423.15 K needs additional supported inputs; this table cannot be continued to
desorption temperatures by clamping, extrapolation or a silent material switch.
