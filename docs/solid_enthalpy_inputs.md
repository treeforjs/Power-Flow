# Explicit solid enthalpy inputs

`mhdlab.solid_enthalpy` evaluates the solid austenitic mass-specific heat capacity
and enthalpy recommended by Mills et al. (2004), printed p. 1663, Eqs. 10–11,
[DOI 10.2355/isijinternational.44.1661](https://doi.org/10.2355/isijinternational.44.1661).
The [source assessment](mills_2004_source_assessment.md) and
[reviewed manifest](../data/materials/mills_2004_source_manifest.json) preserve the
evidence and limits. This is an explicitly selected evaluated family correlation.
It does not characterize the dissertation foil or install a solver default.

The source equations use kelvin and mass-specific SI units:

\[
c_p(T)=472+0.136T-2.82\times10^6/T^2
\quad[\mathrm{J\,kg^{-1}\,K^{-1}}],
\]

\[
H_{\rm printed}(T)=472T+0.068T^2+2.82\times10^6/T-156000
\quad[\mathrm{J\,kg^{-1}}].
\]

The printed enthalpy is labelled relative to H298, but its rounded constant gives
about 157.76 J/kg at 298 K. Table 3 separately prints zero at 298.2 K; the equation
gives about 253.92 J/kg there. These source reference differences remain visible
in metadata. Neither reference is silently substituted for the caller's choice.

## Required selection

```python
from mhdlab.solid_enthalpy import load_mills_2004_solid_enthalpy

solid = load_mills_2004_solid_enthalpy(
    "data/materials/mills_2004_source_manifest.json",
    temperature_min_k=300.0,
    temperature_max_k=1500.0,
    reference_temperature_k=800.0,
    reference_convention="rezeroed_eq11",
    expected_material_family="austenitic stainless steel",
    expected_phase="solid",
    material_qualification="Evaluated austenitic family only; actual foil composition and condition unknown",
    phase_qualification="Selected solid branch over this interval; actual foil phase requires separate verification",
)
values = solid.evaluate([700.0, 800.0, 900.0])
temperature_k = solid.temperature_from_enthalpy_k(values.enthalpy_j_kg)
```

Every keyword is required. The explicit supported interval must satisfy
`300 <= Tmin < Tmax <= 1700 K`; the declared `Tref` must lie within that interval.
This conservative ceiling stays below the source's assumed generic melting
temperature and within its discussed solid heat-capacity uncertainty range. It
does not establish the actual foil's solidus or guarantee austenite under a
particular heating history. Select a narrower interval when the physical phase
qualification requires one. A grade name alone cannot establish that condition.

The factory accepts only the stated family and solid phase. Caller qualification
strings preserve the reason and unresolved applicability; they are declarations,
not independent phase verification. The adapter pins the exact reviewed manifest
SHA256 `ff64f70a2c605cf1024cc8ecfff9ac113d325981feb4efb9cc69dc5340205a8a`,
the official source-PDF metadata SHA256
`bf120dcbdf20d5f9b00494efd5e05643793edc0838deee8ba2c52988f9074836`,
source identity, equation text and evidence pages. A changed assessment needs
review and a revised adapter. Evaluation requires no local PDF or network access;
the factory qualifies the assessed source metadata rather than claiming to
re-download or independently authenticate absent PDF bytes.

## Reference conventions and methods

| Method or attribute | Meaning |
| --- | --- |
| `cp_j_kg_k(T)` | Printed Eq. 10, mass-specific constant-pressure heat capacity |
| `printed_enthalpy_j_kg(T)` | Printed Eq. 11, preserving its rounded constant |
| `delta_enthalpy_from_reference_j_kg(T)` | Stable integral of Cp from the declared Tref |
| `enthalpy_j_kg(T)` | Enthalpy in the explicitly selected reference convention |
| `temperature_from_enthalpy_k(H)` | Bounded inverse of that selected convention |
| `evaluate(T)` | Immutable `SolidEnthalpyValues` with the five quantities named below |
| `supported_temperature_range_k` | Included temperature endpoints |
| `supported_enthalpy_range_j_kg` | Included H endpoints in the selected convention |
| `printed_enthalpy_at_reference_j_kg` | Printed Eq. 11 evaluated at declared Tref |
| `metadata` | Detached JSON-compatible record of the selections and source |
| `source_metadata` | Deeply immutable reviewed source assessment |

`reference_convention="printed_eq11"` preserves the printed absolute expression.
Declaring Tref still records a comparison reference and provides its stable
enthalpy difference; it does not re-zero the printed expression.
`reference_convention="rezeroed_eq11"` sets the declared enthalpy to exactly zero
at the supplied Tref. It evaluates the integral directly:

\[
\Delta H=(T-T_{\rm ref})
\left[472+0.068(T+T_{\rm ref})-
\frac{2.82\times10^6}{TT_{\rm ref}}\right].
\]

This factorization retains adjacent-temperature differences that subtraction of
two large printed enthalpies can lose. Evaluation uses exact rational printed
coefficients and the exact binary values of supplied float64 temperatures,
rounding each final Cp or enthalpy to float64 once. The forward enthalpy endpoints
and inverse bounds therefore share that rounding; no outward tolerance admits
values beyond those bounds. This is numerical evaluation of the printed fit,
not additional accuracy in the underlying physical recommendation.
`evaluate` returns `temperature_k`,
`cp_j_kg_k`, `printed_enthalpy_j_kg`, `delta_enthalpy_from_reference_j_kg` and
`enthalpy_j_kg`; each quantity preserves the input scalar or nonempty array shape.
Scalars are Python floats. Arrays own immutable byte storage and cannot have
their writeable flag re-enabled. Input mutation and mutation of detached metadata
cannot change the model or past results.

All numeric inputs must be finite real scalars or nonempty arrays. Booleans,
strings, complex numbers and masked arrays are rejected, including mixed arrays
that could otherwise coerce booleans into numbers. Every temperature method and
the inverse includes endpoints and rejects every out-of-domain value. There is
no clamping, extrapolation or implicit phase crossing.

Cp is positive and increasing over the admitted interval, so the enthalpy branch
is monotone. The inverse uses safeguarded Newton steps and adjacent-float bracket
checks with the factored integral as its residual. It introduces no smoothing or
interpolated material table. The residual retains exact rational source values
before each Newton update. Float64 precision of a supplied printed-offset H can
limit the recovered temperature resolution; an inverse cannot restore digits
already rounded out of that input. For small differences around Tref, use the
explicitly re-zeroed convention and the factored difference.

## Physical limits

These are evaluated equations based on cited 304/316 heat-capacity relations,
with no certified dissertation-foil composition, oxide, cold work or phase
fraction. The assessment retains the source's distinct uncertainty statements:
about ±2% for solid heat-capacity measurements below 1000 K and ±5% from
1000–1700 K; calculated Cp/enthalpy probably below 5%; and Table 3's separate
footer. The adapter does not turn them into a confidence interval, distribution,
pointwise error bar or nanosecond validation.

Mass-specific `Cp` is not `Cv` or volumetric heat capacity. Enthalpy is not
automatically EOS internal energy. Density, pressure/mechanical closure,
composition, actual phase validity, magnetic effects and pulsed-heating evidence
require separate inputs. No liquid branch, fusion jump, mushy interval,
transport-data splice or Flynn/Wilkes extrapolation is supplied. Using this input
in a physical evolution model requires separately justified equations and phase
coupling; no such coupling is installed by this module.
