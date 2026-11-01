# Mills et al. (2004): thermophysical correlation source assessment

The official full text is accessible and supplies useful **evaluated correlations**, including mass-specific heat capacity and enthalpy for austenitic stainless steels and composition-dependent solid/liquid density expressions. It is a primary source for these recommended equations, rather than a new measurement campaign on a characterized dissertation foil. It does not supply a complete EOS or a calibrated volumetric constant-volume heat capacity.

K. C. Mills, Y. Su, Z. Li and R. F. Brooks, “Equations for the Calculation of the Thermo-physical Properties of Stainless Steel,” *ISIJ International* **44**(10), 1661–1668 (2004), [DOI 10.2355/isijinternational.44.1661](https://doi.org/10.2355/isijinternational.44.1661). [Official article](https://www.jstage.jst.go.jp/article/isijinternational1989/44/10/44_10_1661/_article); [official PDF](https://www.jstage.jst.go.jp/article/isijinternational1989/44/10/44_10_1661/_pdf). The separate [manifest](../data/materials/mills_2004_source_manifest.json) records acquisition, evidence and limits. No executable material adapter, property table or default is installed by this assessment.

## Access and inspected evidence

Both official requests returned HTTP 200 on 2026-10-18. The validated PDF has eight pages, 289,043 bytes and SHA256 `bf120dcbdf20d5f9b00494efd5e05643793edc0838deee8ba2c52988f9074836`. The current official landing page labels the article open access under CC BY-NC-ND 4.0; the PDF bears © 2004 ISIJ. The PDF, landing page, extracted text and inspection images remain in the external source cache, outside the repository. Only this assessment and its source metadata are added here.

All eight PDF pages were read as text. Page images were additionally inspected for printed pp. 1661–1663 and 1667: Table 1 and sample classifications; density notation; Eqs. 10–13; and Table 3's units, phase separation, reference row and footnotes. Acquisition commands, exits, HTTP results, byte counts, digests and compact equation checks are recorded under `/workspace/research_work/mitl-upstream-20261018T1858/mills2004_acquisition/`.

## Material identity, composition and phase

The authors combine cited measurements and recommended datasets for several stainless-steel families, fit equations and calculate recommended tables. Their solid heat-capacity discussion uses 304/316 data from Bogaard et al. and Mills's 2002 compilation; the liquid heat-capacity basis includes cited Chapman/NPL results. The article is not evidence of a known chemistry, oxide, stock condition or thermal history for the dissertation foil.

Table 1, p. 1661, gives **typical**, not certified, compositions in mass percent. Its illustrative 304 entry is C 0.04, Cr 19, Ni 9.3, Mo 0, Mn 1.0 and Si 0.5; a footnote says average C, Si and Mn values are assumed. These numbers must not become a foil composition default. The paper calls the 3-series structure metastable austenite and discusses separate ferritic/martensitic transformations. A grade label does not determine the actual foil's cold-work or phase fractions.

The abstract states an overall 300–1800 K scope. Individual tables and figures extend outside that summary interval: Table 3(a) begins at 298.2 K, while Table 3(b) includes liquid rows through 1900 K. Those entries are explicitly calculated recommendations, not authorization for a universal measured 298–1900 K material law.

Keep two temperature conventions separate. Table 1's illustrative 304 liquidus is **calculated as 1727 K** using Howe's coefficients. The authors **assume 1723 K for the 3-/6-series Table 3 calculations**, and 1773 K for 4-series alloys. Neither value establishes the dissertation foil's solidus, liquidus, superheating, melt fraction or latent-heat distribution. The single fusion jump in the recommended representation does not resolve an alloy's mushy interval.

## Heat capacity, enthalpy and reference check

On p. 1663, Eqs. 10–11 give the solid, austenitic recommendation, with temperature in kelvin:

\[
c_p(T)=472+0.136T-\frac{2.82\times10^6}{T^2}
\quad[\mathrm{J\,kg^{-1}\,K^{-1}}],
\]

\[
H_T-H_{298}=472T+0.068T^2+\frac{2.82\times10^6}{T}-156000
\quad[\mathrm{J\,kg^{-1}}].
\]

These are evaluated fits to experimental heat-capacity relations. The authors explain that ferritic heat capacity is harder to represent near its magnetic transition, roughly 500–1100 K; the smooth austenitic expression cannot be substituted for that peak. Their statement that enthalpy is less sensitive to those transitions does not make the heat capacities interchangeable.

Differentiating the printed Eq. 11 gives Eq. 10 exactly, including its negative inverse-square term. The printed integration constant is rounded: Eq. 11 evaluates to approximately **157.76 J/kg at 298 K**, rather than exact zero. Table 3(a) separately prints zero enthalpy at **298.2 K**; Eq. 11 evaluates to approximately 253.92 J/kg there. These are reference/rounding checks on printed coefficients, not measured offsets or error bars. A future adapter must declare its reference temperature and any explicit re-zeroing, preserve the source expression, and verify derivative and continuity conventions. No correction is installed here.

For liquid, Eq. 12 gives `cp = 800 ± 50 J/(kg K)`. Eq. 13 estimates fusion entropy as `delta_s_fus = delta_h_fus/T_liq = 160 J/(kg K)`. This is a separate fusion contribution; it is not already the integral of the solid heat-capacity polynomial. A phase-aware enthalpy representation must account for it once, define the liquid branch relative to the same reference, and qualify the adopted melting treatment.

**`cp` is mass-specific constant-pressure heat capacity.** It is not `cv`, and multiplying it by an assumed density does not establish the constant-volume volumetric heat capacity used by the diffusion benchmarks. The article does not provide the pressure-dependent thermodynamic closure needed for that conversion. A defensible approximation would need an explicit density/composition/phase convention and a stated mechanical constraint; none is inferred here. Likewise, mass-specific enthalpy is not automatically EOS internal energy.

## Density and unit conventions needing resolution

Eqs. 5–6, p. 1662, give composition-dependent density expressions with separate solid and liquid temperature terms, referenced to 298 K and 1823 K respectively. The solid expression subtracts `0.5*(T-298)`, while the liquid expression subtracts `0.86*(T-1823)`; density is displayed in kg/m³. The authors describe their density mixing assumption as approximate and compare calculated values with cited data. They also discuss a thermal-expansion refinement and slightly low liquid predictions for some 4-series cases.

Two conventions need explicit resolution before an executable density adapter:

- Eq. 4 defines its mixing variable as **mole fraction**; Table 1 lists **mass percent**; Eqs. 5–6 use unqualified `%Fe`, `%Cr`, etc. notation. This pass records that distinction rather than silently selecting a percentage basis or inferring the missing foil chemistry.
- The numerical thermal-expansion values in Eq. 1 correspond to the **micro-per-kelvin scale** shown as `10^6*alpha` in Fig. 1, although the displayed equation labels only K⁻¹. That scale must be reconciled before using its density refinement.

These qualifications do not negate the usefulness of the density correlations as a lead. They prevent treating the generic calculated density column in Table 3 as a measured, composition-certified 304 foil EOS.

## Uncertainties and useful next step

The article contains several distinct estimates: database density measurement uncertainty probably below 2%; solid heat-capacity measurement uncertainty about ±2% below 1000 K and ±5% at 1000–1700 K; calculated heat capacity/enthalpy uncertainty probably below 5%; liquid `cp = 800 ± 50 J/(kg K)`; and Table 3's separate property uncertainty footer. In particular, the liquid ±50 statement and the table's heat-capacity `<5%` footer should be retained as distinct source statements, rather than collapsed into one invented distribution. No coverage factor, confidence level, between-foil variability or nanosecond validation follows from these statements.

This source can support the next **guarded mass-specific enthalpy input contract**: explicitly selected austenitic solid branch, a declared enthalpy reference, bounded temperature/phase support, separately justified fusion treatment and retained evaluated-source provenance. It does not justify automatic stitching to the Flynn plate measurements, a default foil composition, an invented `Cv`, or a runtime curve through melting. The article itself notes that electrical/thermal transport depends on microstructure and that its transport recommendations favor minimum resistivity/maximum conductivity; those values are not a validation of the actual foil's pulsed heating.
