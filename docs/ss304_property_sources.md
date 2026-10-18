# SS304 property source acquisition

The acquired evidence supports bounded comparisons for solid 304 and separately identified 304L. It does **not** yet supply a validated, continuous SS304 material model from room temperature through melting for the dissertation foils. No property curve, extrapolation, material default or solver configuration is installed here.

The companion [`ss304_source_manifest.json`](../data/materials/ss304_source_manifest.json) records source identity, cached-file hashes, access outcomes, evidence locations and material limits. PDFs and inspection images remain outside Git in `/workspace/research_sources/ss304_properties/`. This acquisition was checked on 2026-10-18.

The subsequent [bounded transport-input interface](material_transport_inputs.md)
uses a separately [recovered ten-row Flynn table](../data/materials/flynn_1968_304_transport.csv)
with [source guards and provenance](../data/materials/flynn_1968_304_transport.json).
That interface stops at 423.15 K and supplies no heat capacity, phase transition
or high-temperature continuation. It is not an adopted foil material model.

## Acquired evidence and its bounds

| Source | Material and property evidence | Verified domain and limits |
| --- | --- | --- |
| D. R. Flynn (1968), [NBS Report 9869](https://doi.org/10.6028/NBS.RPT.9869) | Primary measurements of electrical resistivity and thermal conductivity for a sample identified as AISI 304. Full official PDF recovered. | −20 to +150 °C (253.15–423.15 K), solid plate sample. Printed p. 1/PDF p. 5 explicitly says chemical analysis and thermal history were not provided. Smoothed values and estimated errors: printed p. 12/PDF p. 16. |
| K. E. Wilkes et al. (1997), [ORNL/CP-96813](https://doi.org/10.2172/658177), *Thermophysical properties of stainless steel foils* | Primary direct-electrical-heating measurements; the two 304 foils are 0.305 and 0.203 mm nominal thickness (Table I, PDF p. 4). Full official OSTI scan recovered. | Resistivity: room temperature to about 500 °C, followed by cooling (PDF p. 6). Conductivity is reported over a narrower range: source excludes temperatures below about 50 °C and poorly resolved high-temperature gradients (PDF pp. 9–10). Figure 5 has a 0–400 °C axis, with retained 304 measurements occupying only part of it; no numerical curve has been digitized here. |
| C. S. Kim (1975), [ANL-75-55](https://doi.org/10.2172/4152287), *Thermophysical properties of stainless steels* | An evaluation/estimation report for **304L and 316L**, not a new measurement of all listed properties. Full official OSTI PDF recovered. Contains solid enthalpy/cp correlations, estimated liquid cp and fusion enthalpy. | Overall tabulated range 300–3000 K includes extrapolation and constituent-additivity estimates. For 304L, cited solid enthalpy/cp data extend to 1620 K; the report extrapolates toward a stated 1670–1730 K melting range (printed p. 2/PDF p. 11). It uses 1700 K as a representative melting temperature. A consequential liquid-enthalpy table conflict is preserved below. |
| A. Cezairliyan and A. P. Miiller (1980), [DOI 10.1007/BF00506274](https://doi.org/10.1007/BF00506274), *Int. J. Thermophys.* **1**, 83–95 | Primary publisher **abstract** accessible; full article PDF not recovered. The material is explicitly AISI **304L**. Abstract supplies cp/resistivity correlations and estimated inaccuracies from subsecond measurements. | 1400–1700 K for cp/resistivity. Reported melting point 1707 K, estimated inaccuracy 8 K. This is not a measured equilibrium solidus/liquidus pair for the dissertation's SS304. Composition, full experimental qualifications and latent heat remain unverified. Numerical evidence location is the publisher abstract, not an inspected article page. |
| C. Y. Ho and T. K. Chu (1977), *Electrical resistivity and thermal conductivity of nine selected AISI stainless steels*, CINDAS Report 45, DTIC ADA129160, [public scan](https://www.inductor-jmag.ru/files/content/a129160.pdf) | Original report scan recovered from a public mirror after the DTIC PDF endpoint failed. A critically evaluated compilation; its 304 recommendations explicitly pool 304 and 304L measurements (printed p. 13/PDF p. 20). | 300-series recommendations are for a quench-annealed condition. Values above 1500 K are explicitly extrapolated (printed p. 2/PDF p. 9). Table 2's 1672–1727 K melting-range markers are compilation inputs, not independently measured phase boundaries for these foils. Its 304 and 304L tables must not be treated as independent grade-specific evidence. |
| R. S. Graves et al. (1991), [DOI 10.1007/BF00500761](https://doi.org/10.1007/BF00500761), *Int. J. Thermophys.* **12**, 409–415 | Publisher abstract and official [OSTI metadata](https://www.osti.gov/biblio/7253632) recovered; full text unavailable. Identifies **304L** transport measurements. | Abstract states thermal conductivity, resistivity and Seebeck coefficient measurements from 300 to 1000 K and diffusivity from 297 to 423 K. Sample chemistry, numerical tables and full uncertainty budget were not inspected. This is a retrieval lead, not an ingested dataset. |
| J. N. Sweet, E. P. Roth and M. Moss (1987), [DOI 10.1007/BF00503645](https://doi.org/10.1007/BF00503645), *Int. J. Thermophys.* **8**, 593–606 | Publisher abstract recovered; full text unavailable. Identifies **304** comparative/flash-diffusivity measurements, separately from Inconel 718. | Abstract reports 0–700 °C, but the underlying 304 chemistry, cp tables and uncertainty budget were not inspected. The related SAND86-0840 report was identified in official OSTI metadata; no full text was recovered in this pass. |

### Solid 304 numerical anchors, not defaults

Flynn's printed p. 12/PDF p. 16 table was visually inspected. At 20 °C it gives resistivity `7.17e-7 ohm m` and conductivity `14.8 W/(m K)`; at 150 °C, `8.19e-7 ohm m` and `17.9 W/(m K)`. These are smoothed results for that one supplied sample. The author estimates resistivity error no greater than ±0.5% and conductivity uncertainty no greater than 2–3%, chiefly from uncorrected insulation heat loss. The report does not state a modern coverage factor, confidence level or between-heat variability bound. Its unknown composition/history and low upper temperature prevent using these figures as universal foil constants or extending them to 1000 K.

Wilkes's primary foil results are particularly relevant to stock condition. Its 304 resistivity is about 1.5–2% above the CINDAS annealed reference (PDF pp. 6–7); its 304 conductivity lies about 10–15% below that reference at the low end of the measured interval (PDF p. 9). These offsets are **comparisons**, not reported material-property uncertainties. The authors suggest cold work and thickness/roughness effects, but describe the rolled foils as *most likely* cold worked rather than establishing a measured treatment state. No batch composition is supplied with Table I.

The conductivity reduction is sensitive to coupled inputs. Table II, PDF p. 11, examines a **Type 201 specimen**, not 304: an assumed 1% resistivity change produces a 12% change in inferred conductivity, and an assumed 1% emittance change produces 11%. The authors say the chosen perturbations are arbitrary and need not represent actual uncertainties. Neither these sensitivities nor scatter among different specimens can become an SS304 uncertainty distribution. The report supports retaining stock condition and the measurement's correlation structure when a later dataset is recovered.

### Heat capacity, enthalpy and phase provenance

Kim's solid 304L correlation uses the reference `H(T)−H(298.15 K)` in cal/g; printed pp. 2–3/PDF pp. 11–12 give a solid enthalpy polynomial and its cp derivative. The liquid cp `0.190 cal/(g K)` and fusion enthalpy `64.0 cal/g` are **estimated** by constituent additivity using 72% Fe, 19% Cr and 9% Ni by mass. Those proportions are a calculation surrogate, not a measured certificate for a particular 304L heat or the dissertation foil. No quantified cp/fusion uncertainty was found in the inspected discussion. No conversion of these estimates into an adopted SI model is performed here.

There is a verified internal conflict in that report:

- Table 2, printed p. 5/PDF p. 14, prints liquid 304L enthalpy at 1700 K as **226.41 cal/g**.
- Equation 9, printed p. 4/PDF p. 13, is `H(T)−H(298.15 K)=0.1900 T−56.588` cal/g, implying approximately **266.41 cal/g** at 1700 K.
- The table's solid 1700 K value, **202.41 cal/g**, plus the stated **64.0 cal/g** fusion estimate also implies approximately **266.41 cal/g**. The next liquid row, 1800 K, is **285.41 cal/g**.

These exact conflicting entries were visually inspected. They are recorded as source evidence, without silently repairing the table or choosing an ingestion rule. This is not a complete numerical audit of the report's tables.

The Cezairliyan/Miiller publisher abstract supplies the distinct 304L fit `cp=1127−0.7265 T+2.884e-4 T^2 J/(kg K)`, valid over 1400–1700 K, with estimated inaccuracy 3%. The abstract also reports resistivity with estimated inaccuracy 2%; its displayed unit is preserved literally in the manifest pending a full-text unit check. These abstract results are useful evidence to pursue, but no full-page, chemistry or latent-heat verification was possible. Subsecond calorimetry is also not direct validation of a nanosecond heating path, superheating or a two-temperature regime.

Neither the approximately 1670–1730 K compilation range, the CINDAS 1672–1727 K range, nor a single 1707 K pulse-heating melting observation establishes the actual foil's equilibrium solidus, liquidus, latent heat and melt kinetics as a self-consistent set. The sources retain their own grade, reference temperature, material state and estimation method.

## Explicitly excluded substitution

[NBS SP 260-46](https://doi.org/10.6028/NBS.SP.260-46), Hust and Giarratano (1975), is a genuine thermal-conductivity/resistivity SRM publication for austenitic steel, SRMs 735/798, from 4 to 1200 K. It is **not SS304**. Table 1, printed p. 8/PDF p. 20, gives measured melts with approximately 19.90/20.48 wt% Ni and 16.41/15.98 wt% Cr, plus niobium and its own anneal. The attractive SRM uncertainty statements belong to that alloy and preparation. The source is retained to prevent a generic “austenitic stainless” search result from being substituted for grade 304. No 316/316L values from the other sources are used as 304 values either.

## Access, permission and reproduction

Official NIST/OSTI endpoints supplied valid PDFs for Flynn, Wilkes, Kim and SP 260-46. Flynn's PDF p. 3 carries a public-release approval dated October 9, 2015 over its historical restricted-report notice. Wilkes's PDF p. 1 and Kim's front matter state unlimited distribution; Wilkes also prints the government's retained license. These release notices establish the recorded access context, not a newly assigned Creative Commons license. Only citations, bounded factual summaries and a few numerical checks are stored in the repository; source text, figures and complete data tables are not redistributed.

The Ho/Chu scan carries an unclassified/unlimited-distribution report cover (PDF p. 2), but was retrieved from the named mirror: `apps.dtic.mil` failed with curl exit 56, no HTTP response/body. Its identity and title were verified in the scan; it is not represented as an official-host download or a newly open-licensed database.

The Springer `/content/pdf/` requests for Graves and Sweet returned HTTP 200 **HTML subscription previews**, not PDFs. Their OSTI `/servlets/purl/` attempts returned HTTP 404 HTML; available official metadata contains citation links without full-text links. Cezairliyan/Miiller likewise remains publisher-abstract-only after a focused NIST search. These failures do not imply the underlying measurements are invalid; they limit what was verified in this acquisition.

Representative acquisition commands and observed responses are below; the manifest contains exact source paths, bytes and hashes. A `.pdf` suffix was accepted only after file-header/PDF validation.

```sh
curl --max-time 50 -Ls -o /workspace/research_sources/ss304_properties/flynn_1968_nbs9869.pdf -w '%{http_code} %{content_type} %{size_download}\n' https://nvlpubs.nist.gov/nistpubs/Legacy/RPT/nbsreport9869.pdf
# 200 application/pdf 908832
curl --max-time 50 -Ls -o /workspace/research_sources/ss304_properties/wilkes_1997_ornl_cp96813.pdf -w '%{http_code} %{content_type} %{size_download}\n' https://www.osti.gov/servlets/purl/658177
# 200 application/pdf 950045
curl --max-time 45 -Ls -o /workspace/research_sources/ss304_properties/ho_chu_1977_mirror.pdf -w '%{http_code} %{content_type} %{size_download}\n' https://www.inductor-jmag.ru/files/content/a129160.pdf
# 200 application/pdf 2304408
pdftotext -layout /workspace/research_sources/ss304_properties/flynn_1968_nbs9869.pdf /workspace/research_sources/ss304_properties/flynn_1968_nbs9869.txt
```

The Wilkes scan has no usable embedded text; selected pages were rendered and read visually, with local OCR used only for navigation. Numerical applicability claims above were checked against the source page images or explicitly labeled as publisher-abstract evidence. No graph points were inferred or digitized.

## Requirements before a material model is selected

Record the dissertation foil's actual alloy/heat certificate, cold work and anneal, density/reference geometry, surface condition and any phase changes. Keep 304 and 304L records separate even when a compilation pools them. A future material table must distinguish measured, smoothed, estimated and extrapolated intervals and stop outside its authorized domain.

Use resistivity units consistently with electrical conductivity, and distinguish electrical resistivity from mass density. The acquired cp values are mass-specific, isobaric quantities; they are not automatically the volumetric constant-volume heat capacity used in a simplified heating equation or the internal-energy derivative required by an EOS. Enthalpy reference offsets, phase fractions and latent energy need a documented consistent treatment. These sources do not quantify high-field magnetoresistance, pressure/density dependence, or nanosecond non-equilibrium corrections, so no such effects or error bounds are inferred.
