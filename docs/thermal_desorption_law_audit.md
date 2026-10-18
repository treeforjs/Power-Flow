# Thermal desorption law and energy provenance audit

This audit compares the dissertation's printed equations with its cited primary literature. The remaining-inventory, first-order law in [`mhdlab/thermal_release.py`](../mhdlab/thermal_release.py) has support as a **kinetic model choice**; it is not a verbatim implementation of dissertation Eq. 2.30. None of the checked sources establishes a calibrated water release law and wall-energy cost for the dissertation's particular foils. No physical defaults or source code are changed by this audit.

The machine-readable source and claim record is [`literature_manifest.json`](../data/surface_models/literature_manifest.json). The dissertation PDF is cached outside the repository in `/workspace/research_sources/`; the literature PDFs, extracted text and inspection images are in its `surface_desorption/` subdirectory.

## Sources actually inspected

| Source | Access and page convention | Relevant evidence |
| --- | --- | --- |
| Trevor J. Smith, *Characterizing Desorbed Surface Contaminants Leading to Plasma Formation in Power Feeds of Pulsed Power Accelerators* (2025), [validated DeepBlue PDF](https://backend.production.deepblue-documents.lib.umich.edu/server/api/core/bitstreams/24e682a8-896a-410d-85bb-ca4654a9db02/content) | Validated 203-page dissertation. Printed pp. 17–18 are PDF pp. 42–43; bibliography p. 174 is PDF p. 199. | Eqs. 2.30–2.31 and references 28, 30, 32. |
| J. Matthew D. Lane, Kevin Leung, Aidan P. Thompson and Michael E. Cuneo (2018), *J. Phys.: Condens. Matter* **30**, 465002, [DOI 10.1088/1361-648X/aae4af](https://doi.org/10.1088/1361-648X/aae4af), dissertation ref. 28 | Full accepted manuscript from [OSTI 1595030](https://www.osti.gov/servlets/purl/1595030), 32 numbered pages. Citations below refer to these manuscript/PDF page numbers, not the journal's 15-page layout. Manuscript dated October 11, 2018; OSTI identifies it as an accepted manuscript, SAND2019-14427J. | Eq. 1 p. 3; Eq. 2 p. 4; first-order choice p. 22; parameter refits pp. 25–26. Equation and refit pages were visually inspected. |
| N. Bennett et al. (2023), *Electrode plasma formation and melt in z-pinch accelerators*, *Phys. Rev. Accel. Beams* **26**, 040401, [DOI 10.1103/PhysRevAccelBeams.26.040401](https://doi.org/10.1103/PhysRevAccelBeams.26.040401), dissertation ref. 30 | Full 12-page [APS PDF](https://journals.aps.org/prab/pdf/10.1103/PhysRevAccelBeams.26.040401), CC BY 4.0. Printed p. 040401-4 is PDF p. 4. | Eqs. 4–5 and their parameter/species definitions. |
| M. E. Cuneo et al. (1997), *Results of vacuum cleaning techniques on the performance of LiF field-threshold ion sources on extraction applied-B ion diodes at 1–10 TW*, *IEEE Trans. Plasma Sci.* **25**, 229–251, [DOI 10.1109/27.602495](https://doi.org/10.1109/27.602495), dissertation ref. 32 | Bibliographic identity confirmed with Crossref and the [official OSTI 522412 record](https://www.osti.gov/biblio/522412). **Full text unavailable in this retrieval.** | No actual paper pages or numerical surface parameters verified. |

The IOP `/pdf` endpoint returned HTML, not a PDF; the OSTI manuscript supplied the actual Lane full text. For Cuneo 1997, the IEEE landing page returned an empty HTTP 200 response and both attempted publisher PDF routes returned HTTP 418. The OSTI record provided metadata/abstract without a full-text link. This is an access limitation, not evidence that the 1997 paper lacks a relevant result. A separate 1998 review was located during retrieval but is not used to substitute for the cited 1997 paper.

## What the three accessible laws say

Use $N(t)$ for the remaining adsorbed areal population, $N_0=N(0)$, and $N_{\rm ref}$ for an explicitly defined reference areal density. The source symbols are retained in the table to expose their differences.

| Source | Printed rate law | Coverage-dependent rate energy and normalization |
| --- | --- | --- |
| Dissertation p. 18, Eqs. 2.30–2.31 | $\dot n=-\nu_{\rm th}n(0)\exp[-E_i/(k_BT)]$ | $E_i=E_d[1-\alpha n(t)/(f n_{\rm ml,0})]$, with $f n_{\rm ml,0}=n(0)$. |
| Lane manuscript pp. 3–4, Eqs. 1–2 | $\dot\theta=-k_0(T)\theta(t)^x\exp[-E'(\theta)/(k_BT(t))]$ | $E'(\theta)=E_d(1-\alpha\theta)$; $\theta$ is instantaneous coverage relative to a defined monolayer, not relative to each run's initial inventory. |
| Bennett p. 040401-4, Eqs. 4–5 | $\dot n=-\nu_{\rm th}n(t)\exp[-E'(n)/(k_BT(t))]$ | $E'(n)=E_d[1-\alpha n(t)/n_{\rm ML}]$; its hard-coded $n_{\rm ML}$ differs from its initial $n(0)=2n_{\rm ML}$. |

The thesis explicitly calls the prefactor population in Eq. 2.30 the **initial** population in the accompanying text. Its Eq. 2.33 writes $n(0,t)$ at a boundary, which does not resolve the Eq. 2.30 initial-versus-remaining ambiguity. Do not silently describe the first-order implementation as an exact transcription or label the thesis expression a confirmed typographical error.

At constant temperature and with a constant barrier, the thesis's printed initial-population prefactor gives linear depletion until an additional finite-inventory rule is imposed. The first-order law gives $N=N_0\exp(-kt)$ and tends to zero continuously. With a coverage-dependent barrier, the initial-prefactor expression still has coverage dependence through the exponential; calling the entire law simply zero order would obscure this distinction.

Lane lists possible exponents $x=0,1,2$, then explicitly sets $x=1$ throughout its comparisons on p. 22 because molecular water is treated as desorbing as a single entity. The authors note that recombination of separated surface H and OH could instead introduce second-order kinetics, which their nonreactive ClayFF water treatment cannot examine. Thus first-order molecular release is a scoped assumption, not a validation of dissociative water chemistry or bulk-hydrogen release.

## Parameters, surface definitions and validity

**Lane 2018.** The substrates are ideal crystalline hematite Fe₂O₃ facets, principally (0001) and hydroxylated C4 $(1\bar{1}02)$, with Cr₂O₃ comparisons from DFT. They represent possible oxide regions on steel; they are not bare SS304 or a measured oxide composition for the dissertation foils. No Cr₂O₃ classical force-field desorption trajectories are supplied in this work (pp. 2–3, 12, 18).

- On (0001), one monolayer is three water molecules per exposed surface metal ion, about 13 molecules/nm², or approximately $1.3\times10^{19}\,\mathrm{m}^{-2}$ (p. 5). The same areal-density convention is used on $(1\bar{1}02)$, despite a different site geometry; the authors explicitly discuss this normalization limitation on p. 25. The hydroxylated C4 surface retains strongly bound OH groups outside the modeled removable-water inventory (pp. 5, 16).
- The MD heating study uses initial removable-water coverages from 6% to 100% of that monolayer, substrate temperature ramps from 300 to 1300 K, and rates 20–2000 K/ns (durations 50–0.5 ns). Only the oxide substrate is thermostatted; water exchanges heat through the substrate (pp. 7, 18). GCMC at 300 K and $10^{-5}$ Torr predicts about 0.37 monolayers on Fe₂O₃(0001), whereas the C4 facet's removable inventory approaches zero (pp. 15–16, 28). These are particular equilibrium surface predictions, not recovered initial inventories for a shot.
- Eq. 2 supplies the Eyring relation $k_0(T)=\kappa k_BT/h$, with transmission coefficient $\kappa$ often set to one. It quotes $6.2\times10^{12}\,\mathrm{s}^{-1}$ at 300 K, then uses the conventional value $1.0\times10^{13}\,\mathrm{s}^{-1}$ **at 300 K** (p. 4). A constant caller-supplied attempt frequency is a separate implementation choice; do not claim this sentence alone specifies its temperature dependence at every temperature.
- The paper adopts $\alpha=0.5$, citing earlier use; it is not reported as a fit to these dissertation spectra or current traces (pp. 3, 22, 26). Finite-temperature correlated migration produces coverage effects on (0001), whereas the C4 facet is much less coverage dependent (pp. 24–25).
- The unmodified low-coverage ClayFF energies used as $E_d$ are 0.87 eV for (0001) and 0.69 eV for C4 $(1\bar{1}02)$. Fitting dilute single-water MD desorption times at **200 K/ns** instead gives effective $E_d=0.695$ and 0.420 eV, respectively (pp. 22, 26). The paper calls these “effective energies” and says they should depend on ramp rate. Those fitted values are not zero-temperature adsorption energies.
- With the refitted energies and adopted $\alpha=0.5$, agreement is within 10% for the compared (0001) trajectories starting at $\theta=0.33,0.67$ at 200 and 2000 K/ns. The C4 and full-coverage (0001) curves retain shape discrepancies; at 20 K/ns the early full-coverage shape is also wrong (p. 26). This percentage describes agreement with this MD benchmark, not an experimental uncertainty or a general parameter error bound.
- Bulk multilayer evaporation is outside the study. The discussion of a 0.47 eV liquid-water vaporization enthalpy at 273 K is a comparison, not a fitted wall-release parameter for steel (p. 12). The nonreactive force field does not model water dissociation/recombination; DFT energetics inform that approximation (pp. 11–12, 22).

**Bennett 2023.** The simulated neutral surface contaminant is explicitly **H₂**, released from local electrode surface cells, although the introduction discusses both water and hydrogen in real outgassing (pp. 040401-2–4). It uses $\nu_{\rm th}\sim10^{13}\,\mathrm{s}^{-1}$, $E_d=0.8$ eV and $\alpha=0.17$. The authors describe $E_d$ as consistent with MD calculations and $\alpha$ as retained from aggregate simulated-current/data comparisons, citing their ref. 31: Bennett et al., *Phys. Rev. Accel. Beams* **22**, 120401 (2019). That earlier fit was not independently audited here.

Its hard-coded $n_{\rm ML}=8\times10^{15}\,\mathrm{cm}^{-2}=8\times10^{19}\,\mathrm{m}^{-2}$ is a CHICAGO convention. Its initial $n(0)=2n_{\rm ML}$ is described as 16 conventional monolayers. Preserve that definition explicitly; neither its $n_{\rm ML}$ nor its fitted $\alpha$ is interchangeable with Lane's physical monolayer or the thesis's $N_0$ normalization. The simulated H₂ treatment does not establish an H₂O parameter set or state-resolved chemistry.

For the *same* linear barrier represented with a new reference density, algebra requires

$$
\frac{\alpha_{\rm new}}{N_{\rm ref,new}}=
\frac{\alpha_{\rm old}}{N_{\rm ref,old}}.
$$

Changing a denominator while retaining its numerical $\alpha$ changes the model. Rescaling it preserves only this barrier function; it does not repair a different kinetic population prefactor or make a fit transferable between species and surfaces.

**The thesis's approximately 1 eV statement.** Printed p. 18 assigns about 1 eV to H₂O on SS304 and jointly cites refs. 30 and 32. Bennett's accessible equations instead give 0.8 eV in an H₂ release model; the cited 1997 full text was not recovered. Lane supplies facet-dependent oxide energies of roughly this scale, including zero-temperature and effective kinetic quantities that differ, but it is not the direct citation for this particular sentence. The exact numerical origin and experimental meaning of **1 eV for H₂O/SS304 remain unverified**. It must not become a default binding heat, a fitted water barrier, or a stated uncertainty interval.

## Rate barriers and wall-energy accounting are separate inputs

Lane explicitly distinguishes the desorption barrier $E^*(\theta)$ from apparent binding energy $E'(\theta)$ on p. 4. It assumes equality for its water treatment because the potential-energy profile is taken to rise monotonically to a plateau without an additional saddle barrier. That is a stated approximation for those surfaces. The later refit of $E_d$, temperature-dependent configurations and suggested free-energy-barrier calculations (p. 26) prevent interpreting every fitted Arrhenius parameter as a measured surface internal-energy change.

A kinetic activation energy sets a rate. The energy needed to move material from a specified adsorbed state to a specified gas reference sets an energy transfer. A fitted rate can incorporate configurational sampling, entropy, the prefactor convention and non-equilibrium heating. An activation free-energy barrier is also not a surface free energy, and a surface free energy is not automatically its internal energy. Neither source derives a surface internal-energy function for the real foil's inventory, hydroxyl chemistry and oxide composition. Neither is a source for charging a fitted barrier again as an additional gas-state excitation or wall loss.

The independent reference implementation keeps these roles explicit:

- `activation_energy_ev` and `attempt_frequency_s` determine finite-inventory first-order release. Its constant barrier per reservoir omits the cited coverage dependence; multiple supplied reservoirs are not, by themselves, a fitted Temkin law.
- `desorption_energy_ev` is a separately supplied adsorbed-to-specified-gas-reference energy cost. `internal_energy_ev` is the emitted state's energy above that gas reference and is charged once.
- The reported $2k_BT$ translational cost per emitted particle follows the separately declared, thermally accommodating flux-weighted Maxwellian boundary. It is not established by an Arrhenius fit, and differs from assigning every emitted particle the single speed in dissertation Eq. 2.34.
- With prescribed wall temperatures, the wall debit is a diagnostic, not a solved wall-temperature or total-energy feedback. A coupled calculation must debit the same ledger from the wall and credit the appropriate gas components once.

The older [`mhdlab/desorption.py`](../mhdlab/desorption.py) has provisional water-labeled defaults (0.95 eV barrier, 0.15 eV coverage subtraction and $3\times10^{18}\,\mathrm{m}^{-2}$ initial inventory). This audit does not verify those numbers against any of these sources. Its normalization is remaining/initial inventory, and its barrier subtraction is not a directly cited Lane fit. The manufactured configuration [`surface_release_manufactured.json`](../configs/surface_release_manufactured.json) is explicitly a numerical verification case, not a water material model.

### Optional mathematical hypothesis, not a literature result

If independent thermodynamic evidence were to justify a **differential internal-energy binding cost** $b(N)=E_d(1-\alpha N/N_{\rm ref})$, measured in eV per particle, one could define an adsorbed energy per area relative to the empty surface and gas ground state as

$$
U_{\rm ads}(N)=-eE_d\left[N-\frac{\alpha N^2}{2N_{\rm ref}}\right],
\qquad \frac{dU_{\rm ads}}{dN}=-e b(N),
$$

where $e$ is the joule-per-eV conversion and $N$ is in m⁻². A release from $N_a$ to $N_b\leq N_a$ would then cost

$$
\Delta U_{\rm ads}=eE_d\left[(N_a-N_b)-
\frac{\alpha}{2N_{\rm ref}}(N_a^2-N_b^2)\right].
$$

This integral, rather than the initial or final barrier times the entire released inventory, would provide the finite-step binding debit **under that hypothesis**. It requires a defined nonnegative binding-cost range and any relevant temperature/internal-state dependence. The cited fitted kinetic laws do not establish this $U_{\rm ads}$; the expression is not implemented or populated with water defaults here.

## Consequences for the publication model

1. Label the implemented baseline as prescribed-temperature, remaining-inventory first-order release with caller-supplied reservoir parameters. Cite Lane/Bennett for the scoped molecular kinetic form and disclose the difference from thesis Eq. 2.30.
2. Before a physical water fit, record foil alloy, surface oxide/facets/conditioning, removable species and hydroxyl inventory, initial coverage/reference-density convention, temperature history, and whether the fitted law uses constant or Eyring-type prefactors. A current-loss fit alone cannot identify all of these or validate a spectrum-resolved surface composition.
3. Treat coverage-dependent rate parameters and a wall-energy function as separate provenance requirements. If they are equated, document the microscopic assumption and reference states and validate the wall/gas energy balance.
4. Keep the thesis's approximately 1 eV attribution unresolved until the relevant primary measurement or the cited 1997 pages are obtained. Do not use the oxide MD refits, the H₂ current-loss model, or the legacy provisional defaults as an experimentally calibrated water parameter set.

Access was checked on 2026-10-18. PDF hashes and source-specific evidence locations are recorded in the companion manifest; full source texts are not redistributed in this repository.
