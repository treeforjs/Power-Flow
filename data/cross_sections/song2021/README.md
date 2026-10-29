# Recommended electron–water ionic yields from Song et al. (2021)

These are factual numeric transcriptions of Tables 12 and 13 in M.-Y. Song et al., *Cross Sections for Electron Collisions with H2O*, **J. Phys. Chem. Ref. Data 50**, 023103 (2021), [doi:10.1063/5.0035315](https://doi.org/10.1063/5.0035315). The [published PDF at UCL](https://discovery.ucl.ac.uk/id/eprint/10128988/1/5.0035315.pdf) was retrieved on 2026-10-18. Its SHA-256, source pages, uncertainties and per-file coverage are recorded in `metadata.json`.

The ten ionic-yield two-column CSV files described above contain **251 cross-section samples**: six positive-ion yields, total ionization, and three dissociative-attachment anion yields. Energies are in eV; cross sections are in m². The two `table*_original_units.csv` files preserve the printed values and omitted cells. Table 12 uses 10⁻¹⁶ cm² throughout, which converts to m² by multiplying the printed value by 10⁻²⁰. Table 13 uses 10⁻¹⁸, 10⁻¹⁹ and 10⁻²⁰ cm² for H⁻, O⁻ and OH⁻ respectively: the SI multipliers are 10⁻²², 10⁻²³ and 10⁻²⁴.

The PDF was inspected as layout text and rendered at twice its native resolution. Every SI value was then compared independently with the numeric cells extracted by their PDF coordinates; all 251 matched. Table 12 is PDF page 21 / article page 023103-20; Table 13 is PDF page 22 / article page 023103-21. `verification.json` records that check. To repeat it with PyMuPDF installed:

```sh
python tools/verify_song2021_tables.py --pdf /path/to/5.0035315.pdf
```

These tables are usable observables, **not an executable reaction network**. A partial ionic yield does not specify the complete neutral co-products, fragment electronic states, energy partition, or mutually exclusive breakup branch. Those fields remain pending in `metadata.json`; the files are deliberately separate from the production loader manifest. `O_doubly_charged` means atomic O²⁺, not molecular O₂⁺. Several positive-ion recommendations rely on D₂O measurements, as stated in the table caption.

Blank cells, absent rows and samples outside each file's range remain unknown. They have not been replaced by zeros, thresholds or extrapolated values. In particular, no H⁻ sample at 17.0 eV was inserted. Total ionization is comparison data: adding it to the partial yields as another reaction would double count ionization. Table 13 gives no OH⁻ uncertainty; an omitted uncertainty does not mean zero error.

Still required are state-resolved neutral dissociation (the separate Table10 addition below supplies one inclusive OH yield), dissociative excitation including direct production of H(2p), secondary-species reactions, recombination and complete product/state branching. These files alone do not establish full H₂O breakup or discharge-breakdown modeling.

Table10 is now separately transcribed in
`e_H2O_neutral_OH_yield_song2021.csv`: eleven measured neutral OH-production
cross sections over10–300eV, verified against the original PDF words and a
rendered table. Its caption gives40percent uncertainty; the adjacent text
says36percent. Both statements remain in `table10_neutral_OH_yield.json`.
The table does not resolve OH internal populations or recoil energies. Eq3
identifies OH+H+e as the most probable pathway; this does not establish a
complete exclusive reaction network. No values outside the table support,
threshold zeros, or Maxwellian-tail extrapolations are supplied. This
observable remains outside the production mechanism loader until its product
and energy interpretation is bound. A standalone Raikou single-channel
screening uses prescribed monoenergetic electron conditions and retains the
explicitly incomplete energy closure.
