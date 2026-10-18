# H₂ Lyman/Werner radiative line list

This directory retains three original ASCII files from the [Meudon MOLAT archive](https://molat.obspm.fr/index.php?page=pages/Molecules/H2/H2can94.php), attributed there to Abgrall, Roueff, Launay and Roncin, *Canadian Journal of Physics* **72**, 856–865 (1994). `metadata.json` records the exact source URLs, retrieval date, SHA-256 digests and original row counts.

`transitions.csv` contains **28,466** transitions with explicit upper/lower vibrational and rotational quantum numbers:

| Source | Upper state | Lines | Branch coverage |
| --- | --- | ---: | --- |
| `raw/BCan94` | B | 19,585 | Lyman P/R |
| `raw/CmCan94` | C− | 2,986 | Werner Q |
| `raw/CpCan94` | C+ | 5,895 | Werner P/R |

All lower levels belong to X. Source A values are in s⁻¹, and transition wavenumbers are in cm⁻¹. Vacuum wavelengths are derived as `1e7 / wavenumber_cm-1` in nm. Their extra displayed digits do not add precision to the source wavenumbers. C− and C+ preserve the archive's separate parity-component labels; match these to the collisional-data level basis before coupling rates. Maximum upper J is 25.

To regenerate the CSV from the retained raw files, or check it without changing files:

```sh
python tools/curate_molat_h2_lines.py
python tools/curate_molat_h2_lines.py --check
```

The parser verifies each raw checksum, declared row count, six-column format, unique level labels, nonnegative A, positive wavenumber and expected ΔJ selection. Both modes passed on all 28,466 transitions.

These are radiative bound–bound transitions, **not solved level populations**. Electron excitation, cascades from omitted states, dissociation continua, predissociation, ionization, quenching, water-breakup source populations and radiative transfer remain required. The [state-resolved data requirements](../../../docs/rovibrational_data_requirements.md) describe those gaps. No complete CR model or fit to the dissertation spectra is claimed.
