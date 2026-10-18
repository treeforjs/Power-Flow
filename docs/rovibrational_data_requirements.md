# State-resolved hydrogen spectra in the full-water model

The publication model must evolve rovibrational populations rather than distribute a single molecular population over every observed line. The dissertation's §4.4 describes the latter approximation and limits it to general agreement. Its atomic n = 2 population must also be separated into 2s and 2p: Lyα emissivity uses the **2p** population. Metastable 2s has a two-photon decay channel and can be mixed or quenched into 2p by fields/collisions; it is not an automatic Lyα photon. [NBS's electron-impact experiment](https://pmc.ncbi.nlm.nih.gov/articles/PMC6696584/) explicitly controlled 2s quenching when measuring 2p Lyα.

## Explicit levels and radiative transitions

Represent each molecular level by `(electronic_state, v, J, parity)` with a documented ortho/para convention. The minimum VUV system is:

- H₂ X¹Σg⁺(v″, J″): the ground-state distribution feeding excitation and receiving radiative decay.
- H₂ B¹Σu⁺(v′, J′): Lyman B → X transitions.
- H₂ C¹Πu⁺ and C¹Πu⁻(v′, J′): Werner C → X transitions, keeping the parity components distinct.
- Additional electronic levels needed for converged cascades; do not assume B/C constitute a closed system. Track radiative dissociation separately from bound–bound photon production.

The primary [Meudon MOLAT bound–bound archive](https://molat.obspm.fr/index.php?page=pages/Molecules/H2/H2can94.php) provides `BCan94`, `CpCan94` and `CmCan94` for B–X, C⁺–X and C⁻–X. Files label upper/lower v and J, Einstein A in s⁻¹ and transition wavenumber in cm⁻¹; λ[nm] = 10⁷ / wavenumber[cm⁻¹]. C⁻ carries Q transitions and C⁺ carries P/R transitions. The archive also provides B′/D transitions and documents coverage through upper J′ = 25. Its source is Abgrall, Roueff, Launay and Roncin, *Can. J. Phys.* **72**, 856–865 (1994). Preserve the source file, version/hash, level identity and all transitions before applying the instrument bandpass. A line list supplies radiative rates; it does not determine populations.

The three B/C files are now retained with hashes and normalized into **28,466** explicitly labeled radiative transitions in [`data/cr_model/molat_h2`](../data/cr_model/molat_h2/README.md). `tools/curate_molat_h2_lines.py --check` reproduces their numeric conversion and validates raw counts and branch selection. Collision rates and the other processes below remain pending.

The three B/C raw files have now been retained, with checksums and a normalized [28,466-transition table](../data/cr_model/molat_h2/README.md). `python tools/curate_molat_h2_lines.py --check` verifies the source hashes, counts, level labels, selection rules, and derived wavelength conversion. This is the radiative bound–bound subset; it supplies no collision rates or solved populations.

The separate [MOLAT dissociation-continuum archive](https://molat.obspm.fr/index.php?page=pages/Molecules/H2/H2cont97q5.php), based on Abgrall et al., *ApJ* **481**, 557–566 (1997), supplies emission probability densities into X's dissociation continuum. These use s⁻¹/Hartree and require integration over the stated continuum-energy grid; treating them as discrete Einstein A values would be incorrect.

## Population and chemistry inputs that remain necessary

For each included level, document electron-impact excitation/de-excitation from X(v,J), transitions between excited states where relevant, dissociation/ionization, radiative branching, cascade feeding, heavy-particle quenching and wall losses. Initial H₂(v,J) and nascent H₂(v,J) created by water chemistry need their own distributions. Thermal translation, rotation, vibration and the electron distribution are distinct inputs. An assumed Boltzmann rotational distribution must be labeled as a closure, not a solved rovibrational population.

[Bergmayr, Wünderlich and Fantz (2023)](https://doi.org/10.1140/epjd/s10053-023-00706-5) distinguishes electronically resolved CR, vibrational CR and rovibrational corona models. It describes Lyman/Werner rovibrational corona models and identifies missing consistent collision data, particularly heavy-particle data. This is a useful benchmark path, but does not establish that a complete time-dependent rovibrational water-plasma model is available off the shelf. Include missing rates explicitly in a coverage matrix and test level truncation and cascade convergence. A documented corona reduction is permissible only after comparing its relaxation/collision times with the experimental history.

Evaluate electron rates from the adopted EEDF and energy-dependent cross sections. Maxwellian-averaged coefficients alone do not establish rates for arbitrary nonthermal electrons. Couple the resulting stiff population equations to transport and water chemistry, then integrate emissivity and absorption along the actual sightline, gate and instrument response. Width-integrated spectra cannot by themselves localize an emitting source across the foil width.

## Electron energy loss is not chemical or spectral branching

[Song et al. (2021)](https://doi.org/10.1063/5.0035315) is the updated water-collision review; the curated Tables 12/13 in `data/cross_sections/song2021/` are ionic yields only. Separate these from the [Budde et al. (2022) swarm-optimized water set](https://doi.org/10.1088/1361-6463/ac8da3). The latter includes an effective H₂O* energy-loss channel fitted using a Balmer-emission cross section scaled by 100. It must not be assigned a literal molecular population, exclusive product branch or photon yield. Its two-term/isotropic assumptions also need review for magnetized, nonlocal MITL electrons; [Budde et al. (2023)](https://doi.org/10.1088/1361-6463/accaf4) addresses anisotropic water scattering.

Finally, water can produce Lyα directly through dissociative excitation. [The primary electron-impact H₂O VUV study](https://ntrs.nasa.gov/citations/19850035943) measures this channel. A synthetic Lyα signal must therefore distinguish excitation of existing H, molecular breakup into excited H, recombination/cascades and radiative transfer before it is interpreted as a neutral H density.
