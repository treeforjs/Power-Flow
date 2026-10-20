# Rovibrational resolution for all molecules in the water-plasma model

The publication requirement is to evolve explicit rovibrational populations for **every retained molecular species**, including molecular ions, anions and intermediates. This applies to chemistry and energy transfer even when a molecule does not emit in the recorded VUV band. H2 Lyman/Werner synthesis is one part of that model. The dissertation's §4.4 describes unresolved molecular populations and limits the original comparison to general agreement.

The [coverage manifest](../data/cr_model/molecular_state_requirements.yaml) classifies all 18 species in the current mechanism: 12 molecular species requiring rovibrational populations and six atomic species requiring atomic treatment. A new molecular product must receive its own level basis and rate coverage. This is a specification; there is no complete population solver or state-to-state rate library on this branch.

## Species-specific state bases

Use a canonical source-state ID with electronic state, energy, statistical weight, symmetry and interpretable quantum assignments. Approximate assignments can cease to be unique at high excitation; preserve their source definitions and confidence. Do not impose a diatomic `(v,J)` basis on water or the polyatomic ions. [ExoMol's quantum cases](https://www.exomol.com/data/quantum-cases/) documents distinct conventions for these geometries and coupling schemes.

| Species | Required type of assignment |
| --- | --- |
| H2O | Electronic state, three-mode `(v1,v2,v3)`, `J,Ka,Kc`, parity and nuclear-spin symmetry |
| OH | Electronic state, `v,J`, applicable spin-orbit component, parity/Lambda doubling |
| H2 | Electronic state, `v,J`, parity and nuclear-spin species; additional spin/coupling labels for triplet manifolds |
| O2 | Electronic state, `v,N,J`, spin component, parity and nuclear symmetry where applicable |
| H2O+ | Electronic state, three-mode assignment, `N,Ka,Kc,J`, parity and nuclear-spin symmetry; open-shell asymmetric rotor |
| OH+ | Electronic state, `v,N,J`, fine-structure component and parity as applicable |
| H2+ | Electronic state, `v,N,J`, parity and nuclear-spin symmetry |
| H3+ | Degenerate-mode assignment `(v1,v2,l2)`, `J,K`, source-specific combined labels and permutation/nuclear-spin symmetry |
| O2+ | Electronic state, `v,J`, applicable spin-orbit component and parity; molecular O2+, not atomic doubly charged O |
| H3O+ | Multimode vibrational assignment including degenerate-mode angular momenta, `J,K`, inversion parity and permutation/nuclear-spin symmetry |
| OH- | Bound electronic state, `v,J`, parity; detachment continuum and transient resonances treated separately |
| O2- | Bound electronic state, `v,J`, spin-orbit component, parity and nuclear symmetry; detachment/resonance treatment |

The exact labels depend on the electronic manifold and source coupling scheme. For example, [H2O+ spectroscopy](https://doi.org/10.1051/0004-6361/201015087) distinguishes rotational N from spin-coupled J; [H3+ assignments](https://doi.org/10.1093/mnras/stad050) and [H3O+ eXeL](https://doi.org/10.1093/mnras/staa2034) require symmetry-adapted polyatomic labels. The hydronium import must retain its original local-mode assignments as well as any normal-mode mapping.

For temporary dissociative intermediates, use state-dependent resonances/continua and lifetimes instead of fictitious stable vibrational levels. This matters for [electronically excited O2-](https://doi.org/10.1021/jp0357323). Bound-state resolution does not justify omitting dissociation or electron detachment.

## Data needed for every molecule

For each retained state basis, require level energies/weights, allowed radiative processes, electron and relevant heavy-particle state-to-state transfer, formation pumping, state-dependent chemical losses, wall events, and outgoing fragment/internal-state/electron-energy branching. Rates must connect the same isotopologue and level definitions. A list of emission wavelengths or an assumed rotational/vibrational temperature does not satisfy explicit population resolution.

Internal-state energies must enter the reaction energy ledger and the translational product distribution. Initial and surface-desorbed state populations need their own justified distributions; nascent fragments need state and recoil distributions from the producing reaction. Use the [energy and charged-transport requirements](energy_and_magnetized_transport.md) rather than resetting every species to the surface temperature.

Radiative starting sources include the retained H2 MOLAT subset below, [OH MYTHOS](https://www.exomol.com/data/molecules/OH/16O-1H/MYTHOS/), [OH+ MoLLIST](https://www.exomol.com/data/molecules/OH_p/16O-1H_p/MoLLIST/), and the H3+/H3O+ sources above. Their availability does not establish complete collision or chemical branching data. Match isotopes explicitly: the current ExoMol H2+ catalog's HD+ files do not establish homonuclear H2+ coverage. [Posen, Dalgarno and Peek](https://doi.org/10.1016/0092-640X(83)90017-7) is a primary homonuclear H2+ quadrupole-radiation starting source.

Some limited collisional data exist: [LAMDA water](https://home.strw.leidenuniv.nl/~moldata/H2O.html) includes rovibrational H2/electron thermal coefficients over specified level and temperature ranges; [LAMDA OH+](https://home.strw.leidenuniv.nl/~moldata/OH%2B.html) includes restricted rotational transfer. These are candidate inputs to inspect, not already-curated data in this repository. Thermal rate coefficients cannot be applied to an arbitrary magnetized nonthermal EEDF without validating that closure.

## H2 VUV radiative subset

Within H2, the minimum VUV system is:

- H₂ X¹Σg⁺(v″, J″): the ground-state distribution feeding excitation and receiving radiative decay.
- H₂ B¹Σu⁺(v′, J′): Lyman B → X transitions.
- H₂ C¹Πu⁺ and C¹Πu⁻(v′, J′): Werner C → X transitions, keeping the parity components distinct.
- Additional electronic levels needed for converged cascades; do not assume B/C constitute a closed system. Track radiative dissociation separately from bound–bound photon production.

The primary [Meudon MOLAT bound–bound archive](https://molat.obspm.fr/index.php?page=pages/Molecules/H2/H2can94.php) provides `BCan94`, `CpCan94` and `CmCan94` for B–X, C⁺–X and C⁻–X. Files label upper/lower v and J, Einstein A in s⁻¹ and transition wavenumber in cm⁻¹; λ[nm] = 10⁷ / wavenumber[cm⁻¹]. C⁻ carries Q transitions and C⁺ carries P/R transitions. The archive also provides B′/D transitions and documents coverage through upper J′ = 25. Its source is Abgrall, Roueff, Launay and Roncin, *Can. J. Phys.* **72**, 856–865 (1994). Preserve the source file, version/hash, level identity and all transitions before applying the instrument bandpass. A line list supplies radiative rates; it does not determine populations.

The three B/C files are now retained with hashes and normalized into **28,466** explicitly labeled radiative transitions in [`data/cr_model/molat_h2`](../data/cr_model/molat_h2/README.md). `tools/curate_molat_h2_lines.py --check` reproduces their numeric conversion and validates raw counts and branch selection. Collision rates and the other processes below remain pending.

The separate [MOLAT dissociation-continuum archive](https://molat.obspm.fr/index.php?page=pages/Molecules/H2/H2cont97q5.php), based on Abgrall et al., *ApJ* **481**, 557–566 (1997), supplies emission probability densities into X's dissociation continuum. These use s⁻¹/Hartree and require integration over the stated continuum-energy grid; treating them as discrete Einstein A values would be incorrect.

## Population and chemistry inputs that remain necessary

For each included level, document electron-impact excitation/de-excitation from X(v,J), transitions between excited states where relevant, dissociation/ionization, radiative branching, cascade feeding, heavy-particle quenching and wall losses. Initial H₂(v,J) and nascent H₂(v,J) created by water chemistry need their own distributions. Thermal translation, rotation, vibration and the electron distribution are distinct inputs. An assumed Boltzmann rotational distribution must be labeled as a closure, not a solved rovibrational population.

[Bergmayr, Wünderlich and Fantz (2023)](https://doi.org/10.1140/epjd/s10053-023-00706-5) distinguishes electronically resolved CR, vibrational CR and rovibrational corona models. It describes Lyman/Werner rovibrational corona models and identifies missing consistent collision data, particularly heavy-particle data. This is a useful benchmark path, but does not establish that a complete time-dependent rovibrational water-plasma model is available off the shelf. Include missing rates explicitly in a coverage matrix and test level truncation and cascade convergence. A documented corona reduction is permissible only after comparing its relaxation/collision times with the experimental history.

Evaluate electron rates from the adopted EEDF and energy-dependent cross sections. Maxwellian-averaged coefficients alone do not establish rates for arbitrary nonthermal electrons. Couple the resulting stiff population equations to transport and water chemistry, then integrate emissivity and absorption along the actual sightline, gate and instrument response. Width-integrated spectra cannot by themselves localize an emitting source across the foil width.

## Electron energy loss is not chemical or spectral branching

[Song et al. (2021)](https://doi.org/10.1063/5.0035315) is the updated water-collision review; the curated Tables 12/13 in `data/cross_sections/song2021/` are ionic yields only. Separate these from the [Budde et al. (2022) swarm-optimized water set](https://doi.org/10.1088/1361-6463/ac8da3). The latter includes an effective H₂O* energy-loss channel fitted using a Balmer-emission cross section scaled by 100. It must not be assigned a literal molecular population, exclusive product branch or photon yield. Its two-term/isotropic assumptions also need review for magnetized, nonlocal MITL electrons; [Budde et al. (2023)](https://doi.org/10.1088/1361-6463/accaf4) addresses anisotropic water scattering.

Finally, water can produce Lyα directly through dissociative excitation. [The primary electron-impact H₂O VUV study](https://ntrs.nasa.gov/citations/19850035943) measures this channel. A synthetic Lyα signal must therefore distinguish excitation of existing H, molecular breakup into excited H, recombination/cascades and radiative transfer before it is interpreted as a neutral H density.

For atomic H, separate 2s and 2p: Lyα emissivity uses the **2p** population. Metastable 2s has two-photon decay and field/collisional mixing; an n=2 atom is not an automatic Lyα photon. [NBS's electron-impact experiment](https://pmc.ncbi.nlm.nih.gov/articles/PMC6696584/) explicitly controlled 2s quenching when measuring 2p Lyα. Atomic H/O and their atomic ions retain electronic/fine-structure kinetics rather than molecular rovibrational labels.
