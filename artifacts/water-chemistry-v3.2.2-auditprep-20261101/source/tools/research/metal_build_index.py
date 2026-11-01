#!/usr/bin/env python3
"""Write data/sources/METAL-SOURCES.json indexing metal (Cr, Fe, Ni) sources and tables.

Usage:
    metal_build_index.py <runtime_root>

For every data/sources/metal-* directory: lists each raw file with sha256,
size and modification time (UTC), reads retrieved_utc.txt / request_url.txt
when present, and checks saved Crossref JSON against the expected citation
(first-author family name, year, volume, first page or article number).
Indexes every table in data/extracted/metal-*.json. Downloaded files are read
as data only.
"""
import datetime
import glob
import hashlib
import json
import os
import re
import sys
import unicodedata

# fmt: off
SOURCES = [
 {"source_id": "metal-nist-asd-ionization-energies",
  "citation": "A. Kramida, Yu. Ralchenko, J. Reader and NIST ASD Team (2024), NIST Atomic Spectra Database (ver. 5.12), Ionization Energies data, https://physics.nist.gov/asd, accessed 2026-10-31",
  "doi": "10.18434/T4W30F", "doi_registry": "DataCite (Crossref returns 404 as expected; DataCite record saved)",
  "url": "https://physics.nist.gov/cgi-bin/ASD/ie.pl", "tag": "recommended", "access": "open; CSV and HTML output retained",
  "elements": ["Cr", "Fe", "Ni"], "content": "first ionization energies with uncertainties and ground/ionized levels",
  "extracted": ["data/extracted/metal-nist-ionization-energies.json"]},
 {"source_id": "metal-nist-isotopic-compositions",
  "citation": "J. S. Coursey, D. J. Schwab, J. J. Tsai and R. A. Dragoset, Atomic Weights and Isotopic Compositions with Relative Atomic Masses, NIST Physical Measurement Laboratory (landing page credits AME2012 masses, Meija et al. 2013 weights, Berglund and Wieser 2009 compositions), accessed 2026-10-31",
  "doi": None, "doi_registry": "no DOI found on the retrieved pages",
  "url": "https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl", "tag": "recommended", "access": "open; HTML output retained",
  "elements": ["Cr", "Fe", "Ni"], "content": "isotope relative atomic masses, isotopic mole fractions, standard atomic weights",
  "extracted": ["data/extracted/metal-nist-masses.json"]},
 {"source_id": "metal-nist-codata-constants",
  "citation": "P. J. Mohr, D. B. Newell, B. N. Taylor and E. Tiesinga, CODATA recommended values of the fundamental physical constants: 2022, Rev. Mod. Phys. 97, 025002 (2025); values from NIST allascii.txt (2022 CODATA adjustment)",
  "doi": "10.1103/RevModPhys.97.025002", "expect": {"family": "Mohr", "year": 2025, "volume": "97", "page": "025002"},
  "url": "https://physics.nist.gov/cuu/Constants/Table/allascii.txt", "tag": "recommended", "access": "open; constants table retained (paper not downloaded)",
  "elements": [], "content": "electron mass in u, atomic mass constant (used for ion masses)",
  "extracted": ["data/extracted/metal-nist-masses.json"]},
 {"source_id": "metal-lennon1988-jpcrd",
  "citation": "M. A. Lennon, K. L. Bell, H. B. Gilbody, J. G. Hughes, A. E. Kingston, M. J. Murray and F. J. Smith, Recommended data on the electron impact ionization of atoms and ions: fluorine to nickel, J. Phys. Chem. Ref. Data 17, 1285-1363 (1988)",
  "doi": "10.1063/1.555809", "expect": {"family": "Lennon", "year": 1988, "volume": "17", "page": "1285"},
  "url": "https://srd.nist.gov/jpcrdreprint/1.555809.pdf", "tag": "recommended",
  "access": "full text retained (NIST JPCRD reprint, scanned with OCR layer)",
  "elements": ["Cr", "Fe", "Ni"],
  "content": "Table 1 fit parameters (I, A, B_i) of recommended single-ionization cross sections for Cr I, Fe I, Ni I; recommendation for all three follows the McGuire (1977) Born calculation (calculated basis). Estimated errors 40 % (Cr I, Fe I), 200 % (Ni I).",
  "extracted": ["data/extracted/metal-lennon1988-ionization.json"]},
 {"source_id": "metal-mcguire1977-pra",
  "citation": "E. J. McGuire, Electron ionization cross sections in the Born approximation, Phys. Rev. A 16, 62-72 (1977)",
  "doi": "10.1103/PhysRevA.16.62", "expect": {"family": "McGuire", "year": 1977, "volume": "16", "page": "62"},
  "url": "https://doi.org/10.1103/PhysRevA.16.62", "tag": "calculated", "access": "closed; abstract level only, numbers not retained",
  "elements": ["Cr", "Fe", "Ni"],
  "content": "Born approximation (generalized oscillator strengths) ionization cross sections for 20 elements with 19 <= Z <= 54 (abstract). Basis of the Lennon et al. recommendations for Cr I, Fe I and (per Fig. 143) Ni I.",
  "extracted": []},
 {"source_id": "metal-freund1990-pra",
  "citation": "R. S. Freund, R. C. Wetzel, R. J. Shul and T. R. Hayes, Cross-section measurements for electron-impact ionization of atoms, Phys. Rev. A 41, 3575-3595 (1990)",
  "doi": "10.1103/PhysRevA.41.3575", "expect": {"family": "Freund", "year": 1990, "volume": "41", "page": "3575"},
  "url": "https://doi.org/10.1103/PhysRevA.41.3575", "tag": "measured", "access": "closed; abstract level only, numbers not retained",
  "elements": ["Fe"],
  "content": "Absolute single-ionization cross sections 0-200 eV, +-10 %, fast-neutral-beam method, for 16 atoms including Fe (abstract); also double/triple ionization for some atoms. Cr and Ni are not in the abstract's list. Wang et al. (2018) and Atoms 2022 cite these Fe data.",
  "extracted": []},
 {"source_id": "metal-shah1993-jpb",
  "citation": "M. B. Shah, P. McCallion, K. Okuno and H. B. Gilbody, Multiple ionization of iron by electron impact, J. Phys. B 26, 2393-2401 (1993)",
  "doi": "10.1088/0953-4075/26/15/023", "expect": {"family": "Shah", "year": 1993, "volume": "26", "page": "2393"},
  "url": "https://iopscience.iop.org/article/10.1088/0953-4075/26/15/023", "tag": "measured",
  "access": "closed (IOPscience paywall despite OpenAlex 'bronze' flag); landing page with abstract retained; numbers not retained",
  "elements": ["Fe"],
  "content": "Pulsed crossed-beam TOF; relative cross sections for Fe^n+ (n = 1-4) from ground-state Fe, 8-1250 eV, normalized to Freund et al. (1990) (abstract). Note: the published title lists four authors (incl. Okuno), not three.",
  "extracted": []},
 {"source_id": "metal-bartlett2004-adndt",
  "citation": "P. L. Bartlett and A. T. Stelbovics, Electron-impact ionization cross sections for elements Z=1 to Z=54, At. Data Nucl. Data Tables 86, 235-265 (2004)",
  "doi": "10.1016/j.adt.2003.11.006", "expect": {"family": "Bartlett", "year": 2004, "volume": "86", "page": "235"},
  "url": "https://doi.org/10.1016/j.adt.2003.11.006", "tag": "calculated", "access": "closed; title level only (no abstract in Crossref/OpenAlex), numbers not retained",
  "elements": ["Cr", "Fe", "Ni"],
  "content": "Tabulated calculated ionization cross sections for Z = 1-54 (title); Cr, Fe, Ni fall in this range, coverage not verified from the paper. Method presumably that of the companion PRA 2002 paper (Born approximation with Hartree-Fock-Slater orbitals); not verified.",
  "extracted": []},
 {"source_id": "metal-bartlett2002-pra",
  "citation": "P. L. Bartlett and A. T. Stelbovics, Calculation of electron-impact total-ionization cross sections, Phys. Rev. A 66, 012707 (2002)",
  "doi": "10.1103/PhysRevA.66.012707", "expect": {"family": "Bartlett", "year": 2002, "volume": "66", "page": "012707"},
  "url": "https://doi.org/10.1103/PhysRevA.66.012707", "tag": "calculated", "access": "closed; abstract level only",
  "elements": ["Fe"],
  "content": "Analytic Born-approximation ionization amplitude with Hartree-Fock-Slater orbitals and full orthogonalization; results include Fe (abstract).",
  "extracted": []},
 {"source_id": "metal-lotz1967-apjs",
  "citation": "W. Lotz, Electron-impact ionization cross-sections and ionization rate coefficients for atoms and ions, Astrophys. J. Suppl. 14, 207 (1967)",
  "doi": "10.1086/190154", "expect": {"family": "Lotz", "year": 1967, "volume": "14", "page": "207"},
  "url": "https://doi.org/10.1086/190154", "tag": "model_derived", "access": "abstract level (MPG copy behind bot check); not retained",
  "elements": [],
  "content": "Empirical formula and compilation of measured single-ionization data; rate coefficients for the H I - Na I - Na XI triangle (abstract). Does not provide Cr/Fe/Ni constants.",
  "extracted": []},
 {"source_id": "metal-lotz1968-zphys",
  "citation": "W. Lotz, Electron-impact ionization cross-sections and ionization rate coefficients for atoms and ions from hydrogen to calcium, Z. Phys. 216, 241-247 (1968)",
  "doi": "10.1007/BF01392963", "expect": {"family": "Lotz", "year": 1968, "volume": "216", "page": "241"},
  "url": "https://doi.org/10.1007/BF01392963", "tag": "model_derived", "access": "title level; not retained",
  "elements": [], "content": "Covers H to Ca only; not applicable to Cr, Fe, Ni (the Sc-Zn constants are in Lotz 1969).",
  "extracted": []},
 {"source_id": "metal-lotz1969-zphys",
  "citation": "W. Lotz, Electron-impact ionization cross-sections and ionization rate coefficients for atoms and ions from scandium to zinc, Z. Phys. 220, 466-472 (1969); retained text is the report version: W. Lotz, IPP 1/76, Institut fuer Plasmaphysik Garching, February 1968 (37 pp.)",
  "doi": "10.1007/BF01394789", "expect": {"family": "Lotz", "year": 1969, "volume": "220", "page": "466"},
  "url": "https://pure.mpg.de/rest/items/item_2133619_1/component/file_2133618/content", "tag": "model_derived",
  "access": "report version (public file in MPG.PuRe, linked by OpenAlex as submitted version) retained; scanned, no text layer",
  "elements": ["Cr", "Fe", "Ni"],
  "content": "Subshell constants q_i, a_i, b_i, c_i (Table 1) and binding energies P_1 (4s, Table 3), P_2 (3d, Table 4), P_3 (3p, Table 5) for Cr I, Fe I, Ni I. Constants are reasoned estimates (no Sc-Zn cross-section curve was known); stated error +40/-30 %.",
  "extracted": ["data/extracted/metal-lotz-ipp1-76-ionization.json"]},
 {"source_id": "metal-lotz1970-zphys",
  "citation": "W. Lotz, Electron-impact ionization cross-sections for atoms up to Z=108, Z. Phys. 232, 101-107 (1970)",
  "doi": "10.1007/BF01393132", "expect": {"family": "Lotz", "year": 1970, "volume": "232", "page": "101"},
  "url": "https://doi.org/10.1007/BF01393132", "tag": "model_derived", "access": "closed; title level only, numbers not retained",
  "elements": ["Cr", "Fe", "Ni"], "content": "Simplified Lotz formula with subshell binding energies up to Z = 108 (title). Constants not retained; the IPP 1/76 constants are used instead.",
  "extracted": []},
 {"source_id": "metal-dere2007-aa",
  "citation": "K. P. Dere, Ionization rate coefficients for the elements hydrogen through zinc, Astron. Astrophys. 466, 771-792 (2007)",
  "doi": "10.1051/0004-6361:20066728", "expect": {"family": "Dere", "year": 2007, "volume": "466", "page": "771"},
  "url": "https://www.aanda.org/articles/aa/pdf/2007/17/aa6728-06.pdf", "tag": "recommended",
  "access": "article blocked by publisher bot wall (HTTP 403); abstract level; numeric content obtained via CHIANTI files (metal-chianti-diparams)",
  "elements": ["Cr", "Fe", "Ni"],
  "content": "Direct-ionization cross-section fits (FAC calculations and laboratory data) for H-Zn; per-species basis for Cr I, Fe I, Ni I not confirmed.",
  "extracted": ["data/extracted/metal-chianti-dere2007-ionization.json"]},
 {"source_id": "metal-chianti-diparams",
  "citation": "CHIANTI atomic database version 11.0.2, files dbase/{cr,fe,ni}/{el}_1/{el}_1.diparams (Dere 2007, A&A 466, 771); format per ChiantiPy source (github.com/chianti-atomic/ChiantiPy, commit 43773711a592ef01eac100c1e08e515e167f073c: tools/util.py descale_bti, core/Ion.py diCross, tools/io.py diRead)",
  "doi": "10.1051/0004-6361:20066728", "crossref_file_dir": "metal-dere2007-aa", "expect": {"family": "Dere", "year": 2007, "volume": "466", "page": "771"},
  "url": "https://sohoftp.nascom.nasa.gov/solarsoft/packages/chianti/dbase/", "tag": "recommended",
  "access": "open; data files and ChiantiPy source (as .txt, not executed) retained",
  "elements": ["Cr", "Fe", "Ni"],
  "content": "Burgess-Tully-scaled spline knots (6 points) of direct single-ionization cross sections for Cr I, Fe I, Ni I; I = 6.767, 7.902, 7.637 eV.",
  "extracted": ["data/extracted/metal-chianti-dere2007-ionization.json"]},
 {"source_id": "metal-mattleubner2005-ijms",
  "citation": "H. Deutsch, P. Scheier, S. Matt-Leubner, K. Becker and T. D. Maerk, A detailed comparison of calculated and measured electron-impact ionization cross sections of atoms using the Deutsch-Maerk (DM) formalism, Int. J. Mass Spectrom. 243, 215-221 (2005); erratum Int. J. Mass Spectrom. 246, 113 (2005)",
  "doi": "10.1016/j.ijms.2005.03.003", "doi_extra": ["10.1016/j.ijms.2005.08.012"], "expect": {"family": "Deutsch", "year": 2005, "volume": "243", "page": "215"},
  "url": "https://doi.org/10.1016/j.ijms.2005.03.003", "tag": "calculated", "access": "closed; title level only, numbers not retained",
  "elements": ["Cr", "Fe"],
  "content": "DM-formalism ionization cross sections compared with measurements; secondary web listings (not verified from the paper) state Fe is among the compared atoms and Cr is calculated where no data exist. Crossref lists Deutsch as first author.",
  "extracted": []},
 {"source_id": "metal-margreiter1994-ijms",
  "citation": "D. Margreiter, H. Deutsch and T. D. Maerk, A semiclassical approach to the calculation of electron impact ionization cross-sections of atoms: from hydrogen to uranium, Int. J. Mass Spectrom. Ion Processes 139, 127-139 (1994)",
  "doi": "10.1016/0168-1176(94)90024-8", "doi_extra": ["10.1016/0168-1176(94)04078-8"], "expect": {"family": "Margreiter", "year": 1994, "volume": "139", "page": "127"},
  "url": "https://doi.org/10.1016/0168-1176(94)90024-8", "tag": "calculated", "access": "closed; title level only, numbers not retained",
  "elements": ["Cr", "Fe", "Ni"],
  "content": "DM calculations H-U. Cited by Atoms 10, 82 (2022) and arXiv physics/0407030 as the place where measured Ni single-ionization data attributed to 'Koparnski' are presented; that primary measurement was not identified. Two Crossref DOIs exist for the same article.",
  "extracted": []},
 {"source_id": "metal-nelson1976-dtic",
  "citation": "A. N. Nelson, Electron Impact Ionization Cross Sections of Gold, Chromium and Iron, Technical Report AFML-TR-75-198, Massachusetts Institute of Technology (cited as 1975; DTIC accession ADA024208, cover dated 1976)",
  "doi": None, "doi_registry": "no DOI (technical report)",
  "url": "https://apps.dtic.mil/sti/pdfs/ADA024208.pdf", "tag": "measured",
  "access": "not retrieved: DTIC returned HTTP 403 to this host; 403 response retained; numbers not retained",
  "elements": ["Cr", "Fe"],
  "content": "Relative and absolute ionization cross sections for Au, Cr, Fe species (search-engine snippet of the report; content not verified). Only located candidate measurement for Cr.",
  "extracted": []},
 {"source_id": "metal-pavlov1970-jetp",
  "citation": "S. I. Pavlov and G. I. Stotskii, Single and multiple ionization of lead atoms by electrons, Sov. Phys. JETP 31, 61 (1970) [Zh. Eksp. Teor. Fiz. 58, 108 (1970)]",
  "doi": None, "doi_registry": "no DOI (1970 JETP translation); JETP archive record retained",
  "url": "http://www.jetp.ras.ru/cgi-bin/e/index/e/31/1/p61?a=list", "tag": "measured", "access": "record page retained",
  "elements": [], "content": "Concerns lead only; not applicable to Cr or Ni. No Pavlov et al. Cr/Ni measurement was found.",
  "extracted": []},
 {"source_id": "metal-bautista2017-aa",
  "citation": "M. A. Bautista, K. Lind and M. Bergemann, Photoionization and electron impact excitation cross sections for Fe I, Astron. Astrophys. 606, A127 (2017)",
  "doi": "10.1051/0004-6361/201731232", "expect": {"family": "Bautista", "year": 2017, "volume": "606", "page": "A127"},
  "url": "https://www.aanda.org/articles/aa/pdf/2017/10/aa31232-17.pdf", "tag": "calculated",
  "access": "not retrieved: publisher PDF HTTP 403 (bot wall), CDS catalogue host timed out; abstract level only",
  "elements": ["Fe"], "content": "R-matrix electron-impact collision strengths among the lowest 46 LS terms of Fe I (secondary description); data reportedly at CDS.",
  "extracted": []},
 {"source_id": "metal-conroy2020-apj",
  "citation": "A. Conroy, C. P. Ballance, C. A. Ramsbottom and F. P. Keenan, Electron-impact excitation of Fe I, Astrophys. J. 902, 68 (2020)",
  "doi": "10.3847/1538-4357/abb3fd", "expect": {"family": "Conroy", "year": 2020, "volume": "902", "page": "68"},
  "url": "https://doi.org/10.3847/1538-4357/abb3fd", "tag": "calculated", "access": "not retrieved (IOPscience); abstract level only",
  "elements": ["Fe"], "content": "Dirac R-matrix (1000-level CC) fine-structure collision strengths and Maxwell-averaged effective collision strengths for Fe I, 1000-100000 K (abstract).",
  "extracted": []},
 {"source_id": "metal-wang2018-apj",
  "citation": "K. Wang, K. Bartschat and O. Zatsarinny, Electron scattering from neutral Fe and low-energy photodetachment of Fe-, Astrophys. J. 867, 63 (2018)",
  "doi": "10.3847/1538-4357/aaddfa", "expect": {"family": "Wang", "year": 2018, "volume": "867", "page": "63"},
  "url": "https://research-repository.griffith.edu.au/server/api/core/bitstreams/31b7f5e3-8918-473a-bdec-f31acb63ebb6/content", "tag": "calculated",
  "access": "accepted manuscript retained (Griffith repository); machine-readable tables (effective collision strengths) not retrieved",
  "elements": ["Fe"],
  "content": "BSR-221 close-coupling for e + Fe I: effective collision strengths for 24,531 LS transitions (MRT, Maxwell-averaged, 1e2-1e5 K); Fig. 5 shows angle-integrated elastic, elastic+excitation and grand total cross sections (figure only, not digitized).",
  "extracted": []},
 {"source_id": "metal-hassan2020-molphys",
  "citation": "M. M. Hassan, M. Haque, A. K. F. Haque et al., Relativistic study on the scattering of electrons and positrons from atomic iron at energies 1 eV - 10 keV, Mol. Phys. 119, e1849838 (2021)",
  "doi": "10.1080/00268976.2020.1849838", "expect": {"family": "Hassan", "year": 2021, "volume": "119", "page": "e1849838"},
  "url": "https://doi.org/10.1080/00268976.2020.1849838", "tag": "calculated", "access": "closed; abstract level only, numbers not retained",
  "elements": ["Fe"], "content": "Dirac partial-wave complex optical potential: elastic DCS, integral, momentum-transfer, viscosity, absorption and total cross sections for e + Fe, 1 eV-10 keV (abstract).",
  "extracted": []},
 {"source_id": "metal-jakubassa2022-atoms",
  "citation": "D. H. Jakubassa-Amundsen, M. M. Haque, A. K. F. Haque et al., Electron and positron scattering from precious metal atoms in the eV to MeV energy range, Atoms 10, 82 (2022)",
  "doi": "10.3390/atoms10030082", "expect": {"family": "Jakubassa-Amundsen", "year": 2022, "volume": "10", "page": "82"},
  "url": "https://mdpi-res.com/d_attachment/atoms/atoms-10-00082/article_deploy/atoms-10-00082.pdf", "tag": "calculated",
  "access": "open (CC BY); full text retained",
  "elements": ["Ni"],
  "content": "Dirac partial-wave optical-potential calculations for e + Ni (and Cu, Pd, Pt), 1 eV-300 MeV: elastic, momentum-transfer, viscosity, inelastic, total ionization cross sections. Integrated cross sections are given only in Figs. 31-35 (no tables; data availability 'not applicable'); not digitized. Fig. 35a compares Ni ionization with 'Koparnski' data via Margreiter et al.",
  "extracted": []},
 {"source_id": "metal-dolmatov2014-cr-elastic",
  "citation": "V. K. Dolmatov, M. Ya. Amusia and L. V. Chernysheva, Electron elastic scattering off a spin-polarized Cr atom, Phys. Rev. A 90, 032717 (2014); arXiv:1407.4493",
  "doi": "10.1103/PhysRevA.90.032717", "expect": {"family": "Dolmatov", "year": 2014, "volume": "90", "page": "032717"},
  "url": "https://arxiv.org/pdf/1407.4493v1", "tag": "calculated", "access": "arXiv preprint retained",
  "elements": ["Cr"], "content": "Spin-polarized HF and RPAE elastic phase shifts and spin-resolved elastic cross sections for e + Cr up to 15 eV; figures only, no momentum-transfer cross section table; not digitized.",
  "extracted": []},
 {"source_id": "metal-smirnov2000-qe",
  "citation": "Yu. M. Smirnov, Excitation of Fe I and Fe II in electron-atom collisions, Quantum Electron. 30, 1019-1024 (2000)",
  "doi": "10.1070/qe2000v030n11abeh001855", "expect": {"family": "Smirnov", "year": 2000, "volume": "30", "page": "1019"},
  "url": "https://doi.org/10.1070/QE2000v030n11ABEH001855", "tag": "measured", "access": "not retrieved (IOPscience); abstract level only",
  "elements": ["Fe"],
  "content": "Optical crossed-beam measurement of total excitation cross sections of the Fe I z5P3o, x5D4o, z5F5o levels and of Fe II odd levels (abstract). Level-resolved; not transcribed.",
  "extracted": []},
 {"source_id": "metal-smirnov2019-pre",
  "citation": "Yu. M. Smirnov, Excitation of singlet and even triplet levels of the iron atom by electron impact, Plasma Res. Express 1, 015004 (2019)",
  "doi": "10.1088/2516-1067/aadfd3", "expect": {"family": "Smirnov", "year": 2019, "volume": "1", "page": "015004"},
  "url": "https://doi.org/10.1088/2516-1067/aadfd3", "tag": "measured", "access": "not retrieved (IOPscience); abstract level only",
  "elements": ["Fe"],
  "content": "19 singlet and 44 even-triplet Fe I excitation cross sections at 50 eV (abstract). Intercombination, level-resolved; not transcribed.",
  "extracted": []},
 {"source_id": "metal-smirnov2001-ht",
  "citation": "Yu. M. Smirnov, Cross sections of excitation of 3F-levels of nickel atom by electron impact, High Temp. 39, 342-346 (2001)",
  "doi": "10.1023/a:1017538004490", "expect": {"family": "Smirnov", "year": 2001, "volume": "39", "page": "342"},
  "url": "https://doi.org/10.1023/A:1017538004490", "tag": "measured", "access": "closed; title level only",
  "elements": ["Ni"], "content": "Excitation cross sections of Ni I 3F levels (title). Level-resolved; not transcribed.",
  "extracted": []},
 {"source_id": "metal-smirnov2022-ijps",
  "citation": "Yu. M. Smirnov, Excitation of 3G levels of nickel atom by electron impact, Int. J. Phys. Sci. 17, 45-50 (2022)",
  "doi": "10.5897/ijps2022.5001", "expect": {"family": "Yu", "year": 2022, "volume": "17", "page": "45"},
  "url": "https://academicjournals.org/journal/IJPS/article-full-text-pdf/6FC6EBF69682.pdf", "tag": "measured",
  "access": "open (CC BY 4.0); full text retained",
  "elements": ["Ni"],
  "content": "14 excitation cross sections of Ni I 3G° levels at 50 eV and optical excitation functions 0-200 eV for five levels. High-lying levels, not total excitation; not transcribed. Crossref records the author family name as 'Yu' (registry metadata error); title, volume and pages match.",
  "extracted": []},
 {"source_id": "metal-nist-srd64-elastic",
  "citation": "A. Jablonski, F. Salvat and C. J. Powell, NIST Electron Elastic-Scattering Cross-Section Database, NIST SRD 64, Version 5.0 (data content 2002)",
  "doi": "10.18434/T4/1502642", "doi_registry": "DataCite (record saved). The SRD 64 page also prints 10.18434/T4T88K, which DataCite resolves to the XPS database (SRD 20); both DataCite records and the Crossref 404 are saved.",
  "url": "https://srdata.nist.gov/srd64/", "tag": "calculated",
  "access": "not retrieved: tables are produced only through an interactive Blazor session; landing and Fe table pages retained without data",
  "elements": ["Cr", "Fe", "Ni"], "content": "Dirac partial-wave elastic DCS, total elastic and transport (momentum-transfer) cross sections, 50 eV-20 keV (300 keV for DCS), Z = 1-96.",
  "extracted": []},
 {"source_id": "metal-lxcat-species-check",
  "citation": "LXCat open-access database (nl.lxcat.net), species list for all 18 electron scattering cross-section databases, accessed 2026-10-31",
  "doi": None, "doi_registry": "website", "url": "https://nl.lxcat.net/data/set_specB.php", "tag": "recommended",
  "access": "species-list pages retained", "elements": [],
  "content": "No Fe, Cr or Ni species in any LXCat electron database (metals present: Cu, Hg, Mg, Na).",
  "extracted": []},
 {"source_id": "metal-nist-srd107-check",
  "citation": "NIST Electron-Impact Cross Sections for Ionization and Excitation Database (SRD 107), atom index, accessed 2026-10-31",
  "doi": None, "doi_registry": "not checked", "url": "https://physics.nist.gov/PhysRefData/Ionization/atom_index.html", "tag": "calculated",
  "access": "atom index page retained", "elements": [],
  "content": "BEB/BE-scaled atom entries are H, He, Li, B, C, N, O, Al, Ga, In only; no Cr, Fe or Ni.",
  "extracted": []},
]

MISSING = [
 {"element": "Cr", "item": "measured ionization cross sections", "status": "none obtained", "detail": "Only candidate found: Nelson, AFML-TR-75-198 (DTIC ADA024208), blocked (HTTP 403). No Pavlov et al. Cr measurement exists in the located literature (Pavlov and Stotskii 1970 is lead). Lennon et al. (1988) state no experimental data for Cr I."},
 {"element": "Fe", "item": "measured ionization cross sections (Freund et al. 1990; Shah et al. 1993)", "status": "abstract level only", "detail": "Both closed; no tabulated numbers retained. These are the measured Fe single (0-200 eV) and multiple (8-1250 eV) ionization data."},
 {"element": "Ni", "item": "measured ionization cross sections", "status": "none obtained", "detail": "Data attributed to 'Koparnski' appear only via Margreiter et al. (1994) and Atoms 10, 82 (2022) Fig. 35a; primary publication not identified. Lennon et al. Ni I recommendation carries 200 % estimated error."},
 {"element": "Cr, Fe, Ni", "item": "Bartlett and Stelbovics (2004) ADNDT tables; Deutsch et al. (2005) DM values; Lotz (1970) Z<=108 constants; McGuire (1977) tables", "status": "not retained", "detail": "All closed access; only Crossref/OpenAlex metadata saved. McGuire (1977) enters indirectly through the Lennon et al. fits."},
 {"element": "Cr, Fe, Ni", "item": "NIST BEB (SRD 107) tables", "status": "do not exist", "detail": "SRD 107 has no Cr, Fe or Ni entries."},
 {"element": "Cr, Fe, Ni", "item": "LXCat cross-section sets", "status": "do not exist", "detail": "No Fe, Cr, Ni species in LXCat electron databases."},
 {"element": "Cr", "item": "excitation cross sections (total or resonance transitions)", "status": "none obtained", "detail": "No measured or calculated Cr I excitation cross-section tables were located. Hanne et al. (1993) superelastic polarized-electron work on Cr gives orientation parameters, not cross sections (not retrieved)."},
 {"element": "Fe", "item": "excitation cross sections", "status": "not tabulated in retrieved material", "detail": "Wang et al. (2018) BSR-221: effective collision strengths only (machine-readable tables not retrieved); total excitation shown in Fig. 5 only. Bautista et al. (2017) R-matrix: CDS unreachable, article blocked. Conroy et al. (2020): not retrieved, effective collision strengths. Smirnov measurements: level-resolved emission/excitation cross sections, not transcribed."},
 {"element": "Ni", "item": "excitation cross sections", "status": "not obtained in usable form", "detail": "Atoms 10, 82 (2022) gives an inelastic cross section only in figures. Smirnov (2001, 2022) level-resolved measurements (3F, 3G levels) are not total excitation; not transcribed."},
 {"element": "Cr", "item": "elastic and momentum-transfer cross sections", "status": "none tabulated", "detail": "Dolmatov et al. (2014) up to 15 eV, figures only, spin-resolved elastic (no momentum-transfer). NIST SRD 64 (>= 50 eV) not retrievable without interactive session."},
 {"element": "Fe", "item": "elastic and momentum-transfer cross sections", "status": "none tabulated", "detail": "Hassan et al. (2021) closed; Wang et al. (2018) elastic only in Fig. 5; NIST SRD 64 not retrieved."},
 {"element": "Ni", "item": "elastic and momentum-transfer cross sections", "status": "figures only", "detail": "Atoms 10, 82 (2022) Figs. 31-32 (elastic, momentum transfer), no tables; not digitized. NIST SRD 64 not retrieved."},
 {"element": "Cr, Fe, Ni", "item": "Dere (2007) per-species data basis", "status": "unconfirmed", "detail": "A&A article blocked (HTTP 403); CHIANTI file header says 'derived from fits to experimental and theoretical data'."},
 {"element": "Cr, Fe, Ni", "item": "figure digitization", "status": "not performed", "detail": "Tabulated ionization data exist (Lennon, Lotz, CHIANTI), so no ionization figures were digitized. For elastic/excitation only figures exist (Wang 2018 Fig. 5 for Fe; Atoms 2022 Figs. 31-35 for Ni; Dolmatov 2014 for Cr); left undigitized pending a decision."},
]
# fmt: on


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def safe(doi):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", doi)


def norm(s):
    s = unicodedata.normalize("NFKD", s or "")
    return "".join(c for c in s if not unicodedata.combining(c)).lower()


def crossref_check(sdir, doi, expect):
    path = os.path.join(sdir, f"crossref_{safe(doi)}.json")
    if not os.path.exists(path):
        return {"doi": doi, "crossref_file": None, "match": None, "detail": "no Crossref file"}
    try:
        msg = json.load(open(path, encoding="utf-8"))["message"]
    except (ValueError, KeyError, TypeError):
        return {"doi": doi, "crossref_file": os.path.basename(path), "match": False,
                "detail": "Crossref returned no record (DOI not registered at Crossref)"}
    fam = [norm(a.get("family", "")) for a in msg.get("author", [])]
    year = None
    for k in ("published-print", "issued", "published-online"):
        dp = (msg.get(k) or {}).get("date-parts")
        if dp and dp[0] and dp[0][0]:
            year = dp[0][0]
            break
    page = (msg.get("page") or msg.get("article-number") or "").split("-")[0]
    res = {
        "doi": doi,
        "crossref_file": os.path.basename(path),
        "crossref_title": (msg.get("title") or [""])[0],
        "crossref_container": (msg.get("container-title") or [""])[0],
        "crossref_volume": msg.get("volume"),
        "crossref_first_page_or_article": page,
        "crossref_year": year,
        "crossref_first_author": fam[0] if fam else None,
    }
    if expect:
        checks = {
            "first_author": fam[:1] == [norm(expect["family"])],
            "author_listed": norm(expect["family"]) in fam,
            "year": year == expect["year"],
            "volume": str(msg.get("volume")) == str(expect["volume"]),
            "page": page.lower() == str(expect["page"]).lower(),
        }
        res["checks"] = checks
        res["match"] = checks["author_listed"] and checks["volume"] and checks["page"] and checks["year"]
    return res


def utc_mtime(path):
    return datetime.datetime.fromtimestamp(os.path.getmtime(path), datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    root = os.path.abspath(sys.argv[1])
    sroot = os.path.join(root, "data/sources")
    known = {s["source_id"] for s in SOURCES}
    dirs = sorted(d for d in os.listdir(sroot) if d.startswith("metal-") and os.path.isdir(os.path.join(sroot, d)))
    unknown = sorted(set(dirs) - known)
    out_sources = []
    for s in SOURCES:
        sdir = os.path.join(sroot, s["source_id"])
        files = []
        for p in sorted(glob.glob(os.path.join(sdir, "*"))):
            files.append({"file": os.path.relpath(p, root), "bytes": os.path.getsize(p), "sha256": sha256(p),
                          "mtime_utc": utc_mtime(p)})
        rt = os.path.join(sdir, "retrieved_utc.txt")
        retrieved = open(rt).read().split() if os.path.exists(rt) else []
        cdir = os.path.join(sroot, s.get("crossref_file_dir", s["source_id"]))
        cr = []
        if s.get("doi"):
            cr.append(crossref_check(cdir, s["doi"], s.get("expect")))
        for d in s.get("doi_extra", []):
            cr.append(crossref_check(cdir, d, None))
        if s.get("crossref_file_dir"):
            for c in cr:
                if c.get("crossref_file"):
                    c["crossref_file"] = os.path.join("data/sources", s["crossref_file_dir"], c["crossref_file"])
        entry = {k: v for k, v in s.items() if k not in ("expect", "doi_extra", "crossref_file_dir")}
        entry["doi_additional"] = s.get("doi_extra", [])
        entry["crossref"] = cr
        entry["retrieval_utc"] = retrieved if retrieved else sorted({f["mtime_utc"] for f in files})
        entry["files"] = files
        out_sources.append(entry)
    tables = []
    for p in sorted(glob.glob(os.path.join(root, "data/extracted/metal-*.json"))):
        doc = json.load(open(p, encoding="utf-8"))
        for t in doc["tables"]:
            xs = [v for v in t["x"] if isinstance(v, (int, float))]
            tables.append({
                "table_id": t["table_id"], "file": os.path.relpath(p, root), "file_sha256": sha256(p),
                "source_id": t["source_id"], "process": t["process"], "quantity": t["quantity"], "tag": t["tag"],
                "x_name": t["x_name"], "n_points": len(t["x"]),
                "x_range": [min(xs), max(xs)] if xs and t["x_name"] == "energy_eV" else None,
                "threshold_eV": t["threshold_eV"], "locator": t["locator"],
            })
    index = {
        "description": "Index of electron-impact data sources for Cr, Fe, Ni vapour (SS304) gathered for the PIC/MCC "
                       "chemistry component; tags: measured | recommended | calculated | model_derived.",
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "generated_by": "tools/research/metal_build_index.py",
        "provenance_rules": "Numbers transcribed only from retained primary papers or authoritative compilations with "
                            "locators; model_derived tables are explicit evaluations of published formulas; no "
                            "smoothing, filling or figure digitization was performed.",
        "sources": out_sources,
        "tables": tables,
        "missing": MISSING,
        "unindexed_metal_dirs": unknown,
    }
    path = os.path.join(sroot, "METAL-SOURCES.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(index, f, indent=1, ensure_ascii=False)
        f.write("\n")
    bad = [c for s in out_sources for c in s["crossref"] if c.get("match") is False]
    print(f"wrote {os.path.relpath(path, root)}: {len(out_sources)} sources, {len(tables)} tables, "
          f"{len(MISSING)} missing items; unindexed dirs: {unknown}")
    for c in bad:
        print("CROSSREF MISMATCH/ABSENT:", c["doi"], c.get("checks"), c.get("detail"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
