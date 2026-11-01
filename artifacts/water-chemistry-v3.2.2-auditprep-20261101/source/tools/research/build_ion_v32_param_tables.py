#!/usr/bin/env python3
"""Write the v3.2 parameter tables transcribed from page images (Rudd 1992 Table V and Table I(b), Gobet 2004 Table II).

usage: build_ion_v32_param_tables.py
Values were read from the page images (Rudd et al. 1992 PDF page 24 rendered at 220 dpi, page 8 text layer;
Gobet et al. 2004 PDF page 5 rendered at 300 dpi) and checked a second time against the text layer (where the
text layer renders the minus sign of Gobet Table II as the digit 2, the image is authoritative).
"""
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, 'data/extracted')

RUDD_CIT = ('M. E. Rudd, Y.-K. Kim, D. H. Madison, T. J. Gay, Electron production in proton collisions with atoms and molecules: '
            'energy distributions, Rev. Mod. Phys. 64, 441 (1992)')
GOBET_CIT = ('F. Gobet, S. Eden, B. Coupier, J. Tabet, B. Farizon, M. Farizon, M. J. Gaillard, M. Carre, S. Ouaskit, T. D. Mark, '
             'P. Scheier, Ionization of water by (20-150)-keV protons: separation of direct-ionization and electron-capture processes, '
             'Phys. Rev. A 70, 062716 (2004)')

NAMES = ['A1', 'B1', 'C1', 'D1', 'E1', 'A2', 'B2', 'C2', 'D2', 'alpha']
tables = [
    {'table_id': 'ion-rudd1992-table5-sdcs-parameters-h2o', 'source_id': 'ion-rudd-1992-rmp64-proton-sdcs', 'citation': RUDD_CIT,
     'doi': '10.1103/RevModPhys.64.441',
     'locator': "Table V 'Parameters for fitting singly-differential cross sections to the Rudd model, Eqs. (41)-(48)', printed page 463 "
                "(PDF page 24), column H2O",
     'process': 'H+ + H2O -> H+ + e + H2O+ (singly differential cross section in the ejected-electron energy, outer shells)',
     'quantity': 'model_parameter', 'x_name': 'parameter', 'x': NAMES,
     'y': [0.97, 82.0, 0.40, -0.30, 0.38, 1.04, 17.3, 0.76, 0.04, 0.64], 'y_uncertainty_rel': None, 'tag': 'recommended',
     'extraction_method': 'transcribed from page image',
     'notes': 'Applies to the shells whose binding energy is at most twice that of the least tightly bound orbital (text of Sec. VI.D, '
              'printed page 465). Values read from the page image and checked a second time against the text layer.'},
    {'table_id': 'ion-rudd1992-table5-sdcs-parameters-inner-shells', 'source_id': 'ion-rudd-1992-rmp64-proton-sdcs', 'citation': RUDD_CIT,
     'doi': '10.1103/RevModPhys.64.441',
     'locator': "Table V, printed page 463 (PDF page 24), column 'Inner shells'",
     'process': 'proton impact ionization of inner shells (universal set for all targets)', 'quantity': 'model_parameter',
     'x_name': 'parameter', 'x': NAMES, 'y': [1.25, 0.50, 1.00, 1.00, 3.0, 1.10, 1.30, 1.00, 0.0, 0.66], 'y_uncertainty_rel': None,
     'tag': 'recommended', 'extraction_method': 'transcribed from page image',
     'notes': "Inner shells: binding energy exceeding twice the binding energy of the least tightly bound orbital (printed page 465). "
              "For H2O this is 2a1 (32.2 eV) and 1a1 (539.7 eV). The source states the inner-shell parameters are approximate. "
              "Values read from the page image and checked against the text layer; identical to the K-shell set of the Geant4 "
              "G4DNARuddIonisationModel.cc (v11.4.3, lines 726-735) retained in data/sources/ion-geant4-v11.4.3-dna-rudd-ionisation."},
    {'table_id': 'ion-rudd1992-table1b-h2o-orbital-binding-energies', 'source_id': 'ion-rudd-1992-rmp64-proton-sdcs', 'citation': RUDD_CIT,
     'doi': '10.1103/RevModPhys.64.441',
     'locator': "Table I(b) (continued), printed page 447 (PDF page 8), H2O block (MO, binding energy I in eV, kinetic energy U in eV)",
     'process': 'H2O molecular orbitals', 'quantity': 'binding_energy_eV', 'x_name': 'orbital', 'x': ['1b1', '3a1', '1b2', '2a1', '1a1'],
     'y': [12.61, 14.73, 18.55, 32.2, 539.7], 'y_uncertainty_rel': None, 'tag': 'recommended', 'extraction_method': 'transcribed from the text layer',
     'notes': 'Occupation 2 electrons per orbital. Asterisked values in the source (32.2, 539.7) are marked there as from a different '
              'reference set; taken as printed. Checked a second time against the text layer.'},
    {'table_id': 'ion-gobet2004-table2-partial-ionization-fit-parameters', 'source_id': 'ion-gobet-2004-pra70-h2o-direct-vs-capture',
     'citation': GOBET_CIT, 'doi': '10.1103/PhysRevA.70.062716',
     'locator': "Table II 'Values for the parameters in Eq. (1) derived by fitting the partial cross sections in Fig. 5', printed page "
                "062716-4 (PDF page 5); Eq. (1) on the same page",
     'process': 'H+ + H2O -> product ion + at least one free electron (sum of direct ionization and capture with ionization), per product ion',
     'quantity': 'model_parameter', 'x_name': 'product_ion', 'x': ['H2O+', 'OH+', 'O+', 'O2+', 'H+'],
     'y': [[0.74, 4.42, 0.74, 0.57], [1.29, -0.53, 0.52, 0.65], [-0.34, 1.25, 0.08, 0.63], [-0.01, 0.03, 0.01, 2.19], [-0.95, 4.33, 0.27, 0.76]],
     'y_columns': ['A', 'B', 'C', 'D'], 'y_uncertainty_rel': None, 'tag': 'measured', 'extraction_method': 'transcribed from page image',
     'fit_equation': 'sigma(E_p) = 4 pi a0^2 / (x / (A ln(1 + x) + B) + 1 / (C x^D)), x = (E_p / 1836) / 13.6, E_p the proton energy in eV (Eq. (1))',
     'validity_keV': [20.0, 150.0],
     'notes': 'Fits to the measured partial cross sections (Fig. 5) for 20-150 keV protons, normalized to the total sigma- of Rudd et al. '
              '1985 (Phys. Rev. A 31, 492). The text layer renders the minus signs as the digit 2 (e.g. 20.53 for -0.53); signs read from '
              'the page image, checked twice. Outside 20-150 keV the fits are not valid (negative A makes them turn over).'},
]
for t in tables:
    json.dump(t, open(os.path.join(OUT, t['table_id'] + '.json'), 'w'), indent=1)
    print('wrote', t['table_id'])
