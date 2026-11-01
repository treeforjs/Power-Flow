#!/usr/bin/env python3
"""Read-only source metadata extraction. Does not produce a chemistry table.

The raw sources remain in the separate local inventory runtime. This script
records file identities, column labels, finite table extents and counts only.
No cross-section values, interpolation, extrapolation or event law is exported.
"""
import hashlib
import json
import math
from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import ZipFile

HERE = Path(__file__).resolve().parent
REF = Path('/home/treefor/projects/runtimes/water-fragment-source-inventory-20261101/references')
CHEM = Path('/home/treefor/projects/runtimes/pic-water-chemistry-v3.2.2-auditprep-20261101')
NS = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def inspect_oxygen(path):
    summaries = []
    with ZipFile(path) as z:
        strings = [''.join(e.itertext()) for e in ET.fromstring(z.read('xl/sharedStrings.xml'))]
        names = [s.attrib['name'] for s in ET.fromstring(z.read('xl/workbook.xml')).find('m:sheets', NS)]
        relationships = ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
        by_id = {r.attrib['Id']: r.attrib['Target'] for r in relationships}
        sheet_ids = [s.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']
                     for s in ET.fromstring(z.read('xl/workbook.xml')).find('m:sheets', NS)]
        for name, rid in zip(names, sheet_ids):
            target = by_id[rid]
            root = ET.fromstring(z.read('xl/' + target))
            cells, numeric = {}, []
            for row in root.findall('m:sheetData/m:row', NS):
                vals = {}
                for c in row:
                    v = c.find('m:v', NS)
                    if v is None:
                        continue
                    assert c.find('m:f', NS) is None, 'Formula requires separate evaluation'
                    value = strings[int(v.text)] if c.get('t') == 's' else float(v.text)
                    cells[c.attrib['r']] = value
                    if not isinstance(value, str):
                        assert math.isfinite(value)
                        vals[c.attrib['r'].rstrip('0123456789')] = value
                if vals:
                    numeric.append(vals)
            item = {'sheet': name, 'labels': {k: v for k, v in cells.items() if isinstance(v, str)}}
            if name == 'EX':
                item['channels'] = []
                for e_col, s_col in [('A','B'), ('D','E'), ('G','H'), ('J','K'), ('M','N')]:
                    rows = [r for r in numeric if e_col in r and s_col in r]
                    item['channels'].append({
                        'energy_column': e_col, 'sigma_column': s_col,
                        'name': cells[s_col+'1'], 'rows': len(rows),
                        'energy_range_eV': [min(r[e_col] for r in rows), max(r[e_col] for r in rows)],
                        'first_positive_energy_eV': min(r[e_col] for r in rows if r[s_col] > 0),
                        'first_row_zero': rows[0][s_col] == 0,
                        'physical_threshold_verified': False})
            elif name == 'DCS':
                item.update({'rows': len(numeric), 'distinct_energies': len({r['A'] for r in numeric}),
                             'energy_range_eV': [min(r['A'] for r in numeric), max(r['A'] for r in numeric)],
                             'angle_range_deg': [min(r['B'] for r in numeric), max(r['B'] for r in numeric)],
                             'unit_label_requires_resolution': True})
            else:
                item.update({'rows': len(numeric), 'energy_range_eV':
                             [min(r['A'] for r in numeric), max(r['A'] for r in numeric)]})
            summaries.append(item)
    return summaries

def inspect_mccc(path):
    text = path.read_text()
    rows = [tuple(map(float, r.split())) for r in text.splitlines() if r and not r.startswith('#')]
    assert all(len(r) == 2 and all(math.isfinite(x) for x in r) for r in rows)
    assert all(b[0] > a[0] for a, b in zip(rows, rows[1:]))
    assert all(e >= 0 and s >= 0 for e, s in rows)
    return {'file': path.name, 'sha256': digest(path), 'rows': len(rows),
            'energy_range_eV': [rows[0][0], rows[-1][0]],
            'first_positive_energy_eV': next(e for e, s in rows if s > 0),
            'first_row_zero': rows[0][1] == 0,
            'cross_section_unit': 'a0^2', 'source_threshold_is_not_registry_threshold': True}

def main():
    manifest = CHEM / 'data/tables/chemistry-v3.2-manifest.json'
    m = json.loads(manifest.read_text())
    names = [g['name'] for g in m['groups']]
    missing = [p for p in ('e+OH', 'e+O', 'e+H2') if p not in names]
    assert len(missing) == 3
    species = {s['name']: s for s in m['species']}
    q = {}
    for label, products, reactants in [
            ('H2_to_2H', ['H', 'H'], ['H2']), ('H2_to_H2plus', ['H2+'], ['H2']),
            ('OH_to_O_H', ['O', 'H'], ['OH']), ('OH_to_OHplus', ['OH+'], ['OH']),
            ('O_to_Oplus', ['O+'], ['O'])]:
        q[label] = sum(species[s]['ref_eV'] for s in products) - sum(species[s]['ref_eV'] for s in reactants)
    base = {'runtime': str(CHEM), 'files': {}, 'absent_encounter_groups': missing,
            'existing_encounter_groups': names,
            'registry_species': [species[n] for n in ('e','H','H2','H2+','OH','OH+','O','O+')],
            'registry_reference_energy_differences_eV': q,
            'reference_energy_difference_is_not_finite_mass_lab_threshold': True,
            'native_v8_4_dependency_unchanged': '/home/treefor/projects/runtimes/pic-water-chemistry-v3.1.1-20261031',
            'native_v8_4_blob_sha256': '70fcc2adc2fdaa62abed7ed02d68ffc09b1b18554eda8fa99c5451120ae454b8'}
    for name in ['chemistry-v3.2.blob','chemistry-v3.2-spec.json','chemistry-v3.2-manifest.json','chemistry-v3.2-coverage-policy.json']:
        p = CHEM/'data/tables'/name
        base['files'][str(p)] = digest(p)
    assert base['files'][str(CHEM/'data/tables/chemistry-v3.2.blob')] == m['blob_sha256']
    summary = {'schema': 'water-fragment-primary-source-metadata-v1',
               'raw_local_directory': str(REF), 'raw_data_redistributed': False,
               'cross_sections_admitted': False, 'production_qualified': False,
               'network_complete': False, 'physical_model_validated': False,
               'oxygen_atomic_xlsx_sha256': digest(REF/'oxygen-atom-2026.xlsx'),
               'oxygen_sheets': inspect_oxygen(REF/'oxygen-atom-2026.xlsx'),
               'h2_files': [inspect_mccc(REF/n) for n in
                            ['MCCC-el-H2-TICS.X1Sg_vi=0.txt', 'MCCC-el-H2-b3Su_DE.X1Sg_vi=0.txt']]}
    for name, obj in [('BASE-PINS.json', base), ('SOURCE-TABLE-SUMMARY.json', summary)]:
        (HERE/name).write_text(json.dumps(obj, indent=2, sort_keys=True)+'\n')

if __name__ == '__main__':
    main()
