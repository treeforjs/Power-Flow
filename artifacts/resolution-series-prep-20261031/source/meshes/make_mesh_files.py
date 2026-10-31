#!/usr/bin/python -I
"""Write the mesh input files of the resolution series from the mesh report and the Hall admission runs.

  make_mesh_files.py OUTDIR MESH_REPORT.json ADMISSION.json...

Hall corrective v3.1 takes its mesh only as MeshSpec{transverse, axial, air, width_m} (field_contract.hpp) and constructs the
structured Yee complex itself (axial_mesh.hpp); a mesh "file" of the series is therefore that specification together with the
identities and sizes that the sealed constructor produces for it, so that any run can check it is on the intended mesh.
"""
import json, sys
from pathlib import Path

NAMES = {(45, 4, 1): 'mesh-200um', (90, 4, 1): 'mesh-100um', (180, 4, 1): 'mesh-50um', (360, 4, 1): 'mesh-25um',
         (90, 8, 2): 'mesh-100um-axial8', (180, 16, 4): 'mesh-50um-axial16'}


def main(argv):
    out = Path(argv[0]); rep = json.load(open(argv[1]))
    adm = {}
    for f in argv[2:]:
        d = json.load(open(f))
        adm[tuple(d['spec'][:3])] = d
    out.mkdir(parents=True, exist_ok=True)
    index = []
    for m in rep['meshes']:
        s = m['spec']; key = (s['transverse'], s['axial'], s['air'])
        a = adm.get(key, {})
        rec = {
            'schema': 'resolution-series-mesh-input-v1',
            'hall_mesh_spec': {'transverse': s['transverse'], 'axial': s['axial'], 'air': s['air'], 'width_m': s['width_m']},
            'cxx': 'hall_corrective_v3::MeshSpec{%d, %d, %d, .009}' % key,
            'geometry': {'foil_width_m': 0.005, 'foil_thickness_m': 0.0002, 'gap_m': 0.005, 'foil_height_m': 0.025, 'material': 'SS304',
                         'box_m': [0.009, 0.009, m['spacing']['zmax_m']['value'] - m['spacing']['zmin_m']['value']]},
            'spacing': m['spacing'], 'counts': m['counts'], 'cells_across': m['cells_across'], 'feature_positions': m['feature_positions'],
            'hall_mesh_json': m['mesh_json'], 'hall_mesh_json_sha256': m['mesh_json_sha256'],
            'encode_mesh_sha256': m['encode_mesh_sha256'], 'encode_mesh_tag': 'resolution-series-mesh-v1 + hall_corrective_v3::encode_mesh',
            'kinetic_geometry': m['kinetic_geometry'], 'kinetic_geometry_sha256': m['kinetic_geometry_sha256'],
            'device_memory_model': m['device_memory_model'],
            'hall_v31_admission': ({'admitted': a.get('admitted'), 'refusal': a.get('refusal'), 'configuration_sha256': a.get('configuration_sha256'),
                                    'motion_mesh_sha256': a.get('motion_mesh_sha256'), 'model_construct_s': a.get('model_construct_ms', 0) / 1e3,
                                    'dc_solves': [{k: x[k] for k in ('name', 'unknowns', 'cycle_iterations', 'converged')} for x in a['dc_profile']['solves']]}
                                   if a else {'admitted': None, 'note': 'not run'}),
        }
        name = NAMES.get(key, 'mesh-%dx%d+%d' % key)
        (out / (name + '.json')).write_text(json.dumps(rec, indent=1) + '\n')
        index.append({'file': name + '.json', 'spec': key, 'edges': m['counts']['edges'], 'admitted': rec['hall_v31_admission'].get('admitted')})
    (out / 'INDEX.json').write_text(json.dumps(index, indent=1) + '\n')
    print(json.dumps(index, indent=1))


if __name__ == '__main__':
    main(sys.argv[1:])
