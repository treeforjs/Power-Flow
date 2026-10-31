#!/usr/bin/python -I
"""Write the graded-box mesh definitions (transverse node coordinates and cell sizes) of the series.

  make_graded_files.py OUTDIR GRADED_SIZES.json [GRADED_DC.json ...]

The axis is graded::axis (graded/graded_mesh.hpp) in the same binary64 operations: a uniform core of round(2 core / h) cells
of size (2 core)/ncells with node coordinates i dx - core, then cells h q, h q^2, ... capped at hmax up to the box half-width
(the last cell ends on the box; a remainder smaller than half a grown cell is merged into the last cell). The transverse cell
count is checked against the C++ construction (graded-sizes). Axial: 4 foil cells of 6.25 mm and 1 air cell per end, as Hall.
Hall v3.1 cannot run these meshes (uniform MeshSpec only); they define the graded box for a Hall version with per-index
coordinates and record the measured DC port inductance where it was computed.
"""
import json, sys
from pathlib import Path


def axis(h, core, box, q, hmax):
    ncells = round(2 * core / h)
    assert abs(ncells * h - 2 * core) < 1e-12 and box >= core
    dx = (2 * core) / ncells
    grow = []; reach = core; step = dx
    while reach < box - 1e-12:
        step = min(step * q, hmax)
        nxt = reach + step
        if box - nxt < .5 * step:
            grow.append(box - reach); reach = box; break
        grow.append(step); reach = nxt
    nl = len(grow)
    hcell = list(reversed(grow)) + [dx] * ncells + grow
    x = [0.0] * (len(hcell) + 1)
    for i in range(ncells + 1):
        x[nl + i] = i * dx - core
    for k in range(nl):
        x[nl - 1 - k] = x[nl - k] - grow[k]
        x[nl + ncells + 1 + k] = x[nl + ncells + k] + grow[k]
    if nl:
        x[0] = -box; x[-1] = box
    return x, hcell


def main(argv):
    out = Path(argv[0]); sizes = json.load(open(argv[1]))['rows']
    dcs = [json.load(open(f)) for f in argv[2:]]
    out.mkdir(parents=True, exist_ok=True)
    index = []
    for r in sizes:
        if r['box_half_m'] < .01:
            continue
        x, hc = axis(r['h_m'], .0045, r['box_half_m'], r['q'], r['hmax_m'])
        assert len(hc) == r['nx'], (r, len(hc))
        name = 'graded-%dum-box%dmm-q%s' % (round(r['h_m'] * 1e6), round(2e3 * r['box_half_m']), ('%g' % r['q']).replace('.', 'p'))
        L = [d['dc_port_inductance_H'] for d in dcs if abs(d['h_m'] - r['h_m']) < 1e-12 and abs(d['box_half_m'] - r['box_half_m']) < 1e-12
             and abs(d['q'] - r['q']) < 1e-12 and abs(d['hmax_m'] - r['hmax_m']) < 1e-12]
        rec = {'schema': 'resolution-series-graded-mesh-v1', 'core_cell_m': r['h_m'], 'core_box_m': 0.009, 'box_m': 2 * r['box_half_m'],
               'growth_ratio': r['q'], 'max_cell_m': r['hmax_m'], 'axial': {'foil_cells': 4, 'air_cells_per_end': 1, 'dz_m': 0.00625},
               'transverse_node_coordinates_m': x, 'transverse_cell_sizes_m': hc,
               'counts': {k: r[k] for k in ('nx', 'nz', 'nodes', 'edges', 'faces', 'cells', 'metal_cells')},
               'largest_cell_m': r['largest_cell_m'], 'largest_adjacent_ratio': r['largest_adjacent_ratio'],
               'resident_field_bytes': r['resident_field_bytes'], 'pinned_host_bytes': r['pinned_host_bytes'],
               'dc_port_inductance_H': L[0] if L else None,
               'hall_v31': 'not admissible (uniform MeshSpec only); a Hall version with per-index transverse coordinates is required'}
        (out / (name + '.json')).write_text(json.dumps(rec, indent=1) + '\n')
        index.append({'file': name + '.json', 'h_um': r['h_m'] * 1e6, 'box_mm': 2e3 * r['box_half_m'], 'q': r['q'], 'n': r['nx'],
                      'edges': r['edges'], 'resident_GB': r['resident_field_bytes'] / 1e9, 'L_nH': L[0] * 1e9 if L else None})
    (out / 'GRADED-INDEX.json').write_text(json.dumps(index, indent=1) + '\n')
    for i in index:
        print(i)


if __name__ == '__main__':
    main(sys.argv[1:])
