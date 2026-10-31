#!/usr/bin/python -I
"""Sizes and resident-path device memory of Hall axial meshes, including variants too large to construct here.

  size_table.py OUT.json [MESH_REPORT.json]

Counts follow axial::Mesh (axial_mesh.hpp): nz = axial + 2 air, nn = (n+1)^2 (nz+1), ne = 2 n (n+1) (nz+1) + (n+1)^2 nz,
nb = 2 (n+1) n nz + n^2 (nz+1), cells n^2 nz, metal cells 2 x (0.005/h) x (0.0002/h) x axial.
Device bytes: the allocation list of gpu::ResidentField::build() (resident_field.cpp), the same formula as
resident_bytes() in hall-scaling/field_scaling.cpp, which equals the measured own peak device bytes on every mesh run.
"""
import json, sys

MESHES = [
    ('200 um (pilot)', 45, 4, 1), ('100 um', 90, 4, 1), ('50 um', 180, 4, 1), ('25 um', 360, 4, 1),
    ('100 um, axial 8', 90, 8, 2), ('50 um, axial 16', 180, 16, 4), ('25 um, axial 8', 360, 8, 2), ('25 um, axial 16', 360, 16, 4),
    ('25 um, axial 32', 360, 32, 8), ('50 um isotropic-z (500 axial)', 180, 500, 20), ('25 um isotropic-z (1000 axial)', 360, 1000, 40),
]


def counts(n, a, air):
    nz = a + 2 * air
    nn = (n + 1) ** 2 * (nz + 1)
    ne = 2 * n * (n + 1) * (nz + 1) + (n + 1) ** 2 * nz
    nb = 2 * (n + 1) * n * nz + n * n * (nz + 1)
    nc = n * n * nz
    h = 0.009 / n
    metal = 2 * round(0.005 / h) * round(0.0002 / h) * a
    return nn, ne, nb, nc, metal


def resident_bytes(nn, ne, nb, nc, nm):
    r16 = lambda b: max(b, 16)
    t = 0
    t += r16(4 * nb * 4) + r16(4 * nb * 8) + r16(nb * 8) + r16((ne + 1) * 4) + r16(4 * nb * 4) + r16(4 * nb * 8) + r16(4 * nb * 8)
    t += r16(ne * 8) + r16(ne) + r16(ne * 8) + r16(ne * 8)
    t += r16(ne * 4) * 2 + r16(ne * 8)
    t += r16((ne + 1) * 4) + r16(12 * nc * 4) + r16(nc * 8) * 2 + r16(12 * nc * 4) + r16(nc * 4) + r16(nm * 4)
    t += r16((nn + 1) * 4) + r16(2 * ne * 4) * 2
    t += 2 * r16(nn * 8)
    t += r16(ne * 8) * 2 + r16(nn) * 2 + r16(ne * 8) + r16(nn * 8)
    t += 2 * (2 * r16(ne * 8) + 2 * r16(nb * 8) + 8 * r16(nn * 8))
    t += r16(nc * 8) + 2 * r16(ne * 8) + r16(nn * 8)
    t += 10 * r16(nn * 8) + 14 * r16(ne * 8) + 7 * r16(nb * 8) + 10 * r16(nn * 8) + 3 * r16(ne * 8)
    t += r16(nc * 8) + r16(nm * 8) + r16(ne * 8 * 15) + r16(nb * 8 * 1) + r16(nn * 8 * 11) + r16(nc * 8 * 5)
    t += r16(nn * 4) + r16(nc * 4) + 2 * r16(ne * 4) + r16(nb * 4) + 2 * r16(nn * 4) + r16(nc * 4) + 2 * r16(max(ne, nb) * 4)
    t += 32 * 256 * 16 + 3 * 32 * 16
    leaves = 6 * (ne * 8 // 2048 + 2) + 8 * (nn * 8 // 2048 + 2) + 16
    t += leaves * 32 * 3 + 4096 * 4 + 4096 * 2048 // 16
    t += 5 * r16(ne * 8) + r16(nb * 8) + 2 * 512 * 8 + 112 + 3 * 200 * 8 * 2 + 3 * 201 * 8
    return t


def main(argv):
    rows = []
    rep = {}
    if len(argv) > 1:
        for m in json.load(open(argv[1]))['meshes']:
            s = m['spec']; rep[(s['transverse'], s['axial'], s['air'])] = m
    for name, n, a, air in MESHES:
        nn, ne, nb, nc, nm = counts(n, a, air)
        if (n, a, air) in rep:   # constructed meshes: exact counts from the sealed constructor
            c = rep[(n, a, air)]['counts']
            assert (c['nodes'], c['edges'], c['faces'], c['cells'], c['metal_cells']) == (nn, ne, nb, nc, nm), name
        b = resident_bytes(nn, ne, nb, nc, nm)
        int_ok = 2 * n * (n + 1) * (a + 2 * air + 1) + (n + 1) ** 2 * (a + 2 * air) <= 2**31 - 1
        rows.append({'name': name, 'spec': [n, a, air], 'h_um': 9000 / n, 'dz_mm': 25 / a, 'nodes': nn, 'edges': ne, 'faces': nb, 'cells': nc,
                     'metal_cells': nm, 'resident_field_bytes': b, 'bytes_per_edge': b / ne, 'hall_int_index_ok': int_ok,
                     'hall_transverse_limit_ok': 9 <= n <= 500,
                     'fits_3080Ti_12GB_with_1.7GB_used_elsewhere': b + 400 * 2**20 < 12e9 - 1.7e9,
                     'fits_RTXPRO4000_24GB': b + 400 * 2**20 < 24e9 * 0.95,
                     'max_carriers_3080Ti_at_1400B': max(0, int((12e9 - 1.7e9 - 400 * 2**20 - b) / 1400)),
                     'max_carriers_RTXPRO4000_at_1400B': max(0, int((24e9 * 0.95 - 400 * 2**20 - b) / 1400)),
                     'constructed_and_checked': (n, a, air) in rep})
    json.dump({'schema': 'resolution-series-mesh-sizes-v1', 'rows': rows}, open(argv[0], 'w'), indent=1)
    print('%-34s %9s %11s %11s %9s %9s %9s %12s %12s' % ('mesh', 'nodes', 'edges', 'faces', 'cells', 'dev GB', 'B/edge', 'carriers12GB', 'carriers24GB'))
    for r in rows:
        print('%-34s %9d %11d %11d %9d %9.3f %9.0f %12d %12d' % (r['name'], r['nodes'], r['edges'], r['faces'], r['cells'], r['resident_field_bytes'] / 1e9,
                                                              r['bytes_per_edge'], r['max_carriers_3080Ti_at_1400B'], r['max_carriers_RTXPRO4000_at_1400B']))


if __name__ == '__main__':
    main(sys.argv[1:])
