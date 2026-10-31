#!/usr/bin/env python3
"""Write (or check) kernel-binding-v2.json for coverage scope v2.

The binding pins the bytes of the v2 high-energy kernel sources, the outer3
dataset and the sealed bridge/scheduler headers they use, plus the bound
Song-band kernel identity (v1 kernel manifest and binding file). Its
canonical_manifest_sha256 is the SHA-256 of the compact, key-sorted JSON of
these fields; the C++ HighEnergyBands constructor recomputes and checks it.
"""
import hashlib
import json
import sys
from pathlib import Path

R = Path('/home/treefor/projects/runtimes')
P = R / 'pic-water-parent-coverage-v2-20261030'
OUT = P / 'kernel-binding-v2.json'
SOURCES = [
    'pic-water-parent-coverage-v2-20261030/source/beb_rbeb.hpp',
    'pic-water-parent-coverage-v2-20261030/source/high_energy_bands.hpp',
    'pic-water-parent-coverage-v2-20261030/source/beb_sdcs_parent_kernel.hpp',
    'h2o-collision-data-20261029/data/collisions/h2o-gas-20261029-outer3.json',
    'pic-water-full-time-subflow-20261029/inputs/water_parent_bridge.hpp',
    'pic-water-full-time-subflow-20261029/inputs/collision_scheduler.hpp',
]
SONG_KERNEL = '02285c3a28fd7ab554418fc8216cdec143b8bd0869d20d3531de535fc8d9255a'
SONG_BINDING_FILE = '8267201576f7f1bc73c640e7117d67f2149438abb1b68a5b3b3148ce6bd9a123'
KERNEL = 'parent_ground_BEB_SDCS_binary_encounter_v2'


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build() -> str:
    v1_binding = R / 'pic-water-exact-count-subflow-20261029/kernel-binding.json'
    if sha(v1_binding) != SONG_BINDING_FILE:
        raise SystemExit('v1 kernel binding file changed')
    if json.loads(v1_binding.read_text())['canonical_manifest_sha256'] != SONG_KERNEL:
        raise SystemExit('v1 kernel manifest changed')
    manifest = {
        'schema': 'parent-coverage-v2-kernel-binding',
        'sources': {name: sha(R / name) for name in SOURCES},
        'song_band_kernel_manifest_sha256': SONG_KERNEL,
        'song_band_kernel_binding_file_sha256': SONG_BINDING_FILE,
        'high_energy_kernel': KERNEL,
    }
    canonical = json.dumps(manifest, sort_keys=True, separators=(',', ':'), ensure_ascii=False)
    doc = dict(manifest)
    doc['canonical_manifest_sha256'] = hashlib.sha256(canonical.encode()).hexdigest()
    doc['scope'] = ('coverage scope v2 high-energy parent kernel: outer3 BEB rate 1-10 keV, provisional RBEB '
                    'above 10 keV, BEB/RBEB SDCS energy sharing with binary-encounter angles')
    return json.dumps(doc, indent=2, sort_keys=True) + '\n'


def main() -> None:
    text = build()
    if len(sys.argv) > 1 and sys.argv[1] == '--check':
        if not OUT.exists() or OUT.read_text() != text:
            raise SystemExit('kernel-binding-v2.json differs from the current sources')
        print('kernel binding matches sources')
        return
    OUT.write_text(text)
    print(f'wrote {OUT}: {json.loads(text)["canonical_manifest_sha256"]}')


if __name__ == '__main__':
    main()
