"""Read-only exact-base and patch check. Does not apply anything."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    here = Path(__file__).resolve().parent
    manifest = json.loads((here / 'base-sha256.json').read_text())
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('base', nargs='?', default=manifest['base_runtime'])
    root = Path(parser.parse_args().base)
    patch = here / 'swarm-completeness.patch'
    if sha(patch) != manifest['patch_sha256']:
        raise SystemExit('refused: adoption patch identity differs')
    for name, expected in manifest['base_sha256'].items():
        path = root / name
        actual = sha(path) if path.is_file() else None
        if actual != expected or (expected is None and path.exists()):
            raise SystemExit('refused: base path differs: ' + name)
    for name, expected in manifest['patched_sha256'].items():
        if sha(here.parent / 'source' / Path(name).name) != expected:
            raise SystemExit('refused: staged source differs: ' + name)
    subprocess.run(['git', 'apply', '--check', str(patch)], cwd=root, check=True)
    print('Exact base and staged patch verified; no files changed. Apply only in a new reviewed runtime.')


if __name__ == '__main__':
    main()
