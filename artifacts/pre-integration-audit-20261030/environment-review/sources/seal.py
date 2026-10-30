from pathlib import Path
import hashlib
import json
import os

base = Path(__file__).resolve().parent
root = base / 'root'
manifest = {}
for directory, dirs, files in os.walk(root, followlinks=False):
    for name in sorted(dirs + files):
        path = Path(directory) / name
        relative = path.relative_to(root).as_posix()
        if path.is_symlink():
            manifest[relative] = {'symlink': os.readlink(path)}
        elif path.is_file():
            with path.open('rb') as stream:
                digest = hashlib.file_digest(stream, 'sha256').hexdigest()
            manifest[relative] = {'bytes': path.stat().st_size, 'sha256': digest}
destination = base / 'root-files.json'
destination.write_text(json.dumps(manifest, sort_keys=True, separators=(',', ':')) + '\n')
digest = hashlib.sha256(destination.read_bytes()).hexdigest()
(base / 'root-files.sha256').write_text(digest + '  root-files.json\n')
print(json.dumps({'entries': len(manifest), 'manifest_sha256': digest,
                  'regular_file_bytes': sum(row.get('bytes', 0) for row in manifest.values())}))
