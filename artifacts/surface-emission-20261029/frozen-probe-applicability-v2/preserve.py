from pathlib import Path
import hashlib, json, shutil, csv, os
assert os.sched_getaffinity(0)=={7} and os.getpriority(os.PRIO_PROCESS,0)>=15
src=Path('/home/treefor/projects/runtimes/pic-emission-applicability-v2-20261029')
out=Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-emission-applicability-v2-20261029')
backup=out/'source-backup-applicability-v2'
backup.mkdir(exist_ok=False)
for p in src.iterdir():
    if p.is_file():shutil.copy2(p,backup/p.name)
(backup/'evidence').mkdir()
for p in out.iterdir():
    if p.is_file():shutil.copy2(p,backup/'evidence'/p.name)
inputs={}
for name,root in [('source_v1',Path('/home/treefor/projects/runtimes/pic-thermal-field-emission-20261029/snapshots/usable-api-v1')),('guard_v1',Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-magnetized-emission-20261029/source-backup-guard-map-v1'))]:
    manifest=root/'SOURCE-SHA256.json'
    for rel,sha in json.loads(manifest.read_text()).items():assert hashlib.sha256((root/rel).read_bytes()).hexdigest()==sha
    inputs[name]={'root':str(root),'manifest_sha256':hashlib.sha256(manifest.read_bytes()).hexdigest(),'all_manifest_files_verified':True}
(backup/'input-seals.json').write_text(json.dumps(inputs,indent=2)+'\n')
# Preserve exact CSV record order in bounded, header-bearing publication shards.
shards=backup/'evidence'/'bound-probes';shards.mkdir()
with (out/'frozen-probes-with-historical-temperature.csv').open(newline='') as f:
    reader=csv.reader(f);header=next(reader);rows=list(reader)
for start in range(0,len(rows),64):
    with (shards/f'probes-{start}.csv').open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(header);writer.writerows(rows[start:start+64])
(backup/'PUBLICATION.txt').write_text('Complete isolated source and evidence, including all 87,760 native predicate calls and their diagnostic phase draws, are retained in this outside-project backup. Git publishes all source, exact bound-probe CSV shards, 1,056-row maps, checks, refinement records, contract and input hashes. The large, deterministic per-draw phase CSVs and executable are omitted from Git; source regenerates them and their hashes remain in SOURCE-SHA256.json. Original field/thermal binaries remain under their established immutable input seals. No production or GPU action.\n')
manifest={str(p.relative_to(backup)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(backup.rglob('*')) if p.is_file()}
(backup/'SOURCE-SHA256.json').write_text(json.dumps(manifest,indent=2)+'\n')
seal={'schema':'frozen-source-original-v1-applicability-v2-seal','backup':str(backup),'manifest_sha256':hashlib.sha256((backup/'SOURCE-SHA256.json').read_bytes()).hexdigest(),'checks':88816,'GPU_executed':False,'production_resumed':False,'native_physical_source_admitted':False}
(backup/'SEALED.json').write_text(json.dumps(seal,indent=2)+'\n');(out/'SEALED.json').write_text(json.dumps(seal,indent=2)+'\n')
print(json.dumps(seal))
