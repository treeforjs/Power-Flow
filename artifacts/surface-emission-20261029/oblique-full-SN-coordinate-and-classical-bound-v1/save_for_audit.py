from pathlib import Path
import hashlib,json,shutil,os
assert os.sched_getaffinity(0)=={7} and os.getpriority(os.PRIO_PROCESS,0)>=15
src=Path('/home/treefor/projects/runtimes/pic-oblique-phase-map-20261029')
out=Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-oblique-phase-map-20261029')
backup=out/'source-backup-full-SN-coordinates-and-classical-bound-v1';backup.mkdir(exist_ok=False)
for p in src.iterdir():
    if p.is_file():shutil.copy2(p,backup/p.name)
(backup/'evidence').mkdir();(backup/'dependencies').mkdir()
for p in out.iterdir():
    if p.is_file():shutil.copy2(p,backup/'evidence'/p.name)
dep=Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-oblique-emission-20261029/source-backup-local-saddle-diagnostic-v1/saddle_modes.py')
shutil.copy2(dep,backup/'dependencies'/dep.name)
provenance={'previous_modes_file':str(dep),'previous_modes_sha256':hashlib.sha256(dep.read_bytes()).hexdigest(),'source_v1_manifest':'9bd0fa5771b8cba405cc0b63b99c0e32a24da8d735eef240606c36b1fb49f0ca','guard_map_manifest':'d4727a2f8904db8bf98e88000b30778b599f5dbee66c625715ab01f96b02322e','frozen_probe_map_commit':'22f32422dbefc2d35710a1d2f8afb3e44cfa75ad','bound_probes_sha256':hashlib.sha256(Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-emission-applicability-v2-20261029/frozen-probes-with-historical-temperature.csv').read_bytes()).hexdigest(),'scientific_status':'exact coordinate identity and conditional classical bound only; no deployable full oblique quantum source','tests_completed_before_parent_audit_only_switch':True,'new_full_coupled_test_run':False,'GPU_work':False,'production_action':False}
(backup/'PROVENANCE.json').write_text(json.dumps(provenance,indent=2)+'\n')
(backup/'PUBLICATION.txt').write_text('Outside-project backup retains all source, binary and full CSV evidence. Git publication may omit the regenerated executable and the 4224-row coordinate CSV, retaining their hashes and regeneration source. Rejected-window and classical bound tables remain diagnostic, not physical emission corrections.\n')
manifest={str(p.relative_to(backup)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(backup.rglob('*')) if p.is_file()}
(backup/'SOURCE-SHA256.json').write_text(json.dumps(manifest,indent=2)+'\n')
seal={'schema':'saved-oblique-full-SN-coordinate-and-conditional-classical-bound-v1','backup':str(backup),'manifest_sha256':hashlib.sha256((backup/'SOURCE-SHA256.json').read_bytes()).hexdigest(),'coordinate_kinematic_checks':236545,'rejected_window_partition_checks':3168,'classical_bound_crossing_checks':24826,'full_quantum_source_admitted':False,'new_full_coupled_test_run':False,'GPU_work':False,'production_action':False}
(backup/'SEALED.json').write_text(json.dumps(seal,indent=2)+'\n');(out/'SEALED.json').write_text(json.dumps(seal,indent=2)+'\n');print(json.dumps(seal))
