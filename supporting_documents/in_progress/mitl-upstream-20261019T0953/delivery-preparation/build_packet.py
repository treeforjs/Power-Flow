"""Assemble only the declared small evidence and exact public patch."""
import hashlib
import json
from pathlib import Path
import re
import shutil
import tempfile
import tarfile

from packet_common import WORK, REPO, BASE, OUT, PUBLIC, HELPERS, contract, clean_head, records, sha, write_new, prefix, exists

PROJECTIONS = (
    (bytes.fromhex('2f6f70742f636f6465782f72756e74696d65732f636f6465782d7072696d6172792d72756e74696d652f646570656e64656e636965732f707974686f6e'), b'[runtime Python directory]'),
    (bytes.fromhex('2f6574632f636f6465782f'), b'[runtime policy directory]/'),
)


def main():
    value = contract()
    clean_head(value)
    records(value)
    assert not exists(OUT) and not exists(OUT.with_suffix('.tar.gz'))
    stage = Path(tempfile.mkdtemp(prefix='.'+OUT.name+'-', dir=OUT.parent))
    temporary_archive = stage.with_suffix('.tar.gz')
    copied = []
    commands = None
    try:
        for name in PUBLIC:
            target = stage/'source'/name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(REPO/name, target)
        names = set(value['evidence_sha256'])|{'commands.jsonl', 'lease.json', 'status-before-lease.json',
            'record_command.py', 'delivery-preparation/contract.json'}|{'delivery-preparation/'+name for name in HELPERS}
        for name in sorted(names):
            source = WORK/name
            assert source.is_file() and not source.is_symlink()
            data = source.read_bytes()
            if name == 'commands.jsonl':
                commands = data
            raw_sha = hashlib.sha256(data).hexdigest()
            projections = []
            for original, replacement in PROJECTIONS:
                if original in data:
                    assert replacement not in data
                    projections.append({'original_hex': original.hex(), 'replacement': replacement.decode(), 'count': data.count(original)})
                    data = data.replace(original, replacement)
            target = stage/'evidence'/name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(data)
            copied.append({'source': str(source), 'packet_path': 'evidence/'+name,
                'original_sha256': raw_sha, 'packet_sha256': hashlib.sha256(data).hexdigest(), 'projections': projections})
        assert commands is not None
        from packet_common import git
        (stage/'npl-processability.patch').write_bytes(git('diff', '--binary', BASE, value['approved_head']))
        write_new(stage/'packet.json', {'schema': 'mitl-npl-processability-evidence-packet-v1',
            'contract': value, 'copied_inputs': copied, 'command_snapshot': prefix(commands, commands),
            'omitted_raw_inputs': value['omitted_raw_inputs'], 'raw_local_files_unchanged': True,
            'source_acquisition_only': True, 'scientific_replay_performed': False,
            'physical_parameters_adopted': 0, 'host_budget_accepted': False,
            'limits': 'Citation identity, access and alloy-coverage evidence only. Source PDFs/cache bodies remain external. No runtime model, actual-foil qualification or receiver acceptance follows.'})
        (stage/'systems-mitl-npl-processability-source-10192026.md').write_text(
            '# NPL processability report source assessment\n\n'
            + value['handoff_summary']+'\n\n'
            + 'Two public files and an exact patch are supplied over baseline '+BASE+'. '
            + 'Source/cache PDFs are omitted with hashes. Request outcomes and any corrections retain their original records. '
            + 'This packet adopts no physical parameters and executes no model or prior source acquisition. '
            + 'Receiver, actual-foil and physical-input dependencies remain pending.\n')
        forbidden = re.compile(r'\b(?:'+'|'.join(bytes.fromhex(s).decode() for s in ['6f70656e6169','63686174677074','636c61756465','636f646578'])+r')\b', re.I)
        manifest = {}
        for path in stage.rglob('*'):
            assert not path.is_symlink()
            if path.is_file():
                assert path.suffix.lower() not in {'.pdf', '.png', '.jpg', '.jpeg', '.tif', '.tiff', '.pyc'}
                assert not forbidden.search(path.read_text())
                manifest[str(path.relative_to(stage))] = sha(path)
        write_new(stage/'manifest.sha256.json', manifest)
        clean_head(value)
        records(value)
        for item in copied:
            if item['packet_path'] == 'evidence/commands.jsonl':
                prefix(commands, Path(item['source']).read_bytes())
            else:
                assert sha(Path(item['source'])) == item['original_sha256']
        with tarfile.open(temporary_archive, 'w:gz') as tar:
            tar.add(stage, arcname=OUT.name)
        clean_head(value)
        records(value)
        for item in copied:
            if item['packet_path'] == 'evidence/commands.jsonl':
                prefix(commands, Path(item['source']).read_bytes())
            else:
                assert sha(Path(item['source'])) == item['original_sha256']
        assert not exists(OUT) and not exists(OUT.with_suffix('.tar.gz'))
        stage.rename(OUT)
        temporary_archive.rename(OUT.with_suffix('.tar.gz'))
        archive = OUT.with_suffix('.tar.gz')
        write_new(WORK/'frozen_delivery.json', {'folder': str(OUT), 'archive': str(archive),
            'source_commit': value['approved_head'], 'base_commit': BASE, 'sha256': sha(archive),
            'bytes': archive.stat().st_size, 'manifest_entries': len(manifest),
            'manifest_sha256': sha(OUT/'manifest.sha256.json'), 'public_files_byte_identical': 2})
        print((WORK/'frozen_delivery.json').read_text())
    finally:
        if stage.exists():
            shutil.rmtree(stage)
        if temporary_archive.exists():
            temporary_archive.unlink()


if __name__ == '__main__':
    main()
