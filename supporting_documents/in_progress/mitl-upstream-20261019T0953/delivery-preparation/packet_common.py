"""Small checksum and archive contract for a bounded source-acquisition packet."""
import datetime
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import subprocess
import tarfile

WORK = Path(__file__).resolve().parents[1]
REPO = Path('/workspace/Power-Flow')
BASE = '5019b86aabf0472037b0fa3262922c3d0dcd2e95'
OUT = Path('/workspace/handoffs/mitl-npl-processability-20261019T0953')
PUBLIC = ('data/materials/npl_processability_1996_source_manifest.json',
          'docs/npl_processability_1996_source_assessment.md')
HELPERS = ('packet_common.py', 'build_packet.py', 'verify_packet.py', 'complete_pass.py')


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def pairs(items):
    result = {}
    for key, value in items:
        if key in result:
            raise ValueError('duplicate JSON key')
        result[key] = value
    return result


def decode(data):
    def finite(value):
        result = float(value)
        if not math.isfinite(result):
            raise ValueError('nonfinite JSON number')
        return result
    return json.loads(data, object_pairs_hook=pairs,
                      parse_float=finite,
                      parse_constant=lambda value: (_ for _ in ()).throw(ValueError(value)))


def load(path):
    return decode(Path(path).read_text())


def encoded(value):
    return (json.dumps(value, indent=2, allow_nan=False)+'\n').encode()


def write_new(path, value):
    with Path(path).open('xb') as stream:
        stream.write(encoded(value))


def exists(path):
    return os.path.lexists(path)


def regular(path):
    path = Path(path)
    assert path.is_absolute()
    assert all(not parent.is_symlink() for parent in [path, *path.parents])
    assert path.is_file()
    return path


def relative(name):
    path = PurePosixPath(name)
    assert type(name) is str and not path.is_absolute() and path.parts
    assert '..' not in path.parts and '.' not in path.parts and path.as_posix() == name
    return name


def git(*args):
    return subprocess.check_output(['git', *args], cwd=REPO)


def contract():
    value = load(WORK/'delivery-preparation/contract.json')
    assert value['schema'] == 'mitl-npl-processability-packet-contract-v1'
    assert value['ready_for_root_execution'] is True and value['pending'] == []
    assert value['baseline_commit'] == BASE and value['approved_head'] != BASE
    assert set(value['public_sha256']) == set(PUBLIC)
    assert set(value['helper_sha256']) == set(HELPERS)
    assert value['evidence_sha256'] and value['record_requirements']
    assert isinstance(value['handoff_summary'], str) and value['handoff_summary'].strip()
    assert isinstance(value['completion_milestone'], str) and value['completion_milestone'].strip()
    assert value['physical_parameters_adopted'] == 0
    assert value['cache_PDFs_distributed'] is False and value['scientific_replay_performed'] is False
    assert value['next_task_selected'] is None
    return value


def clean_head(value):
    assert git('rev-parse', 'HEAD').decode().strip() == value['approved_head']
    assert git('branch', '--show-current').decode().strip() == 'research/flash-publication-foundation'
    assert not git('status', '--porcelain').strip()
    subprocess.run(['git', 'merge-base', '--is-ancestor', BASE, value['approved_head']], cwd=REPO, check=True)
    assert set(git('diff', '--name-only', BASE, value['approved_head']).decode().splitlines()) == set(PUBLIC)
    assert all(sha(regular(REPO/name)) == digest for name, digest in value['public_sha256'].items())
    assert all(sha(regular(WORK/'delivery-preparation'/name)) == digest for name, digest in value['helper_sha256'].items())


def prefix(snapshot, live):
    assert snapshot.endswith(b'\n') and live.startswith(snapshot)
    lines = snapshot.splitlines()
    assert lines and all(isinstance(decode(line), dict) for line in lines)
    return {'bytes': len(snapshot), 'complete_records': len(lines),
            'sha256': hashlib.sha256(snapshot).hexdigest()}


def records(value):
    for name, digest in value['evidence_sha256'].items():
        relative(name)
        assert sha(regular(WORK/name)) == digest
    for name, requirements in value['record_requirements'].items():
        assert name in value['evidence_sha256']
        record = load(WORK/name)
        for row in requirements:
            selected = record
            for key in row['path']:
                selected = selected[key]
            assert json.dumps(selected, sort_keys=True, allow_nan=False) == json.dumps(row['equals'], sort_keys=True, allow_nan=False)
    for item in value['omitted_raw_inputs']:
        assert sha(regular(Path(item['path']))) == item['sha256']
        assert Path(item['path']).stat().st_size == item['bytes']


def verify_archive(value):
    """Read actual compressed bytes; never import acquisition or scientific code."""
    delivery = load(WORK/'frozen_delivery.json')
    archive = OUT.with_suffix('.tar.gz')
    assert delivery['folder'] == str(OUT) and delivery['archive'] == str(archive)
    assert delivery['sha256'] == sha(regular(archive)) and delivery['bytes'] == archive.stat().st_size
    assert delivery['source_commit'] == value['approved_head'] and delivery['base_commit'] == BASE
    with tarfile.open(archive, 'r:gz') as tar:
        files, seen = {}, set()
        for member in tar.getmembers():
            path = PurePosixPath(member.name)
            assert not path.is_absolute() and '..' not in path.parts and path.parts[0] == OUT.name
            assert path.as_posix() == member.name
            assert member.name not in seen and not member.issym() and not member.islnk()
            seen.add(member.name)
            if member.isdir():
                continue
            assert member.isfile()
            name = PurePosixPath(*path.parts[1:]).as_posix()
            assert name != '.'
            files[name] = tar.extractfile(member).read()
    manifest = decode(files['manifest.sha256.json'])
    assert set(files) == set(manifest)|{'manifest.sha256.json'}
    assert all(hashlib.sha256(files[name]).hexdigest() == digest for name, digest in manifest.items())
    assert sha(OUT/'manifest.sha256.json') == delivery['manifest_sha256']
    assert len(manifest) == delivery['manifest_entries']
    assert {str(path.relative_to(OUT)) for path in OUT.rglob('*') if path.is_file()} == set(files)
    assert all(not path.is_symlink() for path in OUT.rglob('*'))
    assert all((OUT/name).read_bytes() == data for name, data in files.items())
    assert all(not name.lower().endswith(('.pdf', '.png', '.jpg', '.jpeg', '.tif', '.tiff', '.pyc')) for name in files)
    assert decode(files['packet.json'])['contract'] == value
    for name, digest in value['public_sha256'].items():
        assert hashlib.sha256(files['source/'+name]).hexdigest() == digest == sha(REPO/name)
    assert files['npl-processability.patch'] == git('diff', '--binary', BASE, value['approved_head'])
    packet = decode(files['packet.json'])
    for item in packet['copied_inputs']:
        data = files[item['packet_path']]
        assert hashlib.sha256(data).hexdigest() == item['packet_sha256']
        for projection in reversed(item['projections']):
            replacement = projection['replacement'].encode()
            assert data.count(replacement) == projection['count']
            data = data.replace(replacement, bytes.fromhex(projection['original_hex']))
        assert hashlib.sha256(data).hexdigest() == item['original_sha256']
        original = regular(Path(item['source']))
        if item['packet_path'] == 'evidence/commands.jsonl':
            assert prefix(data, original.read_bytes()) == packet['command_snapshot']
        else:
            assert sha(original) == item['original_sha256']
    clean_head(value)
    records(value)
    return {**delivery, 'status': 'passed', 'actual_archive_regular_files': len(files),
            'unsafe_or_duplicate_members': 0, 'public_files_byte_identical': 2,
            'raw_PDFs_distributed': False, 'source_requests_or_scientific_checks_repeated': False,
            'command_snapshot_complete_record_prefix': True}
