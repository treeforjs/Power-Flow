from pathlib import Path
import os, json, hashlib, subprocess, time, resource, tarfile

P = Path(__file__).resolve().parent
O = Path('/home/treefor/Documents/Codex/2026-10-28/we/outputs/pic-water-full-time-subflow-20261029')
assert P.name == 'pic-water-full-time-subflow-20261029'
if not (P/'native-handoff.json').is_file():
    raise SystemExit('Native full-time implementation handoff is not installed; draft declarations are insufficient.')
assert os.sched_getaffinity(0) == {7}
assert os.getpriority(os.PRIO_PROCESS, 0) == 15
resource.setrlimit(resource.RLIMIT_AS, (3*1024**3, 3*1024**3))
F = '/home/treefor/projects/runtimes/arch-base-20261028/run-frozen'
M = Path('/home/treefor/projects/runtimes/mfem-4.10-pic-20261029')
for name, want in json.loads((P/'kernel-binding.json').read_text())['sources'].items():
    assert hashlib.sha256((P/name).read_bytes()).hexdigest() == want, name
handoff = json.loads((P/'native-handoff.json').read_text())
identity = handoff['source_identity']
archive = P/'native-source.tar.gz'
assert hashlib.sha256(archive.read_bytes()).hexdigest() == handoff['source_archive_sha256']
with tarfile.open(archive) as tar:
    for entry in tar.getmembers():
        if entry.isfile():
            want=hashlib.sha256(tar.extractfile(entry).read()).hexdigest()
            assert hashlib.sha256((P/'native'/entry.name).read_bytes()).hexdigest() == want, entry.name
includes = ['-I'+str(P/'source'), '-I'+str(P/'native/source/include'),
            '-isystem', str(P/'inputs'), '-isystem', str(P/'inputs/include'),
            '-isystem', str(P/'native/inputs'),
            '-isystem', str(P/'native/vendor')]
for suffix in ['install-mpi/include', 'tpl/hypre/include', 'tpl/metis/include']:
    includes += ['-isystem', str(M/suffix)]
cc = [F, '/usr/bin/mpic++', '-std=c++20', '-O2', '-g', '-ffunction-sections',
      '-fdata-sections', '-DPIC_REACTION_NATIVE_SOURCE_SHA256="'+identity+'"'] + includes
commands=[]
def run(command, name, timeout=120):
    commands.append(command)
    (O/'commands.json').write_text(json.dumps(commands, indent=2)+'\n')
    with (O/(name+'.log')).open('w') as log:
        return subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=timeout)
objects=[]
for path, name, strict in [(P/'native/source/src/mfem_gap_pic.cpp', 'native', False),
                           (P/'native/source/src/mfem_pic_arena.cpp', 'arena', True),
                           (P/'tests/test_native_parent_subflow.cpp', 'subflow', True)]:
    obj=P/'build'/(name+'.o'); objects.append(str(obj))
    run(cc + (['-Wall','-Wextra','-Wpedantic','-Werror'] if strict else []) +
        ['-c',str(path),'-o',str(obj)], 'build-'+name)
libs=[str(M/x) for x in ['install-mpi/lib/libmfem.a','tpl/hypre/lib/libHYPRE.a',
                         'tpl/metis/lib/libmetis.a','tpl/GKlib/lib/libGKlib.a']]+['-lcrypto']
exe=P/'build/test_native_parent_subflow'
run([F,'/usr/bin/mpic++','-Wl,--gc-sections']+objects+libs+['-o',str(exe)],'link')
start=time.monotonic()
run([F,str(exe),str(P/'inputs/h2o-gas.json'),str(O/'run-binding.json'),str(O),identity],'checks',180)
result=json.loads((O/'result.json').read_text())
assert result['status']=='passed'
source={str(f.relative_to(P)):hashlib.sha256(f.read_bytes()).hexdigest()
        for d in ['source','tests','inputs'] for f in sorted((P/d).rglob('*')) if f.is_file()}
(O/'build-receipt.json').write_text(json.dumps({'status':'passed','source_files':source,
    'native_source_identity':identity,'native_archive_sha256':handoff['source_archive_sha256'],
    'executable_sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),'elapsed_checks_s':time.monotonic()-start,
    'cpu':7,'nice':15,'threads':1,'gpu_used':False,'frozen_launcher':F,
    'checks':result['check_count'],'physical_run':False,
    'scope':'legacy regression only; two-stage full-time integration is not yet tested'},indent=2)+'\n')
print(json.dumps({'checks':result['check_count'],'first_events':result['first_events'],
                  'elapsed_checks_s':time.monotonic()-start}))
