from pathlib import Path
import os,json,hashlib,subprocess,resource,tarfile,gzip,re
P=Path(__file__).resolve().parent
O=Path(os.environ.get('PIC_STREAMING_WRAPPER_OUTPUT','/home/treefor/sim-workspace/powerflow/outputs/pic-generic-water-streaming-wrapper-20261030'))
assert os.sched_getaffinity(0)=={7} and os.getpriority(os.PRIO_PROCESS,0)==15
resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
resource.setrlimit(resource.RLIMIT_CPU,(180,180))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
summary=json.loads((O/'summary.json').read_text());assert summary['status']=='passed'
build=json.loads((O/'build-receipt.json').read_text())
assert build['wrapper_source_sha256']==sha(P/'source/streaming_water_bridge.hpp')
assert build['observed_builder_sha256']==sha(P/'source/observed_candidates.hpp')
files=[p for p in sorted(P.rglob('*')) if p.is_file() and p.relative_to(P).parts[0] not in ['build','__pycache__'] and p.name!='source-manifest.json']
manifest={str(p.relative_to(P)):sha(p) for p in files}
manifest_path=P/'source-manifest.json';manifest_path.write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
source_id=sha(manifest_path);(O/'source-manifest.json').write_bytes(manifest_path.read_bytes())
archive=O/'generic-water-streaming-wrapper-source-v1.tar.gz'
with archive.open('wb') as raw,gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as compressed,tarfile.open(fileobj=compressed,mode='w') as tar:
 for p in sorted(files+[manifest_path]):
  info=tar.gettarinfo(str(p),arcname=str(p.relative_to(P)));info.mtime=0;info.uid=info.gid=0;info.uname=info.gname=''
  with p.open('rb') as f:tar.addfile(info,f)
ldd=subprocess.check_output(['/usr/bin/ldd',str(P/'build/test_wrapper')],text=True)
runtime={path:sha(Path(path)) for path in sorted(set(re.findall(r'/[^\s()]+',ldd))) if Path(path).is_file()}
toolchain={'compiler_version':subprocess.check_output(['/usr/bin/g++','--version'],text=True),'compiler_sha256':sha(Path('/usr/bin/g++')),'linked_library_hashes':runtime,'ldd':ldd,'no_dependencies_installed':True}
(O/'frozen-toolchain.json').write_text(json.dumps(toolchain,indent=2)+'\n')
proposal={'schema':'generic-v6-streaming-wrapper-integration-proposal-v1','entrypoint':'raikou::streaming_water::Bridge::candidates','physical_reference_unchanged':True,'native_or_joint_admission':False,'complete_mechanism':False,'old_generic_v6_immutable':True,'new_evidence_schema':'generic-v6-streaming-wrapper-evidence-v1','preserve':['one competing candidate set','all eight physical callbacks and original ordering','original target-rest energy/flux/unequal-weight rate/remainder bytes','original query/candidate budgets and exact unsupported-kernel refusals','all raw unapplied positive/unknown diagnostic records','physical product and scheduler RNG operations'],'host_requirements':['bind exact replayable packet/cell input, physical positions, authoritative common stage, geometry, provider/kernel/model/missing-rate contracts','keep batch/evidence private through any later accepted transaction','use a separately versioned native completion adapter; this evidence is not a completion receipt','retain explicit numerical-domain/binding/work refusals; no lossy moments or discarded unsupported positives','do not run the old 8.9s full-provider diagnostic every production step; no production budget qualification here'],'evidence_weight_measure':'projectile-carrier weight per evaluated electron-H2O pair; repeated pairs are not a unique population census','remaining_limits':['full conditional chemistry audit unchanged','native multi-channel/product-family/completion admission remains separate','exact-site attachment policy unchanged','1keV model seam and missing positive rates unchanged','finite moment/weight totals and preserved subnormal/rounding environment required','complete per-query metadata hashing still costs work; no every-step production timing claim']}
(O/'integration-proposal.json').write_text(json.dumps(proposal,indent=2)+'\n')
(O/'report.txt').write_bytes((P/'README.txt').read_bytes())
handoff={'status':'frozen_bounded_wrapper_qualification_passed','source_identity_sha256':source_id,'source_archive':str(archive),'source_archive_sha256':sha(archive),'wrapper_source_sha256':build['wrapper_source_sha256'],'observed_builder_source_sha256':build['observed_builder_sha256'],'stream_source_sha256':build['stream_header_sha256'],'original_generic_source_sha256':manifest['provider/source/generic_water_bridge.hpp'],'original_scheduler_source_sha256':manifest['provider/inputs/collision_scheduler.hpp'],'assertions':summary['assertions'],'independent_process_case_receipts_identical':summary['independent_process_case_receipts_identical'],'complete_provider_projections':4139,'physical_error_parity_cases':14,'new_evidence_numeric_domain_examples':2,'native_field_joint_execution':False,'complete_mechanism':False,'summary':str(O/'summary.json'),'report':str(O/'report.txt'),'integration_proposal':str(O/'integration-proposal.json'),'cpu':7,'nice':15,'threads':1}
(O/'handoff.json').write_text(json.dumps(handoff,indent=2)+'\n')
outputs={str(p.relative_to(O)):sha(p) for p in sorted(O.rglob('*')) if p.is_file() and p.name not in ['output-manifest.json','READY']}
(O/'output-manifest.json').write_text(json.dumps(outputs,indent=2,sort_keys=True)+'\n')
(O/'READY').write_text('Frozen wrapper/evidence component only; no native/field/joint execution or complete mechanism.\n')
print(json.dumps(handoff,indent=2))
