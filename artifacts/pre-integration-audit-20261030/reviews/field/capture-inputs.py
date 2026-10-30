from pathlib import Path
import hashlib,json,shutil,subprocess

here=Path(__file__).resolve().parent
repo=Path('/home/treefor/projects/runtimes/hall-coupled-transport-20261029/source')
commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip()
assert commit=='c2f81ab87ddc295e7bdad808d33135fa781b3942'
assert not subprocess.check_output(['git','status','--porcelain'],cwd=repo,text=True)
root=here/'inputs';root.mkdir(exist_ok=False)
entries=[]
def copy(source,relative):
 data=source.read_bytes();target=root/relative;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data)
 assert source.read_bytes()==data
 entries.append({'original':str(source),'path':str(target.relative_to(here)),
                 'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()})
for name in ['hall-canonical-driver-20261029','hall-finite-foil-em-20261029','hall-coupled-transport-20261029']:
 for path in sorted((repo/'development'/name).rglob('*')):
  if path.is_file():copy(path,Path('hall-source')/path.relative_to(repo))
for name in ['json.h','json.hpp']:
 copy(repo/'include/raikou'/name,Path('hall-source/include/raikou')/name)
native=Path('/home/treefor/projects/runtimes/pic-reaction-transactions-20261029/native-full-time')
for relative in ['MOTION-RECEIPT-SCHEMA.txt','source/include/raikou/mfem_pic_full_time.h',
                 'source/include/raikou/mfem_pic_arena.h','source/src/mfem_pic_arena.cpp',
                 'source/src/mfem_pic_full_time.inc']:
 copy(native/relative,Path('native-owner-draft')/relative)
parent=Path('/home/treefor/projects/runtimes/pic-water-full-time-subflow-20261029/source')
for name in ['native_parent_subflow.hpp','closed_parent_joint_step.hpp']:
 copy(parent/name,Path('parent-owner-draft')/name)
data=Path('/home/treefor/projects/Power-Flow/data/dissertation')
for name in ['shot_17133_current.json','shot_17133_current.csv','figure_4_28_simulation_drive.json',
             'figure_4_29_simulation_drive.json','figure_4_28_magnetic_gap.json',
             'figure_4_29_magnetic_gap.json','figure_4_30_electromagnetic.json',
             'pilot_comparison_manifest.json','target_optics_geometry.json']:
 copy(data/name,Path('thesis-evidence')/name)
source_review=here.parent/'source-v3-review-001'
for path in sorted(source_review.rglob('*')):
 if path.is_file():copy(path,Path('source-v3-review')/path.relative_to(source_review))
(here/'input-manifest.json').write_text(json.dumps({'schema':'hall-full-audit-inputs-v1',
 'hall_commit':commit,'hall_git_clean':True,'source_v3_source_snapshot':'source-v3-review-001 verified 29-file owner seal',
 'native_and_parent_status':'draft APIs captured immutably for this audit; no full-link qualification inferred',
 'files':entries},indent=2)+'\n')
print(json.dumps({'hall_commit':commit,'captured_files':len(entries)}))
