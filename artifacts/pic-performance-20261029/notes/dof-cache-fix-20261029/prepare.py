import pathlib,shutil,json,hashlib,datetime,difflib,os
out=pathlib.Path(__file__).parent
shared=pathlib.Path("/home/treefor/projects/runtimes/pic-cpu-candidate-20261029/source")
source=out/"source"
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
before={str(p.relative_to(shared)):sha(p) for p in shared.rglob("*") if p.is_file()}
(out/"shared-source-before.json").write_text(json.dumps(before,indent=2)+"\n")
shutil.copytree(shared,source)
p=source/"include/raikou/q1_scalar_gather.h"
s=p.read_text()
s=s.replace("// use the ordinary GridFunction operations.\n", "// use the ordinary GridFunction operations. Connectivity must stay unchanged\n// while cached: call invalidate() before unversioned DOF renumbering (for\n// example, FiniteElementSpace::ReorderElementToDofTable()).\n",1)
s=s.replace("    bool available() const {\n", "    // Permanently fall back to MFEM. Construct a new gather to cache again.\n    void invalidate() { dofs_.clear(); }\n\n    bool available() const {\n",1)
p.write_text(s)
p=source/"tests/test_q1_scalar_gather.cpp"
s=p.read_text()
needle="void compare(mfem::Mesh& mesh,int order,bool supported) {\n"
addition="""void dof_renumbering_invalidation() {
    auto mesh=mfem::Mesh::MakeCartesian2D(3,2,mfem::Element::QUADRILATERAL,true,2.,3.,false);
    mfem::H1_FECollection collection(1,2);
    mfem::FiniteElementSpace space(&mesh,&collection);
    mfem::GridFunction field(&space);
    for(int i=0;i<field.Size();++i) field[i]=i+1.;
    Q1ScalarGather gather(field);
    check(gather.available(),\"Q1 gather initially cached\");
    const long mesh_sequence=mesh.GetSequence(),space_sequence=space.GetSequence();
    // MFEM does not advance either sequence when it renumbers element DOFs.
    gather.invalidate();
    space.ReorderElementToDofTable();
    check(mesh.GetSequence()==mesh_sequence && space.GetSequence()==space_sequence,
          \"DOF renumbering regression must exercise unchanged sequence counters\");
    check(!gather.available(),\"explicit invalidation rejects renumbered connectivity\");
    gather.invalidate(); // Idempotent, including an already-disabled cache.
    Q1ScalarGather rebuilt(field);
    check(rebuilt.available(),\"a new gather caches the renumbered connectivity\");
    for(int e=0;e<mesh.GetNE();++e) {
        for(double x:{0.,.23,1.}) for(double y:{0.,.61,1.}) {
            mfem::IntegrationPoint point;point.Set2(x,y);
            const double value=field.GetValue(e,point);
            identical(gather.value(e,point),value);
            identical(rebuilt.value(e,point),value);
            auto* transformation=mesh.GetElementTransformation(e);
            transformation->SetIntPoint(&point);
            mfem::Vector reference(2);field.GetGradient(*transformation,reference);
            const auto fallback=gather.gradient(e,point),cached=rebuilt.gradient(e,point);
            for(int i=0;i<2;++i) {
                identical(fallback[i],reference[i]);identical(cached[i],reference[i]);
            }
        }
    }
}
"""
assert needle in s
s=s.replace(needle,addition+needle,1)
s=s.replace("    try {\n", "    try {\n        dof_renumbering_invalidation();\n",1)
s=s.replace("Q1 scalar gathers and fallbacks match MFEM bit for bit", "Q1 scalar gathers, explicit DOF invalidation and fallbacks match MFEM bit for bit")
p.write_text(s)
patch=""
changes=[]
for rel,digest in before.items():
    p=source/rel
    if sha(p)!=digest:
        changes.append({"path":rel,"shared_sha256":digest,"fixed_sha256":sha(p)})
        patch+="".join(difflib.unified_diff((shared/rel).read_text().splitlines(keepends=True),p.read_text().splitlines(keepends=True),fromfile="a/"+rel,tofile="b/"+rel))
assert {c["path"] for c in changes}=={"include/raikou/q1_scalar_gather.h","tests/test_q1_scalar_gather.cpp"}
(out/"dof-cache-invalidation.patch").write_text(patch)
(out/"source-changes.json").write_text(json.dumps(changes,indent=2)+"\n")
(out/"progress.json").write_text(json.dumps({"status":"isolated_fix_prepared","created_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"cpu":6,"nice":15,"source":str(source),"shared_candidate_modified":False,"gpu_work":False,"fix":"explicit invalidation for unversioned connectivity mutations; no hot-path changes"},indent=2)+"\n")
print(json.dumps({"isolated_source":str(source),"changed_files":changes,"affinity":sorted(os.sched_getaffinity(0)),"nice":os.getpriority(os.PRIO_PROCESS,0)},indent=2))
