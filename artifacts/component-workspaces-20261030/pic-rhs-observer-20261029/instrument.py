import difflib
import hashlib
import json
import pathlib
import re

base = pathlib.Path(__file__).parent
meta = json.loads((base / "evidence/original-sources.json").read_text())
source = base / "source"
names = ["src/bench_profile.hpp", "src/mfem_gap_pic.cpp",
         "src/mfem_pic_distributed.hpp", "examples/mfem_heating_history_pic.cpp"]
before = {name: (source / name).read_text() for name in names}
after = dict(before)


def marked(text):
    return text + " // RHS_OBSERVATION\n"


def once(text, needle, replacement):
    assert text.count(needle) == 1, (needle, text.count(needle))
    return text.replace(needle, replacement)


after["src/bench_profile.hpp"] += marked('#include "rhs_observer.hpp"')
name = "src/mfem_gap_pic.cpp"
text = after[name]
needle = '    form.FormLinearSystem(essential,trial,rhs,matrix,x,b);\n'
text = once(text, needle, marked(
    '    const auto observed_call=rhs_observer::before_form(trial.GetData(),full_rhs.GetData(),trial.Size(),essential.GetData(),essential.Size(),reinterpret_cast<std::uintptr_t>(poisson_.get()));') + needle)
needle = '    if (!sparse) throw std::runtime_error("Gap PIC expected sparse Poisson matrix");\n'
text = once(text, needle, needle + marked(
    '    rhs_observer::after_form(observed_call,sparse->GetI(),sparse->GetJ(),sparse->GetData(),sparse->Height(),sparse->GetI()[sparse->Height()],x.GetData(),b.GetData(),reinterpret_cast<std::uintptr_t>(sparse));'))
needle = '        throw std::invalid_argument("Gap PIC advance requires positive finite dt");\n'
text = once(text, needle, needle + marked(
    '    rhs_observer::AdvanceScope observed_advance(!particles_.empty(),particles_.size(),ledger_.time_s,dt,options_.lower_potential_V,options_.upper_potential_V);'))
begin = text.index('void GapPic::advance(double dt) {')
end = text.index('void GapPic::advance_empty', begin)
segment = text[begin:end]
needle = '    drift(.5*dt); solve_fields();\n'
assert segment.count(needle) == 2
for half in [1, 2]:
    segment = segment.replace(needle, '    drift(.5*dt);\n' +
                              marked(f'    rhs_observer::set_half({half});') +
                              '    solve_fields();\n', 1)
text = text[:begin] + segment + text[end:]
after[name] = text
name = "src/mfem_pic_distributed.hpp"
needle = '      MFEM_VERIFY(std::isfinite(residual),"Nonfinite solver residual");\n'
after[name] = once(after[name], needle, needle + marked(
    '      rhs_observer::accepted(rhs.GetData(),answer.GetData(),error.GetData(),rhs.Size(),last_residual_norm,norm,id,reinterpret_cast<std::uintptr_t>(this));'))
name = "examples/mfem_heating_history_pic.cpp"
needle = '            if(!injected && pic.particles().empty())pic.advance_empty(dt);else pic.advance(dt);++step;\n'
after[name] = once(after[name], needle, marked(
    '            rhs_observer::driver_step(step+1,t,dt);') + needle)
needle = 'pic_bench::resources((out/"process_resources.json").string());'
after[name] = once(after[name], needle, needle + '\n' +
                   marked('        rhs_observer::finish();') + '        ')
diffs, proof = [], []
for name, text in after.items():
    stripped = ''.join(line for line in text.splitlines(keepends=True)
                       if not line.rstrip().endswith('// RHS_OBSERVATION'))
    assert re.sub(r'\s+', '', stripped) == re.sub(r'\s+', '', before[name]), name
    (source / name).write_text(text)
    proof.append({
        "path": name,
        "original_sha256": hashlib.sha256(before[name].encode()).hexdigest(),
        "observed_sha256": hashlib.sha256(text.encode()).hexdigest(),
        "scientific_source_identical_after_removing_observation_scaffolding": True})
    diffs.extend(difflib.unified_diff(
        before[name].splitlines(True), text.splitlines(True),
        fromfile=meta["original_source"] + "/" + name, tofile=str(source / name)))
(base / "evidence/source.diff").write_text(''.join(diffs))
(base / "evidence/source-removal-proof.json").write_text(json.dumps({
    "original_source": meta["original_source"], "files": proof,
    "all_original_physics_tokens_and_order_preserved": True,
    "CudssPoisson_source_unmodified": True}, indent=2) + "\n")
print("Observation hooks inserted; removal restores all original source tokens and ordering.")
