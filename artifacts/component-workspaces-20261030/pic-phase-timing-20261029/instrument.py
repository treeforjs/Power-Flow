import difflib, hashlib, json, pathlib, re, shutil

base = pathlib.Path(__file__).parent
originals = json.loads((base / "evidence/original-sources.json").read_text())
summary = {}


def marked(text):
    return text + " // PHASE_TIMING\n"


def once(text, needle, replacement):
    assert text.count(needle) == 1, (needle, text.count(needle))
    return text.replace(needle, replacement)


for variant, meta in originals.items():
    src = base / variant / "source"
    shutil.copy2(base / "phase_profile.hpp", src / "src/phase_profile.hpp")
    names = ["src/bench_profile.hpp", "src/cudss_poisson.hpp",
             "src/mfem_pic_distributed.hpp", "src/mfem_gap_pic.cpp",
             "examples/mfem_heating_history_pic.cpp"]
    files = {name: (src / name).read_text() for name in names}
    changed = dict(files)
    changed["src/bench_profile.hpp"] += marked('#include "phase_profile.hpp"')
    name = "src/cudss_poisson.hpp"
    text = changed[name]
    needle = '  pic_bench::Timer timer("cuda_factor_setup");\n'
    text = once(text, needle, needle + marked(
        '  pic_phase::Scope phase_factor(pic_phase::Phase::FactorSetup);'))
    needle = '  pic_bench::Timer timer("cuda_solve_including_transfers");answer.SetSize(n);\n'
    text = once(text, needle, needle + marked(
        '  pic_phase::Scope phase_cuda(pic_phase::solve_phase(pic_phase::Phase::CudaFirst,pic_phase::Phase::CudaSecond,pic_phase::Phase::CudaOther));'))
    changed[name] = text
    name = "src/mfem_pic_distributed.hpp"
    text = changed[name]
    needle = '      Vector error(rhs.Size());matrix.Mult(answer,error);error-=rhs;\n'
    text = once(text, needle, marked(
        '      pic_phase::Scope phase_residual(pic_phase::solve_phase(pic_phase::Phase::ResidualFirst,pic_phase::Phase::ResidualSecond,pic_phase::Phase::ResidualOther));') + needle)
    changed[name] = text
    name = "src/mfem_gap_pic.cpp"
    text = changed[name]
    needle = ' pic_bench::Timer timer("advance_total");\n'
    text = once(text, needle, needle + marked(
        ' pic_phase::Scope phase_advance(pic_phase::Phase::Advance);'))
    begin = text.index('void GapPic::advance(double dt) {')
    end = text.index('void GapPic::advance_empty', begin)
    segment = text[begin:end]
    needle = '    drift(.5*dt); solve_fields();\n'
    assert segment.count(needle) == 2
    for which in ["First", "Second"]:
        suffix = which.lower()
        replacement = (
            marked(f'    pic_phase::Scope phase_drift_{suffix}(pic_phase::Phase::Drift{which});') +
            '    drift(.5*dt);\n' +
            marked(f'    phase_drift_{suffix}.stop();') +
            marked(f'    pic_phase::Scope phase_fields_{suffix}(pic_phase::Phase::SolveFields{which});') +
            '    solve_fields();\n' +
            marked(f'    phase_fields_{suffix}.stop();'))
        segment = segment.replace(needle, replacement, 1)
    needle = '    for (auto& p : particles_) {\n'
    segment = once(segment, needle, marked(
        '    pic_phase::Scope phase_force(pic_phase::Phase::ForceLoop);') + needle)
    needle = '    pic_phase::Scope phase_drift_second'
    segment = segment.replace(needle, marked('    phase_force.stop();') + needle, 1)
    text = text[:begin] + segment + text[end:]
    needle = 'void GapPic::advance_empty(double dt) {\n'
    text = once(text, needle, needle + marked(
        '    pic_phase::Scope phase_empty(pic_phase::Phase::AdvanceEmpty);'))
    needle = 'GapPicDiagnostics GapPic::diagnostics() const {\n'
    text = once(text, needle, needle + marked(
        '    pic_phase::Scope phase_diagnostics(pic_phase::Phase::Diagnostics);'))
    changed[name] = text
    name = "examples/mfem_heating_history_pic.cpp"
    text = changed[name]
    needle = '        while(start+pic.diagnostics().time_s<end-1e-22) {\n'
    replacement = (
        '        while((pic_phase::begin_step(step+1), start+pic.diagnostics().time_s<end-1e-22)) {\n' +
        marked('            pic_phase::Scope phase_selection(pic_phase::Phase::DriverSelection);'))
    text = once(text, needle, replacement)
    needle = '            auto T=thermal.at(t+.5*dt);double jmax=0;\n'
    text = once(text, needle, marked('            phase_selection.stop();') +
                marked('            pic_phase::Scope phase_source(pic_phase::Phase::ThermalSource);') + needle)
    needle = '            if(!injected && pic.particles().empty())pic.advance_empty(dt);else pic.advance(dt);++step;\n'
    text = once(text, needle, marked('            phase_source.stop();') +
                marked('            pic_phase::interval(t,dt);') + needle)
    needle = '            if(start+d.time_s>=next_frame*1e-9-1e-21)'
    text = once(text, needle, marked(
        '            pic_phase::Scope phase_output(pic_phase::Phase::PeriodicOutput);') + needle)
    needle = '            if(std::chrono::duration<double>(clock()-saved).count()>60)save();\n        }\n'
    text = once(text, needle,
                '            if(std::chrono::duration<double>(clock()-saved).count()>60)save();\n' +
                marked('            phase_output.stop();') +
                marked('            pic_phase::end_step();') + '        }\n' +
                marked('        pic_phase::end_terminal_condition();'))
    needle = 'pic_bench::resources((out/"process_resources.json").string());'
    text = once(text, needle, needle + '\n' +
                marked('        pic_phase::write((out/"phase_profile.json").string());'))
    changed[name] = text
    diffs, proof = [], []
    for name, text in changed.items():
        before = files[name]
        stripped = ''.join(line for line in text.splitlines(keepends=True)
                           if not line.rstrip().endswith('// PHASE_TIMING'))
        stripped = stripped.replace(
            'while((pic_phase::begin_step(step+1), start+pic.diagnostics().time_s<end-1e-22)) {',
            'while(start+pic.diagnostics().time_s<end-1e-22) {')
        assert re.sub(r'\s+', '', before) == re.sub(r'\s+', '', stripped), name
        (src / name).write_text(text)
        diffs.extend(difflib.unified_diff(
            before.splitlines(True), text.splitlines(True),
            fromfile=meta['source'] + '/' + name, tofile=str(src / name)))
        proof.append({
            'path': name, 'only_timing_scaffolding_after_removal': True,
            'original_sha256': hashlib.sha256(before.encode()).hexdigest(),
            'instrumented_sha256': hashlib.sha256(text.encode()).hexdigest()})
    (base / 'evidence' / f'{variant}.diff').write_text(''.join(diffs))
    summary[variant] = {
        'original_source': meta['source'], 'changed_files': proof,
        'new_header': 'src/phase_profile.hpp',
        'physics_token_stream_identical_after_scaffolding_removal': True}
(base / 'evidence/instrumentation-validation.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps(summary, indent=2))
