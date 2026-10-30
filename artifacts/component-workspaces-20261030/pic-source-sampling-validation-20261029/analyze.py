#!/usr/bin/env python3
"""Single-threaded standard-library analysis of deterministic native source replay."""
import csv
import hashlib
import json
import math
import os
import statistics as st
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = Path('/home/treefor/sim-workspace/powerflow/outputs/pic-profile/source-sampling-validation')
N = 32
T95 = 2.0395134463964077  # Student t, 31 degrees of freedom, two-sided 95%.
CHI95 = (17.538738216309, 48.231889640522)  # chi square, 31 df, 2.5/97.5%.
GROUPS = ['back upper', 'back lower', 'gap upper', 'gap lower', 'left upper', 'left lower', 'right upper', 'right lower', 'all']

def read_csv(name):
    with (OUT/name).open() as stream:
        return [{k:float(v) for k,v in row.items()} for row in csv.DictReader(stream)]

def describe(values, expected=None, predicted_variance=None):
    mean = st.mean(values)
    sd = st.stdev(values)
    half = T95*sd/math.sqrt(len(values))
    result = dict(mean=mean, sd=sd, ci95_mean=[mean-half,mean+half])
    if expected is not None:
        result.update(expected=expected, error=mean-expected,
                      relative_error=(mean/expected-1) if expected else None)
    if predicted_variance is not None and predicted_variance>0:
        result.update(predicted_sd=math.sqrt(predicted_variance),
                      mean_z=(mean-expected)/math.sqrt(predicted_variance/len(values)),
                      empirical_to_predicted_variance=sd*sd/predicted_variance,
                      ci95_sd_normal_approx=[sd*math.sqrt((N-1)/CHI95[1]),sd*math.sqrt((N-1)/CHI95[0])])
    return result

if os.sched_getaffinity(0)!={6} or os.getpriority(os.PRIO_PROCESS,0)!=15:
    raise RuntimeError('Analysis requires CPU6 nice15')
run = json.loads((OUT/'run-summary.json').read_text())
assert run['complete'] and run['source_only'] and not run['coupled_convergence']
rows = read_csv('seed-metrics.csv')
expected = {(int(r['M']),int(r['group'])):r for r in read_csv('expected-metrics.csv')}
actual = {(int(r['M']),int(r['group'])):[] for r in rows}
for row in rows:
    actual[int(row['M']),int(row['group'])].append(row)
for values in actual.values():
    values.sort(key=lambda r:r['seed'])
    assert [r['seed'] for r in values]==list(range(N))
report = dict(source_only=True,coupled_convergence=False,run=run,counts={},paired_comparisons={},confidence='32 deterministic seed offsets; independent-seed Student t intervals; variance intervals are normal approximations')
maximum_cumulative_charge_error = 0
for M in [64,32,16,8]:
    groups = {}
    for g,name in enumerate(GROUPS):
        values=actual[M,g]; exp=expected[M,g]
        max_charge=max(abs(r['physical_electrons']/exp['physical_electrons']-1) for r in values)
        maximum_cumulative_charge_error=max(maximum_cumulative_charge_error,max_charge)
        groups[name] = dict(
            events=int(exp['events']),samples=int(values[0]['samples']),
            expected_physical_electrons=exp['physical_electrons'],max_charge_relative_error=max_charge,
            nonrel_source_energy_eV=describe([r['nonrel_source_energy_eV'] for r in values],exp['source_energy_eV'],exp['variance_source_energy_eV2']),
            native_source_energy_eV=describe([r['source_energy_eV'] for r in values],exp['source_energy_eV'],exp['variance_source_energy_eV2']),
            slab_physical_electrons=describe([r['slab_physical_electrons'] for r in values],exp['slab_physical_electrons'],exp['variance_slab_electrons2']),
            patch_squared_error=describe([r['patch_squared_error'] for r in values],exp['expected_patch_squared_error']),
            missed_probability_mass=describe([r['missed_probability_mass'] for r in values],exp['expected_missed_probability_mass']),
            weight_effective_samples=st.mean([r['weight_effective_samples'] for r in values]),
            mean_native_relativistic_correction_eV=st.mean([r['source_energy_eV']-r['nonrel_source_energy_eV'] for r in values]))
    report['counts'][M]=groups
for M,ref in [(32,64),(16,32),(8,32),(16,64),(8,64)]:
    groups={}
    for g,name in enumerate(GROUPS):
        a,b=actual[M,g],actual[ref,g]
        exp=expected[M,g]
        energy=describe([(x['source_energy_eV']-y['source_energy_eV'])/exp['source_energy_eV'] for x,y in zip(a,b)])
        slab=describe([(x['slab_physical_electrons']-y['slab_physical_electrons'])/exp['slab_physical_electrons'] for x,y in zip(a,b)])
        energy['ci95_within_5pct']=max(abs(v) for v in energy['ci95_mean'])<=.05
        slab['ci95_within_5pct']=max(abs(v) for v in slab['ci95_mean'])<=.05
        groups[name]=dict(source_energy_relative_difference=energy,slab_charge_relative_difference=slab)
    report['paired_comparisons'][f'{M}_vs_{ref}']=groups
patch_rows=read_csv('patch-metrics.csv')
report['patch_ensemble']={}
def count_tail_bound(k,lam):
    # Poisson-binomial Chernoff bounds; valid for independent nonidentical Bernoulli trials.
    # The upper factorial-moment bound is useful when a rare cell is hit once.
    if lam==0:
        return 1. if k==0 else 0.
    if k==0:
        return math.exp(-lam)
    log_chernoff=-lam+k-k*math.log(k/lam)
    if k>lam:
        log_factorial=k*math.log(lam)-math.lgamma(k+1)
        return min(1.,math.exp(min(log_chernoff,log_factorial)))
    return min(1.,math.exp(log_chernoff))
family_size=len(patch_rows)
for M in [64,32,16,8]:
    values=[r for r in patch_rows if r['M']==M and r['predicted_variance']>0]
    regular=[r for r in values if r['ensemble_expected_hits']>=25 and r['ensemble_variance_hits']>=25]
    max_z=max(abs(r['z_mean']) for r in values)
    min_tail=min(count_tail_bound(int(r['ensemble_hits']),r['ensemble_expected_hits']) for r in values)
    charge_weighted_var=sum(r['ensemble_variance'] for r in values)/sum(r['predicted_variance'] for r in values)
    report['patch_ensemble'][M]=dict(nonzero_probability_patches=len(values),max_abs_mean_z=max_z,
        normal_eligible_patches=len(regular),max_abs_normal_eligible_mean_z=max(abs(r['z_mean']) for r in regular),
        summed_empirical_to_predicted_variance=charge_weighted_var,mean_z_over_six=sum(abs(r['z_mean'])>6 for r in values),
        minimum_count_tail_upper_bound=min_tail,bonferroni_family_size=family_size,
        any_count_tail_bound_below_family_001=min_tail<.01/family_size)
assert maximum_cumulative_charge_error<1e-12
assert all(not r['any_count_tail_bound_below_family_001'] for r in report['patch_ensemble'].values())
# These are calibration guards, not equivalence declarations for later coupled tests.
assert all(abs(report['counts'][M]['all']['nonrel_source_energy_eV']['mean_z'])<4 for M in [64,32,16,8])
assert all(abs(report['counts'][M]['all']['slab_physical_electrons']['mean_z'])<4 for M in [64,32,16,8])
report['maximum_cumulative_charge_relative_error']=maximum_cumulative_charge_error
report['common_injected_charge_C']=-1.602176634e-19*expected[32,8]['physical_electrons']
report['numerical_and_ensemble_guards_passed']=True
(OUT/'analysis.json').write_text(json.dumps(report,indent=2)+'\n')

lines=[
    'Source-sampling validation, 2026-10-29', '',
    'SOURCE ONLY. This is not coupled trajectory, field, density, timestep, mesh, or physics convergence.', '',
    'Completed CPU-only validation at M=64,32,16,8 with 32 paired deterministic seeds.',
    f"Captured {run['events']} electrode events over 100 ps starting at {run['absolute_start_s']*1e9:.12g} ns.",
    f"The common replay emits {expected[32,8]['physical_electrons']:.12g} physical electrons, independent of M within the 1e-12 guard.",
    f"Common physical injected charge: {report['common_injected_charge_C']:.12g} C (negative electron charge).",
    f"Maximum face/event relative charge error: {run['max_event_charge_relative_error']:.6g}; maximum cumulative error: {maximum_cumulative_charge_error:.6g}.",
    f"Exact position/velocity prefix checks: {run['prefix_checks']:,}; scaled-weight difference: {run['max_prefix_scaled_weight_relative_error']:.3g}.",
    f"Zero-current patch draws: {run['zero_current_patch_draws']}; additional synthetic face/orientation draws: {run['synthetic_orientation_draws']:,}.",
    f"Replay wall time including generation and refinement check: {run['wall_s']:.3f} s; CPU6, nice15, one thread.",
    'The final build, replay and analysis hold the shared testing-core.lock lease. Earlier development runs used CPU0',
    'under the previous assignment; the final results supersede those measurements after the CPU6 allocation correction.', '',
    'Native-code equivalence', '',
    'extract_native.py copies eight byte-identical source slices from the frozen source and checks their SHA256 digests.',
    'The sampled methods are the native inject_surface_profile and inject_thermal_electrons, including the same mt19937_64,',
    'libstdc++ discrete/uniform/normal distributions, geometry arithmetic, proper-velocity conversion, finite/weight checks,',
    'kinetic-energy routine, constants, particle/options/ledger types, and the native thermal interpolation and reservoir loop.',
    'A storage-only host replaces the MFEM-owning class. Every field read or written by those methods has its native type.',
    'Options are copied from the audited checkpoint. Particle IDs and ledgers start from their native defaults for each',
    'event; their offsets do not feed the random draws, positions, velocities or weights. The comparison concerns emitted',
    'source particles and ledger increments, not continuation of checkpoint IDs or its cumulative ledger. No native emission',
    'statement is rewritten; no GPU/MPI/MFEM library is linked or initialized.',
    'Compiler: /usr/sbin/c++, same -O3 -DNDEBUG floating-point settings, no fast-math. native-equivalence.json records line',
    'bounds and slice/file hashes. This establishes source/math equivalence for source sampling on this compiler/stdlib;',
    'it is not a comparison of whole-executable machine code or a cross-toolchain bitwise guarantee.', '',
    'Common physical source and clock', '',
    'Input: frozen heating-history-full.bin and the prior audited copied checkpoint_400040.bin/restart.meta.',
    'The original copied options and per-patch integral/temperature reservoirs are loaded without modifying either input.',
    'The native interpolation, Richardson law (phi=4.3 eV, A=1201732.289489034), one-electron face/electrode threshold,',
    'current integrals, weighted temperatures, four faces, lower-side reversal, reservoir reset and native base seed formula',
    'are preserved. The event sequence is captured once with a prescribed 10 fs midpoint quadrature; every M replays the',
    'same event time, span, electrode, surface, patch temperatures and current integrals. events.bin stores full T/J arrays;',
    'events.csv is the human-readable index. seeds.csv records all paired offsets. Seed s uses native_seed +',
    's*0x9e3779b97f4a7c15 modulo 2^64; seed zero retains the captured native seed. Count choices share nested prefixes.',
    'The prescribed clock is not the active adaptive trajectory clock. No field, displacement, plasma or residual guard',
    'was evaluated or changed. No live output/checkpoint, sealed executable, source, dependency, queue or reservation was edited.',
    f"Initial pending: {run['initial_pending_electrons']:.12g}; terminal pending (10 fs): {run['terminal_pending_electrons']:.12g} electrons.",
    f"Total emitted+pending (10 fs): {run['coarse_emitted_plus_pending_electrons']:.12g}; 5 fs: {run['fine_emitted_plus_pending_electrons']:.12g} electrons.",
    f"Relative supplied-charge difference: {run['quadrature_relative_difference']:.6g}; event counts: {run['events']} versus {run['fine_quadrature_events']}.",
    f"Terminal pending (5 fs): {run['fine_terminal_pending_electrons']:.12g}; emitted (5 fs): {run['fine_emitted_electrons']:.12g} electrons.",
    'The refinement difference in flush counts is caused by threshold timing/overshoot and pending source. Equality of',
    'integrated supply does not imply equality of event times or terminal emitted charge. The 5 fs sequence is a source',
    'quadrature sensitivity check only and is never mixed into the M comparison.', '',
    'Effusive energy and variance model', '',
    'At fixed T the two tangential Gaussian velocities contribute an exponential kinetic energy with mean kT.',
    'The flux-weighted normal speed contributes another exponential with mean kT. Thus nonrelativistic K~Gamma(2,kT),',
    'E[K]=2kT and Var[K]=2(kT)^2. For a patch mixture p_i=J_i/sum J, mean=2 sum p_i kT_i and',
    'variance=6 sum p_i(kT_i)^2 - (2 sum p_i kT_i)^2. Event total energy variance=N_event^2 Var[K]/M.',
    'Independent captured events add variances. Actual native relativistic kinetic energy is reported separately:',
    'the sampled nonrelativistic energy tests the exact Gamma moments; the proper-velocity correction is small but retained.',
    'Patch charge variance is N_event^2 p_i(1-p_i)/M. Summing over patches gives the predicted squared-error trace.',
    'The expected missed source-probability mass is sum p_i(1-p_i)^M, averaged equally over captured events.',
    'The 1 mm midplane slab is 4% of the fixed 25 mm height: E[N_slab]=0.04 N, Var=sum N_event^2*0.04*0.96/M.', '',
    'All-event ensemble results (one 100 ps source replay per seed)', '',
    'M | macros/seed | source energy eV mean +/- SD [95% mean CI] | predicted SD | slab electrons mean +/- SD [95% mean CI] | predicted SD',
]
for M in [64,32,16,8]:
    r=report['counts'][M]['all']; en=r['native_source_energy_eV']; slab=r['slab_physical_electrons']
    lines.append(f"{M} | {r['samples']} | {en['mean']:.6f} +/- {en['sd']:.6f} [{en['ci95_mean'][0]:.6f},{en['ci95_mean'][1]:.6f}] | {en['predicted_sd']:.6f} | {slab['mean']:.6f} +/- {slab['sd']:.6f} [{slab['ci95_mean'][0]:.6f},{slab['ci95_mean'][1]:.6f}] | {slab['predicted_sd']:.6f}")
lines += ['',f"Expected total nonrelativistic energy: {expected[32,8]['source_energy_eV']:.12g} eV; expected mean per represented electron: {expected[32,8]['source_energy_eV']/expected[32,8]['physical_electrons']:.12g} eV.",
          f"Expected slab source: {expected[32,8]['slab_physical_electrons']:.12g} physical electrons.", '',
          'M | event-mean missed probability mass observed / predicted | summed patch variance observed / predicted | maximum eligible patch mean |z| | smallest count-tail upper bound']
for M in [64,32,16,8]:
    miss=report['counts'][M]['all']['missed_probability_mass']; p=report['patch_ensemble'][M]
    lines.append(f"{M} | {miss['mean']:.6%} / {miss['expected']:.6%} | {p['summed_empirical_to_predicted_variance']:.6f} | {p['max_abs_normal_eligible_mean_z']:.6f} | {p['minimum_count_tail_upper_bound']:.6g}")
lines += ['', 'Sampling SD scales as 1/sqrt(M): relative to32,16 increases SD by sqrt(2),8 by2, and64 reduces it by sqrt(2).',
          'Native geometry subtraction gives tiny patch-dependent roundoff in weights, most visible on side faces; all charge',
          'comparisons meet 1e-12. Positive-weight and finite checks run inside the unchanged native method.', '',
          'Paired all-event differences, normalized by the exact source expectation (mean and 95% CI)', '',
          'pair | source energy | slab source charge | both CIs within +/-5%']
for pair in ['32_vs_64','16_vs_32','8_vs_32','16_vs_64','8_vs_64']:
    r=report['paired_comparisons'][pair]['all']; e=r['source_energy_relative_difference']; s=r['slab_charge_relative_difference']
    lines.append(f"{pair} | {e['mean']:.4%} [{e['ci95_mean'][0]:.4%},{e['ci95_mean'][1]:.4%}] | {s['mean']:.4%} [{s['ci95_mean'][0]:.4%},{s['ci95_mean'][1]:.4%}] | {e['ci95_within_5pct'] and s['ci95_within_5pct']}")
lines += ['', 'The pooled all-source intervals fit +/-5% for these aggregate source observables. Per-face conclusions are',
          'weaker. Each back electrode emits only one event in this window, so its source-slab intervals are too wide',
          'for a relative5% test even at32 versus64. At8 versus32, five additional face/electrode slab intervals extend',
          'outside +/-5%: gap lower, left upper/lower and right upper/lower. At16 versus32 those five intervals fit',
          'inside +/-5%, while the sparse back faces remain unresolved on that relative scale. Near-zero regions need',
          'a predeclared absolute tolerance. Aggregate agreement must not be read as all-region equivalence.', '',
          '8 versus32 per-face slab source relative difference: mean [95% CI]']
for name,r in report['paired_comparisons']['8_vs_32'].items():
    if name=='all':
        continue
    s=r['slab_charge_relative_difference']
    lines.append(f"{name}: {s['mean']:.4%} [{s['ci95_mean'][0]:.4%},{s['ci95_mean'][1]:.4%}]")
lines += ['', 'Intervals describe the mean over randomized numerical source samples, not a physical shot-noise model.',
          'They use Student t with31 degrees of freedom and no multiplicity correction; per-face results in analysis.json are',
          'diagnostic, not family-wise equivalence claims. Raw patch z scores are retained but are misleading for rare cells:',
          'the first Gaussian all-patch check rejected a single hit where the ensemble expectation was only about3.7e-5.',
          'The final analysis uses Poisson-binomial Chernoff/factorial-moment count-tail upper bounds with a Bonferroni1%',
          'family threshold across all36864 patch/count tests. No bound falls below that threshold. Passing this conservative',
          'screen establishes no detected inconsistency, not proof of unbiased sampling for each rare patch. Gaussian z',
          'summaries are limited to patches with expected hit count and hit-count variance>=25. SD intervals in analysis.json assume',
          'approximately normal seed totals; sparse back-face/slab totals can be discrete and those intervals are approximate.',
          'Full per-seed/per-face results, expected moments and per-patch ensemble statistics are retained as CSV.', '',
          'Criteria for a later parent-controlled coupled test', '',
          'These are proposed screening gates, not a recommendation to change the active count or permission to launch a branch.',
          '1. First compare64 versus32 and16 versus32 with common checkpoint particles/reservoirs, the same accepted physical',
          '   times, replayed events and a common schedule satisfying every unchanged stability bound. Include8 only after16',
          '   is adequate. Preserve source current, one-electron threshold, thermal law, geometry and density/slab volumes.',
          '2. Measure fresh emitted and retained charge. If fresh retained charge is less than20% of terminal active charge,',
          '   classify the prospective count test inconclusive or extend under a newly declared cap. A late restart retains',
          '   the original32-sample population and cannot establish whole-history convergence or immediate proportional speedup.',
          '3. Use fixed coarse physical bins/regions and declared time averages for weighted full-height and1mm-slab density',
          '   integrals. Never use a single25um pixel maximum or screenshot similarity as the primary criterion. At4% slab',
          '   volume, per-event expected macro occupancy is2.56/1.28/0.64/0.32 at64/32/16/8; spatial subdivision is sparser.',
          '4. Predeclare engineering tolerances: paired95% CIs wholly inside +/-5% for non-negligible integrated active/source/',
          '   sink charge, energy moments/quantiles and regional densities; weighted energy-CDF distance<=0.05; coarse density',
          '   L1 difference<=10% of total represented population. For near-zero sinks use an absolute scale, e.g.1% of emitted',
          '   charge. An overlapping wide CI does not establish equivalence. Check consistency of32 with64 and increase seeds',
          '   when intervals are too wide. The32 seeds here apply to the source-only replay and do not replace coupled seeds.',
          '5. Preserve finite/positive states and all frozen solver/continuity/residual/conservation guards. Verify mapped',
          '   physical charge, injected charge to1e-12, and charge-balance residual<=1e-10 of cumulative emitted charge with',
          '   an electron-charge denominator floor. Report potential/field norms against predeclared physical scales.',
          '6. After isolating the count with replay, confirm a promising count with the original adaptive source/step rules.',
          '   Report supply, event-clock and pending-reservoir changes separately. Do not label their changes count-only.',
          '7. Record wall time per accepted physical ps, steps, fresh/total macro count, weighted effective sample sizes,',
          '   solver/particle phases, caps and incomplete status. Run scheduling, integration and deployment remain outside',
          '   this validation. No runtime speedup is established here.', '',
          'Limits', '',
          'No trajectories, electric/magnetic fields, Poisson solve, collective response, absorption or escape were computed.',
          'Source patch expectation and effusive moments validate sampling only. Lower sampling can amplify field/orbit noise,',
          'and source-level aggregate agreement cannot establish transported slab-density or energy-distribution convergence.',
          'Fractional macro weights describe ensemble density; physical electron shot-noise fidelity is not established.',
          'No complete boundary/external-work ledger is available, so no total-energy conservation claim is made.',
          'The existing3D-trajectory/height-averaged2D-Poisson, conditional voltage, mirrored lower thermal state, prescribed',
          'magnetic fields and omitted collisions/gas/ions/phase-change limits remain. Full60ns source history was not tested.', '',
          'Artifacts', '',
          f'Runtime and reproducible scripts: {ROOT}',f'Results and provenance: {OUT}',
          'run-summary.json, analysis.json, native-equivalence.json, provenance.json, seed-metrics.csv, expected-metrics.csv,',
          'patch-metrics.csv, seeds.csv, events.csv, events.bin, build.log, run.log and source-sampling-report.txt.', '']
(OUT/'source-sampling-report.txt').write_text('\n'.join(lines))
with (OUT/'seeds.csv').open('w') as stream:
    stream.write('seed_index,offset_decimal,offset_hex\n')
    for s in range(N):
        v=(s*0x9e3779b97f4a7c15)&((1<<64)-1)
        stream.write(f'{s},{v},0x{v:016x}\n')
paths=list(ROOT.iterdir())+[OUT/n for n in ['run-summary.json','analysis.json','native-equivalence.json','seed-metrics.csv','expected-metrics.csv','patch-metrics.csv','seeds.csv','events.csv','events.bin','build.log','run.log','source-sampling-report.txt']]
inputs=[Path('/home/treefor/sim-workspace/powerflow/outputs/heating-history-full.bin'),
        Path('/home/treefor/sim-workspace/powerflow/outputs/pic-profile/particle-audit/copied-checkpoint/restart.meta'),
        Path('/home/treefor/sim-workspace/powerflow/outputs/pic-profile/particle-audit/copied-checkpoint/checkpoint_400040.bin'),
        Path('/home/treefor/sim-workspace/powerflow/outputs/pic-profile/particle-audit/audit-report.txt'),
        Path('/home/treefor/sim-workspace/powerflow/outputs/pic-profile/particle-audit/inventory.json')]
def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
provenance=dict(compiler=subprocess.check_output(['/usr/sbin/c++','--version'],text=True).splitlines()[0],
    flags='-std=c++20 -O3 -DNDEBUG -Wall -Wextra -Werror=return-type',cpu=6,nice=15,threads=1,
    testing_lease='/home/treefor/sim-workspace/powerflow/outputs/pic-profile/testing-core.lock',
    linked_libraries=subprocess.check_output(['ldd',str(ROOT/'validate')],text=True),
    inputs_sha256={str(p):digest(p) for p in inputs},files_sha256={str(p):digest(p) for p in paths if p.is_file()})
equivalence=json.loads((OUT/'native-equivalence.json').read_text())
assert provenance['inputs_sha256']==equivalence['inputs_sha256_before']
assert all(digest(Path(path))==value for path,value in equivalence['frozen_sha256'].items())
sealed={r['path']:r['sha256'] for r in json.loads(Path('/home/treefor/projects/runtimes/mfem-gpu-pic-20261029/sealed-files.json').read_text())}
assert all(sealed[path]==value for path,value in equivalence['frozen_sha256'].items())
provenance['input_and_frozen_source_before_after_hashes_equal']=True
provenance['native_source_hashes_match_sealed_manifest']=True
assert not any(s in provenance['linked_libraries'].lower() for s in ['cuda','cudss','mpi','mfem','gomp','openblas'])
(OUT/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
print('Analysis passed numerical and ensemble guards; report:',OUT/'source-sampling-report.txt')
