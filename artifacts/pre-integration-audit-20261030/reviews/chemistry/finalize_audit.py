from pathlib import Path
import os,json,hashlib,csv,math
assert os.sched_getaffinity(0)=={7} and os.getpriority(os.PRIO_PROCESS,0)==15
O=Path(__file__).parent;P=Path('/home/treefor/projects/runtimes/pic-generic-water-bridge-20261029');W=Path('/home/treefor/Documents/Codex/2026-10-28/we')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(P/'source/generic_water_bridge.hpp')=='b16d517d0d2cb9aa5da70dde60ddf7a660804ca62018058a3295d8f99189b191'
assert sha(P/'source/generic_water_native_adapter.hpp')=='307244054950c827c57ca138ced3ed19fae9017f9c7e31177de079c91045e838'
D=P/'products/data/collisions/h2o-gas-20261029-elastic.json';data=json.loads(D.read_text());assert sha(D)=='6ea91f535096beafc19b3d69029597addb10909a60c1fae5064417f0df38b140'
G=str(P/'source/generic_water_bridge.hpp');A=str(P/'source/generic_water_native_adapter.hpp');N=str(P/'native/source/src/mfem_pic_arena.cpp');S=str(P/'inputs/collision_scheduler.hpp');B=str(P/'products/tools/build_h2o_outer_valence_extensions.py')
findings=[]
def add(id,severity,kind,title,detail,remedy,refs):findings.append(dict(id=id,severity=severity,classification=kind,title=title,finding=detail,remediation=remedy,evidence=refs))
add('PHY-01','P1','missing_physics_coverage','Cold/secondary-electron transport is incomplete','No elastic event kernel below1keV, no rotation/vibrational event channels, one selected electronic state only. At20eV the unapplied elastic integral is1.36e-19m2 versus selected parent+A1B1 sum5.466014492753623e-21m2. Scalar/momentum/energy-loss observables are available but are not exclusive differential events.','Bind supported low-energy angle/internal-state kernels and consistent state/reference inventories, or quantify declared omitted-process bounds. Do not use momentum-transfer/effective/stopping cross sections as extra event rates.',[{'path':G,'lines':[88,95,119,131]},{'path':str(D),'json_pointer':'channels[id=H2O:elastic_integral_song2021]'}])
add('PHY-02','P1','missing_physics_coverage','Attachment and fragment-ion channels are not complete','Hminus is applied only6.25..6.75eV although its source yield is positive at8.5/12eV; Ominus is never an event. At8.5eV source Hminus/Ominus yields are1.26e-22/8e-24m2; at12eV8e-24/2.74e-23m2. Below1keV only parent ionization is sampled, while OH+/H+/O+/other partial yields exist without qualified exclusive completion.','Add source-qualified event multiplicities, all neutral/electron products and state/angle/internal closures. Treat yield curves as observables until that contract exists; do not synthesize branches by arbitrary yield normalization.',[{'path':G,'lines':[88,95]},{'path':str(D),'json_pointer':'channels[id=H2O:attachment_yield:Hminus,Ominus;H2O:ion_yield:*]'}, {'url':'https://arxiv.org/abs/1007.5144'}])
add('PHY-03','P1','missing_physics_coverage','High-energy coverage remains partial','At8084.470604111418eV2a1 contributes1.5808384508973464e-22m2 without v6 products; core/multiple/Auger and fragment-species continuation remain unknown. Mapped fraction0.8898794893 is only of the original four-valence BEB model, not all-shell physical coverage.','Retain positive/unknown ledgers and refuse full-mechanism admission. Add source-qualified missing rates/products and continuation, or make a quantitatively bounded explicit subset claim. Published v7 has only a conditional27% mapping of2a1 and still73%positive unmapped; it is not integrated or audited here.',[{'path':G,'lines':[115,135]},{'path':B,'lines':[18,40]},{'path':str(W/'outputs/h2o-collision-data-20261029/coordination-handoff-v7.json')}])
add('PHY-04','P1','unimplemented_integration','Native multichannel and OHv family admission are absent','Generic native proofs use one fixed outcome-specific registry under forced_kernel_integration_only. Current single-channel completion requires every SSA child to match one channel/model; Hminus needs five allowed product alternatives under one base rate. Latest full-time candidate remains single-channel and unqualified.','Introduce a separately versioned immutable channel-set/product-family schema, preserve old single-channel checks, bind every child data/kernel/model/clock/stage/prefix, and qualify commit/replay/rollback. Do not relabel forced tests or duplicate full rates.',[{'path':A,'lines':[9,14,27,35]},{'path':N,'lines':[885,888,1438,1444,1488,1508]}])
add('PHY-05','P1','unimplemented_spatial_physics','Spatially separated capture cannot be admitted','DEA rates remain in the race, but sample requires exact3D coincidence. Native also refuses net charge transfer between distinct reactant sites. Independent macro e/H2O pairs normally need not coincide; selected capture will abort, not provide a complete reacting history.','Specify/qualify a physical collision-site and compatible reaction-current/work path before nonlocal capture. Retain the refusal meanwhile. Do not drop positive capture pairs, skip selected events, or move ions to the electron site merely to pass.',[{'path':G,'lines':[158,166]},{'path':N,'lines':[867,891]},{'path':str(W/'outputs/pic-generic-water-bridge-20261029/native-summary.json')}])
add('PHY-06','P1_for_physical_EEDF','unqualified_product_models','Conservation does not validate the EEDF/angles','Parent/outer3 prescribe equal-energy back-to-back CM electron pairs; A1B1 uses fixed7.49eV loss and ground/isotropic breakup; DEA factorizes a narrow6.5eV fit; OHminus uses ground/isotropic endpoint; IAM is an incoherent O+2H free-atom approximation.','Retain conditional labels, compare differential/secondary/internal observables to primary data and bound model sensitivity. Replace spectra only with an explicitly sourced, conservatively admitted kernel. Do not infer thermalization/avalanche from scalar rate or closure tests.',[{'path':str(P/'products/include/raikou/h2o_conditional_parent.h'),'lines':[38,51]},{'path':B,'lines':[31,36]},{'path':str(D),'json_pointer':'conditional_event_models[*].assumptions'},{'url':'https://physics.nist.gov/PhysRefData/Ionization/intro.html'},{'url':'https://nvlpubs.nist.gov/nistpubs/NSRDS/NIST.NSRDS.64-2023.pdf'}])
add('PHY-07','P2','confirmed_model_join_artifact','Selected collision sum jumps at1keV','Adjacent binary64 energies around1000eV give selected sums4.32e-21 then1.3102231042042155e-20m2, ratio3.0329238523. Parent alone jumps4.32e-21 to5.476053113888497e-21; IAM/fragment branches activate. No overlap double count was found, but continuity/physical transport is unqualified.','Obtain compatible supported low-energy/angular coverage and assess model matching/sensitivity. Do not hide the seam by arbitrary rescaling, endpoint holds or smoothing.',[{'path':G,'lines':[45,62,88,95]},{'path':str(O/'cap-probes.json'),'json_pointer':'one_keV_join'}])
add('PHY-08','P2','weighted_macro_model_validation_gap','Physical first moments do not establish stochastic convergence','With k=sigma*g/V, cohort hazard=max(Wp,Wn)k and event weight=min(Wp,Wn) preserve WpWn*k. Frozen-population count variance scales as WpWn*min(Wp,Wn)k, so weights/decomposition alter correlations and noise even at fixed mean physical source.','Converge macro weights/target carrier representation, cell size and interval with fixed physical source/neutral inventory. Report density, free-electron EEDF, species charge/current, ionization/attachment counts and variance; do not claim physical breakdown-delay fluctuations from first-moment tests.',[{'path':S,'lines':[202,218,289,306]},{'path':str(W/'outputs/pic-generic-water-bridge-20261029/clock-ensembles.json')}])
add('PHY-09','P2','confirmed_component_resource_limit','Coverage evidence cap is below a broad kinetic population','Default256 distinct relative energies refuses257; hard maximum1024 admits1024 and refuses1025. The old pre-kick fixture has46726 distinct stationary-H2O relative energies. This is an honest bounded refusal, not a physical zero or current physics-rate bug.','Fork a bounded streaming/hash/weighted-coverage evidence path that retains every supported/positive/unknown contribution, without changing rates or omitting pairs; qualify resource equivalence. Do not simply suppress the cap.',[{'path':G,'lines':[137,156]},{'path':str(O/'cap-probes.json')}])
add('PHY-10','P2','experimental_calibration_and_source_model_limit','RGA does not fix local H2O density/composition','The current recovered50-channel figure is published vector geometry, not raw calibrated neutral partial pressures. At assumed5e-5Torr and300K ideal-gas total density is1.6094166e12cm^-3; water is x_H2O times that unknown fraction. Desorbed oxide-water inventory/ground-state effusion is a declared surrogate scenario, not measured SS304 coverage.','Bind p,T,x species and sensitivity as explicit background scenarios or recover calibration/cracking matrix; use finite transported desorption inventory and gas internal-state assumptions. Keep microscopic source/binding/translation energy debits distinct and avoid normalizing RGA peaks into neutral fractions.',[{'path':str(W/'outputs/rga-channel-recovery/figure_2_3_rga_recovery.json'),'json_pointer':'limitations;physical_neutral_composition_ready'},{'path':str(W/'outputs/pic-electrode-heat-feedback-20261029/source-snapshot-v2/water-scenario-v2.json')}])
add('PHY-11','P2','model_resolution_and_unit_admission','Cell volume, density and field resolution remain host-owned','Rates use m2,m/s,m3 once. Generic fixtures use authored1e-12m3, not a verified full foil collision grid. Input positions/cell volume/native epoch are not authenticated by the sampler. Height-averaged field charge and a midplane diagnostic slab are not alternate gas populations.','Bind actual physical3D vacuum cell volumes and native private census/phase; keep all species number densities cm^-3, net charge separately, and converge collision-cell/Debye/plasma/gyro/coupling scales. Do not substitute slab volume or apply height twice.',[{'path':S,'lines':[202,208]},{'path':G,'lines':[71,72,137,167]},{'path':str(W/'outputs/water-breakdown-tonight-20261029/coupling-port-requirements.json'),'json_pointer':'dimensions;screening_observables'}])
(O/'findings.json').write_text(json.dumps(findings,indent=2)+'\n')
with (O/'findings.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=['id','severity','classification','title','finding','remediation']);w.writeheader();w.writerows({k:r[k] for k in w.fieldnames} for r in findings)
contract=json.loads((W/'outputs/pic-generic-water-bridge-20261029/mechanism-contract.json').read_text())
with (O/'channel-coverage.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['channel','model','conditional_support_eV','primary_capture','product_alternatives','electron_net_gain','scope'])
 for ch in contract['channels']:
  gain=-1 if ch['captures_primary'] else 0 if ch['model'] in ['A1B1_fixed_loss_dipole_ground_v1','IAM_NIST64_piecewise_mu_elastic_v1'] else 1
  w.writerow([ch['channel'],ch['model'],str(ch['data_support_eV']),ch['captures_primary'],json.dumps(ch['allowed_product_multisets']),gain,'conditional subset only; parent excludes exact1keV'])
primary=[{'title':'NIST Electron-Impact Ionization Database, Introduction','url':'https://physics.nist.gov/PhysRefData/Ionization/intro.html','fresh_web_read':True,'supports':'BEB supplies orbital/total ionization; it does not supply H2O exclusive branching or a H2O secondary spectrum; stated typical accuracy through about1keV is not an8keV guarantee.'},{'title':'Murakami et al., Single and multiple electron removal processes in proton–water-molecule collisions (2012)','url':'https://www.yorku.ca/marko/IonAtom/PRA85_052704.pdf','doi':'10.1103/PhysRevA.85.052704','fresh_web_read':True,'supports':'Eqs27–30 high-impact-energy single-vacancy branching; proton collision accuracy and applicability do not calibrate the present electron-impact composite model.'},{'title':'Ram, Prabhudesai and Krishnakumar, Dissociative Electron Attachment to Polyatomic Molecules–I: Water','url':'https://arxiv.org/abs/1007.5144','fresh_web_read':True,'supports':'Primary DEA observations include6.5,8.5,12eV resonances and Hminus/Ominus; a6.25..6.75eV Hminus fit is narrower than the overall attachment physics.'},{'title':'NIST Electron Elastic-Scattering Cross-Section Database v5 User Guide (2023)','url':'https://nvlpubs.nist.gov/nistpubs/NSRDS/NIST.NSRDS.64-2023.pdf','doi':'10.6028/NIST.NSRDS.64-2023','fresh_web_read':True,'supports':'Relativistic free-atom DCS source, not molecular-water differential validation.'}]
for name in ['song2021','munoz2007','hargreaves2012']:
 r=data['sources'][name];primary.append({'title':r.get('title',name),'url':r['url'],'doi':r.get('doi'),'fresh_web_read':False,'basis':'pinned v6 source/model provenance; fresh full-text fetch failed or timed out; not claimed newly reviewed'})
(O/'primary-sources.json').write_text(json.dumps(primary,indent=2)+'\n')
pressure=5e-5*133.32236842105263;density=pressure/(1.380649e-23*300)/1e6
(O/'density-scenario.json').write_text(json.dumps({'pressure_Torr_assumed':5e-5,'temperature_K_assumed':300,'ideal_gas_total_density_cm_minus3':density,'H2O_density_cm_minus3':'x_H2O * total_density; x_H2O not calibrated by recovered figure','local_desorption_density_inferred':False},indent=2)+'\n')
inputs=[P/'source/generic_water_bridge.hpp',P/'source/generic_water_native_adapter.hpp',P/'inputs/collision_scheduler.hpp',P/'inputs/water_parent_bridge.hpp',P/'products-source.tar.gz',D,P/'products/data/collisions/h2o-missing-rate-ledger-v6.json',P/'native/source/src/mfem_pic_arena.cpp',P/'inputs/native-source.tar.gz',W/'outputs/pic-generic-water-bridge-20261029/handoff.json',W/'outputs/pic-generic-water-bridge-20261029/validation.json',W/'outputs/pic-generic-water-bridge-20261029/coverage.json',W/'outputs/pic-generic-water-bridge-20261029/mechanism-contract.json',W/'outputs/pic-profile/em-field-consistency-audit/alignment/particles-with-fields.csv',W/'outputs/pic-profile/em-field-consistency-audit/alignment/decode-validation.json',W/'outputs/rga-channel-recovery/figure_2_3_rga_recovery.json',W/'outputs/pic-electrode-heat-feedback-20261029/source-snapshot-v2/water-scenario-v2.json',W/'outputs/h2o-collision-data-20261029/coordination-handoff-v7.json',W/'outputs/pic-profile/reaction-transactions-20261029/native-full-time/SELF-AUDIT.txt']
(O/'input-manifest.json').write_text(json.dumps({str(p):{'sha256':sha(p),'bytes':p.stat().st_size} for p in inputs},indent=2)+'\n')
(O/'physics-audit.txt').write_text('''Pre-integration physics/reaction audit — 2026-10-29, US/Mountain

Full physical H2O breakdown/EEDF integration is not ready. The frozen generic-v6
bridge is an internally conservative conditional subset with explicit refusals;
its tests do not establish a complete reaction mechanism, calibrated product
spectra, spatial capture law, or converged kinetic history. No additional wrong-
sign, duplicate-debit or reaction arithmetic error was confirmed in the eight
models' declared domain. Important omissions and integration/resource limits
below remain. Other native full-time defects are on a separate documented hold.

Review boundary and evidence
The generic bridge is unchanged: source SHA b16d517d0d2cb9aa5da70dde60ddf7a660804ca62018058a3295d8f99189b191,
native adapter307244054950c827c57ca138ced3ed19fae9017f9c7e31177de079c91045e838,
v6 product dataset6ea91f535096beafc19b3d69029597addb10909a60c1fae5064417f0df38b140.
Runtime /home/treefor/projects/runtimes/pic-generic-water-bridge-20261029.
findings.json gives exact absolute source paths/line references and remediation;
findings.csv is the mergeable severity list; channel-coverage.csv is the model map.
input-manifest.json binds every local source/evidence input, including the old
fixture, RGA metadata and source scenario. Native forced admissions use the
immutable v6.2 source5fa3a59f…/archived14f1e2…; their48 records are24 accepted
restart/rollback controls plus24 required nonlocal-capture refusals. No new native,
field, coupled, source or physical-time solve was run for this audit.

The only new execution is the bounded read-only cap_probe.cpp: direct scalar
queries at adjacent1keV energies, small indexed candidate/evidence-cap fixtures,
and relative-energy inventory of an already qualified CSV using its named
pre-kick u columns. CPU7/nice15/one thread/frozen userspace/lease; no GPU or
production edits. It is not a new water EEDF. Existing512 actual-kernel outcomes
and140-digit conservation evidence are reused, not presented as new validation.

P1: coverage gaps block a complete reaction/transport claim (PHY-01..03)
Eight models do not span the available observables. Below1keV there is no elastic
angle kernel, rotation or vibrational state evolution. Only one dissociative
excited molecular state is active. At20eV the compiled/recommended elastic
integral1.36e-19m2 is24.88 times the selected parent+A1B1 sum5.466e-21m2. This is a
cross-section comparison, not an exclusive-event probability or a claim about
absolute collisionality at an unknown density. Positive low-energy scalar data
must not be silently treated as absence of scattering. Momentum-transfer,
effective-momentum and stopping moments cannot simply be added as extra events.

The Hminus sampler is a6.25..6.75eV comparison even though the source Hminus yield
is positive at8.5 and12eV; Ominus has no active product kernel. At8.5eV the source
Hminus/Ominus yields are1.26e-22/8e-24m2; at12eV8e-24/2.74e-23m2. OHminus alone is
not the complete attachment set. Below1keV the OH+/H+/O+/H2+/atomic O(2+) partial-yield
observables have no corresponding qualified exclusive completion in this bridge.
Do not normalize those yields into event branches without neutral completion and
multiplicity evidence. The primary DEA study covers several resonances and both
Hminus/Ominus; the narrow selected fit does not replace them:
https://arxiv.org/abs/1007.5144 .

At8084.470604111418eV the v6 positive unmapped2a1 rate is1.5808384508973464e-22m2.
The88.99% mapped fraction is only of the four-valence BEB model. Core/multiple/
Auger products and reactions of OH/H/ions/anions are not complete; the target
validator accepts H2O only. Fragment transport as a species is not its chemical
continuation. Assess ion-neutral/charge-exchange, electron detachment/recombination
and field-dependent state destruction on the chosen density/time/field scales;
this audit does not invent rates or declare every omitted process dominant.
Published v7 reports two conditional2a1 branches mapping27% of that hole while73%
remains positive/unmapped; it is a separate release, not integrated or numerically
reviewed here, and does not resolve the v6/full-mechanism readiness flag.

BEB is an orbital/total ionization model, not exclusive branching or a H2O
secondary spectrum. The official introduction's typical accuracy statement is
through roughly1keV for many targets, not an8keV guarantee:
https://physics.nist.gov/PhysRefData/Ionization/intro.html .
The adopted high-energy branching lineage appears in Eqs27–30 of a proton-water
study. Its asymptotic single-vacancy argument is a conditional basis for the
composite electron model; the proton collision accuracy/range is not an electron
calibration. The same paper discusses failure of fixed branching when multiple
removal matters. No proton impact rate was imported:
https://www.yorku.ca/marko/IonAtom/PRA85_052704.pdf .

P1: native mechanism and spatial admission are incomplete (PHY-04..05)
The24 successful native comparisons register the already selected product
multiset and use forced_kernel_integration_only. They do not establish one
native multichannel SSA exposure. The existing completion binds one channel,
model and fixed multiset. Hminus needs5 explicit OHv alternatives under ONE base
rate, not5 full-rate aliases or sequential full-duration clocks. A new immutable
channel-set schema must preserve all current per-child data/kernel/model/stage/
clock/prefix and strictly ordered local-time validation. Keep the single-channel
schema intact. New validator feature work is paused at its saved design snapshot.

Capture at separated3D sites has no qualified current/work law. Such rates stay
in the candidate race; selection/sample refuses rather than relocating ions or
silently deleting events. A continuous independently sampled e/H2O population
will not generally be exactly co-located. A physical collision-site/current path
must be specified before this becomes a usable complete capture operator. Full
3D coincidence is required even when height-averaged2D nodal charge could match.
The native guard permits other reactions only when charge is conserved at each
reactant site, or when both sites truly coincide. The generic placement follows
that rule: scattered primary at primary, other products at target. It does not
invent a current. Native Q1 assembly differences up to4.8148e-35C are within the
separately derived floating-operation bound; this is not unexplained physical
charge transport, and no native guard was widened.

P1 for physical EEDF interpretation: prescribed spectra (PHY-06)
The parent and outer3 kernels impose equal-energy opposite CM electron pairs;
this fixes energy sharing rather than sampling measured H2O joint spectra.
A1B1 substitutes a fixed7.49eV loss and ground/isotropic breakup for a broad
continuum. Hminus factorizes a6.5eV fit; OHminus uses a ground/isotropic endpoint.
Conservation is necessary but does not validate the resulting EEDF, pitch angles,
secondary multiplication or cooling. IAM elastic uses incoherent O+2H free-atom
DCS with exact molecular recoil. The official source describes free atoms;
molecular interference/bonding and transport validation remain separate:
https://nvlpubs.nist.gov/nistpubs/NSRDS/NIST.NSRDS.64-2023.pdf .
Do not interpret the numerical angular-integral consistency or mesh difference
(up to1.91% from rounded atomic totals) as measured molecular accuracy. Conditional
product-model sensitivities and differential/secondary data are needed before a
physical kinetic/EEDF claim. Fresh fetches of the frozen Hargreaves/Munoz/Song
PDF URLs timed out; their precise conditional assumptions here are verified from
the pinned model metadata, not asserted as a new full-text review.

P2 confirmed model seam: no overlap, but a3.03-fold jump (PHY-07)
cap-probes.json records selected sums at the adjacent binary64 energies:
999.9999999999999eV:4.3199999999999996e-21m2;
1000eV:1.3102231042042155e-20m2;
1000.0000000000001eV:1.3102231042042152e-20m2.
The exact1keV policy prevents simultaneous old-parent/outer3 double counting.
Nevertheless, parent alone changes4.32e-21→5.4760531e-21m2 and fragment/IAM
channels activate, giving a3.03292385 subset-sum ratio. This seam follows the
incomplete support/model selection, not a physical threshold. It can affect
transport as particles cross1keV. Add compatible supported coverage and model
matching/sensitivity; do not conceal it with unsourced rescaling or smoothing.

No overlap was found within the declared high-energy three-branch construction:
orbital weights are parent[0,.08,1,1],OHp_H[0,.70,0,0],Hp_OH[0,.22,0,0], and
unmapped[1,0,0,0]. Their sum is[1,1,1,1]; the positive unmapped term is not sampled.
Do not also add full BEB, old partial-ion/parent curves, aggregate inelastic or
neutral-OH yield on top. A1B1 is one selected molecular event; the neutral yield
is not an independent extra dissociation. There is no duplicated Hminus base rate.

P2 weighting/convergence and resource scope (PHY-08..09)
For k=sigma*g/V, hazard lambda=max(Wp,Wn)*k and cohort w=min(Wp,Wn) give mean
lambda*w=Wp*Wn*k. This is the tested physical first moment. In a frozen-population
Poisson-clock comparison the weighted count variance rate is lambda*w^2=
Wp*Wn*min(Wp,Wn)*k. Thus changing macro weights or exact-count carrier decomposition
changes noise/correlations. It is not an unweighted physical shot-noise model;
weights below1 can smooth, above1 can amplify. No macro/cell/coupling convergence
or physical breakdown-time distribution has been established. Keep source,
inventory and physical volumes fixed while converging these quantities.

The evidence ledger also limits broad populations: default256 distinct relative
energies refuses257; maximum1024 admits1024 and refuses1025. The frozen old
pre-kick fixture has46726 distinct stationary-H2O relative energies (maximum
8100.0814eV), so the current small component cannot simply take that whole kinetic
population. This is an explicit refusal, not a rate bias. A bounded streaming,
weighted-coverage/hash evidence path needs equivalence/resource qualification
before broad-population use; merely dropping records or enlarging caps is not a
physics/resource solution. No solver or new trajectory was used to show this.

P2 density/source/reference scope (PHY-10..11)
The recovered current RGA artifact is a50-channel published vector plot, not raw
acquisition or calibrated neutral partial pressure/composition. Ionizer cracking,
sensitivities, electron energy, timing and temperature were not recovered from
that figure. A peak at a fragment mass is not its neutral mole fraction or a
preexisting plasma ion count. Under the explicit illustrative p=5e-5Torr,T=300K
ideal-gas assumption, total density is1.6094e12cm^-3 and nH2O=xH2O*n_total with
xH2O unknown. A100%water case may be an explicit scenario, not a recovered fact.
Local nonequilibrium desorption density cannot be read from the RGA figure.

The pinned water-source scenario uses an Fe2O3-like oxide surrogate on nominal
SS304, assumed coverage0.33 and finite patch inventory, with ground-state effusive
water. It is not a measured shot coverage, oxide composition or gas rovibrational
population. Retain the declared scenario/range sensitivities and joint surface/
work-function assumption. Its gas H2O reference0 matches this reaction catalog;
adsorption binding is paid by the source reservoir and is not paid again as a
fragment or ionization cost. A1B1's7.49eV intermediate loss is not an extra heat
debit: its5.10eV final chemical reference and remaining fragment motion already
close the selected model. Ionization references remain in products; recombination
or collection must own their eventual energy transfer explicitly.

No mass/reference or unit/sign inconsistency was found inside this reduced
catalog: negative ions include the captured electron mass, positive ions remove
it, OH vibrational energy changes internal/reference bookkeeping and state mass
consistently. Implied0.75/1.83eV affinities are source-rounded model conventions,
not precision species masses. The140-digit review confirms K+internal+reference
and momentum in its admitted numerical domain; extreme bulk conditioning may
correctly refuse. Cross sections arem2, fluxm/s, cellsm3, counts represented
physical particles. Display number densities ascm^-3 (m^-3 divided by1e6).
The host must bind actual3D vacuum cell volume, private census/epoch and field
geometry. Do not use a plotted midplane slab as an additional gas population or
apply the25mm invariant height twice. Free-electron density/EEDF must be separate
from negative-ion density and net charge. Plasma shielding/residence requires
species current, field and resolution evidence, not electron-cloud appearance.

Current native full-time status is a separate hold
The latest owner-provided audit-candidate-002 source is eccf8960… (archiveaf0824d8…),
still one channel/fixed product set, zero qualified joint steps. Its SELF-AUDIT.txt
records standalone ownership/receipt/restart/configuration issues and history
bounds; those are not independently re-executed here or silently attributed to
the older qualified forced-reaction source. The owner also reports only368
legacy+39 manufactured storage controls. Generic reaction proofs cannot qualify
that new owner, field/heat/source/collection integration or a60ns history. Root's
broader audit must resolve that hold before joint execution. No result here
identifies water ionization as the cause or cure of the old density halo.

The mechanism feature remains paused. The bridge, sealed dependencies, native
sources and production inputs remain unchanged. Remediation proposals are future
work and do not authorize dropping positive missing rates or relaxing guards.
''')
summary={'status':'preintegration_physics_audit_complete_full_mechanism_not_ready','findings':len(findings),'confirmed_additional_reaction_arithmetic_bugs':0,'confirmed_artifacts':['one_keV_selected_subset_rate_seam','broad_energy_coverage_evidence_cap'],'new_execution':'bounded scalar/candidate/read-only CSV probes only','no_new_native_field_or_joint_solve':True,'generic_bridge_unchanged':True,'native_sources_unchanged':True,'mechanism_feature_paused':True,'full_water_mechanism_ready':False,'physical_EEDF_validated':False,'latest_native_full_time_candidate_qualified':False}
(O/'audit-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
h={'status':summary['status'],'report':str(O/'physics-audit.txt'),'findings':str(O/'findings.json'),'findings_csv':str(O/'findings.csv'),'coverage_table':str(O/'channel-coverage.csv'),'primary_sources':str(O/'primary-sources.json'),'narrow_probes':str(O/'cap-probes.json'),'paused_feature':str(W/'outputs/pic-water-mechanism-validator-20261029/PAUSED-DESIGN.json'),'generic_source_sha256':'b16d517d0d2cb9aa5da70dde60ddf7a660804ca62018058a3295d8f99189b191','generic_source_archive_sha256':'90ccb83772920f19121985097684535075eed8a8eaded54d5324e962928a07eb','full_mechanism_ready':False,'CPU7_lease_released_on_completion':True}
(O/'handoff.json').write_text(json.dumps(h,indent=2)+'\n')
(O/'output-manifest.json').write_text(json.dumps({f.name:{'sha256':sha(f),'bytes':f.stat().st_size} for f in sorted(O.iterdir()) if f.is_file() and f.name!='output-manifest.json'},indent=2)+'\n')
print(json.dumps(h))
