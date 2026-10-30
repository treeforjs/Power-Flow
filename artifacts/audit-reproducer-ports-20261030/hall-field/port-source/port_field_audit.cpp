// Port of the pre-integration field audit counterexamples (findings A09 and A10)
// to the corrected Hall field owner, namespace hall_corrective, commit 8dddec0.
//
// Every input is taken from the audit probe
//   Power-Flow/artifacts/pre-integration-audit-20261030/reviews/field/probe_operators.cpp
// default axial::Mesh, sigma = 1/(71.9e-8) S/m in foil cells, t = 7.4978805549199997e-10 s,
// dt = 1e-14 s, I = 6012.0713604700004 A, manufactured -e*1e9 point charge at field
// position (.00021,-.00017,.001) m, B[bx(22,22,kshort+1)] += 0.01 T, zero plasma current.
// The unchanged axial::Fields operator is re-executed first (original behavior), then the
// identical state, charge and perturbation are passed through the corrected API.
//
// Scope: no native PreparedArena, SSA, source population, ID allocation, joint native
// publication or physical run. Census rows below are constructed by this port; they are
// the corrected API's census type (native::ArenaRateEntry), not rows produced by a native arena.
#include "root_api.hpp"
#include "../hall-finite-foil-em-20261029/axial_particles.hpp"
#include <iostream>
#include <optional>

namespace h=hall_corrective;namespace n=raikou::mfem_backend;namespace hd=hall_driver;
using h::Json;using axial::Vec;

namespace {
constexpr double t0=7.4978805549199997e-10,dt=1e-14,I0=6012.0713604700004;
constexpr double audit_weight=1e9;
const axial::V3 audit_point{.00021,-.00017,.001};
constexpr long double machine_eps=std::numeric_limits<double>::epsilon();

template<class F>Json attempt(F&&f){
    try{f();return Json{{"accepted",true}};}
    catch(const std::exception&e){return Json{{"accepted",false},{"rejection",e.what()}};}
}
// Same option construction as the corrective owner test.
n::GapPicOptions geometry(const axial::Mesh&g){n::GapPicOptions o;o.nx=o.ny=g.nx;o.width_m=o.gap_m=g.side;o.depth_m=g.depth;
    o.finite_foils=o.finite_axial_extent=true;o.short_z_m=-g.depth/2;o.feed_z_m=g.depth/2;
    o.lower_foil={g.side/2-.0025,g.side/2+.0025,g.side/2-.0027,g.side/2-.0025};
    o.upper_foil={g.side/2-.0025,g.side/2+.0025,g.side/2+.0025,g.side/2+.0027};return o;}
bool conductor_node(const axial::Mesh&g,int k){return g.node_body[k]>=0||g.short_node[k];}
h::FreeCharge no_free_charge(const axial::Mesh&g){return {h::WordVec(g.nn),Vec(g.nn)};}
// Declaration rule used by the corrective tests: conductor/PEC and exterior banks carry
// the supplied field's own charge on their supported nodes. Free charge is supplied.
h::ChargeLedger dc_declaration(const axial::Mesh&g,const axial::State&s,h::FreeCharge supplied){
    h::ChargeLedger q;q.particles=std::move(supplied);q.conductor_C.assign(g.nn,h::Words());q.exterior_C.assign(g.nn,h::Words());
    q.conductor_declaration=h::sha256("audit port: DC induced conductor/PEC charge declared from the supplied field");
    q.exterior_declaration=h::sha256("audit port: DC exterior termination charge declared from the supplied field");
    const auto Q=h::field_charge(g,s);
    for(int k=0;k<g.nn;++k){if(conductor_node(g,k))q.conductor_C[k]=Q[k];else if(g.outer_node[k])q.exterior_C[k]=Q[k];}
    return q;}
// Most permissive declaration the API admits: every supported conductor/exterior node
// absorbs field minus free charge. Only vacuum interior nodes remain constrained.
h::ChargeLedger absorbing_declaration(const axial::Mesh&g,const axial::State&s,h::FreeCharge supplied){
    auto q=dc_declaration(g,s,std::move(supplied));const auto Q=h::field_charge(g,s);
    q.conductor_declaration=h::sha256("audit port: conductor bank absorbs field minus free charge");
    q.exterior_declaration=h::sha256("audit port: exterior bank absorbs field minus free charge");
    for(int k=0;k<g.nn;++k){if(conductor_node(g,k))q.conductor_C[k]=Q[k]-q.particles.net_C[k];
        else if(g.outer_node[k])q.exterior_C[k]=Q[k]-q.particles.net_C[k];}
    return q;}
const char*node_class(const axial::Mesh&g,int k){return conductor_node(g,k)?"conductor_or_short":g.outer_node[k]?"em_exterior":"vacuum_interior";}
// Independent per-node recomputation of the absolute Gauss predicate in admit_constraints.
Json gauss_scan(const axial::Mesh&g,const axial::State&s,const h::ChargeLedger&q,double solver_bound){
    const auto E=h::words(s.E,s.E_low);std::vector<long double>activity(g.nn);
    for(int e=0;e<g.ne;++e){const auto&edge=g.edges[e];const long double a=h::magnitude(E[e])*g.M[e]/edge.length;activity[edge.a]+=a;activity[edge.b]+=a;}
    const auto Q=h::field_charge(g,s);long double max_defect=0,bound_at_max=0,max_ratio=0;int node_max=-1,node_ratio=-1,first=-1;
    std::size_t violations=0,vacuum_violations=0;double first_defect=0,first_bound=0;
    for(int k=0;k<g.nn;++k){const auto d=Q[k]-q.particles.net_C[k]-q.conductor_C[k]-q.exterior_C[k];
        const long double bound=solver_bound+256*machine_eps*(activity[k]+q.particles.gross_C[k]+h::magnitude(q.conductor_C[k])+h::magnitude(q.exterior_C[k]));
        const long double ad=std::abs(d.wide());
        if(!(ad<=bound)){++violations;if(!conductor_node(g,k)&&!g.outer_node[k])++vacuum_violations;
            if(first<0){first=k;first_defect=double(d.wide());first_bound=double(bound);}}
        if(ad>max_defect){max_defect=ad;bound_at_max=bound;node_max=k;}
        if(ad/bound>max_ratio){max_ratio=ad/bound;node_ratio=k;}}
    return {{"max_abs_defect_C",double(max_defect)},{"bound_at_max_defect_C",double(bound_at_max)},{"node_of_max_defect",node_max},
        {"node_of_max_defect_class",node_max>=0?node_class(g,node_max):"none"},{"max_defect_over_bound",double(max_ratio)},{"node_of_max_ratio",node_ratio},
        {"nodes_violating_bound",violations},{"vacuum_interior_nodes_violating_bound",vacuum_violations},
        {"first_violating_node",first},{"first_violation_defect_C",first_defect},{"first_violation_bound_C",first_bound}};}
// Independent per-cell recomputation of the absolute divB predicate in admit_constraints.
Json divb_scan(const axial::Mesh&g,const axial::State&s){
    const auto B=h::words(s.B,s.B_low);long double max_div=0,bound_at_max=0,max_ratio=0;std::size_t violations=0;int cell_max=-1;
    for(int k=0;k<g.nz;++k)for(int j=0;j<g.ny;++j)for(int i=0;i<g.nx;++i){h::Words d;long double scale=0;
        const std::array<std::tuple<int,int,double>,3>faces{{{g.bx(i+1,j,k),g.bx(i,j,k),g.dx},{g.by(i,j+1,k),g.by(i,j,k),g.dy},{g.bz(i,j,k+1),g.bz(i,j,k),g.dz}}};
        for(const auto&[a,b,hh]:faces){d=d+(B[a]-B[b])/hh;scale+=(h::magnitude(B[a])+h::magnitude(B[b]))/hh;}
        const long double bound=128*machine_eps*scale,ad=std::abs(d.wide());if(!(ad<=bound))++violations;
        if(ad>max_div){max_div=ad;bound_at_max=bound;cell_max=(k*g.ny+j)*g.nx+i;}if(bound>0&&ad/bound>max_ratio)max_ratio=ad/bound;}
    return {{"max_abs_divB_T_per_m",double(max_div)},{"bound_at_max_T_per_m",double(bound_at_max)},{"cell_of_max",cell_max},
        {"max_divB_over_bound",double(max_ratio)},{"cells_violating_bound",violations}};}
std::string mesh_digest(const axial::Mesh&g){return h::sha256(h::mesh_record(g).dump());}
hd::MotionBinding binding_for(const h::Owner&owner,const n::GapPicOptions&o,std::span<const n::ArenaRateEntry>rows){
    hd::MotionBinding b;b.phase=hd::ParticlePhase::full_time_position_and_proper_velocity;b.accepted_step=owner.revision();
    b.staged_arena_revision=1;b.begin_s=owner.time_s();b.dt_s=dt;b.end_s=owner.time_s()+dt;
    b.staged_arena_sha256=hd::canonical_census_sha256(rows);b.mesh_sha256=owner.motion_mesh_sha256();
    b.material_sha256=hd::material_binding_sha256();b.species_sha256=hd::text_sha256("audit port fixture: electron and proton carriers");
    b.common_clock_sha256=hd::text_sha256(Json{{"begin",b.begin_s},{"dt",dt},{"end",b.end_s},{"step",b.accepted_step},{"fixture","audit port"}}.dump());
    b.kinetic_policy_sha256=hd::kinetic_policy_sha256(hd::declared_kinetic_policy(o));return b;}
n::ArenaRateEntry carrier(std::uint64_t id,const n::GapPicOptions&o,const axial::V3&field_x,n::PicVector u,double charge,double mass,double weight,const std::string&tag){
    n::ArenaRateEntry row;row.packet_id="audit-port/packet/"+tag;row.state="audit-port/resident/"+tag;row.native.id=id;
    row.native.position_m={field_x[0]+o.width_m/2,field_x[1]+o.gap_m/2,field_x[2]};row.native.proper_velocity_m_s=u;
    row.native.charge_C=charge;row.native.mass_kg=mass;row.native.weight=weight;return row;}
Vec compensated_average(const Vec&a,const Vec&alo,const Vec&b,const Vec&blo){Vec out(a.size());
    for(std::size_t i=0;i<a.size();++i)out[i]=((h::Words(a[i],alo[i])+h::Words(b[i],blo[i]))*.5).rounded();return out;}
Json high_word_comparison(const Vec&old_hi,const Vec&new_hi,const Vec&authoritative){std::size_t differing=0;double max_difference=0;
    for(std::size_t i=0;i<authoritative.size();++i){const double hw=.5*(old_hi[i]+new_hi[i]);
        if(hw!=authoritative[i])++differing;max_difference=std::max(max_difference,std::abs(hw-authoritative[i]));}
    return {{"differing",differing},{"components",authoritative.size()},{"max_abs_difference",max_difference}};}
} // namespace

int main(int argc,char**argv){MPI_Init(&argc,&argv);int rc=0;try{
    h::need(argc==3,"usage: port_field_audit <shot_17133_current.csv> <output directory>");
    const std::filesystem::path trace=argv[1],out=argv[2];h::need(std::filesystem::is_directory(out),"output directory must exist");
    Json result{{"schema","hall-field-audit-reproducer-port-v1"},{"hall_api",h::api},{"trace",trace.string()}};

    // ---------------------------------------------------------------- original operator (probe inputs)
    axial::Mesh g;Vec sigma(g.cells.size());for(std::size_t c=0;c<sigma.size();++c)if(g.cells[c].body>=0)sigma[c]=1/(71.9e-8);
    axial::Fields fields(g,sigma,I0);const auto initial=fields.start(t0);
    auto divergence=[&](const axial::State&s){auto d=g.magnetic_divergence(s.B),low=g.magnetic_divergence(s.B_low);
        for(std::size_t i=0;i<d.size();++i)d[i]+=low[i];return axial::maximum(d);};
    const auto qfield=g.field_charge(initial.E);Vec qauthored(g.nn);axial::point_deposit(g,qauthored,audit_point,-pic_e*audit_weight);
    double initial_vacuum=0,authored_mismatch=0,authored_max=0;
    for(int k=0;k<g.nn;++k)if(g.node_body[k]<0&&!g.short_node[k]&&!g.outer_node[k]){
        initial_vacuum=std::max(initial_vacuum,std::abs(qfield[k]));
        authored_mismatch=std::max(authored_mismatch,std::abs(qfield[k]-qauthored[k]));
        authored_max=std::max(authored_max,std::abs(qauthored[k]));}
    const auto step=fields.prepare(initial,Vec(g.ne),sigma,dt,I0+.008);
    double mismatch_after=0;{auto hi=g.field_charge(step.first.E),lo=g.field_charge(step.first.E_low);
        for(int k=0;k<g.nn;++k)if(g.node_body[k]<0&&!g.short_node[k]&&!g.outer_node[k])mismatch_after=std::max(mismatch_after,std::abs(hi[k]+lo[k]-qauthored[k]));}
    std::size_t differing=0;double max_difference=0;
    for(int e=0;e<g.ne;++e){const double midpoint=.5*(initial.E[e]+step.first.E[e]);
        if(midpoint!=step.second.mid_E[e])++differing;max_difference=std::max(max_difference,std::abs(midpoint-step.second.mid_E[e]));}
    const int face=g.bx(22,22,g.kshort+1);auto perturbed=initial;perturbed.B[face]+=.01;
    Json original_divb{{"divB_before_T_per_m",divergence(perturbed)}};
    try{auto trial=fields.prepare(perturbed,Vec(g.ne),sigma,dt,I0);original_divb["accepted_by_Fields_prepare"]=true;
        original_divb["divB_after_T_per_m"]=divergence(trial.first);original_divb["gauss_increment_C"]=trial.second.gauss_increment_C;
        original_divb["energy_balance_J"]=trial.second.balance_J;}
    catch(const std::exception&e){original_divb["accepted_by_Fields_prepare"]=false;original_divb["rejection"]=e.what();}
    result["original_operator"]={{"operator","axial::Fields::prepare (unchanged since c2f81ab)"},
        {"initial_divB_T_per_m",divergence(initial)},{"initial_vacuum_field_charge_max_C",initial_vacuum},
        {"manufactured_authored_charge_max_C",authored_max},{"initial_field_vs_manufactured_authored_charge_max_C",authored_mismatch},
        {"A09_gauss",{{"api_operand_for_particle_charge",false},{"Fields_prepare_from_inconsistent_state_accepted",true},
            {"gauss_increment_C",step.second.gauss_increment_C},{"energy_balance_J",step.second.balance_J},
            {"field_vs_manufactured_charge_max_after_one_step_C",mismatch_after}}},
        {"A09_divB",original_divb},
        {"A10_midpoint",{{"current_A",I0+.008},{"different_edges",differing},{"edges",g.ne},{"max_abs_difference_V_per_m",max_difference},
            {"comparison","high-word average .5*(old.E+new.E) vs Maxwell receipt.mid_E"}}}};

    // ---------------------------------------------------------------- corrected model, same inputs
    auto model=std::make_shared<const h::Model>(h::MeshSpec{},trace);const auto&m=model->mesh();const auto o=geometry(m);
    const auto seed=model->dc_seed();const double solver_bound=model->initial_gauss_solve_bound_C();
    const bool same_mesh=mesh_digest(g)==mesh_digest(m);h::need(same_mesh,"corrected model mesh differs from the probe mesh; indices not portable");
    Json identity{{"probe_mesh_equals_model_mesh",same_mesh},{"probe_sigma_equals_model_initial_sigma",sigma==model->initial_conductivity()},
        {"probe_initial_equals_model_dc_seed",initial.E==seed.E&&initial.B==seed.B&&initial.time_s==seed.time_s},
        {"probe_time_equals_C19_start",t0==model->drive().start()},{"C19_start_current_A",model->drive().current(model->drive().start())},
        {"probe_current_A",I0},{"configuration_sha256",model->configuration_sha256()},{"initial_gauss_solver_bound_C",solver_bound}};
    {double sigma_diff=0;for(std::size_t c=0;c<sigma.size();++c)sigma_diff=std::max(sigma_diff,std::abs(sigma[c]-model->initial_conductivity()[c]));identity["max_abs_sigma_difference_S_per_m"]=sigma_diff;}
    {double e_diff=0;for(int e=0;e<m.ne;++e)e_diff=std::max(e_diff,std::abs(initial.E[e]-seed.E[e]));identity["max_abs_E_probe_minus_seed_V_per_m"]=e_diff;}
    result["fixture_identity"]=identity;
    const auto control=dc_declaration(m,initial,no_free_charge(m));

    // ---------------------------------------------------------------- A09 Gauss: inconsistent initial charge/field
    {
        Json r;const std::vector<n::ArenaRateEntry>census{carrier(1,o,audit_point,{0,0,0},-pic_e,pic_me,audit_weight,"A09-manufactured-electron")};
        const auto census_charge=h::charge_of_census(m,o,census);double deposit_parity=0;Json support=Json::array();
        for(int k=0;k<m.nn;++k){deposit_parity=std::max(deposit_parity,double(std::abs(census_charge.net_C[k].wide()-(long double)qauthored[k])));
            if(census_charge.net_C[k].hi!=0)support.push_back({{"node",k},{"class",node_class(m,k)},{"net_C",census_charge.net_C[k].hi}});}
        r["actual_census_rows"]=census.size();r["census_sha256"]=hd::canonical_census_sha256(census);
        r["census_deposit_vs_probe_point_deposit_max_abs_C"]=deposit_parity;r["deposit_support"]=support;
        const auto honest=dc_declaration(m,initial,census_charge),absorbing=absorbing_declaration(m,initial,census_charge);
        r["control_admit_field_without_free_charge"]=attempt([&]{model->admit(initial,control);});
        r["control_owner_empty_census"]=attempt([&]{h::Owner owner(model,initial,control,o,{});});
        r["Model_admit_census_charge_dc_declaration"]=attempt([&]{model->admit(initial,honest);});
        r["Model_admit_census_charge_absorbing_declaration"]=attempt([&]{model->admit(initial,absorbing);});
        r["Model_prepare_from_inconsistent_state"]=attempt([&]{(void)model->prepare(initial,honest,Vec(m.ne),honest.particles,sigma,dt);});
        r["Owner_actual_census_with_matching_ledger"]=attempt([&]{h::Owner owner(model,initial,honest,o,census);});
        r["Owner_actual_census_with_absorbing_ledger"]=attempt([&]{h::Owner owner(model,initial,absorbing,o,census);});
        r["Owner_actual_census_with_field_only_ledger"]=attempt([&]{h::Owner owner(model,initial,control,o,census);});
        r["Owner_empty_census_with_charged_ledger"]=attempt([&]{h::Owner owner(model,initial,honest,o,{});});
        r["gauss_scan_control"]=gauss_scan(m,initial,control,solver_bound);
        r["gauss_scan_dc_declaration"]=gauss_scan(m,initial,honest,solver_bound);
        r["gauss_scan_absorbing_declaration"]=gauss_scan(m,initial,absorbing,solver_bound);
        r["solver_bound_C"]=solver_bound;r["audit_mismatch_over_solver_bound"]=authored_mismatch/solver_bound;
        r["elementary_charge_over_solver_bound"]=pic_e/solver_bound;
        // Supplementary looseness probe: a co-located neutral electron/proton pair of weight W
        // raises gross charge and therefore the 256-epsilon arithmetic term at the same nodes.
        Json inflation=Json::array();
        for(double W:{1e21,1e22}){std::vector<n::ArenaRateEntry>rows{census[0],
                carrier(2,o,audit_point,{0,0,0},-pic_e,pic_me,W,"inflation-electron"),carrier(3,o,audit_point,{0,0,0},pic_e,hall::proton_mass,W,"inflation-proton")};
            Json item{{"pair_weight_each",W}};const auto q=dc_declaration(m,initial,h::charge_of_census(m,o,rows));
            item["gauss_scan"]=gauss_scan(m,initial,q,solver_bound);
            std::unique_ptr<h::Owner>owner;item["Owner_bootstrap"]=attempt([&]{owner=std::make_unique<h::Owner>(model,initial,q,o,rows);});
            if(owner)item["closed_proposal_from_admitted_state"]=attempt([&]{auto b=binding_for(*owner,o,rows);
                (void)h::prepare_closed_proposal(*owner,o,rows,b,h::sha256("audit port: unchanged closed source identity"));});
            inflation.push_back(item);}
        r["supplementary_gross_inflation_probe"]=inflation;
        result["A09_gauss"]=r;
    }

    // ---------------------------------------------------------------- A09 divB: 0.01 T single-face perturbation
    {
        Json r;auto bad=initial;bad.B[face]+=.01;r["face"]=face;r["divB_scan_unperturbed"]=divb_scan(m,initial);r["divB_scan_perturbed"]=divb_scan(m,bad);
        r["Model_admit"]=attempt([&]{model->admit(bad,control);});
        r["Model_prepare_from_perturbed_state"]=attempt([&]{(void)model->prepare(bad,control,Vec(m.ne),control.particles,sigma,dt);});
        r["Owner_bootstrap_empty_census"]=attempt([&]{h::Owner owner(model,bad,control,o,{});});
        {   const auto saved=out/"A09-divB-owner-checkpoint",edited=out/"A09-divB-owner-checkpoint-perturbed";
            h::Owner live(model,initial,control,o,{});live.save(saved);
            r["restore_unmodified_checkpoint_control"]=attempt([&]{(void)h::Owner::restore(model,o,{},saved);});
            std::filesystem::copy(saved,edited,std::filesystem::copy_options::recursive);Json j;{std::ifstream in(saved/"checkpoint.json");in>>j;}
            j["field"]["B"][face]=j["field"]["B"][face].get<double>()+.01;
            j["field_sha256"]=model->state_sha256(h::checkpoint_detail::field(j.at("field"),m),h::checkpoint_detail::charge(j.at("charge"),m));
            {std::ofstream file(edited/"checkpoint.json",std::ios::binary);file<<j.dump()<<'\n';}
            r["restore_resealed_perturbed_checkpoint"]=attempt([&]{(void)h::Owner::restore(model,o,{},edited);});}
        Json sweep=Json::array();
        for(int p=2;p<=16;++p){const double delta=std::pow(10.,-p);auto s=initial;s.B[face]+=delta;
            sweep.push_back({{"delta_T",delta},{"divB_scan",divb_scan(m,s)},{"Model_admit",attempt([&]{model->admit(s,control);})}});}
        r["supplementary_detection_sweep"]=sweep;
        result["A09_divB"]=r;
    }

    // ---------------------------------------------------------------- A10 midpoint, field-only step (probe step)
    {
        Json r;auto f=model->prepare(initial,control,Vec(m.ne),control.particles,sigma,dt);const auto&s1=f.field();
        r["current_A"]=f.receipt().current_A;r["probe_current_A"]=I0+.008;r["current_minus_probe_current_A"]=f.receipt().current_A-(I0+.008);
        r["E_high_word_vs_authoritative"]=high_word_comparison(initial.E,s1.E,f.midpoint().E());
        r["B_high_word_vs_authoritative"]=high_word_comparison(initial.B,s1.B,f.midpoint().B());
        r["E_authoritative_equals_compensated_average"]=f.midpoint().E()==compensated_average(initial.E,initial.E_low,s1.E,s1.E_low);
        r["B_authoritative_equals_compensated_average"]=f.midpoint().B()==compensated_average(initial.B,initial.B_low,s1.B,s1.B_low);
        r["Maxwell_receipt_mid_E_equals_authoritative"]=f.receipt().mid_E==f.midpoint().E();
        r["joule_cell_power_from_authoritative"]=f.receipt().cell_joule_W==m.cell_joule_W(sigma,f.midpoint().E());
        std::size_t nonzero_E_low=0;for(double v:s1.E_low)nonzero_E_low+=v!=0;r["candidate_nonzero_E_low_words"]=nonzero_E_low;
        r["candidate_constraints"]={{"max_gauss_C",f.constraints().max_gauss_C},{"max_gauss_bound_C",f.constraints().max_gauss_bound_C},
            {"max_divB_T_per_m",f.constraints().max_divB_T_m},{"max_divB_bound_T_per_m",f.constraints().max_divB_bound_T_m}};
        r["energy_balance_J"]=f.receipt().balance_J;
        // Original operator at the corrected model's current, for a same-current comparison.
        const auto same=fields.prepare(initial,Vec(g.ne),sigma,dt,f.receipt().current_A);
        r["original_operator_same_current_high_word_vs_receipt_mid_E"]=high_word_comparison(initial.E,same.first.E,same.second.mid_E);
        {std::size_t d=0;double mx=0;for(int e=0;e<m.ne;++e){if(same.second.mid_E[e]!=f.midpoint().E()[e])++d;mx=std::max(mx,std::abs(same.second.mid_E[e]-f.midpoint().E()[e]));}
            r["original_receipt_mid_E_vs_corrected_authoritative"]={{"differing",d},{"max_abs_difference",mx}};}
        result["A10_field_step"]=r;
    }

    // ---------------------------------------------------------------- A10 consumers through prepare_closed_proposal
    {
        Json cases=Json::array();const auto source=h::sha256("audit port: unchanged closed source identity");
        auto closed=[&](const std::string&label,const std::vector<n::ArenaRateEntry>&rows){
            Json r{{"case",label},{"census_rows",rows.size()}};const auto ledger=dc_declaration(m,seed,h::charge_of_census(m,o,rows));
            std::unique_ptr<h::Owner>owner;r["Owner_bootstrap"]=attempt([&]{owner=std::make_unique<h::Owner>(model,seed,ledger,o,rows);});
            if(!owner){cases.push_back(r);return;}
            const auto b=binding_for(*owner,o,rows);std::optional<h::JointProposal>p;
            r["prepare_closed_proposal"]=attempt([&]{p.emplace(h::prepare_closed_proposal(*owner,o,rows,b,source));});
            if(!p){cases.push_back(r);return;}
            const auto&mo=p->motion();const auto&mid=p->midpoint();
            r["iterations"]=p->iterations();r["owner_preflight"]=owner->preflight(*p,rows);
            r["motion_binding_hash_equals_authoritative_midpoint"]=mo.binding().trial_field_sha256==hd::total_field_sha256(mid.E(),mid.B());
            {auto replay=b;replay.trial_field_sha256=hd::total_field_sha256(mid.E(),mid.B());
                auto again=hd::prepare_canonical_motion(m,o,rows,replay,{},mid.E(),mid.B());
                r["motion_replay_from_authoritative_arrays_identical_seal"]=again.seal_sha256()==mo.seal_sha256();}
            r["motion_field_work_uses_authoritative_E"]=mo.budget().field_work_J==dt*axial::dot(mid.E(),mo.current_A_m());
            r["kinetic_change_J"]=mo.budget().kinetic_change_J;r["motion_field_work_J"]=mo.budget().field_work_J;
            std::size_t moving=0;for(double v:mo.current_A_m())moving+=v!=0;r["nonzero_current_edges"]=moving;
            auto f=model->prepare(owner->field(),owner->charge(),mo.current_A_m(),h::charge_of_endpoints(m,o,mo.endpoints()),model->initial_conductivity(),dt);
            const bool same_candidate=f.after_sha256()==p->candidate_field_sha256();r["parallel_candidate_identical_to_proposal"]=same_candidate;
            if(same_candidate){const auto&s0=owner->field();const auto&s1=f.field();
                r["midpoint_identical_to_parallel_candidate"]=f.midpoint().E()==mid.E()&&f.midpoint().B()==mid.B()&&f.midpoint().seal()==mid.seal();
                r["E_authoritative_equals_compensated_average"]=mid.E()==compensated_average(s0.E,s0.E_low,s1.E,s1.E_low);
                r["B_authoritative_equals_compensated_average"]=mid.B()==compensated_average(s0.B,s0.B_low,s1.B,s1.B_low);
                r["E_high_word_vs_authoritative"]=high_word_comparison(s0.E,s1.E,mid.E());
                r["B_high_word_vs_authoritative"]=high_word_comparison(s0.B,s1.B,mid.B());
                r["Maxwell_plasma_work_J"]=f.receipt().plasma_work_J;r["Maxwell_receipt_mid_E_equals_authoritative"]=f.receipt().mid_E==mid.E();
                r["joule_cell_power_from_authoritative"]=f.receipt().cell_joule_W==m.cell_joule_W(model->initial_conductivity(),mid.E());
                const double work=mo.budget().kinetic_change_J-f.receipt().plasma_work_J;
                r["particle_minus_Maxwell_work_J"]=work;
                r["work_guard_J"]=1e-20+1e-10*(std::abs(mo.budget().kinetic_change_J)+std::abs(f.receipt().plasma_work_J));
                r["candidate_constraints"]={{"max_gauss_C",f.constraints().max_gauss_C},{"max_gauss_bound_C",f.constraints().max_gauss_bound_C},
                    {"max_divB_T_per_m",f.constraints().max_divB_T_m}};
                Vec hwE(m.ne),hwB(m.nb);for(int e=0;e<m.ne;++e)hwE[e]=.5*(s0.E[e]+s1.E[e]);for(int k=0;k<m.nb;++k)hwB[k]=.5*(s0.B[k]+s1.B[k]);
                auto counter=b;counter.trial_field_sha256=hd::total_field_sha256(hwE,hwB);
                Json cf;try{auto mh=hd::prepare_canonical_motion(m,o,rows,counter,{},hwE,hwB);
                    cf={{"executed",true},{"seal_differs_from_authoritative",mh.seal_sha256()!=mo.seal_sha256()},
                        {"current_differs",mh.current_A_m()!=mo.current_A_m()},{"binding_hash_differs",counter.trial_field_sha256!=mo.binding().trial_field_sha256}};}
                catch(const std::exception&e){cf={{"executed",false},{"rejection",e.what()}};}
                r["counterfactual_high_word_motion"]=cf;}
            cases.push_back(r);};
        closed("empty census",{});
        for(double W:{1e9,1e11,1e12})closed("co-located electron/proton pair, weight "+Json(W).dump(),
            {carrier(1,o,audit_point,{2e6,-1e6,5e5},-pic_e,pic_me,W,"pair-electron"),carrier(2,o,audit_point,{0,0,0},pic_e,hall::proton_mass,W,"pair-proton")});
        result["A10_closed_consumers"]=cases;
        result["A10_source_consumer"]={{"portable",false},{"reason","hall_corrective closed port exposes no source sampling or source population API; prepare_closed_proposal requires an unchanged source identity digest only"}};
    }

    // ---------------------------------------------------------------- Gauss bound: independent solver-residual reconstruction
    {
        Json r;const auto ini=axial::dc_initial(m,model->initial_conductivity(),model->drive().current(model->drive().start()));
        r["recomputed_dc_initial_equals_seed"]=ini.E==seed.E&&ini.B==seed.B;
        std::vector<bool>fixed(m.nn);for(int k=0;k<m.nn;++k)fixed[k]=conductor_node(m,k);
        Vec known(m.nn),vac(m.nn);for(int k=0;k<m.nn;++k)(fixed[k]?known[k]:vac[k])=ini.potential[k];
        auto flux=m.gradient(known);for(int e=0;e<m.ne;++e)flux[e]*=m.VE[e];auto b=m.gradient_transpose(flux);
        for(int k=0;k<m.nn;++k)b[k]=fixed[k]?0:-b[k];
        auto apply=[&](Vec v){for(int k=0;k<m.nn;++k)if(fixed[k])v[k]=0;auto d=m.gradient(v);for(int e=0;e<m.ne;++e)d[e]*=m.VE[e];
            auto a=m.gradient_transpose(d);for(int k=0;k<m.nn;++k)if(fixed[k])a[k]=0;return a;};
        const auto residual=axial::subtract(b,apply(vac));const double nb=axial::norm(b),target=std::max(1e-30,2e-13*nb);
        r["vacuum_rhs_L2"]=nb;r["solver_target_L2"]=target;r["achieved_residual_L2"]=axial::norm(residual);r["achieved_residual_Linf"]=axial::maximum(residual);
        r["eps0_times_target_C"]=axial::eps*target;r["model_solver_bound_C"]=solver_bound;r["model_bound_equals_eps0_times_target"]=axial::eps*target==solver_bound;
        r["eps0_times_achieved_residual_L2_C"]=axial::eps*axial::norm(residual);r["eps0_times_achieved_residual_Linf_C"]=axial::eps*axial::maximum(residual);
        const auto Q=h::field_charge(m,seed);const auto E=h::words(seed.E,seed.E_low);std::vector<long double>activity(m.nn);
        for(int e=0;e<m.ne;++e){const auto&edge=m.edges[e];const long double a=h::magnitude(E[e])*m.M[e]/edge.length;activity[edge.a]+=a;activity[edge.b]+=a;}
        long double max_q=0,max_gap=0,max_old_ratio=0;std::size_t old_violations=0,vacuum_nodes=0;int node_old=-1;
        for(int k=0;k<m.nn;++k){if(fixed[k]||m.outer_node[k])continue;++vacuum_nodes;const long double q=Q[k].wide();
            max_q=std::max(max_q,std::abs(q));max_gap=std::max(max_gap,std::abs(q+(long double)axial::eps*residual[k]));
            const long double old_bound=1e-12L*activity[k];if(!(std::abs(q)<=old_bound)){++old_violations;}
            if(old_bound>0&&std::abs(q)/old_bound>max_old_ratio){max_old_ratio=std::abs(q)/old_bound;node_old=k;}}
        r["vacuum_interior_nodes"]=vacuum_nodes;r["max_vacuum_field_charge_C"]=double(max_q);
        r["max_abs_field_charge_plus_eps0_residual_C"]=double(max_gap);
        r["bound_over_max_vacuum_defect"]=double(solver_bound/max_q);
        r["prototype_local_bound_1e-12_activity"]={{"vacuum_nodes_violating",old_violations},{"max_defect_over_local_bound",double(max_old_ratio)},{"node_of_max_ratio",node_old}};
        const int k=2163;if(k<m.nn)r["node_2163"]={{"class",node_class(m,k)},{"field_charge_C",double(Q[k].wide())},{"minus_eps0_residual_C",-axial::eps*residual[k]},
            {"field_activity_C",double(activity[k])},{"prototype_local_bound_C",double(1e-12L*activity[k])},
            {"corrected_bound_C",double(solver_bound+256*machine_eps*activity[k])}};
        result["gauss_bound_reconstruction"]=r;
    }

    result["status"]="executed";std::ofstream(out/"port-observations.json")<<result.dump(2)<<'\n';std::cout<<result.dump(2)<<'\n';
}catch(const std::exception&e){std::cerr<<"port aborted: "<<e.what()<<'\n';rc=1;}MPI_Finalize();return rc;}
