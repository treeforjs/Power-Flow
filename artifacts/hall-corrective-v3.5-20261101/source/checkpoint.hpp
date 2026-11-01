#pragma once
// Hall corrective owner v3 checkpoint export (human-readable JSON plus the
// unchanged heat slices). Derived from v2 / 8dddec0 checkpoint.hpp; the stored
// digests are the v3 binary identities; v3 adds the cumulative kinetic sink,
// surface birth and surface load ledgers. Not on the per-step path.
// v3.2: the admitted Gauss ledger is stored under "admitted_gauss_C" (word pairs; empty = zero)
// and is covered by field_sha256 (state digest).
#include "owner.hpp"

namespace hall_corrective_v3 {
namespace checkpoint_detail {
inline double number(const Json&j){need(j.is_number(),"checkpoint numeric type");double x=j.get<double>();need(std::isfinite(x),"checkpoint nonfinite scalar");return x;}
inline std::uint64_t revision(const Json&j){need(j.is_number_unsigned()||(j.is_number_integer()&&j.get<std::int64_t>()>=0),"checkpoint revision must be an exact nonnegative integer");
    auto n=j.get<std::uint64_t>();need(n<UINT64_MAX,"checkpoint revision exhausted");return n;}
inline void keys(const Json&j,std::initializer_list<const char*>required){need(j.is_object()&&j.size()==required.size(),"checkpoint schema shape");
    for(auto key:required)need(j.contains(key),"checkpoint schema key");}
inline Vec vector(const Json&j,std::size_t n){need(j.is_array()&&j.size()==n,"checkpoint vector dimensions");Vec out;out.reserve(n);for(const auto&v:j)out.push_back(number(v));return out;}
inline WordVec word_vector(const Json&j,std::size_t n){need(j.is_array()&&j.size()==n,"checkpoint word-vector dimensions");WordVec out;out.reserve(n);
    for(const auto&v:j){need(v.is_array()&&v.size()==2,"checkpoint word-pair dimensions");Words x;x.hi=number(v[0]);x.lo=number(v[1]);
        need(canonical(x),"checkpoint noncanonical words");out.push_back(x);}return out;}
inline WordVec ledger_words(const Json&j,std::size_t n){need(j.is_array()&&(j.empty()||j.size()==n),"checkpoint admitted Gauss ledger dimensions");
    return j.empty()?WordVec():word_vector(j,n);}
inline Json ledger_json(const WordVec&v){return v.empty()?Json::array():record(v);}
// v3.4: divB gross ledger (one non-negative value per cell; empty = zero)
inline Vec divb_gross(const Json&j,std::size_t cells){need(j.is_array()&&(j.empty()||j.size()==cells),"checkpoint divB gross ledger dimensions");
    Vec out;out.reserve(j.size());for(const auto&x:j){const double v=number(x);need(v>=0,"checkpoint divB gross ledger sign");out.push_back(v);}return out;}
inline std::string hash(const Json&j){need(j.is_string(),"checkpoint hash type");auto s=j.get<std::string>();need(digest(s),"checkpoint hash format");return s;}
inline axial::State field(const Json&j,const Mesh&g){keys(j,{"E","B","E_low","B_low","time_s","revision"});axial::State s;
    s.E=vector(j.at("E"),g.ne);s.B=vector(j.at("B"),g.nb);s.E_low=vector(j.at("E_low"),g.ne);s.B_low=vector(j.at("B_low"),g.nb);
    s.time_s=number(j.at("time_s"));s.revision=revision(j.at("revision"));return s;}
inline ChargeLedger charge(const Json&j,const Mesh&g){keys(j,{"free_C","gross_free_C","conductor_C","exterior_C","conductor_declaration","exterior_declaration"});ChargeLedger q;
    q.particles={word_vector(j.at("free_C"),g.nn),vector(j.at("gross_free_C"),g.nn)};q.conductor_C=word_vector(j.at("conductor_C"),g.nn);
    q.exterior_C=word_vector(j.at("exterior_C"),g.nn);q.conductor_declaration=hash(j.at("conductor_declaration"));q.exterior_declaration=hash(j.at("exterior_declaration"));return q;}
inline Json words_json(Words w){Json a=Json::array();a.push_back(w.hi);a.push_back(w.lo);return a;}
inline Words words_of(const Json&j){need(j.is_array()&&j.size()==2,"checkpoint word pair");Words w;w.hi=number(j[0]);w.lo=number(j[1]);need(canonical(w),"checkpoint noncanonical words");return w;}
inline Json tally_json(const ChargeTally&t){Json j=Json::object();j["count"]=t.count;j["non_elementary"]=t.non_elementary;j["weight"]=exact_json(t.weight);j["kinetic_exact"]=exact_json(t.kinetic_exact);
    j["positive_e"]=exact_json(t.positive_e);j["negative_e"]=exact_json(t.negative_e);j["charge_C"]=words_json(t.charge_C);j["kinetic_J"]=words_json(t.kinetic_J);return j;}
inline ExactCount exact_of(const Json&j){need(j.is_array()&&j.size()<=ExactCount::maximum_carriers,"checkpoint exact count");std::vector<double>w;for(const auto&x:j)w.push_back(number(x));
    auto c=ExactCount::from_positive_carriers(w);need(c.positive_carriers()==w,"checkpoint noncanonical exact count");return c;}
inline ChargeTally tally_of(const Json&j){keys(j,{"count","non_elementary","weight","kinetic_exact","positive_e","negative_e","charge_C","kinetic_J"});ChargeTally t;t.kinetic_exact=exact_of(j.at("kinetic_exact"));
    t.count=revision(j.at("count"));t.non_elementary=revision(j.at("non_elementary"));t.weight=exact_of(j.at("weight"));t.positive_e=exact_of(j.at("positive_e"));
    t.negative_e=exact_of(j.at("negative_e"));t.charge_C=words_of(j.at("charge_C"));t.kinetic_J=words_of(j.at("kinetic_J"));return t;}
inline Json ledgers_json(const KineticLedgers&l){Json j=Json::object();Json s=Json::array(),b=Json::array();
    for(const auto&t:l.sinks.boundary)s.push_back(tally_json(t));for(const auto&t:l.births.origin)b.push_back(tally_json(t));
    j["sinks"]=std::move(s);j["births"]=std::move(b);
    Json c=Json::object();for(const auto&[k,v]:l.load_channel_J)c[k]=words_json(v);j["load_channel_J"]=std::move(c);j["loads"]=l.loads;
    j["wall_heat_J"]=words_json(l.wall_heat_J);j["external_J"]=words_json(l.external_J);j["short_bridge_J"]=words_json(l.short_bridge_J);j["escaped_J"]=words_json(l.escaped_J);return j;}
inline KineticLedgers ledgers_of(const Json&j){
    keys(j,{"sinks","births","load_channel_J","loads","wall_heat_J","external_J","short_bridge_J","escaped_J"});KineticLedgers l;
    need(j.at("sinks").is_array()&&j.at("sinks").size()==5&&j.at("births").is_array()&&j.at("births").size()==3,"checkpoint ledger dimensions");
    for(int b=0;b<5;++b)l.sinks.boundary[b]=tally_of(j.at("sinks")[b]);for(int c=0;c<3;++c)l.births.origin[c]=tally_of(j.at("births")[c]);
    need(j.at("load_channel_J").is_object(),"checkpoint load channels");for(const auto&[k,v]:j.at("load_channel_J").items())l.load_channel_J[k]=words_of(v);
    l.loads=revision(j.at("loads"));l.wall_heat_J=words_of(j.at("wall_heat_J"));l.external_J=words_of(j.at("external_J"));
    l.short_bridge_J=words_of(j.at("short_bridge_J"));l.escaped_J=words_of(j.at("escaped_J"));return l;}
} // namespace checkpoint_detail
namespace genesis_detail {
inline Json provenance_json(const GenesisProvenance& p) {
    if (p.empty()) return Json(nullptr);
    Json j = Json::object();
    j["checkpoint_schema"] = p.checkpoint_schema; j["configuration_sha256"] = p.configuration_sha256; j["field_sha256"] = p.field_sha256;
    j["heat_sha256"] = p.heat_sha256; j["kinetic_ledgers_sha256"] = p.kinetic_ledgers_sha256; j["motion_mesh_sha256"] = p.motion_mesh_sha256;
    j["revision"] = p.revision; j["time_s"] = p.time_s; j["transfer_text"] = p.transfer_text; j["sha256"] = p.sha256();
    return j;
}
inline GenesisProvenance provenance_of(const Json& j) {
    using namespace checkpoint_detail;
    GenesisProvenance p;
    if (j.is_null()) return p;
    keys(j, {"checkpoint_schema", "configuration_sha256", "field_sha256", "heat_sha256", "kinetic_ledgers_sha256", "motion_mesh_sha256",
             "revision", "time_s", "transfer_text", "sha256"});
    need(j.at("checkpoint_schema").is_string() && !j.at("checkpoint_schema").get<std::string>().empty() && j.at("transfer_text").is_string(),
         "checkpoint genesis provenance types");
    p.checkpoint_schema = j.at("checkpoint_schema").get<std::string>(); p.configuration_sha256 = hash(j.at("configuration_sha256"));
    p.field_sha256 = hash(j.at("field_sha256")); p.heat_sha256 = hash(j.at("heat_sha256"));
    p.kinetic_ledgers_sha256 = hash(j.at("kinetic_ledgers_sha256")); p.motion_mesh_sha256 = hash(j.at("motion_mesh_sha256"));
    p.revision = checkpoint_detail::revision(j.at("revision")); p.time_s = number(j.at("time_s")); p.transfer_text = j.at("transfer_text").get<std::string>();
    need(hash(j.at("sha256")) == p.sha256(), "checkpoint genesis provenance digest mismatch");
    return p;
}
inline std::string limits_text(const CarrierLimits& l) { return std::to_string(l.maximum_carriers) + "/" + std::to_string(l.memory_budget_bytes); }
inline std::string declaration_text(const ThermalDeclaration& d) { return d.empty() ? std::string("none") : d.sha256(); }
inline std::string g17(double x) { char b[40]; std::snprintf(b, sizeof b, "%.17g", x); return b; }
inline bool same_limits(const CarrierLimits& a, const CarrierLimits& b) {
    return a.maximum_carriers == b.maximum_carriers && a.memory_budget_bytes == b.memory_budget_bytes;
}
// Same heat-only patch (the fields bound into the declaration digest).
inline bool same_patch(const eh::Patch& a, const eh::Patch& b) {
    return a.id == b.id && a.face == b.face && a.begin == b.begin && a.end == b.end && a.desorption_enabled == b.desorption_enabled &&
           a.initial_particles == b.initial_particles && a.initial_low_particles == b.initial_low_particles && a.source_hash == b.source_hash;
}
// `now` only adds heat-only patches to `recorded`: every recorded slice keeps its source hash and its patches as the
// leading patches of the same slice; added patches and slices pass the thermal owner's heat-only rule (ThermalV3::checked
// at the Model's construction: no desorption, zero inventory, patch source hash = slice source hash).
inline bool declaration_extends(const ThermalDeclaration& recorded, const ThermalDeclaration& now) {
    for (const auto& [key, r] : recorded.slices) {
        const auto it = now.slices.find(key);
        if (it == now.slices.end() || it->second.source_hash != r.source_hash || it->second.patches.size() < r.patches.size()) return false;
        for (std::size_t i = 0; i < r.patches.size(); ++i) if (!same_patch(r.patches[i], it->second.patches[i])) return false;
    }
    for (const auto& [key, s] : now.slices)
        for (const auto& p : s.patches)
            if (p.desorption_enabled || p.initial_particles != 0 || p.initial_low_particles != 0 || p.coverage_temkin || p.source_hash != s.source_hash)
                return false;
    return true;
}
inline bool vacuum_charge(const ChargeLedger& q) {
    for (const auto& w : q.particles.net_C) if (w.hi != 0 || w.lo != 0) return false;
    for (double g : q.particles.gross_C) if (g != 0) return false;
    return q.residue_C.empty();
}
} // namespace genesis_detail

inline void Owner::save(const std::filesystem::path&dir)const{
    clock_guard();model_->admit(field_,charge_);need(!std::filesystem::exists(dir),"checkpoint destination already exists");
    std::filesystem::create_directory(dir);thermal_->save(dir/"heat");
    Json j=Json::object();
    j["schema"]=api;j["configuration_sha256"]=configuration_sha256();j["motion_mesh_sha256"]=motion_.mesh_sha256;
    j["field"]=state_record(field_);j["charge"]=record(charge_);j["field_sha256"]=field_sha256_;j["heat_sha256"]=heat_sha256_;
    j["kinetic_ledgers"]=checkpoint_detail::ledgers_json(ledgers_);j["kinetic_ledgers_sha256"]=ledgers_sha256_;j["bank_residue_C"]=charge_.residue_C;
    j["admitted_gauss_C"]=checkpoint_detail::ledger_json(charge_.gauss_admitted_C);
    j["divb_gross_T_m"]=charge_.divb_gross_T_m;
    j["genesis_provenance"]=genesis_detail::provenance_json(provenance_);   // v3.5
    std::string text=j.dump();drain(j);
    std::ofstream out(dir/"checkpoint.json",std::ios::binary);out<<text<<'\n';out.close();need(bool(out),"checkpoint metadata write failed");
}
namespace genesis_detail {
// Host Owner checkpoint (v3.3: no divB gross ledger; v3.4: + divb_gross_T_m; v3.5: + genesis_provenance).
struct HostCheckpoint {
    std::string recorded, motion_mesh_sha256, field_sha256, heat_sha256, ledgers_sha256;
    axial::State state; ChargeLedger ledger; KineticLedgers ledgers; GenesisProvenance provenance; bool v33 = false, v35 = false;
};
inline HostCheckpoint parse_host_checkpoint(const Mesh& g, const std::filesystem::path& dir, bool allow_v33) {
    using namespace checkpoint_detail;
    const auto file = dir / "checkpoint.json";
    const auto limit = std::uintmax_t(128) * (2 * std::uintmax_t(g.ne) + 2 * std::uintmax_t(g.nb) + 9 * std::uintmax_t(g.nn) + g.cells.size()) + 1048576;
    need(std::filesystem::is_regular_file(file) && std::filesystem::file_size(file) <= limit, "checkpoint byte cap");
    std::ifstream in(file, std::ios::binary); Json j = Json::parse(in);
    struct Guard { Json& j; ~Guard() { drain(j); } } guard{j};
    HostCheckpoint c;
    need(j.is_object(), "checkpoint schema shape");
    if (j.contains("genesis_provenance"))
        keys(j, {"schema", "configuration_sha256", "motion_mesh_sha256", "field", "charge", "field_sha256", "heat_sha256", "kinetic_ledgers",
                 "kinetic_ledgers_sha256", "bank_residue_C", "admitted_gauss_C", "divb_gross_T_m", "genesis_provenance"});
    else if (j.contains("divb_gross_T_m"))
        keys(j, {"schema", "configuration_sha256", "motion_mesh_sha256", "field", "charge", "field_sha256", "heat_sha256", "kinetic_ledgers",
                 "kinetic_ledgers_sha256", "bank_residue_C", "admitted_gauss_C", "divb_gross_T_m"});
    else {
        need(allow_v33, "checkpoint schema shape");
        keys(j, {"schema", "configuration_sha256", "motion_mesh_sha256", "field", "charge", "field_sha256", "heat_sha256", "kinetic_ledgers",
                 "kinetic_ledgers_sha256", "bank_residue_C", "admitted_gauss_C"});
        c.v33 = true;
    }
    need(j.at("schema") == api, "checkpoint schema");
    c.recorded = hash(j.at("configuration_sha256")); c.motion_mesh_sha256 = hash(j.at("motion_mesh_sha256"));
    c.state = checkpoint_detail::field(j.at("field"), g); c.ledger = checkpoint_detail::charge(j.at("charge"), g);
    {   const auto& rj = j.at("bank_residue_C"); need(rj.is_array() && (rj.empty() || rj.size() == std::size_t(g.nn)), "checkpoint bank residue dimensions");
        for (const auto& x : rj) { const double v = number(x); need(v >= 0, "checkpoint bank residue sign"); c.ledger.residue_C.push_back(v); } }
    c.ledger.gauss_admitted_C = ledger_words(j.at("admitted_gauss_C"), g.nn);
    if (!c.v33) c.ledger.divb_gross_T_m = divb_gross(j.at("divb_gross_T_m"), g.cells.size());
    if (j.contains("genesis_provenance")) { c.v35 = true; c.provenance = provenance_of(j.at("genesis_provenance")); }
    c.field_sha256 = hash(j.at("field_sha256")); c.heat_sha256 = hash(j.at("heat_sha256"));
    c.ledgers = ledgers_of(j.at("kinetic_ledgers")); c.ledgers_sha256 = hash(j.at("kinetic_ledgers_sha256"));
    return c;
}
} // namespace genesis_detail
inline std::unique_ptr<Owner> Owner::restore(std::shared_ptr<const Model>model,const native::GapPicOptions&o,
    std::span<const native::ArenaRateEntry>census,const std::filesystem::path&dir){
    model=checked_model(std::move(model));const auto&g=model->mesh();
    auto c=genesis_detail::parse_host_checkpoint(g,dir,false);
    need(c.recorded==model->configuration_sha256()&&c.motion_mesh_sha256==model->motion_mesh_sha256(o),"checkpoint field/forcing/mesh configuration changed");
    if(!c.provenance.empty())c.ledger.genesis_provenance=c.provenance.sha256();   // v3.5: a rebased owner keeps its provenance
    need(c.field_sha256==model->state_sha256(c.state,c.ledger),"checkpoint field/charge seal mismatch");
    auto candidate=std::unique_ptr<Owner>(new Owner(model,std::move(c.state),std::move(c.ledger),o,census,RestoreTag{true}));
    candidate->provenance_=c.provenance;
    // Every heat slice loads into a wholly private aggregate; a late corrupt slice
    // destroys this candidate without changing any accepted owner.
    candidate->thermal_->load(dir/"heat");candidate->clock_guard();
    candidate->heat_sha256_=heat_sha256_of(candidate->heat_identity_,*candidate->thermal_);
    need(c.heat_sha256==candidate->heat_sha256_,"checkpoint thermal seal mismatch");
    candidate->ledgers_=std::move(c.ledgers);candidate->ledgers_sha256_=candidate->ledgers_.sha256();
    need(c.ledgers_sha256==candidate->ledgers_sha256_,"checkpoint kinetic ledger seal mismatch");
    return candidate;
}
// v3.5 genesis of the host Owner: the checks a-f of ResidentOwner::restore_as_genesis (resident_owner.hpp) on the host
// checkpoint format; the recorded configuration is this Model's composition with its own dt_max, solve backend and
// motion forms and the checkpoint's declaration, carrier limits and rule form (v3.3 when the checkpoint has no divB
// gross ledger).
inline std::unique_ptr<Owner> Owner::restore_as_genesis(std::shared_ptr<const Model> model, const native::GapPicOptions& o,
                                                        std::span<const native::ArenaRateEntry> census, const std::filesystem::path& dir,
                                                        GenesisTransfer transfer) {
    using namespace genesis_detail;
    model = checked_model(std::move(model)); const auto& g = model->mesh();
    auto c = parse_host_checkpoint(g, dir, true);
    const ConfigurationForm form = c.v33 ? ConfigurationForm::v33 : ConfigurationForm::v34;
    std::vector<ThermalDeclaration> decls;
    if (transfer.recorded_declaration) decls.push_back(*transfer.recorded_declaration);
    else { decls.push_back(model->thermal_declaration()); if (!model->thermal_declaration().empty()) decls.push_back(ThermalDeclaration{}); }
    std::vector<CarrierLimits> lims{model->limits()};
    for (const auto& l : {transfer.recorded_limits, std::optional<CarrierLimits>(CarrierLimits{})})
        if (l && std::none_of(lims.begin(), lims.end(), [&](const CarrierLimits& x) { return same_limits(x, *l); })) lims.push_back(*l);
    const ThermalDeclaration* dr = nullptr; const CarrierLimits* lr = nullptr;
    for (const auto& d : decls) {
        for (const auto& l : lims)
            if (model->configuration_transfer(model->dt_max(), model->solver().name(), model->solver().descriptor_sha256(), model->motion_forms(), d, l, form) ==
                c.recorded) { dr = &d; lr = &l; break; }
        if (dr) break;
    }
    need(dr != nullptr, "genesis checkpoint configuration is not a transfer of the restoring Model (mesh, drive, material law, seed, rules, guards or operands differ)");
    // (a) seals under the recorded configuration
    const std::string prior = c.provenance.empty() ? std::string() : c.provenance.sha256();
    c.ledger.genesis_provenance = prior;
    need(model->state_sha256_under(c.recorded, c.state, c.ledger, form) == c.field_sha256, "checkpoint field/charge seal mismatch");
    ThermalV3 recorded_heat(g, model->drive().start(), *dr, model->dt_max());
    recorded_heat.load(dir / "heat");
    need(recorded_heat.time_s() == c.state.time_s && recorded_heat.revision() == c.state.revision, "accepted field/heat clock mismatch");
    need(heat_sha256_of(genesis_heat_identity(Model::material_sha256_of(*dr), prior), recorded_heat) == c.heat_sha256, "checkpoint thermal seal mismatch");
    need(c.ledgers.sha256() == c.ledgers_sha256, "checkpoint kinetic ledger seal mismatch");
    // (b) transfer, (c) vacuum
    const std::string current_mesh = model->motion_mesh_sha256(o);
    const bool same_declaration = dr->sha256() == model->thermal_declaration().sha256() && dr->empty() == model->thermal_declaration().empty();
    need(same_declaration || (transfer.thermal_declaration && declaration_extends(*dr, model->thermal_declaration())),
         "genesis thermal declaration may only add heat-only source patches");
    const bool vacuum = vacuum_charge(c.ledger) && c.ledgers.sha256() == KineticLedgers{}.sha256();
    need(same_limits(*lr, model->limits()) || (transfer.carrier_limits && vacuum), "genesis carrier limits differ");
    need(c.motion_mesh_sha256 == current_mesh || (transfer.native_geometry && vacuum), "genesis native geometry differs");
    need(vacuum, "genesis checkpoint is not a vacuum state (free charge, bank residue or kinetic ledgers)");
    // (d) admission in the restoring Model (census charge, host admission)
    auto owner = std::unique_ptr<Owner>(new Owner(model, c.state, c.ledger, o, census, RestoreTag{true}));
    // (e) rebase
    GenesisProvenance p;
    p.checkpoint_schema = c.v33 ? "hall-v3.3-host-checkpoint" : c.v35 ? "hall-v3.5-host-checkpoint" : "hall-v3.4-host-checkpoint";
    p.configuration_sha256 = c.recorded; p.field_sha256 = c.field_sha256; p.heat_sha256 = c.heat_sha256;
    p.kinetic_ledgers_sha256 = c.ledgers_sha256; p.motion_mesh_sha256 = c.motion_mesh_sha256; p.revision = c.state.revision; p.time_s = c.state.time_s;
    p.transfer_text = std::string("hall-v3.5-genesis-transfer-v1") + "; rule_form=" + (c.v33 ? "v3.3" : "v3.4") +
        "; thermal_declaration=" + declaration_text(*dr) + "->" + declaration_text(model->thermal_declaration()) +
        "; carrier_limits=" + limits_text(*lr) + "->" + limits_text(model->limits()) +
        "; motion_mesh=" + c.motion_mesh_sha256 + "->" + current_mesh +
        "; divb_gross_ledger=" + (c.v33 ? "zero (v3.3 checkpoint)" : "kept") +
        "; prior_provenance=" + (prior.empty() ? std::string("none") : prior);
    owner->thermal_->install_genesis(recorded_heat);
    owner->field_.revision = 0;
    owner->ledgers_ = KineticLedgers{}; owner->ledgers_sha256_ = owner->ledgers_.sha256();
    // (f) identity
    owner->provenance_ = std::move(p); owner->charge_.genesis_provenance = owner->provenance_.sha256();
    owner->heat_identity_ = genesis_heat_identity(model->material_sha256(), owner->charge_.genesis_provenance);
    owner->field_sha256_ = model->state_sha256(owner->field_, owner->charge_);
    owner->heat_sha256_ = heat_sha256_of(owner->heat_identity_, *owner->thermal_);
    owner->accepted_verified_generation_ = UINT64_MAX;
    owner->clock_guard();
    return owner;
}
} // namespace hall_corrective_v3
