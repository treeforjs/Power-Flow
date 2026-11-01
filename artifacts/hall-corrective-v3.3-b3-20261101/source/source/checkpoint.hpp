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
inline void Owner::save(const std::filesystem::path&dir)const{
    clock_guard();model_->admit(field_,charge_);need(!std::filesystem::exists(dir),"checkpoint destination already exists");
    std::filesystem::create_directory(dir);thermal_->save(dir/"heat");
    Json j=Json::object();
    j["schema"]=api;j["configuration_sha256"]=configuration_sha256();j["motion_mesh_sha256"]=motion_.mesh_sha256;
    j["field"]=state_record(field_);j["charge"]=record(charge_);j["field_sha256"]=field_sha256_;j["heat_sha256"]=heat_sha256_;
    j["kinetic_ledgers"]=checkpoint_detail::ledgers_json(ledgers_);j["kinetic_ledgers_sha256"]=ledgers_sha256_;j["bank_residue_C"]=charge_.residue_C;
    j["admitted_gauss_C"]=checkpoint_detail::ledger_json(charge_.gauss_admitted_C);
    std::string text=j.dump();drain(j);
    std::ofstream out(dir/"checkpoint.json",std::ios::binary);out<<text<<'\n';out.close();need(bool(out),"checkpoint metadata write failed");
}
inline std::unique_ptr<Owner> Owner::restore(std::shared_ptr<const Model>model,const native::GapPicOptions&o,
    std::span<const native::ArenaRateEntry>census,const std::filesystem::path&dir){
    model=checked_model(std::move(model));const auto&g=model->mesh();const auto file=dir/"checkpoint.json";
    const auto limit=std::uintmax_t(128)*(2*std::uintmax_t(g.ne)+2*std::uintmax_t(g.nb)+9*std::uintmax_t(g.nn))+1048576;
    need(std::filesystem::is_regular_file(file)&&std::filesystem::file_size(file)<=limit,"checkpoint byte cap");
    std::ifstream in(file,std::ios::binary);Json j=Json::parse(in);using namespace checkpoint_detail;
    struct Guard {Json&j;~Guard(){drain(j);}} guard{j};
    keys(j,{"schema","configuration_sha256","motion_mesh_sha256","field","charge","field_sha256","heat_sha256","kinetic_ledgers","kinetic_ledgers_sha256","bank_residue_C",
            "admitted_gauss_C"});
    need(j.at("schema")==api&&hash(j.at("configuration_sha256"))==model->configuration_sha256()&&
        hash(j.at("motion_mesh_sha256"))==model->motion_mesh_sha256(o),"checkpoint field/forcing/mesh configuration changed");
    auto state=checkpoint_detail::field(j.at("field"),g);auto ledger=checkpoint_detail::charge(j.at("charge"),g);
    {const auto&rj=j.at("bank_residue_C");need(rj.is_array()&&(rj.empty()||rj.size()==std::size_t(g.nn)),"checkpoint bank residue dimensions");
     for(const auto&x:rj){const double v=number(x);need(v>=0,"checkpoint bank residue sign");ledger.residue_C.push_back(v);}}
    ledger.gauss_admitted_C=ledger_words(j.at("admitted_gauss_C"),g.nn);
    need(hash(j.at("field_sha256"))==model->state_sha256(state,ledger),"checkpoint field/charge seal mismatch");
    auto candidate=std::unique_ptr<Owner>(new Owner(model,std::move(state),std::move(ledger),o,census,RestoreTag{true}));
    // Every heat slice loads into a wholly private aggregate; a late corrupt slice
    // destroys this candidate without changing any accepted owner.
    candidate->thermal_->load(dir/"heat");candidate->clock_guard();
    candidate->heat_sha256_=heat_sha256_of(model->material_sha256(),*candidate->thermal_);
    need(hash(j.at("heat_sha256"))==candidate->heat_sha256_,"checkpoint thermal seal mismatch");
    candidate->ledgers_=checkpoint_detail::ledgers_of(j.at("kinetic_ledgers"));candidate->ledgers_sha256_=candidate->ledgers_.sha256();
    need(hash(j.at("kinetic_ledgers_sha256"))==candidate->ledgers_sha256_,"checkpoint kinetic ledger seal mismatch");
    return candidate;
}
} // namespace hall_corrective_v3
