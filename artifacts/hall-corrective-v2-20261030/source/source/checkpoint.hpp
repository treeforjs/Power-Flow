#pragma once
// Hall corrective owner v2 checkpoint export (human-readable JSON plus the
// unchanged v2 heat slices). Derived from 8dddec0 checkpoint.hpp; the stored
// digests are the v2 binary identities. Not on the per-step path.
#include "owner.hpp"

namespace hall_corrective_v2 {
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
inline std::string hash(const Json&j){need(j.is_string(),"checkpoint hash type");auto s=j.get<std::string>();need(digest(s),"checkpoint hash format");return s;}
inline axial::State field(const Json&j,const Mesh&g){keys(j,{"E","B","E_low","B_low","time_s","revision"});axial::State s;
    s.E=vector(j.at("E"),g.ne);s.B=vector(j.at("B"),g.nb);s.E_low=vector(j.at("E_low"),g.ne);s.B_low=vector(j.at("B_low"),g.nb);
    s.time_s=number(j.at("time_s"));s.revision=revision(j.at("revision"));return s;}
inline ChargeLedger charge(const Json&j,const Mesh&g){keys(j,{"free_C","gross_free_C","conductor_C","exterior_C","conductor_declaration","exterior_declaration"});ChargeLedger q;
    q.particles={word_vector(j.at("free_C"),g.nn),vector(j.at("gross_free_C"),g.nn)};q.conductor_C=word_vector(j.at("conductor_C"),g.nn);
    q.exterior_C=word_vector(j.at("exterior_C"),g.nn);q.conductor_declaration=hash(j.at("conductor_declaration"));q.exterior_declaration=hash(j.at("exterior_declaration"));return q;}
} // namespace checkpoint_detail
inline void Owner::save(const std::filesystem::path&dir)const{
    clock_guard();model_->admit(field_,charge_);need(!std::filesystem::exists(dir),"checkpoint destination already exists");
    std::filesystem::create_directory(dir);thermal_->save(dir/"heat");
    Json j=Json::object();
    j["schema"]=api;j["configuration_sha256"]=configuration_sha256();j["motion_mesh_sha256"]=motion_.mesh_sha256;
    j["field"]=state_record(field_);j["charge"]=record(charge_);j["field_sha256"]=field_sha256_;j["heat_sha256"]=heat_sha256_;
    std::string text=j.dump();drain(j);
    std::ofstream out(dir/"checkpoint.json",std::ios::binary);out<<text<<'\n';out.close();need(bool(out),"checkpoint metadata write failed");
}
inline std::unique_ptr<Owner> Owner::restore(std::shared_ptr<const Model>model,const native::GapPicOptions&o,
    std::span<const native::ArenaRateEntry>census,const std::filesystem::path&dir){
    model=checked_model(std::move(model));const auto&g=model->mesh();const auto file=dir/"checkpoint.json";
    const auto limit=std::uintmax_t(128)*(2*std::uintmax_t(g.ne)+2*std::uintmax_t(g.nb)+7*std::uintmax_t(g.nn))+8192;
    need(std::filesystem::is_regular_file(file)&&std::filesystem::file_size(file)<=limit,"checkpoint byte cap");
    std::ifstream in(file,std::ios::binary);Json j=Json::parse(in);using namespace checkpoint_detail;
    struct Guard {Json&j;~Guard(){drain(j);}} guard{j};
    keys(j,{"schema","configuration_sha256","motion_mesh_sha256","field","charge","field_sha256","heat_sha256"});
    need(j.at("schema")==api&&hash(j.at("configuration_sha256"))==model->configuration_sha256()&&
        hash(j.at("motion_mesh_sha256"))==model->motion_mesh_sha256(o),"checkpoint field/forcing/mesh configuration changed");
    auto state=checkpoint_detail::field(j.at("field"),g);auto ledger=checkpoint_detail::charge(j.at("charge"),g);
    need(hash(j.at("field_sha256"))==model->state_sha256(state,ledger),"checkpoint field/charge seal mismatch");
    auto candidate=std::unique_ptr<Owner>(new Owner(model,std::move(state),std::move(ledger),o,census,RestoreTag{true}));
    // Every heat slice loads into a wholly private aggregate; a late corrupt slice
    // destroys this candidate without changing any accepted owner.
    candidate->thermal_->load(dir/"heat");candidate->clock_guard();
    candidate->heat_sha256_=heat_sha256_of(model->material_sha256(),*candidate->thermal_);
    need(hash(j.at("heat_sha256"))==candidate->heat_sha256_,"checkpoint thermal seal mismatch");
    return candidate;
}
} // namespace hall_corrective_v2
