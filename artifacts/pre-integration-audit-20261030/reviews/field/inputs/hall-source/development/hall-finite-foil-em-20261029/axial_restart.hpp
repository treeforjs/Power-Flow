#pragma once
#include "axial_system.hpp"

namespace axial {
struct Ledger {
    std::array<double,8> work{}; // dU, dK, Joule, circuit, outer loss, outer supply, dH, thermal numerical residual
    Vec last_current_A_m;double last_feed_A=0,last_terminal_V=0;std::uint64_t next_id=1;
    void observe(const JointTrial&t){const auto&r=t.receipt;
        std::array<double,8>inc{r.field_delta_J,t.transport.kinetic_change_J,r.joule_J,r.circuit_work_J,
            r.outer_loss_J,r.outer_supply_J,t.heat.thermal_delta_J,t.heat.numerical_residual_J};
        for(int i=0;i<8;++i)work[i]+=inc[i];last_current_A_m=t.transport.current_load_A_m;
        last_feed_A=r.current_A;last_terminal_V=r.terminal_V;
    }
};
inline void save_checkpoint(const std::filesystem::path&target,const Mesh&g,const State&field,
    std::span<const GapParticle>arena,const Thermal&thermal,const Ledger&ledger){
    hall::require(!std::filesystem::exists(target)&&field.time_s==thermal.time_s()&&field.revision==thermal.revision(),"checkpoint destination or common clock");
    auto staging=target;staging+=".staging";hall::require(!std::filesystem::exists(staging),"checkpoint staging exists");std::filesystem::create_directory(staging);
    Json particles=Json::array();for(const auto&p:arena)particles.push_back({{"id",p.id},{"x",p.position_m},{"u",p.proper_velocity_m_s},{"q",p.charge_C},{"m",p.mass_kg},{"weight",p.weight}});
    Json state{{"E",field.E},{"B",field.B},{"E_low",field.E_low},{"B_low",field.B_low},{"time_s",field.time_s},{"revision",field.revision},{"arena",particles},
        {"work",ledger.work},{"last_current_A_m",ledger.last_current_A_m},{"last_feed_A",ledger.last_feed_A},
        {"last_terminal_V",ledger.last_terminal_V},{"next_id",ledger.next_id}};
    std::ofstream data(staging/"state.json");data<<state.dump()<<'\n';data.close();hall::require(bool(data),"checkpoint data write");
    thermal.save(staging/"thermal");Json hashes=Json::object();
    for(const auto&entry:std::filesystem::recursive_directory_iterator(staging))if(entry.is_regular_file())
        hashes[std::filesystem::relative(entry.path(),staging).string()]=hall::file_checksum(entry.path());
    Json manifest{{"schema","axial-total-field-native-thermal-checkpoint-v2"},{"mesh",g.json()},
        {"material_seal",ss304_source_seal_sha256},{"rho_ref_kg_m3",8000.},{"mode","evaluated_solid_1500"},
        {"field_owner","Maxwell total only"},{"inventory_owner","caller canonical arena"},{"particle_phase","full_time_work_conjugate_v1"},{"executable_fnv1a64",hall::file_checksum("/proc/self/exe")},
        {"file_fnv1a64",hashes},{"event_source_scope","Joule-only pilot; no reaction/emission event journal"}};
    std::ofstream meta(staging/"manifest.json");meta<<manifest.dump(2)<<'\n';meta.close();hall::require(bool(meta),"checkpoint manifest write");
    std::filesystem::rename(staging,target);
}
inline void load_checkpoint(const std::filesystem::path&source,const Mesh&g,State&field,
    std::vector<GapParticle>&arena,Thermal&thermal,Ledger&ledger){
    hall::require(field.revision==0&&thermal.revision()==0&&arena.empty(),"restart requires fresh component owners");
    std::ifstream m(source/"manifest.json");Json manifest;m>>manifest;
    hall::require(manifest.at("schema")=="axial-total-field-native-thermal-checkpoint-v2"&&manifest.at("mesh")==g.json()&&
        manifest.at("material_seal")==ss304_source_seal_sha256&&manifest.at("rho_ref_kg_m3")==8000.&&
        manifest.at("mode")=="evaluated_solid_1500"&&manifest.at("particle_phase")=="full_time_work_conjugate_v1"&&manifest.at("executable_fnv1a64")==hall::file_checksum("/proc/self/exe"),"restart code/mesh/material/phase mismatch");
    for(auto it=manifest.at("file_fnv1a64").begin();it!=manifest.at("file_fnv1a64").end();++it){
        std::filesystem::path rel=it.key();hall::require(!rel.is_absolute()&&rel.string().find("..") == std::string::npos,"invalid checkpoint file path");
        hall::require(hall::file_checksum(source/rel)==it.value().get<std::uint64_t>(),"checkpoint file checksum mismatch");}
    std::ifstream f(source/"state.json");Json state;f>>state;State candidate;
    candidate.E=state.at("E").get<Vec>();candidate.B=state.at("B").get<Vec>();candidate.E_low=state.at("E_low").get<Vec>();candidate.B_low=state.at("B_low").get<Vec>();candidate.time_s=state.at("time_s");candidate.revision=state.at("revision");
    hall::require(candidate.E.size()==std::size_t(g.ne)&&candidate.B.size()==std::size_t(g.nb)&&candidate.E_low.size()==std::size_t(g.ne)&&candidate.B_low.size()==std::size_t(g.nb)&&
        finite(candidate.E)&&finite(candidate.B)&&finite(candidate.E_low)&&finite(candidate.B_low)&&
        std::isfinite(candidate.time_s)&&candidate.time_s>=field.time_s,"invalid restored field");
    std::vector<GapParticle>next;std::set<std::uint64_t>ids;std::uint64_t largest=0;
    for(const auto&v:state.at("arena")){GapParticle p;p.id=v.at("id");p.position_m=v.at("x").get<V3>();p.proper_velocity_m_s=v.at("u").get<V3>();
        p.charge_C=v.at("q");p.mass_kg=v.at("m");p.weight=v.at("weight");
        hall::require(p.id>0&&ids.insert(p.id).second&&p.mass_kg>0&&p.weight>0&&p.charge_C!=0&&
            std::isfinite(p.mass_kg)&&std::isfinite(p.weight)&&std::isfinite(p.charge_C),"invalid restored particle");
        for(auto values:{p.position_m,p.proper_velocity_m_s})for(double x:values)hall::require(std::isfinite(x),"nonfinite restored particle");
        largest=std::max(largest,p.id);next.push_back(p);}
    Ledger l;l.work=state.at("work").get<std::array<double,8>>();l.last_current_A_m=state.at("last_current_A_m").get<Vec>();
    l.last_feed_A=state.at("last_feed_A");l.last_terminal_V=state.at("last_terminal_V");l.next_id=state.at("next_id");
    hall::require(l.next_id>largest&&l.last_current_A_m.size()==std::size_t(g.ne)&&finite(l.last_current_A_m)&&
        std::isfinite(l.last_feed_A)&&std::isfinite(l.last_terminal_V),"invalid restored current/ID ledger");
    for(double x:l.work)hall::require(std::isfinite(x),"nonfinite restored work");
    thermal.load(source/"thermal");hall::require(thermal.time_s()==candidate.time_s&&thermal.revision()==candidate.revision,"restored owner clocks mismatch");
    field=std::move(candidate);arena=std::move(next);ledger=std::move(l);
}
} // namespace axial
