// Hall motion provider v2. Derived from the sealed hall-canonical-driver-20261029/
// native_arena_motion.cpp (commit 8dddec0). The CPU arithmetic, refusal order and
// texts, census/mesh/field digests and the sealed record (schema
// canonical-full-time-total-field-motion-v2, seal over "motion-receipt-record-v1"
// + record JSON) are unchanged. Changes:
//  * the per-call mesh identity comes from a precomputed MotionContext (the
//    sealed function built a JSON value on every call); its value is identical;
//  * a device-computed motion can be assembled into the same immutable receipt
//    (MotionAssembly::assemble), sharing the host admission and the seal;
//  * every JSON value built here is cleared bottom-up before destruction, so its
//    destructor never allocates (the frozen JSON destructor allocates a
//    flattening stack under noexcept and terminates on allocation failure).
#include "motion_receipt.hpp"
#include "hall-finite-foil-em-20261029/axial_particles.hpp"
#include <bit>
#include <chrono>
#include <map>
#include <openssl/sha.h>

namespace hall_driver {
namespace {
struct Encoding {
    std::string bytes;
    void integer(std::uint64_t v){for(int i=0;i<8;++i)bytes.push_back(char((v>>(8*i))&255));}
    void scalar(double v){integer(std::bit_cast<std::uint64_t>(v));}
    void text(const std::string&s){integer(s.size());bytes.append(s);}
    void vector(const std::vector<double>&v){integer(v.size());for(double x:v)scalar(x);}
    void particle(const native::GapParticle&p){integer(p.id);for(double v:p.position_m)scalar(v);
        for(double v:p.proper_velocity_m_s)scalar(v);scalar(p.charge_C);scalar(p.mass_kg);scalar(p.weight);}
    std::string finish()const{unsigned char d[SHA256_DIGEST_LENGTH];
        SHA256(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size(),d);
        static constexpr char hex[]="0123456789abcdef";std::string s(64,'0');
        for(int i=0;i<32;++i){s[2*i]=hex[d[i]>>4];s[2*i+1]=hex[d[i]&15];}return s;}
};
bool digest(const std::string&s){return s.size()==64&&s.find_first_not_of("0123456789abcdef")==std::string::npos;}
// Clears nested JSON bottom-up; afterwards the destructor releases without allocating.
void drain(native::ArenaJson&j)noexcept{
    if(j.is_object()){auto&o=j.get_ref<native::ArenaJson::object_t&>();for(auto&item:o)drain(item.second);o.clear();}
    else if(j.is_array()){auto&a=j.get_ref<native::ArenaJson::array_t&>();for(auto&item:a)drain(item);a.clear();}
}
struct Drained {native::ArenaJson&j;~Drained(){drain(j);}};
std::vector<const native::ArenaRateEntry*> ordered(std::span<const native::ArenaRateEntry> census,std::size_t cap=100000){
    hall::require(census.size()<=cap,cap==100000?"motion census exceeds sealed native carrier cap":"motion census exceeds the configured Hall carrier cap");
    std::vector<const native::ArenaRateEntry*> rows;rows.reserve(census.size());for(const auto&row:census)rows.push_back(&row);
    std::sort(rows.begin(),rows.end(),[](auto a,auto b){return a->native.id<b->native.id;});
    std::uint64_t prior=0;for(auto row:rows){hall::require(row->native.id>prior&&!row->packet_id.empty()&&!row->state.empty(),"motion missing/duplicate native/packet identity");prior=row->native.id;}
    return rows;
}
void geometry_guard(const axial::Mesh&g,const native::GapPicOptions&o){
    hall::require(o.finite_foils&&o.finite_axial_extent&&!o.distributed_fields&&
        o.width_m<=g.side&&o.gap_m<=g.side&&o.depth_m==g.depth&&
        o.short_z_m==-g.depth/2&&o.feed_z_m==g.depth/2,"motion native/field geometry mismatch");
    (void)declared_kinetic_policy(o);
    const double roundoff=8*std::numeric_limits<double>::epsilon()*g.side;
    auto equal=[&](double a,double b){return std::isfinite(a)&&std::abs(a-b)<=roundoff;};
    for(int body=0;body<2;++body){const auto&f=body?o.upper_foil:o.lower_foil;
        hall::require(equal(f.left_m,o.width_m/2-.0025)&&equal(f.right_m,o.width_m/2+.0025)&&
            equal(f.bottom_m,o.gap_m/2+(body?.0025:-.0027))&&
            equal(f.top_m,o.gap_m/2+(body?.0027:-.0025)),"motion physical foil/frame mismatch");}
}
// The record is built with assignments and moves only (no initializer-list
// temporaries holding nested values), so no nonempty JSON temporary is ever
// destroyed on the normal path. Object keys are sorted (std::map), so the dump
// is byte-identical to the sealed initializer-list construction.
using J=native::ArenaJson;
J array_of(std::initializer_list<double> values){J a=J::array();for(double x:values)a.push_back(x);return a;}
J array_of(const native::PicVector&v){J a=J::array();for(double x:v)a.push_back(x);return a;}
J binding_json(const MotionBinding&b){J j=J::object();
    j["phase"]="full_time_position_and_proper_velocity";j["accepted_step"]=b.accepted_step;
    j["staged_arena_revision"]=b.staged_arena_revision;j["trial_stage"]=b.trial_stage;
    j["begin_s"]=b.begin_s;j["end_s"]=b.end_s;j["dt_s"]=b.dt_s;j["staged_arena_sha256"]=b.staged_arena_sha256;
    j["trial_field_sha256"]=b.trial_field_sha256;j["mesh_sha256"]=b.mesh_sha256;j["material_sha256"]=b.material_sha256;
    j["species_sha256"]=b.species_sha256;j["common_clock_sha256"]=b.common_clock_sha256;j["kinetic_policy_sha256"]=b.kinetic_policy_sha256;
    return j;}
std::string array_hash(const std::vector<double>&v){Encoding e;e.text("motion-array-binary64-v1");e.vector(v);return e.finish();}
J record_without_seal(const MotionReceipt&r){
    J endpoints=J::array(),paths=J::array();
    endpoints.get_ref<J::array_t&>().reserve(r.endpoints().size());paths.get_ref<J::array_t&>().reserve(r.paths().size());
    for(const auto&p:r.endpoints()){J e=J::object();e["native_id"]=p.id;e["position_m"]=array_of(p.position_m);
        e["proper_velocity_m_s"]=array_of(p.proper_velocity_m_s);e["charge_C"]=p.charge_C;e["mass_kg"]=p.mass_kg;e["weight"]=p.weight;
        endpoints.push_back(std::move(e));}
    for(const auto&p:r.paths()){J e=J::object();e["native_id"]=p.native_id;e["field_start_m"]=array_of(p.field_start_m);
        e["field_displacement_m"]=array_of(p.field_displacement_m);e["flight_s"]=p.flight_s;e["motion_owner"]="ordinary_total_field";
        paths.push_back(std::move(e));}
    const auto&b=r.budget();
    const auto&k=r.kinetic_policy();
    J policy=J::object();policy["schema"]=kinetic_boundary_api;policy["width_m"]=k.width_m;policy["gap_m"]=k.gap_m;
    policy["short_z_m"]=k.short_z_m;policy["feed_z_m"]=k.feed_z_m;
    policy["lower_foil"]=array_of({k.lower.left_m,k.lower.right_m,k.lower.bottom_m,k.lower.top_m});
    policy["upper_foil"]=array_of({k.upper.left_m,k.upper.right_m,k.upper.bottom_m,k.upper.top_m});
    policy["short_bridge"]="lower foil x span and inner y gap only";policy["xy_policy"]="absorbing escape";
    policy["feed_policy"]="axial escape";policy["tie_order"]="x y lower upper axial; later exact tie wins";
    J budget=J::object();budget["kinetic_change_J"]=b.kinetic_change_J;budget["field_work_J"]=b.field_work_J;
    budget["work_defect_J"]=b.work_defect_J;budget["continuity_C"]=b.continuity_C;budget["plasma_dt"]=b.plasma_dt;
    budget["maximum_gyro"]=b.maximum_gyro;budget["maximum_displacement_cells"]=b.maximum_displacement_cells;
    J j=J::object();j["schema"]=motion_api;j["binding"]=binding_json(r.binding());j["census_sha256"]=r.census_sha256();
    j["kinetic_policy"]=std::move(policy);j["native_to_field_m"]=array_of(r.frame().native_to_field_m);
    j["endpoints"]=std::move(endpoints);j["paths"]=std::move(paths);
    j["current_A_m_sha256"]=array_hash(r.current_A_m());j["initial_charge_C_sha256"]=array_hash(r.initial_charge_C());
    j["endpoint_charge_C_sha256"]=array_hash(r.endpoint_charge_C());j["budget"]=std::move(budget);
    j["coverage"]="all_native_carriers_no_sources_no_boundaries";j["neutral_motion_owner"]="same_full_time_q0_drift";
    j["native_Poisson_called"]=false;j["native_DKD_called"]=false;
    return j;
}
// Seal of the record: identical bytes to the sealed provider. Every partially
// or fully built JSON value is drained before its destructor runs.
std::string seal_of(const MotionReceipt&r){
    const auto t=std::chrono::steady_clock::now();
    native::ArenaJson record;Drained guard{record};
    record=record_without_seal(r);
    Encoding seal;seal.text("motion-receipt-record-v1");seal.text(record.dump());auto out=seal.finish();
    auto&T=hall_corrective_v2::motion_host_timers();T.record_seal_s+=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();++T.record_seals;
    return out;
}
std::string mesh_digest(const axial::Mesh&g,const native::GapPicOptions&o,const std::string&mesh_json_dump){
    geometry_guard(g,o);Encoding e;e.text("canonical-motion-geometry-v2");e.text(mesh_json_dump);
    e.integer(o.nx);e.integer(o.ny);for(double x:{o.width_m,o.gap_m,o.depth_m,o.short_z_m,o.feed_z_m})e.scalar(x);
    for(const auto&f:{o.lower_foil,o.upper_foil})for(double x:{f.left_m,f.right_m,f.bottom_m,f.top_m})e.scalar(x);
    return e.finish();
}
std::string mesh_json_of(const axial::Mesh&g){native::ArenaJson j;Drained guard{j};j=g.json();return j.dump();}
} // namespace
namespace {
std::string census_digest(const std::vector<const native::ArenaRateEntry*>&rows){
    Encoding e;e.bytes.reserve(16+rows.size()*96);e.text("canonical-motion-census-v1");e.integer(rows.size());
    for(const auto*row:rows){e.text(row->packet_id);e.text(row->state);e.particle(row->native);}return e.finish();
}
}
std::string canonical_census_sha256(std::span<const native::ArenaRateEntry> census){return census_digest(ordered(census));}
std::string motion_mesh_sha256(const axial::Mesh&g,const native::GapPicOptions&o){return mesh_digest(g,o,mesh_json_of(g));}
std::string total_field_sha256(const std::vector<double>&E,const std::vector<double>&B){
    Encoding e;e.text("total-midpoint-E-B-binary64-v1");e.vector(E);e.vector(B);return e.finish();
}
native::ArenaJson motion_record(const MotionReceipt&r){
    hall::require(digest(r.seal_sha256()),"motion record is not sealed for native admission (MotionAssembly::seal_v2)");
    auto j=record_without_seal(r);j["seal_sha256"]=r.seal_sha256();return j;}
MotionReceipt prepare_canonical_motion(const axial::Mesh&g,const native::GapPicOptions&o,
    std::span<const native::ArenaRateEntry>census,const MotionBinding&binding,
    std::span<const std::uint64_t> active_layer_ids,const std::vector<double>&E,const std::vector<double>&B){
    // Same check sequence as the sealed provider up to its mesh identity check;
    // the context (whose construction runs the geometry guard) is built after it.
    hall::require(binding.phase==ParticlePhase::full_time_position_and_proper_velocity,"motion requires owner-certified full-time phase");
    hall::require(std::isfinite(binding.begin_s)&&binding.begin_s>=0&&std::isfinite(binding.end_s)&&
        std::isfinite(binding.dt_s)&&binding.dt_s>0&&binding.dt_s<=1e-14&&binding.end_s>binding.begin_s&&
        binding.end_s==binding.begin_s+binding.dt_s,"motion common interval/10 fs guard");
    for(const auto*s:{&binding.staged_arena_sha256,&binding.trial_field_sha256,&binding.mesh_sha256,
        &binding.material_sha256,&binding.species_sha256,&binding.common_clock_sha256,&binding.kinetic_policy_sha256})hall::require(digest(*s),"motion requires sealed stage identities");
    hall::require(active_layer_ids.empty(),"motion requires source owner for every active layer; no implicit exclusion/freeze");
    hall::require(E.size()==std::size_t(g.ne)&&B.size()==std::size_t(g.nb)&&axial::finite(E)&&axial::finite(B),"motion total-field dimensions/finiteness");
    hall::require(binding.trial_field_sha256==total_field_sha256(E,B)&&binding.mesh_sha256==motion_mesh_sha256(g,o),"motion changed trial field or mesh binding");
    const auto context=hall_corrective_v2::make_motion_context(g,o,mesh_json_of(g));
    auto r=hall_corrective_v2::MotionAssembly::cpu(context,census,binding,active_layer_ids,E,B);
    hall_corrective_v2::MotionAssembly::seal_v2(r);return r;
}
} // namespace hall_driver

namespace hall_driver {
namespace {
// Binary motion v3 record: the complete content of the v2 record in canonical
// binary encoding (u64 LE, f64 bit patterns, length-prefixed text, counted
// sequences), fields in the order of the v2 record keys (sorted JSON keys).
std::string record_v3_digest(const MotionReceipt&r){
    const auto t=std::chrono::steady_clock::now();
    Encoding e;e.bytes.reserve(4096+r.endpoints().size()*80+r.paths().size()*64);
    e.text("canonical-full-time-total-field-motion-v3-binary");e.text(motion_api);
    const auto&b=r.binding();
    e.text("binding");e.text("full_time_position_and_proper_velocity");
    e.integer(b.accepted_step);e.integer(b.staged_arena_revision);e.integer(b.trial_stage);
    e.scalar(b.begin_s);e.scalar(b.end_s);e.scalar(b.dt_s);
    for(const std::string*x:{&b.staged_arena_sha256,&b.trial_field_sha256,&b.mesh_sha256,&b.material_sha256,
        &b.species_sha256,&b.common_clock_sha256,&b.kinetic_policy_sha256})e.text(*x);
    const auto&bu=r.budget();
    e.text("budget");for(double x:{bu.kinetic_change_J,bu.field_work_J,bu.work_defect_J,bu.continuity_C,bu.plasma_dt,bu.maximum_gyro,bu.maximum_displacement_cells})e.scalar(x);
    e.text("census_sha256");e.text(r.census_sha256());
    e.text("coverage");e.text("all_native_carriers_no_sources_no_boundaries");
    e.text("current_A_m_sha256");e.text(array_hash(r.current_A_m()));
    e.text("endpoint_charge_C_sha256");e.text(array_hash(r.endpoint_charge_C()));
    e.text("endpoints");e.integer(r.endpoints().size());
    for(const auto&p:r.endpoints()){e.integer(p.id);for(double v:p.position_m)e.scalar(v);for(double v:p.proper_velocity_m_s)e.scalar(v);
        e.scalar(p.charge_C);e.scalar(p.mass_kg);e.scalar(p.weight);}
    e.text("initial_charge_C_sha256");e.text(array_hash(r.initial_charge_C()));
    const auto&k=r.kinetic_policy();
    e.text("kinetic_policy");e.text(kinetic_boundary_api);
    for(double x:{k.width_m,k.gap_m,k.short_z_m,k.feed_z_m,k.lower.left_m,k.lower.right_m,k.lower.bottom_m,k.lower.top_m,
        k.upper.left_m,k.upper.right_m,k.upper.bottom_m,k.upper.top_m})e.scalar(x);
    for(const char*x:{"lower foil x span and inner y gap only","absorbing escape","axial escape","x y lower upper axial; later exact tie wins"})e.text(x);
    e.text("native_DKD_called");e.integer(0);e.text("native_Poisson_called");e.integer(0);
    e.text("native_to_field_m");for(double v:r.frame().native_to_field_m)e.scalar(v);
    e.text("neutral_motion_owner");e.text("same_full_time_q0_drift");
    e.text("paths");e.integer(r.paths().size());
    for(const auto&p:r.paths()){e.integer(p.native_id);for(double v:p.field_start_m)e.scalar(v);for(double v:p.field_displacement_m)e.scalar(v);e.scalar(p.flight_s);}
    e.text("ordinary_total_field");
    auto out=e.finish();
    auto&T=hall_corrective_v2::motion_host_timers();T.record_v3_s+=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();++T.record_v3;
    return out;
}
}
}
namespace hall_corrective_v2 {
using namespace hall_driver;
MotionHostTimers& motion_host_timers()noexcept{static MotionHostTimers timers;return timers;}
std::string canonical_census_sha256(std::span<const native::ArenaRateEntry>census,std::size_t cap){
    return hall_driver::census_digest(hall_driver::ordered(census,cap));}
void MotionAssembly::seal_v2(MotionReceipt&r){if(r.seal_sha256_.empty())r.seal_sha256_=hall_driver::seal_of(r);}
std::string MotionAssembly::record_v3(const MotionReceipt&r){return hall_driver::record_v3_digest(r);}
MotionContext make_motion_context(const axial::Mesh&g,const native::GapPicOptions&o,const std::string&mesh_json_dump,std::size_t cap){
    hall::require(cap>0,"Hall carrier cap must be positive");
    MotionContext c;c.mesh=&g;c.maximum_carriers=cap;c.options=o;c.mesh_sha256=hall_driver::mesh_digest(g,o,mesh_json_dump);
    c.policy=declared_kinetic_policy(o);c.policy_sha256=kinetic_policy_sha256(c.policy);return c;
}
void require_census_identity(std::span<const native::ArenaRateEntry>census,std::size_t cap){(void)hall_driver::ordered(census,cap);}
std::string motion_mesh_sha256(const axial::Mesh&g,const native::GapPicOptions&o,const std::string&mesh_json_dump){
    return hall_driver::mesh_digest(g,o,mesh_json_dump);
}
AdmittedCensus MotionAssembly::admit(const MotionContext&c,std::span<const native::ArenaRateEntry>census,
    const MotionBinding&binding,std::span<const std::uint64_t> active_layer_ids,const std::vector<double>&E,const std::vector<double>&B){
    const auto&g=*c.mesh;
    hall::require(binding.phase==ParticlePhase::full_time_position_and_proper_velocity,"motion requires owner-certified full-time phase");
    hall::require(std::isfinite(binding.begin_s)&&binding.begin_s>=0&&std::isfinite(binding.end_s)&&
        std::isfinite(binding.dt_s)&&binding.dt_s>0&&binding.dt_s<=1e-14&&binding.end_s>binding.begin_s&&
        binding.end_s==binding.begin_s+binding.dt_s,"motion common interval/10 fs guard");
    for(const auto*s:{&binding.staged_arena_sha256,&binding.trial_field_sha256,&binding.mesh_sha256,
        &binding.material_sha256,&binding.species_sha256,&binding.common_clock_sha256,&binding.kinetic_policy_sha256})hall::require(hall_driver::digest(*s),"motion requires sealed stage identities");
    hall::require(active_layer_ids.empty(),"motion requires source owner for every active layer; no implicit exclusion/freeze");
    hall::require(E.size()==std::size_t(g.ne)&&B.size()==std::size_t(g.nb)&&axial::finite(E)&&axial::finite(B),"motion total-field dimensions/finiteness");
    hall::require(binding.trial_field_sha256==total_field_sha256(E,B)&&binding.mesh_sha256==c.mesh_sha256,"motion changed trial field or mesh binding");
    hall::require(binding.kinetic_policy_sha256==c.policy_sha256,"motion kinetic policy changed independently of field mesh");
    const auto t=std::chrono::steady_clock::now();
    AdmittedCensus a;a.rows=hall_driver::ordered(census,c.maximum_carriers);a.census_sha256=hall_driver::census_digest(a.rows);
    auto&T=motion_host_timers();T.census_hash_s+=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();++T.census_hashes;
    return a;
}
MotionReceipt MotionAssembly::cpu(const MotionContext&c,std::span<const native::ArenaRateEntry>census,
    const MotionBinding&binding,std::span<const std::uint64_t> active_layer_ids,const std::vector<double>&E,const std::vector<double>&B){
    const auto&g=*c.mesh;const auto&o=c.options;
    auto admitted=admit(c,census,binding,active_layer_ids,E,B);
    MotionReceipt r;r.binding_=binding;r.frame_.native_to_field_m={-o.width_m/2,-o.gap_m/2,0};
    r.kinetic_policy_=c.policy;
    r.census_sha256_=std::move(admitted.census_sha256);const auto&rows=admitted.rows;
    r.current_A_m_.resize(g.ne);r.initial_charge_C_.resize(g.nn);r.endpoint_charge_C_.resize(g.nn);
    axial::Vec omega(g.nn),gross(g.nn);std::map<std::string,const native::ArenaRateEntry*> packets;
    r.endpoints_.reserve(rows.size());r.paths_.reserve(rows.size());
    for(const auto*row:rows){const auto&p=row->native;
        hall::require(std::isfinite(p.weight)&&p.weight>0&&std::isfinite(p.mass_kg)&&p.mass_kg>0&&std::isfinite(p.charge_C),"motion invalid carrier");
        for(const auto&v:{p.position_m,p.proper_velocity_m_s})for(double x:v)hall::require(std::isfinite(x),"motion nonfinite kinematics");
        require_kinetic_interior(r.kinetic_policy_,p.position_m);
        auto [where,inserted]=packets.emplace(row->packet_id,row);
        if(!inserted){const auto&q=where->second->native;hall::require(row->state==where->second->state&&p.position_m==q.position_m&&
            p.proper_velocity_m_s==q.proper_velocity_m_s&&p.charge_C==q.charge_C&&p.mass_kg==q.mass_kg,"motion divergent Count2 packet carriers");}
        auto field_x=p.position_m;for(int d=0;d<3;++d)field_x[d]+=r.frame_.native_to_field_m[d];
        double load=p.weight*p.charge_C*p.charge_C/(axial::eps*p.mass_kg);
        hall::require(std::isfinite(load)&&(p.charge_C==0||load>0),"motion unresolved plasma operand");
        axial::point_deposit(g,omega,field_x,load);
    }
    for(int n=0;n<g.nn;++n)r.budget_.plasma_dt=std::max(r.budget_.plasma_dt,binding.dt_s*std::sqrt(omega[n]/g.VN[n]));
    hall::require(std::isfinite(r.budget_.plasma_dt)&&r.budget_.plasma_dt<=.1,"motion full-census plasma guard");
    long double kinetic=0;
    for(const auto*row:rows){const auto&p=row->native;auto field_p=p;
        for(int d=0;d<3;++d)field_p.position_m[d]+=r.frame_.native_to_field_m[d];
        axial::Push moved{};
        if(p.charge_C==0){auto stream=prepare_neutral_stream(r.kinetic_policy_,p,binding.begin_s,binding.end_s,binding.dt_s);
            hall::require(stream.hit.sink==KineticSink::none,"neutral first hit requires native sink owner");
            moved.particle=field_p;moved.hit={1,-1};moved.flight_s=binding.dt_s;
            for(int d=0;d<3;++d)moved.D[d]=stream.endpoint.position_m[d]-p.position_m[d];
            moved.displacement_cells=std::sqrt(hall::dot(moved.D,moved.D))/std::min({g.dx,g.dy,g.dz});
            hall::require(moved.displacement_cells<=.2,"neutral full-interval displacement guard");
        }else moved=axial::push(g,field_p,E,B,binding.dt_s);
        hall::require(moved.hit.surface<0&&moved.flight_s==binding.dt_s,"motion needs boundary owner before collection");
        hall::require(first_kinetic_hit(r.kinetic_policy_,p.position_m,moved.D).sink==KineticSink::none,"charged first hit requires native sink owner");
        auto endpoint=p;endpoint.proper_velocity_m_s=moved.particle.proper_velocity_m_s;
        for(int d=0;d<3;++d)endpoint.position_m[d]=p.position_m[d]+moved.D[d];
        require_kinetic_interior(r.kinetic_policy_,endpoint.position_m);
        native::PicVector D{};for(int d=0;d<3;++d)D[d]=(endpoint.position_m[d]+r.frame_.native_to_field_m[d])-field_p.position_m[d];
        hall::require(axial::intersection(g,field_p.position_m,D).surface<0,"rounded canonical endpoint crosses field boundary");
        axial::Path path(g,field_p.position_m,D);double charge=p.weight*p.charge_C;
        hall::require(std::isfinite(charge)&&(p.charge_C==0||charge!=0),"motion unresolved macrocharge");
        path.deposit(r.current_A_m_,charge,D,binding.dt_s);
        axial::point_deposit(g,r.initial_charge_C_,field_p.position_m,charge);
        auto field_end=endpoint.position_m;for(int d=0;d<3;++d)field_end[d]+=r.frame_.native_to_field_m[d];
        axial::point_deposit(g,r.endpoint_charge_C_,field_end,charge);
        axial::point_deposit(g,gross,field_p.position_m,std::abs(charge));
        kinetic+=static_cast<long double>(p.weight)*(pic_kinetic_energy_J(endpoint)-pic_kinetic_energy_J(p));
        r.budget_.maximum_gyro=std::max(r.budget_.maximum_gyro,moved.gyro);
        r.budget_.maximum_displacement_cells=std::max(r.budget_.maximum_displacement_cells,moved.displacement_cells);
        r.endpoints_.push_back(endpoint);r.paths_.push_back({p.id,field_p.position_m,D,binding.dt_s});
    }
    r.budget_.kinetic_change_J=double(kinetic);r.budget_.field_work_J=binding.dt_s*axial::dot(E,r.current_A_m_);
    r.budget_.work_defect_J=r.budget_.kinetic_change_J-r.budget_.field_work_J;
    hall::require(std::isfinite(r.budget_.work_defect_J)&&std::abs(r.budget_.work_defect_J)<=
        1e-20+1e-10*(std::abs(r.budget_.kinetic_change_J)+std::abs(r.budget_.field_work_J)),"canonical motion work mismatch");
    auto div=g.gradient_transpose(r.current_A_m_);
    for(int n=0;n<g.nn;++n)r.budget_.continuity_C=std::max(r.budget_.continuity_C,
        std::abs(r.endpoint_charge_C_[n]-r.initial_charge_C_[n]-binding.dt_s*div[n]));
    hall::require(axial::finite(r.current_A_m_)&&std::isfinite(r.budget_.continuity_C)&&
        r.budget_.continuity_C<=1e-12*std::max(axial::maximum(gross),1e-30),"canonical motion charge continuity mismatch");
    r.record_v3_sha256_=hall_driver::record_v3_digest(r);
    return r;
}
MotionReceipt MotionAssembly::assemble(const MotionContext&c,const AdmittedCensus&a,const MotionBinding&binding,MotionArrays&&m){
    const auto&g=*c.mesh;const auto&o=c.options;const std::size_t n=a.rows.size();
    hall::require(m.endpoints.size()==n&&m.paths.size()==n&&m.current_A_m.size()==std::size_t(g.ne)&&
        m.initial_charge_C.size()==std::size_t(g.nn)&&m.endpoint_charge_C.size()==std::size_t(g.nn),"device motion array dimensions");
    const native::PicVector shift{-o.width_m/2,-o.gap_m/2,0};
    for(std::size_t i=0;i<n;++i){const auto&p=a.rows[i]->native;const auto&e=m.endpoints[i];const auto&path=m.paths[i];
        // The device returns kinematics only; identity, species and count stay the census values.
        hall::require(e.id==p.id&&e.charge_C==p.charge_C&&e.mass_kg==p.mass_kg&&e.weight==p.weight&&path.native_id==p.id&&
            path.flight_s==binding.dt_s,"device motion altered carrier identity, species, count or flight");
        for(int d=0;d<3;++d)hall::require(path.field_start_m[d]==p.position_m[d]+shift[d]&&
            path.field_displacement_m[d]==(e.position_m[d]+shift[d])-path.field_start_m[d],"device canonical path start/endpoint operands");
    }
    const auto&b=m.budget;
    hall::require(std::isfinite(b.work_defect_J)&&b.work_defect_J==b.kinetic_change_J-b.field_work_J&&std::abs(b.work_defect_J)<=
        1e-20+1e-10*(std::abs(b.kinetic_change_J)+std::abs(b.field_work_J))&&b.plasma_dt>=0&&b.plasma_dt<=.1&&
        b.maximum_gyro>=0&&b.maximum_gyro<=.05&&b.maximum_displacement_cells>=0&&b.maximum_displacement_cells<=.2&&
        std::isfinite(b.continuity_C)&&b.continuity_C>=0,"device motion budget outside the sealed guards");
    MotionReceipt r;r.binding_=binding;r.frame_.native_to_field_m=shift;r.kinetic_policy_=c.policy;r.census_sha256_=a.census_sha256;
    r.budget_=b;r.endpoints_=std::move(m.endpoints);r.paths_=std::move(m.paths);r.current_A_m_=std::move(m.current_A_m);
    r.initial_charge_C_=std::move(m.initial_charge_C);r.endpoint_charge_C_=std::move(m.endpoint_charge_C);
    r.record_v3_sha256_=hall_driver::record_v3_digest(r);
    return r;
}
} // namespace hall_corrective_v2
