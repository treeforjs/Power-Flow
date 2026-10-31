// Hall motion provider v3 (first-hit kinetic sinks and surface births on top of v2).
// v3: a charged carrier whose represented path hits a field surface or the declared
// kinetic boundary within the step is truncated at its first kinetic hit by
// sink_core.hpp (shared with the GPU provider); every other carrier follows the v2
// (sealed) arithmetic unchanged. Without sinks and births the receipt, record,
// seal and v3 digest are those of Hall v2.
// v2 description: derived from the sealed hall-canonical-driver-20261029/
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
#include "sink_core.hpp"
#include "sink_ledger.hpp"
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
    if(!r.v2_identical())j=hall_corrective_v3::MotionAssembly::extend_json_record(r,std::move(j));
    return j;
}
// Seal of the record: identical bytes to the sealed provider. Every partially
// or fully built JSON value is drained before its destructor runs.
std::string seal_of(const MotionReceipt&r){
    const auto t=std::chrono::steady_clock::now();
    native::ArenaJson record;Drained guard{record};
    record=record_without_seal(r);
    Encoding seal;seal.text("motion-receipt-record-v1");seal.text(record.dump());auto out=seal.finish();
    auto&T=hall_corrective_v3::motion_host_timers();T.record_seal_s+=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();++T.record_seals;
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
    const auto context=hall_corrective_v3::make_motion_context(g,o,mesh_json_of(g));
    auto r=hall_corrective_v3::MotionAssembly::cpu(context,census,binding,active_layer_ids,E,B);
    hall_corrective_v3::MotionAssembly::seal_v2(r);return r;
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
    const bool v2=r.v2_identical();
    e.text(v2?"canonical-full-time-total-field-motion-v3-binary":hall_corrective_v3::motion_v3_flows_schema);
    e.text(motion_api);
    const auto&b=r.binding();
    e.text("binding");e.text("full_time_position_and_proper_velocity");
    e.integer(b.accepted_step);e.integer(b.staged_arena_revision);e.integer(b.trial_stage);
    e.scalar(b.begin_s);e.scalar(b.end_s);e.scalar(b.dt_s);
    for(const std::string*x:{&b.staged_arena_sha256,&b.trial_field_sha256,&b.mesh_sha256,&b.material_sha256,
        &b.species_sha256,&b.common_clock_sha256,&b.kinetic_policy_sha256})e.text(*x);
    const auto&bu=r.budget();
    e.text("budget");for(double x:{bu.kinetic_change_J,bu.field_work_J,bu.work_defect_J,bu.continuity_C,bu.plasma_dt,bu.maximum_gyro,bu.maximum_displacement_cells})e.scalar(x);
    e.text("census_sha256");e.text(r.census_sha256());
    e.text("coverage");e.text(v2?"all_native_carriers_no_sources_no_boundaries":hall_corrective_v3::flows_coverage);
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
    if(!v2){
        // native v6 flows section (API-006.txt section 4): sinks, then births
        e.text("sinks");e.integer(r.sinks().size());
        for(const auto&k:r.sinks()){e.integer(k.native_id);e.integer(std::uint64_t(k.kind));e.integer(std::uint64_t(k.boundary));e.integer(std::uint64_t(k.axis));
            e.scalar(k.plane_m);e.scalar(k.fraction);for(const auto*v:{&k.segment_m,&k.hit_point_m,&k.incident_proper_velocity_m_s})for(double x:*v)e.scalar(x);
            e.scalar(k.flight_s);e.scalar(k.discarded_flight_s);}
        e.text("births");e.integer(r.births().size());
        for(const auto&b:r.births()){e.integer(b.native_id);e.scalar(b.birth_time_s);e.scalar(b.flight_s);}
    }
    auto out=e.finish();
    auto&T=hall_corrective_v3::motion_host_timers();T.record_v3_s+=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();++T.record_v3;
    return out;
}
}
}
namespace hall_corrective_v3 {
using namespace hall_driver;
MotionHostTimers& motion_host_timers()noexcept{static MotionHostTimers timers;return timers;}
std::string canonical_census_sha256(std::span<const native::ArenaRateEntry>census,std::size_t cap){
    return hall_driver::census_digest(hall_driver::ordered(census,cap));}
void MotionAssembly::seal_v2(MotionReceipt&r){if(r.seal_sha256_.empty())r.seal_sha256_=hall_driver::seal_of(r);}
std::string MotionAssembly::record_v3(const MotionReceipt&r){return hall_driver::record_v3_digest(r);}
namespace {
hgp::Mesh sink_mesh(const axial::Mesh&g){return {g.nx,g.ny,g.nz,g.nex,g.ney,g.nbx,g.nby,g.nn,g.ne,g.nb,g.side,g.depth,g.dx,g.dy,g.dz,g.zmin};}
hgp::V3 v3of(const native::PicVector&a){return {{a[0],a[1],a[2]}};}
native::PicVector picof(const hgp::V3&a){return {a.x[0],a.x[1],a.x[2]};}
}
SinkRecord make_sink_record(const native::ArenaRateEntry&row,const native::GapParticle&endpoint,const hgp::SinkOut&so,double K0,double K1){
    SinkRecord k;const auto&p=row.native;k.native_id=p.id;k.packet_id=row.packet_id;k.state=row.state;k.boundary=so.sink;
    k.collected=so.sink==0||so.sink==1||so.sink==3;k.kind=k.collected?1:2;k.conductor=so.sink==3?2:(k.collected?so.sink:-1);k.axis=so.axis;k.plane_m=so.plane;
    k.fraction=so.fraction;k.segment_m=picof(so.segment);k.hit_point_m=endpoint.position_m;k.field_point_m=picof(so.field_point);
    k.incident_proper_velocity_m_s=endpoint.proper_velocity_m_s;k.flight_s=so.flight;k.discarded_flight_s=so.discarded;
    k.charge_C=p.charge_C;k.mass_kg=p.mass_kg;k.weight=p.weight;k.kinetic_J=K1;k.weighted_kinetic_J=p.weight*K1;k.initial_weighted_kinetic_J=p.weight*K0;
    k.gyro=so.gyro;k.displacement_cells=so.dcells;return k;
}
namespace {
hgp::SinkGeometry sink_geometry(const MotionContext&c){
    hgp::SinkGeometry s{};const auto&k=c.policy;
    s.p={k.width_m,k.gap_m,k.short_z_m,k.feed_z_m,{k.lower.left_m,k.lower.right_m,k.lower.bottom_m,k.lower.top_m},
         {k.upper.left_m,k.upper.right_m,k.upper.bottom_m,k.upper.top_m}};
    s.off[0]=-c.options.width_m/2;s.off[1]=-c.options.gap_m/2;s.off[2]=0;
    for(int a=0;a<5;++a){for(int d=0;d<3;++d){s.lo[a][d]=c.sink_lo[a][d];s.hi[a][d]=c.sink_hi[a][d];}s.supported[a]=c.sink_supported[a];}
    s.aligned=c.sink_support_aligned?1:0;return s;
}
// Node-index support of every sink kind; a kind is supported when every node of
// its range may hold the bank charge (conductor: node_body>=0 or short node;
// exterior: EM exterior node or target-end plane).
void fill_sink_support(MotionContext&c){
    const auto&g=*c.mesh;const auto m=sink_mesh(g);const auto&k=c.policy;
    const double off[3]={-c.options.width_m/2,-c.options.gap_m/2,0};
    bool aligned=true;
    auto plane=[&](int d,double v){int p=0;if(hgp::nearest_plane(m,d,v,off[d],p)!=hgp::ok){aligned=false;return 0;}return p;};
    const int xl=plane(0,k.lower.left_m),xr=plane(0,k.lower.right_m),x0=plane(0,0.),xw=plane(0,k.width_m);
    const int yb0=plane(1,k.lower.bottom_m),yt0=plane(1,k.lower.top_m),yb1=plane(1,k.upper.bottom_m),yt1=plane(1,k.upper.top_m),y0=plane(1,0.),yg=plane(1,k.gap_m);
    const int zs=plane(2,k.short_z_m),zf=plane(2,k.feed_z_m);
    const int ranges[5][3][2]={{{xl,xr},{yb0,yt0},{zs,zf}},{{xl,xr},{yb1,yt1},{zs,zf}},{{x0,xw},{y0,yg},{zs,zf}},{{xl,xr},{yt0,yb1},{zs,zs}},
                               {{x0,xw},{y0,yg},{zs,zf}}};
    c.sink_support_aligned=aligned&&plane(0,k.upper.left_m)==xl&&plane(0,k.upper.right_m)==xr;
    std::ostringstream text;text<<"aligned="<<c.sink_support_aligned;
    for(int a=0;a<5;++a){bool ok=c.sink_support_aligned;
        for(int d=0;d<3;++d){c.sink_lo[a][d]=ranges[a][d][0];c.sink_hi[a][d]=ranges[a][d][1];ok=ok&&ranges[a][d][0]<=ranges[a][d][1];}
        const bool collect=a==0||a==1||a==3;
        if(ok){// The hit plane of a kind: foils any face; xy escape the four outer planes; bridge/axial the end planes.
            for(int kk=c.sink_lo[a][2];ok&&kk<=c.sink_hi[a][2];++kk)for(int j=c.sink_lo[a][1];ok&&j<=c.sink_hi[a][1];++j)for(int i=c.sink_lo[a][0];ok&&i<=c.sink_hi[a][0];++i){
                const bool face=a<2||a==3?true:a==2?(i==x0||i==xw||j==y0||j==yg):(kk==zs||kk==zf);
                if(!face)continue;const int n=g.node(i,j,kk);
                ok=collect?(g.node_body[n]>=0||g.short_node[n]):(g.outer_node[n]||kk==g.kshort||kk==g.kfeed);}}
        c.sink_supported[a]=ok?1:0;text<<";"<<sink_boundary_names[a]<<"="<<ok<<":";
        for(int d=0;d<3;++d)text<<"["<<c.sink_lo[a][d]<<","<<c.sink_hi[a][d]<<"]";}
    c.sink_support_text=text.str();
}
}
MotionContext make_motion_context(const axial::Mesh&g,const native::GapPicOptions&o,const std::string&mesh_json_dump,std::size_t cap){
    hall::require(cap>0,"Hall carrier cap must be positive");
    MotionContext c;c.mesh=&g;c.maximum_carriers=cap;c.options=o;c.mesh_sha256=hall_driver::mesh_digest(g,o,mesh_json_dump);
    c.policy=declared_kinetic_policy(o);c.policy_sha256=kinetic_policy_sha256(c.policy);fill_sink_support(c);return c;
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
MotionAssembly::Flights MotionAssembly::flights(const MotionContext&c,const AdmittedCensus&a,const MotionBinding&binding,std::span<const Newborn>newborns){
    Flights f;const std::size_t n=a.rows.size();f.flight.assign(n,binding.dt_s);f.tau.assign(n,binding.begin_s);
    f.face_body.assign(n,-1);f.face_axis.assign(n,-1);f.newborn.assign(n,0);
    if(newborns.empty())return f;
    std::map<std::uint64_t,double> tau;
    for(const auto&b:newborns){hall::require(tau.emplace(b.native_id,b.birth_time_s).second,"newborn listed twice");
        hall::require(std::isfinite(b.birth_time_s)&&b.birth_time_s>=binding.begin_s&&b.birth_time_s<binding.end_s,"newborn birth time outside the step or without remaining flight");}
    const auto&k=c.policy;const hgp::Policy pol{k.width_m,k.gap_m,k.short_z_m,k.feed_z_m,{k.lower.left_m,k.lower.right_m,k.lower.bottom_m,k.lower.top_m},
                                                {k.upper.left_m,k.upper.right_m,k.upper.bottom_m,k.upper.top_m}};
    for(std::size_t i=0;i<n;++i){const auto it=tau.find(a.rows[i]->native.id);if(it==tau.end())continue;
        const auto&p=a.rows[i]->native;f.newborn[i]=1;f.tau[i]=it->second;f.flight[i]=binding.end_s-it->second;++f.count;
        hall::require(f.flight[i]>0,"newborn birth time outside the step or without remaining flight");
        int body=-1,axis=-1;const int code=hgp::newborn_start_rule(pol,v3of(p.position_m),v3of(p.proper_velocity_m_s),body,axis);
        hall::require(code==hgp::ok,hgp::sink_code_message(code));f.face_body[i]=body;f.face_axis[i]=axis;}
    hall::require(f.count==tau.size(),"newborn native ID absent from the census");
    return f;
}
MotionReceipt MotionAssembly::cpu(const MotionContext&c,std::span<const native::ArenaRateEntry>census,
    const MotionBinding&binding,std::span<const std::uint64_t> active_layer_ids,const std::vector<double>&E,const std::vector<double>&B,
    std::span<const Newborn> newborns,bool birth_scope){
    const auto&g=*c.mesh;const auto&o=c.options;
    auto admitted=admit(c,census,binding,active_layer_ids,E,B);
    const auto fl=flights(c,admitted,binding,newborns);
    MotionReceipt r;r.binding_=binding;r.frame_.native_to_field_m={-o.width_m/2,-o.gap_m/2,0};
    r.kinetic_policy_=c.policy;r.birth_scope_=birth_scope||!newborns.empty();
    r.census_sha256_=std::move(admitted.census_sha256);const auto&rows=admitted.rows;
    r.current_A_m_.resize(g.ne);r.initial_charge_C_.resize(g.nn);r.endpoint_charge_C_.resize(g.nn);
    axial::Vec omega(g.nn),gross(g.nn);std::map<std::string,const native::ArenaRateEntry*> packets;
    const auto hm=sink_mesh(g);const auto sg=sink_geometry(c);
    r.endpoints_.reserve(rows.size());r.paths_.reserve(rows.size());
    for(std::size_t i=0;i<rows.size();++i){const auto*row=rows[i];const auto&p=row->native;
        hall::require(std::isfinite(p.weight)&&p.weight>0&&std::isfinite(p.mass_kg)&&p.mass_kg>0&&std::isfinite(p.charge_C),"motion invalid carrier");
        for(const auto&v:{p.position_m,p.proper_velocity_m_s})for(double x:v)hall::require(std::isfinite(x),"motion nonfinite kinematics");
        if(fl.face_body[i]<0)require_kinetic_interior(r.kinetic_policy_,p.position_m);   // a newborn face start was admitted by flights()
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
    // Truncated (sink) carrier: path to the hit point, kinetic change to the incident velocity.
    auto sink_out=[&](const native::ArenaRateEntry&row,const native::GapParticle&field_p,const hgp::SinkOut&so){
        const auto&p=row.native;
        auto endpoint=p;endpoint.position_m=picof(so.point);endpoint.proper_velocity_m_s=picof(so.u);
        native::PicVector D{};for(int d=0;d<3;++d)D[d]=(endpoint.position_m[d]+r.frame_.native_to_field_m[d])-field_p.position_m[d];
        axial::Path path(g,field_p.position_m,D);double charge=p.weight*p.charge_C;
        hall::require(std::isfinite(charge)&&(p.charge_C==0||charge!=0),"motion unresolved macrocharge");
        path.deposit(r.current_A_m_,charge,D,binding.dt_s);
        axial::point_deposit(g,r.initial_charge_C_,field_p.position_m,charge);
        auto field_end=endpoint.position_m;for(int d=0;d<3;++d)field_end[d]+=r.frame_.native_to_field_m[d];
        axial::point_deposit(g,r.endpoint_charge_C_,field_end,charge);
        axial::point_deposit(g,gross,field_p.position_m,std::abs(charge));
        const double K1=pic_kinetic_energy_J(endpoint),K0=pic_kinetic_energy_J(p);
        kinetic+=static_cast<long double>(p.weight)*(K1-K0);
        r.budget_.maximum_gyro=std::max(r.budget_.maximum_gyro,so.gyro);
        r.budget_.maximum_displacement_cells=std::max(r.budget_.maximum_displacement_cells,so.dcells);
        r.paths_.push_back({p.id,field_p.position_m,D,so.flight});
        r.sinks_.push_back(hall_corrective_v3::make_sink_record(row,endpoint,so,K0,K1));
    };
    for(std::size_t i=0;i<rows.size();++i){const auto*row=rows[i];const auto&p=row->native;auto field_p=p;
        for(int d=0;d<3;++d)field_p.position_m[d]+=r.frame_.native_to_field_m[d];
        const double F=fl.flight[i];const int face=fl.face_body[i];
        if(fl.newborn[i]){
            NewbornRecord nb;nb.native_id=p.id;nb.birth_time_s=fl.tau[i];nb.flight_s=F;nb.origin=face;nb.axis=fl.face_axis[i];nb.start_m=p.position_m;
            nb.packet_id=row->packet_id;nb.state=row->state;nb.charge_C=p.charge_C;nb.weight=p.weight;nb.birth_weighted_kinetic_J=p.weight*pic_kinetic_energy_J(p);
            r.births_.push_back(std::move(nb));
        }
        axial::Push moved{};hgp::SinkOut so{};bool sink=false;
        if(fl.newborn[i]){
            // Newborn over [tau, end] from its birth position (ports of the sealed arithmetic, face-aware first hit).
            if(p.charge_C==0){hgp::V3 end{};double dcells=0;
                const int code=hgp::neutral_drift(hm,sg,v3of(p.position_m),v3of(p.proper_velocity_m_s),F,face,so,end,dcells);
                hall::require(code==hgp::ok,hgp::sink_code_message(code));
                if(so.sink!=hgp::sink_none)sink=true;
                else{moved.particle=field_p;moved.hit={1,-1};moved.flight_s=F;for(int d=0;d<3;++d)moved.D[d]=end.x[d]-p.position_m[d];
                    moved.displacement_cells=dcells;hall::require(dcells<=.2,"neutral full-interval displacement guard");}
            }else{hgp::PushResult pr{};
                int code=hgp::push(hm,v3of(field_p.position_m),v3of(p.proper_velocity_m_s),p.charge_C,p.mass_kg,E.data(),B.data(),F,pr);
                hall::require(code==hgp::ok,hgp::sink_code_message(code));
                hgp::KineticHit kh{};code=hgp::first_hit_from(sg.p,v3of(p.position_m),pr.D,face,kh);hall::require(code==hgp::ok,hgp::sink_code_message(code));
                if(!(pr.hit.surface<0&&pr.flight_s==F)||kh.sink!=hgp::sink_none){
                    code=hgp::sink_truncate(hm,sg,v3of(p.position_m),v3of(field_p.position_m),v3of(p.proper_velocity_m_s),p.charge_C,p.mass_kg,E.data(),B.data(),F,so,face);
                    hall::require(code==hgp::ok,hgp::sink_code_message(code));sink=true;}
                else{moved.particle=field_p;moved.particle.proper_velocity_m_s=picof(pr.u);moved.D=picof(pr.D);moved.hit={1,-1};moved.flight_s=F;
                    moved.gyro=pr.gyro;moved.displacement_cells=pr.displacement_cells;}
            }
        }else if(p.charge_C==0){auto stream=prepare_neutral_stream(r.kinetic_policy_,p,binding.begin_s,binding.end_s,binding.dt_s);
            if(stream.hit.sink!=KineticSink::none){   // v3: neutral first hit -> sink (v2 refused)
                hgp::V3 end{};double dcells=0;const int code=hgp::neutral_drift(hm,sg,v3of(p.position_m),v3of(p.proper_velocity_m_s),F,-1,so,end,dcells);
                hall::require(code==hgp::ok&&so.sink!=hgp::sink_none,hgp::sink_code_message(code==hgp::ok?hgp::lost_kinetic_surface:code));sink=true;}
            else{moved.particle=field_p;moved.hit={1,-1};moved.flight_s=binding.dt_s;
                for(int d=0;d<3;++d)moved.D[d]=stream.endpoint.position_m[d]-p.position_m[d];
                moved.displacement_cells=std::sqrt(hall::dot(moved.D,moved.D))/std::min({g.dx,g.dy,g.dz});
                hall::require(moved.displacement_cells<=.2,"neutral full-interval displacement guard");}
        }else{moved=axial::push(g,field_p,E,B,binding.dt_s);
            // v3: first-hit truncation of a charged carrier whose represented path hits
            // a field surface or the kinetic boundary (sink_core.hpp, shared with the GPU).
            if(!(moved.hit.surface<0&&moved.flight_s==binding.dt_s)||first_kinetic_hit(r.kinetic_policy_,p.position_m,moved.D).sink!=KineticSink::none){
                const int code=hgp::sink_truncate(hm,sg,v3of(p.position_m),v3of(field_p.position_m),v3of(p.proper_velocity_m_s),p.charge_C,p.mass_kg,
                                                  E.data(),B.data(),binding.dt_s,so);
                hall::require(code==hgp::ok,hgp::sink_code_message(code));sink=true;}
        }
        if(sink){sink_out(*row,field_p,so);continue;}
        hall::require(moved.hit.surface<0&&moved.flight_s==F,"motion needs boundary owner before collection");
        if(!fl.newborn[i])hall::require(first_kinetic_hit(r.kinetic_policy_,p.position_m,moved.D).sink==KineticSink::none,"charged first hit requires native sink owner");
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
        r.endpoints_.push_back(endpoint);r.paths_.push_back({p.id,field_p.position_m,D,F});
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
    std::map<std::uint64_t,double> newborn_flight;for(const auto&b:m.births)newborn_flight.emplace(b.native_id,b.flight_s);
    hall::require(m.endpoints.size()+m.sinks.size()==n&&m.paths.size()==n&&m.current_A_m.size()==std::size_t(g.ne)&&
        m.initial_charge_C.size()==std::size_t(g.nn)&&m.endpoint_charge_C.size()==std::size_t(g.nn),"device motion array dimensions");
    const native::PicVector shift{-o.width_m/2,-o.gap_m/2,0};
    std::size_t ie=0,is=0;
    for(std::size_t i=0;i<n;++i){const auto&p=a.rows[i]->native;const auto&path=m.paths[i];
        const bool sink=is<m.sinks.size()&&m.sinks[is].native_id==p.id;
        native::PicVector end{};
        if(sink){const auto&k=m.sinks[is++];
            hall::require(k.packet_id==a.rows[i]->packet_id&&k.state==a.rows[i]->state&&k.charge_C==p.charge_C&&k.mass_kg==p.mass_kg&&k.weight==p.weight&&
                path.native_id==p.id&&path.flight_s==k.flight_s,"device sink altered carrier identity, species, count or flight");
            const auto nb=newborn_flight.find(p.id);const double F=nb==newborn_flight.end()?binding.dt_s:nb->second;
            hall::require(k.flight_s==F*k.fraction&&k.discarded_flight_s==F-k.flight_s,"device sink flight differs from F*fraction");
            end=k.hit_point_m;}
        else{hall::require(ie<m.endpoints.size(),"device motion endpoint order");const auto&e=m.endpoints[ie++];
            // The device returns kinematics only; identity, species and count stay the census values.
            const auto nb=newborn_flight.find(p.id);
            hall::require(e.id==p.id&&e.charge_C==p.charge_C&&e.mass_kg==p.mass_kg&&e.weight==p.weight&&path.native_id==p.id&&
                path.flight_s==(nb==newborn_flight.end()?binding.dt_s:nb->second),"device motion altered carrier identity, species, count or flight");
            end=e.position_m;}
        for(int d=0;d<3;++d)hall::require(path.field_start_m[d]==p.position_m[d]+shift[d]&&
            path.field_displacement_m[d]==(end[d]+shift[d])-path.field_start_m[d],"device canonical path start/endpoint operands");
    }
    hall::require(ie==m.endpoints.size()&&is==m.sinks.size(),"device motion sink/endpoint partition");
    const auto&b=m.budget;
    hall::require(std::isfinite(b.work_defect_J)&&b.work_defect_J==b.kinetic_change_J-b.field_work_J&&std::abs(b.work_defect_J)<=
        1e-20+1e-10*(std::abs(b.kinetic_change_J)+std::abs(b.field_work_J))&&b.plasma_dt>=0&&b.plasma_dt<=.1&&
        b.maximum_gyro>=0&&b.maximum_gyro<=.05&&b.maximum_displacement_cells>=0&&b.maximum_displacement_cells<=.2&&
        std::isfinite(b.continuity_C)&&b.continuity_C>=0,"device motion budget outside the sealed guards");
    MotionReceipt r;r.binding_=binding;r.frame_.native_to_field_m=shift;r.kinetic_policy_=c.policy;r.census_sha256_=a.census_sha256;
    r.budget_=b;r.endpoints_=std::move(m.endpoints);r.paths_=std::move(m.paths);r.current_A_m_=std::move(m.current_A_m);
    r.initial_charge_C_=std::move(m.initial_charge_C);r.endpoint_charge_C_=std::move(m.endpoint_charge_C);r.sinks_=std::move(m.sinks);r.births_=std::move(m.births);
    r.birth_scope_=m.birth_scope||!r.births_.empty();
    r.record_v3_sha256_=hall_driver::record_v3_digest(r);
    return r;
}

// ---------------------------------------------------------------- v3 record sections
namespace {
J vec(const native::PicVector&v){return array_of(v);}
J sink_json(const SinkRecord&k){J j=J::object();
    j["record"]=sink_record_api;j["native_id"]=k.native_id;j["kind"]=k.kind;j["kind_name"]=k.kind==1?"collected_at_conductor":"escaped_exterior";
    j["surface"]=k.boundary;j["surface_name"]=sink_boundary_names[k.boundary];j["axis"]=k.axis;j["plane_m"]=k.plane_m;j["fraction"]=k.fraction;
    j["segment_m"]=vec(k.segment_m);j["hit_point_m"]=vec(k.hit_point_m);j["incident_proper_velocity_m_s"]=vec(k.incident_proper_velocity_m_s);
    j["flight_s"]=k.flight_s;j["discarded_flight_s"]=k.discarded_flight_s;j["packet_id"]=k.packet_id;j["state"]=k.state;
    j["conductor"]=k.collected?std::string(conductor_names[k.conductor]):std::string("none");j["field_point_m"]=vec(k.field_point_m);
    j["charge_C"]=k.charge_C;j["mass_kg"]=k.mass_kg;j["weight"]=k.weight;j["kinetic_J"]=k.kinetic_J;j["weighted_kinetic_J"]=k.weighted_kinetic_J;
    j["initial_weighted_kinetic_J"]=k.initial_weighted_kinetic_J;return j;}
J birth_json(const NewbornRecord&b){J j=J::object();
    j["native_id"]=b.native_id;j["birth_time_s"]=b.birth_time_s;j["flight_s"]=b.flight_s;j["origin"]=b.origin<0?std::string("interior"):std::string(conductor_names[b.origin])+"-face";
    j["axis"]=b.axis;j["start_m"]=vec(b.start_m);j["packet_id"]=b.packet_id;j["state"]=b.state;j["charge_C"]=b.charge_C;j["weight"]=b.weight;
    j["birth_weighted_kinetic_J"]=b.birth_weighted_kinetic_J;return j;}
}
native::ArenaJson MotionAssembly::extend_json_record(const MotionReceipt&r,native::ArenaJson j){
    j["schema"]=motion_api_sinks;j["coverage"]=coverage_sinks;
    J sinks=J::array(),births=J::array();for(const auto&k:r.sinks())sinks.push_back(sink_json(k));for(const auto&b:r.births())births.push_back(birth_json(b));
    j["sinks"]=std::move(sinks);j["births"]=std::move(births);
    {J l=SinkLedger::of(r.sinks()).json();j["sink_ledger"]=std::move(l);}
    {J l=BirthLedger::of(r.births()).json();j["birth_ledger"]=std::move(l);}
    return j;
}
std::string MotionAssembly::sink_birth_digest(const MotionReceipt&r){
    Digest d("hall-v3-motion-sink-birth-section-v1");d.text(sink_record_api);d.u64(r.sinks().size());
    for(const auto&k:r.sinks()){d.u64(k.native_id);d.i64(k.kind);d.i64(k.boundary);d.i64(k.axis);d.f64(k.plane_m);d.f64(k.fraction);
        for(const auto*v:{&k.segment_m,&k.hit_point_m,&k.incident_proper_velocity_m_s,&k.field_point_m})for(double x:*v)d.f64(x);
        d.text(k.packet_id);d.text(k.state);d.boolean(k.collected);d.i64(k.conductor);
        for(double x:{k.flight_s,k.discarded_flight_s,k.charge_C,k.mass_kg,k.weight,k.kinetic_J,k.weighted_kinetic_J,
                      k.initial_weighted_kinetic_J,k.gyro,k.displacement_cells})d.f64(x);}
    d.u64(r.births().size());
    for(const auto&b:r.births()){d.u64(b.native_id);d.f64(b.birth_time_s);d.f64(b.flight_s);d.i64(b.origin);d.i64(b.axis);for(double x:b.start_m)d.f64(x);
        d.text(b.packet_id);d.text(b.state);d.f64(b.charge_C);d.f64(b.weight);d.f64(b.birth_weighted_kinetic_J);}
    SinkLedger::of(r.sinks()).encode(d);BirthLedger::of(r.births()).encode(d);
    return d.hex();
}
} // namespace hall_corrective_v3
