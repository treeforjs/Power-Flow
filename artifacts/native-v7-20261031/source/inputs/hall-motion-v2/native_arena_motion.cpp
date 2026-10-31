#include "motion_receipt.hpp"
#include "../hall-finite-foil-em-20261029/axial_particles.hpp"
#include <bit>
#include <iomanip>
#include <openssl/sha.h>
#include <sstream>

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
        std::ostringstream out;for(auto b:d)out<<std::hex<<std::setfill('0')<<std::setw(2)<<unsigned(b);return out.str();}
};
bool digest(const std::string&s){return s.size()==64&&s.find_first_not_of("0123456789abcdef")==std::string::npos;}
std::vector<const native::ArenaRateEntry*> ordered(std::span<const native::ArenaRateEntry> census){
    hall::require(census.size()<=100000,"motion census exceeds sealed native carrier cap");
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
native::ArenaJson binding_json(const MotionBinding&b){return {{"phase","full_time_position_and_proper_velocity"},
    {"accepted_step",b.accepted_step},{"staged_arena_revision",b.staged_arena_revision},{"trial_stage",b.trial_stage},
    {"begin_s",b.begin_s},{"end_s",b.end_s},{"dt_s",b.dt_s},{"staged_arena_sha256",b.staged_arena_sha256},
    {"trial_field_sha256",b.trial_field_sha256},{"mesh_sha256",b.mesh_sha256},{"material_sha256",b.material_sha256},
    {"species_sha256",b.species_sha256},{"common_clock_sha256",b.common_clock_sha256},{"kinetic_policy_sha256",b.kinetic_policy_sha256}};}
std::string array_hash(const std::vector<double>&v){Encoding e;e.text("motion-array-binary64-v1");e.vector(v);return e.finish();}
native::ArenaJson record_without_seal(const MotionReceipt&r){
    native::ArenaJson endpoints=native::ArenaJson::array(),paths=native::ArenaJson::array();
    for(const auto&p:r.endpoints())endpoints.push_back({{"native_id",p.id},{"position_m",p.position_m},
        {"proper_velocity_m_s",p.proper_velocity_m_s},{"charge_C",p.charge_C},{"mass_kg",p.mass_kg},{"weight",p.weight}});
    for(const auto&p:r.paths())paths.push_back({{"native_id",p.native_id},{"field_start_m",p.field_start_m},
        {"field_displacement_m",p.field_displacement_m},{"flight_s",p.flight_s},{"motion_owner","ordinary_total_field"}});
    const auto&b=r.budget();
    const auto&k=r.kinetic_policy();
    return {{"schema",motion_api},{"binding",binding_json(r.binding())},{"census_sha256",r.census_sha256()},
        {"kinetic_policy",{{"schema",kinetic_boundary_api},{"width_m",k.width_m},{"gap_m",k.gap_m},
            {"short_z_m",k.short_z_m},{"feed_z_m",k.feed_z_m},
            {"lower_foil",{k.lower.left_m,k.lower.right_m,k.lower.bottom_m,k.lower.top_m}},
            {"upper_foil",{k.upper.left_m,k.upper.right_m,k.upper.bottom_m,k.upper.top_m}},
            {"short_bridge","lower foil x span and inner y gap only"},{"xy_policy","absorbing escape"},
            {"feed_policy","axial escape"},{"tie_order","x y lower upper axial; later exact tie wins"}}},
        {"native_to_field_m",r.frame().native_to_field_m},{"endpoints",endpoints},{"paths",paths},
        {"current_A_m_sha256",array_hash(r.current_A_m())},{"initial_charge_C_sha256",array_hash(r.initial_charge_C())},
        {"endpoint_charge_C_sha256",array_hash(r.endpoint_charge_C())},
        {"budget",{{"kinetic_change_J",b.kinetic_change_J},{"field_work_J",b.field_work_J},
            {"work_defect_J",b.work_defect_J},{"continuity_C",b.continuity_C},{"plasma_dt",b.plasma_dt},
            {"maximum_gyro",b.maximum_gyro},{"maximum_displacement_cells",b.maximum_displacement_cells}}},
        {"coverage","all_native_carriers_no_sources_no_boundaries"},{"neutral_motion_owner","same_full_time_q0_drift"},
        {"native_Poisson_called",false},{"native_DKD_called",false}};
}
} // namespace
std::string canonical_census_sha256(std::span<const native::ArenaRateEntry> census){
    auto rows=ordered(census);Encoding e;e.text("canonical-motion-census-v1");e.integer(rows.size());
    for(const auto*row:rows){e.text(row->packet_id);e.text(row->state);e.particle(row->native);}return e.finish();
}
std::string motion_mesh_sha256(const axial::Mesh&g,const native::GapPicOptions&o){
    geometry_guard(g,o);Encoding e;e.text("canonical-motion-geometry-v2");e.text(g.json().dump());
    e.integer(o.nx);e.integer(o.ny);for(double x:{o.width_m,o.gap_m,o.depth_m,o.short_z_m,o.feed_z_m})e.scalar(x);
    for(const auto&f:{o.lower_foil,o.upper_foil})for(double x:{f.left_m,f.right_m,f.bottom_m,f.top_m})e.scalar(x);
    // Native prescribed fields and potentials are deliberately absent: they are
    // not sampled by this total-field operator. Native arena owns their options.
    return e.finish();
}
std::string total_field_sha256(const std::vector<double>&E,const std::vector<double>&B){
    Encoding e;e.text("total-midpoint-E-B-binary64-v1");e.vector(E);e.vector(B);return e.finish();
}
native::ArenaJson motion_record(const MotionReceipt&r){auto j=record_without_seal(r);j["seal_sha256"]=r.seal_sha256();return j;}
MotionReceipt prepare_canonical_motion(const axial::Mesh&g,const native::GapPicOptions&o,
    std::span<const native::ArenaRateEntry>census,const MotionBinding&binding,
    std::span<const std::uint64_t> active_layer_ids,const std::vector<double>&E,const std::vector<double>&B){
    hall::require(binding.phase==ParticlePhase::full_time_position_and_proper_velocity,"motion requires owner-certified full-time phase");
    hall::require(std::isfinite(binding.begin_s)&&binding.begin_s>=0&&std::isfinite(binding.end_s)&&
        std::isfinite(binding.dt_s)&&binding.dt_s>0&&binding.dt_s<=1e-14&&binding.end_s>binding.begin_s&&
        binding.end_s==binding.begin_s+binding.dt_s,"motion common interval/10 fs guard");
    for(const auto*s:{&binding.staged_arena_sha256,&binding.trial_field_sha256,&binding.mesh_sha256,
        &binding.material_sha256,&binding.species_sha256,&binding.common_clock_sha256,&binding.kinetic_policy_sha256})hall::require(digest(*s),"motion requires sealed stage identities");
    hall::require(active_layer_ids.empty(),"motion requires source owner for every active layer; no implicit exclusion/freeze");
    hall::require(E.size()==std::size_t(g.ne)&&B.size()==std::size_t(g.nb)&&axial::finite(E)&&axial::finite(B),"motion total-field dimensions/finiteness");
    hall::require(binding.trial_field_sha256==total_field_sha256(E,B)&&binding.mesh_sha256==motion_mesh_sha256(g,o),"motion changed trial field or mesh binding");
    MotionReceipt r;r.binding_=binding;r.frame_.native_to_field_m={-o.width_m/2,-o.gap_m/2,0};
    r.kinetic_policy_=declared_kinetic_policy(o);
    hall::require(binding.kinetic_policy_sha256==kinetic_policy_sha256(r.kinetic_policy_),"motion kinetic policy changed independently of field mesh");
    r.census_sha256_=canonical_census_sha256(census);const auto rows=ordered(census);
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
    Encoding seal;seal.text("motion-receipt-record-v1");seal.text(record_without_seal(r).dump());r.seal_sha256_=seal.finish();
    return r;
}
} // namespace hall_driver
