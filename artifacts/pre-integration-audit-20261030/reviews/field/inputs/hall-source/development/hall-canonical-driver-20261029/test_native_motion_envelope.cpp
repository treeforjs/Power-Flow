#include "native_joint_binding.hpp"
#include <iostream>

namespace h=hall_driver;namespace n=raikou::mfem_backend;
Json checks=Json::array();
void check(bool pass,const std::string&name){checks.push_back({{"name",name},{"passed",pass}});if(!pass)throw std::runtime_error(name);}
template<class F>void reject(F f,const std::string&name){bool caught=false;try{f();}catch(const std::exception&){caught=true;}check(caught,name);}
int main(int argc,char**argv){MPI_Init(&argc,&argv);int code=0;try{
    hall::require(argc==2,"output directory required");std::filesystem::path out=argv[1];axial::Mesh mesh;
    n::GapPicOptions geometry;geometry.nx=geometry.ny=45;geometry.width_m=geometry.gap_m=.009;geometry.depth_m=.025;
    geometry.finite_foils=geometry.finite_axial_extent=true;geometry.short_z_m=-.0125;geometry.feed_z_m=.0125;
    geometry.lower_foil={.0045-.0025,.0045+.0025,.0045-.0027,.0045-.0025};
    geometry.upper_foil={.0045-.0025,.0045+.0025,.0045+.0025,.0045+.0027};
    n::GapParticle particle;particle.id=1;particle.position_m={.00421,.00433,.001};particle.proper_velocity_m_s={1e5,-2e5,3e5};particle.weight=1e5;
    std::vector<n::ArenaRateEntry>census{{"envelope-fixture/electron","electron",particle}};
    axial::Vec E(mesh.ne,1e4),B(mesh.nb,.2);
    h::MotionBinding binding;binding.phase=h::ParticlePhase::full_time_position_and_proper_velocity;
    binding.accepted_step=7;binding.staged_arena_revision=29;binding.begin_s=7.4978805549199997e-10;
    binding.dt_s=1e-14;binding.end_s=binding.begin_s+binding.dt_s;binding.staged_arena_sha256=h::canonical_census_sha256(census);
    binding.trial_field_sha256=h::total_field_sha256(E,B);binding.mesh_sha256=h::motion_mesh_sha256(mesh,geometry);
    binding.material_sha256=h::material_binding_sha256();binding.species_sha256=h::text_sha256("envelope-fixture/electron");
    binding.common_clock_sha256=h::text_sha256("envelope-fixture/common-clock");binding.kinetic_policy_sha256=h::kinetic_policy_sha256(h::declared_kinetic_policy(geometry));
    auto motion=h::prepare_canonical_motion(mesh,geometry,census,binding,{},E,B);
    n::FullTimeMotionBinding admission;auto&common=admission.common;common.accepted_step=binding.accepted_step;common.accepted_arena_revision=22;
    common.begin_s=binding.begin_s;common.end_s=binding.end_s;common.dt_s=binding.dt_s;common.common_clock_sha256=binding.common_clock_sha256;
    common.mesh_sha256=binding.mesh_sha256;common.material_sha256=binding.material_sha256;
    common.accepted_arena_sha256=h::text_sha256("accepted-arena-before-first-exposure");
    common.accepted_field_sha256=h::text_sha256("accepted-field-fixture");
    admission.transport_trial_sha256=motion.seal_sha256();admission.field_before_sha256=common.accepted_field_sha256;
    admission.field_after_sha256=h::text_sha256("candidate-field-fixture");admission.heat_after_sha256=h::text_sha256("candidate-heat-fixture");
    admission.source_after_sha256=h::text_sha256("unchanged-source-fixture");
    auto envelope=h::native_motion_envelope(admission,motion);std::set<std::string>keys;
    for(auto it=envelope.begin();it!=envelope.end();++it)keys.insert(it.key());
    check(keys==std::set<std::string>{"schema","transport_trial_sha256","common_clock_sha256","mesh_sha256","field_before_sha256","field_after_sha256","begin_s","end_s","dt_s","owner_audit"},"exact native envelope keys without extra fields");
    check(envelope.at("schema")=="arena_total_field_motion_receipt_v1","agreed native envelope version");
    check(envelope.at("owner_audit")==h::motion_record(motion)&&!envelope.at("owner_audit").empty(),"actual provider motion record retained unchanged as owner audit");
    check(envelope.at("transport_trial_sha256")==motion.seal_sha256()&&motion.seal_sha256()!=binding.staged_arena_sha256,
        "transport identity is provider proposal seal and distinct from arena prefix");
    check(envelope.at("common_clock_sha256")==common.common_clock_sha256&&envelope.at("mesh_sha256")==common.mesh_sha256&&
        envelope.at("field_before_sha256")==admission.field_before_sha256&&envelope.at("field_after_sha256")==admission.field_after_sha256&&
        envelope.at("begin_s")==common.begin_s&&envelope.at("end_s")==common.end_s&&envelope.at("dt_s")==common.dt_s,
        "exact native common and field identity operands");
    auto wrong=admission;wrong.transport_trial_sha256=binding.staged_arena_sha256;
    reject([&]{h::native_motion_envelope(wrong,motion);},"arena-prefix substitution for provider seal rejected");
    wrong=admission;wrong.common.dt_s=std::nextafter(wrong.common.dt_s,INFINITY);
    reject([&]{h::native_motion_envelope(wrong,motion);},"changed authoritative dt rejected");
    wrong=admission;wrong.common.common_clock_sha256=h::text_sha256("another-common-clock");
    reject([&]{h::native_motion_envelope(wrong,motion);},"foreign common clock rejected");
    std::ofstream(out/"native-envelope.json")<<envelope.dump(2)<<'\n';
    std::ofstream(out/"validation.json")<<Json{{"status","passed"},{"count",checks.size()},{"checks",checks},
        {"scope","caller envelope only; actual motion receipt with illustrative external field/heat identities"},
        {"native_owner_linked",false},{"native_admission_exercised",false}}.dump(2)<<'\n';
    std::cout<<checks.size()<<" envelope checks passed\n";
    }catch(const std::exception&e){std::cerr<<e.what()<<'\n';code=1;}MPI_Finalize();return code;}
