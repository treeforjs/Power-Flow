#include "axial_joint_v2.hpp"
#include <iostream>
namespace h=hall_driver;namespace n=raikou::mfem_backend;
using Rows=std::vector<n::ArenaRateEntry>;Json checks=Json::array();int forbidden_calls=0;
extern "C" void __wrap__ZN6raikou12mfem_backend6GapPic12solve_fieldsEv(void*){++forbidden_calls;throw std::runtime_error("native Poisson forbidden");}
extern "C" void __wrap__ZN6raikou12mfem_backend6GapPic7advanceEd(void*,double){++forbidden_calls;throw std::runtime_error("native D-K-D forbidden");}
void check(bool ok,const std::string&name,Json detail=Json::object()){
    checks.push_back({{"name",name},{"passed",ok},{"detail",detail}});if(!ok)throw std::runtime_error(name);}
template<class F>void reject(F f,const std::string&name){bool caught=false;try{f();}catch(const std::exception&){caught=true;}check(caught,name);}
n::GapPicOptions options(const axial::Mesh&g){n::GapPicOptions o;o.nx=o.ny=g.nx;o.width_m=o.gap_m=g.side;o.depth_m=g.depth;
    o.finite_foils=o.finite_axial_extent=true;o.short_z_m=-g.depth/2;o.feed_z_m=g.depth/2;
    o.lower_foil={g.side/2-.0025,g.side/2+.0025,g.side/2-.0027,g.side/2-.0025};
    o.upper_foil={g.side/2-.0025,g.side/2+.0025,g.side/2+.0025,g.side/2+.0027};return o;}
Rows authored(){Rows rows;std::uint64_t id=1;for(int species=0;species<3;++species){
    for(double weight:n::PicArena::positive_carriers({1e9,std::ldexp(1.,-30)})){n::GapParticle p;p.id=id++;
        p.position_m={.00421,.00433,.001};p.proper_velocity_m_s={1e5,-2e5,3e5};p.weight=weight;
        p.charge_C=species==0?-pic_e:species==1?pic_e:0;p.mass_kg=species==0?pic_me:species==1?hall::proton_mass:2.9915e-26;
        rows.push_back({"manufactured/packet/"+std::to_string(species),"manufactured/state/"+std::to_string(species),p});}}return rows;}
h::MotionBinding binding(const axial::Mesh&g,const n::GapPicOptions&o,const Rows&rows,const axial::State&s,double dt){
    h::MotionBinding b;b.phase=h::ParticlePhase::full_time_position_and_proper_velocity;b.accepted_step=s.revision;
    b.staged_arena_revision=100+2*s.revision;b.begin_s=s.time_s;b.dt_s=dt;b.end_s=s.time_s+dt;
    b.staged_arena_sha256=h::canonical_census_sha256(rows);b.mesh_sha256=h::motion_mesh_sha256(g,o);
    b.material_sha256=h::material_binding_sha256();b.species_sha256=h::text_sha256("manufactured/electron-proton-neutral");
    b.common_clock_sha256=h::text_sha256(Json{{"begin",b.begin_s},{"dt",dt},{"end",b.end_s},{"step",b.accepted_step}}.dump());
    b.kinetic_policy_sha256=h::kinetic_policy_sha256(h::declared_kinetic_policy(o));return b;}
Rows proposed_rows(const Rows&rows,const h::MotionReceipt&r){auto next=rows;std::map<std::uint64_t,n::GapParticle> map;
    for(auto p:r.endpoints())map.emplace(p.id,p);for(auto&row:next)row.native=map.at(row.native.id);return next;}
int main(int argc,char**argv){MPI_Init(&argc,&argv);int rc=0;try{
    hall::require(argc==2,"output required");std::filesystem::path out=argv[1];axial::Mesh g;auto geometry=options(g);
    h::CurrentTrace trace((std::filesystem::path(__FILE__).parent_path()/"../hall-finite-foil-em-20261029/data/shot_17133_current.csv").string());
    h::Thermal heat(g,trace.start_s());auto sigma0=heat.conductivity();axial::Fields fields(g,sigma0,trace.current(trace.start_s()));auto field=fields.start(trace.start_s());auto rows=authored();
    double max_work=0,max_balance=0,total_joule=0,max_T=300;unsigned max_iterations=0;Json history=Json::array();
    const auto original_field=h::field_state_sha256(field),original_heat=h::thermal_state_sha256(heat),original_rows=h::canonical_census_sha256(rows);
    auto first=h::prepare_axial_joint(g,fields,field,geometry,rows,binding(g,geometry,rows,field,1e-14),heat,trace);
    check(h::field_state_sha256(field)==original_field&&h::thermal_state_sha256(heat)==original_heat&&h::canonical_census_sha256(rows)==original_rows,
        "C19 joint preparation changes no field thermal or census input");
    auto stale=rows;stale.front().native.proper_velocity_m_s[0]+=1;
    check(!h::preflight_axial_joint(field,heat,stale,first),"changed staged census rejects joint preflight");
    auto altered=field;altered.E.front()+=1;check(!h::preflight_axial_joint(altered,heat,rows,first),"changed accepted field rejects joint preflight");
    h::Thermal foreign(g,trace.start_s());check(!h::preflight_axial_joint(field,foreign,rows,first),"foreign thermal owner rejects joint preflight");
    auto wrong=binding(g,geometry,rows,field,1e-14);wrong.material_sha256=std::string(64,'a');
    reject([&]{h::prepare_axial_joint(g,fields,field,geometry,rows,wrong,heat,trace);},"unselected material scenario rejected");
    std::vector<std::uint64_t>layer{1};reject([&]{h::prepare_axial_joint(g,fields,field,geometry,rows,binding(g,geometry,rows,field,1e-14),heat,trace,layer);},"joint refuses active source layer before any owner mutation");
    check(h::field_state_sha256(field)==original_field&&h::thermal_state_sha256(heat)==original_heat&&h::canonical_census_sha256(rows)==original_rows,
        "joint failures and stale preflight leave every input unchanged");
    for(unsigned i=0;i<12;++i){double dt=std::array<double,4>{1e-14,5e-15,2.5e-15,7.5e-15}[i%4];
        auto trial=h::prepare_axial_joint(g,fields,field,geometry,rows,binding(g,geometry,rows,field,dt),heat,trace);
        check(h::preflight_axial_joint(field,heat,rows,trial),"C19 joint preflight "+std::to_string(i));
        axial::Vec E(g.ne),B(g.nb);for(int e=0;e<g.ne;++e)E[e]=.5*(field.E[e]+trial.field.E[e]);
        for(int b=0;b<g.nb;++b)B[b]=.5*(field.B[b]+trial.field.B[b]);
        check(trial.motion->binding().trial_field_sha256==h::total_field_sha256(E,B),"identical accepted midpoint field seal "+std::to_string(i));
        max_iterations=std::max(max_iterations,trial.iterations);max_work=std::max(max_work,std::abs(trial.particle_work_defect_J));
        max_balance=std::max(max_balance,std::abs(trial.combined_balance_J));total_joule+=trial.field_receipt.joule_J;
        auto record=trial.field_receipt.json();record["thermal_delta_J"]=trial.heat.thermal_delta_J;record["thermal_residual_J"]=trial.heat.numerical_residual_J;
        record["particle_work_defect_J"]=trial.particle_work_defect_J;record["combined_balance_J"]=trial.combined_balance_J;
        record["motion_seal_sha256"]=trial.motion->seal_sha256();record["dt_s"]=dt;record["iterations"]=trial.iterations;
        // Test-only transient kinematics recurrence, not a native arena commit.
        const auto candidate_heat_sha=h::thermal_trial_state_sha256(trial.heat);
        auto next=proposed_rows(rows,*trial.motion);h::publish_field_heat_prevalidated(field,heat,std::move(trial));rows.swap(next);
        check(h::thermal_state_sha256(heat)==candidate_heat_sha,"proposed v2 heat state seal equals published owner state "+std::to_string(i));
        max_T=std::max(max_T,heat.max_temperature_K());record["time_s"]=field.time_s;history.push_back(record);
    }
    check(field.revision==12&&heat.revision()==12&&field.time_s==heat.time_s(),"C19 variable-step field/v2 heat clocks agree");
    check(!h::preflight_axial_joint(field,heat,rows,first),"old whole-joint token cannot be reused");
    check(max_T>300&&total_joule>0&&axial::maximum(axial::subtract(heat.conductivity(),sigma0))>0,"actual C19 Joule changes v2 temperature and next conductivity");
    check(max_work<=1e-20,"C19 full-time native motion and Maxwell work close",{{"max_defect_J",max_work},{"maximum_iterations",max_iterations},{"max_combined_balance_J",max_balance}});
    check(forbidden_calls==0,"C19 joint called neither native Poisson nor native D-K-D");
    std::ofstream(out/"joint-history.json")<<Json{{"material",h::Thermal::material_binding()},{"history",history},{"maximum_T_K",max_T},
        {"actual_SSA",false},{"accepted_native_arena_commit",false}}.dump(2)<<'\n';
    std::ofstream(out/"joint-validation.json")<<Json{{"status","passed"},{"count",checks.size()},{"checks",checks},
        {"actual_SSA",false},{"accepted_native_arena_commit",false}}.dump(2)<<'\n';
    std::cout<<checks.size()<<" C19 joint port checks passed\n";
    }catch(const std::exception&e){std::cerr<<e.what()<<'\n';rc=1;}MPI_Finalize();return rc;}
