#include "axial_system.hpp"
#include "axial_surface_path_port.hpp"
#include "axial_restart.hpp"
#include "axial_refinement.hpp"
#include "axial_split_stage.hpp"
#include <iostream>
#include <chrono>
using namespace axial;
Json checks=Json::array();
void check(const std::string&name,bool pass,Json detail=Json::object()){
    checks.push_back({{"name",name},{"passed",pass},{"detail",detail}});
    std::cout<<checks.back().dump()<<std::endl;
}
void operators(){
    Mesh g;Vec phi(g.nn),E(g.ne),B(g.nb);
    for(int n=0;n<g.nn;++n)phi[n]=std::sin(.13*n);
    for(int e=0;e<g.ne;++e)E[e]=std::cos(.31*e);
    for(int b=0;b<g.nb;++b)B[b]=std::sin(.23*b);
    auto grad=g.gradient(phi),cg=g.curl(grad),ce=g.curl(E),dc=g.magnetic_divergence(ce);
    double a=maximum(cg)/(maximum(grad)/g.dx),b=maximum(dc)/(maximum(ce)/g.dx);
    check("3D scalar curl cancellation",a<=2e-14,{{"relative",a}});
    check("3D div curl identity",b<=2e-14,{{"relative",b}});
    check("3D curl adjoint",std::abs(dot(ce,B)-dot(E,g.curl_transpose(B)))<=1e-12*std::abs(dot(ce,B)));
    check("feed port current normalization",std::abs(g.upper_region_inflow(g.feed_shape(),g.kshort+1)+1)<=2e-14);
    Vec sigma(g.cells.size());for(std::size_t c=0;c<sigma.size();++c)if(g.cells[c].body>=0)sigma[c]=1/(71.9e-8);
    auto S=g.conduction(sigma),cell=g.cell_joule_W(sigma,E);double heat=std::accumulate(cell.begin(),cell.end(),0.);
    double work=0;for(int e=0;e<g.ne;++e)work+=S[e]*E[e]*E[e];
    check("solid Joule edge and cell quadrature identical",std::abs(work-heat)<=2e-14*heat);
    Fields solver(g,sigma,6012.0713604700004);auto state=solver.start(7.4978805549199997e-10);auto old=state;
    std::cout<<Json{{"geometry",g.json()},{"initial_current_A",solver.initial.current_A},
        {"initial_terminal_V",solver.initial.terminal_V},{"dc_continuity_A",solver.initial.dc_continuity_A},
        {"interior_Ampere_residual_A_m",solver.initial.magnetic_residual_A_m},{"initial_energy_J",solver.energy(state)}}.dump()<<std::endl;
    check("DC current closes feed foils and short",solver.initial.dc_continuity_A<=2e-10*6012.0713604700004);
    check("initial scalar E has no fictitious curl",maximum(g.curl(state.E))<=2e-10*maximum(state.E)/g.dx);
    auto step=solver.prepare(state,Vec(g.ne),sigma,1e-14,6012.0713604700004);
    double er=maximum(subtract(step.first.E,state.E))/maximum(state.E),br=maximum(subtract(step.first.B,state.B))/maximum(state.B);
    check("supported DC state is stationary",er<=1e-10&&br<=1e-10,{{"relative_E",er},{"relative_B",br}});
    check("predicted DC voltage matches resistive foil pair",std::abs(solver.initial.terminal_V-2*.025*6012.0713604700004*71.9e-8/(.005*.0002))<=1e-8,
          {{"voltage_V",solver.initial.terminal_V}});
    double worst=0;Json history=Json::array();
    for(int k=0;k<16;++k){Vec J(g.ne);J[g.ex(22,22,g.kshort+2)]=1e-6;
        J[g.ey(22,22,g.kshort+2)]=-2e-6;J[g.ez(22,22,g.kshort+1)]=3e-6;
        auto advance=solver.prepare(state,J,sigma,1e-14,6012.0713604700004+(k+.5)*.016);
        const auto&r=advance.second;double sum=r.region_conduction_A+r.region_plasma_A+r.region_displacement_A+r.region_boundary_A;
        worst=std::max(worst,std::abs(sum-r.current_A));history.push_back(r.json());state=std::move(advance.first);
    }
    check("total feed includes conduction plasma displacement and boundary",worst<=2e-7,{{"max_error_A",worst}});
    check("all three authored plasma currents feed Maxwell",maximum(subtract(state.E,old.E))>0&&maximum(subtract(state.B,old.B))>0);
    std::cout<<Json{{"family","axial_current_control"},{"history",history}}.dump()<<std::endl;
}
void particles_and_heat(){
    Mesh g;Vec E(g.ne),B(g.nb);for(int e=0;e<g.ne;++e)E[e]=std::cos(.21*e);for(int b=0;b<g.nb;++b)B[b]=std::sin(.11*b);
    double continuity=0,adjoint=0;
    for(double sign:{-1.,1.})for(V3 D:std::vector<V3>{{.3*g.dx,-.2*g.dy,.1*g.dz},{2.7*g.dx,1.4*g.dy,-1.2*g.dz},{0,0,.2*g.dz}}){
        V3 x{.123*g.dx,-.45*g.dy,.2*g.dz},y=x;for(int d=0;d<3;++d)y[d]+=D[d];
        Path path(g,x,D);Vec J(g.ne),q0(g.nn),q1(g.nn);path.deposit(J,sign,D,1e-14);point_deposit(g,q0,x,sign);point_deposit(g,q1,y,sign);
        auto div=g.gradient_transpose(J);for(int n=0;n<g.nn;++n)continuity=std::max(continuity,std::abs(q1[n]-q0[n]-1e-14*div[n]));
        adjoint=std::max(adjoint,std::abs(1e-14*dot(E,J)-sign*hall::dot(path.gather_E(E),D)));
    }
    check("3D path charge continuity with axial cell crossings",continuity<=1e-12,{{"relative",continuity}});
    check("3D path current and gather work are conjugate",adjoint<=1e-12,{{"defect",adjoint}});
    CurrentTrace trace((std::filesystem::path(__FILE__).parent_path()/"data/shot_17133_current.csv").string());
    Thermal thermal(g,trace.start_s());auto sigma0=thermal.conductivity();Fields solver(g,sigma0,trace.current(trace.start_s()));auto state=solver.start(trace.start_s());
    std::vector<GapParticle>arena;std::uint64_t id=1;
    for(int k=0;k<2;++k)for(int j=0;j<8;++j)for(int i=0;i<8;++i)for(int sp=0;sp<2;++sp){
        GapParticle p;p.id=id++;p.position_m={-.002+(i+.5)*.004/8,-.002+(j+.5)*.004/8,-.00625+(k+.5)*.0125/2};
        p.proper_velocity_m_s={sp?0.:1e4,sp?0.:2e4,sp?0.:1e5};p.charge_C=sp?pic_e:-pic_e;p.mass_kg=sp?hall::proton_mass:pic_me;
        p.weight=1e18*.004*.004*.0125/128;arena.push_back(p);
    }
    auto old_field=state;auto old_arena=arena;double max_work=0,max_energy=0,total_joule=0,total_heat=0;
    Json history=Json::array();
    for(int step=0;step<32;++step){auto trial=prepare_joint(g,solver,state,arena,thermal,trace,1e-14);
        max_work=std::max(max_work,std::abs(trial.particle_work_defect_J));max_energy=std::max(max_energy,std::abs(trial.total_balance_J));
        total_joule+=trial.receipt.joule_J;total_heat+=trial.heat.thermal_delta_J;
        auto record=trial.receipt.json();record["thermal_delta_J"]=trial.heat.thermal_delta_J;record["thermal_numerical_residual_J"]=trial.heat.numerical_residual_J;
        record["particle_work_defect_J"]=trial.particle_work_defect_J;record["total_balance_J"]=trial.total_balance_J;record["time_s"]=trial.field.time_s;
        record["iterations"]=trial.iterations;record["path_continuity_C"]=trial.transport.continuity_C;
        check("common field arena thermal commit "+std::to_string(step),commit_joint(state,arena,thermal,std::move(trial)));
        record["max_T_K"]=thermal.max_temperature_K();history.push_back(record);
    }
    auto sigma1=thermal.conductivity();double sigma_change=maximum(subtract(sigma1,sigma0));
    check("actual Joule work heats the sealed enthalpy owner",thermal.max_temperature_K()>300&&total_joule>0&&thermal.energy_J()>0,
        {{"temperature_K",thermal.max_temperature_K()},{"joule_J",total_joule},{"thermal_delta_J",total_heat},{"residual_J",total_heat-total_joule}});
    check("accepted thermal state changes next-stage solid conductivity",sigma_change>0,{{"max_sigma_change_S_m",sigma_change}});
    auto changed=solver.prepare(state,Vec(g.ne),sigma1,1e-14,trace.current(state.time_s+5e-15));
    auto frozen=solver.prepare(state,Vec(g.ne),sigma0,1e-14,trace.current(state.time_s+5e-15));
    check("conductivity feedback changes actual next Maxwell solve",maximum(subtract(changed.first.E,frozen.first.E))>0,
        {{"max_E_feedback_V_m",maximum(subtract(changed.first.E,frozen.first.E))}});
    check("actual native particles change all current components",arena[0].position_m!=old_arena[0].position_m&&maximum(subtract(state.E,old_field.E))>0);
    check("actual current-driven magnetic field evolves",maximum(subtract(state.B,old_field.B))>0,{{"max_delta_B_T",maximum(subtract(state.B,old_field.B))}});
    check("all accepted component clocks agree",state.time_s==thermal.time_s()&&state.revision==thermal.revision()&&state.revision==32);
    check("particle and total-field work close",max_work<=1e-20,{{"max_work_defect_J",max_work},{"max_combined_balance_J",max_energy}});
    auto trial=prepare_joint(g,solver,state,arena,thermal,trace,1e-14);auto stale_arena=arena;stale_arena[0].proper_velocity_m_s[0]+=1;
    double beforeT=thermal.max_temperature_K();auto before=field_fingerprint(state);auto beforeArena=arena_fingerprint(arena);
    check("stale canonical arena rejects whole common commit",!commit_joint(state,stale_arena,thermal,std::move(trial))&&thermal.max_temperature_K()==beforeT&&field_fingerprint(state)==before);
    auto excessive=changed.second;for(auto&p:excessive.cell_joule_W)p*=1e20;excessive.joule_J*=1e20;
    bool rejected=false;try{thermal.prepare(excessive,1e-14,field_fingerprint(state));}catch(const std::exception&){rejected=true;}
    check("out-of-domain thermal trial retains accepted EM heat and arena",rejected&&thermal.max_temperature_K()==beforeT&&field_fingerprint(state)==before&&arena_fingerprint(arena)==beforeArena);
    std::cout<<Json{{"family","axial_native_field_enthalpy"},{"density_cm3",1e12},{"particles",arena.size()},
        {"history",history},{"phase","evaluated solid"},{"physical_water_breakdown",false}}.dump()<<std::endl;
}
void surface_and_rejections(){
    Mesh g;Vec E(g.ne);for(int e=0;e<g.ne;++e)E[e]=g.edges[e].axis==0?2e6:g.edges[e].axis==1?-3e6:1e5;
    foil::CommonInterval c{4,1,2,1e-9,1e-9+1e-14,1e-9+5e-15,1e-14,"field","arena","mesh","material"};
    V3 origin{0,.0025,0};std::vector<foil::SurfacePathSI>paths;
    auto add=[&](V3 a,V3 b,double t0,double t1){V3 d{};for(int k=0;k<3;++k)d[k]=b[k]-a[k];
        paths.push_back({a,b,t0,t1,-1e-12,-1e-12*hall::dot(V3{2e6,-3e6,1e5},d)});};
    add({0,0,0},{0,-2e-9,0},c.begin_s,c.begin_s);add({0,-2e-9,0},{1e-9,-3e-9,2e-9},c.begin_s,c.end_s);
    auto r=import_surface_paths(g,c,c,origin,17,{17},paths,E);
    auto global=paths;for(auto&p:global)for(int d=0;d<3;++d){p.from_m[d]+=origin[d];p.to_m[d]+=origin[d];}
    auto absolute=import_global_surface_paths(g,c,c,17,{17},global,E);
    check("source v2 global paths are not translated twice",absolute.current_load_C_m_s==r.current_load_C_m_s&&absolute.endpoint_charge_impulse_C==r.endpoint_charge_impulse_C);
    check("3D surface bridge deposits integrated q dx over accepted dt",r.zero_duration_bridges==1&&r.continuity_C<=1e-24&&std::abs(r.work_defect_J)<=1e-20,
        {{"continuity_C",r.continuity_C},{"work_defect_J",r.work_defect_J}});
    bool rejected=false;try{import_surface_paths(g,c,c,origin,17,{},paths,E);}catch(const std::exception&){rejected=true;}
    check("3D surface ordinary-motion ownership cannot duplicate",rejected);
    auto other_stage=c;other_stage.field_hash="different field at same clock";rejected=false;
    try{import_surface_paths(g,c,other_stage,origin,17,{17},paths,E);}catch(const std::exception&){rejected=true;}
    check("3D source path rejects different trial field at same clock",rejected);
    CurrentTrace trace((std::filesystem::path(__FILE__).parent_path()/"data/shot_17133_current.csv").string());
    rejected=false;try{trace.current(std::nextafter(trace.start_s(),-INFINITY));}catch(const std::exception&){rejected=true;}
    check("C19 prehistory is rejected without padding",rejected);
    rejected=false;try{trace.current(std::nextafter(trace.end_s(),INFINITY));}catch(const std::exception&){rejected=true;}
    check("C19 post-support extrapolation is rejected",rejected);
}
void restart_run(const std::filesystem::path&out,const std::string&mode){
    Mesh g;CurrentTrace trace((std::filesystem::path(__FILE__).parent_path()/"data/shot_17133_current.csv").string());
    Thermal thermal(g,trace.start_s());Fields fields(g,thermal.conductivity(),trace.current(trace.start_s()));auto state=fields.start(trace.start_s());
    std::vector<GapParticle>arena;Ledger ledger;ledger.last_current_A_m.resize(g.ne);
    if(mode=="resume")load_checkpoint(out/"split",g,state,arena,thermal,ledger);
    else for(int sp=0;sp<2;++sp){GapParticle p;p.id=17+sp;p.position_m={.00013,-.00021,.003};p.proper_velocity_m_s={1e5,2e5,3e5};
        p.charge_C=sp?pic_e:-pic_e;p.mass_kg=sp?hall::proton_mass:pic_me;p.weight=1e7;arena.push_back(p);ledger.next_id=19;}
    int end=mode=="split"?5:12;
    while(state.revision<std::uint64_t(end)){auto trial=prepare_joint(g,fields,state,arena,thermal,trace,1e-14);
        ledger.observe(trial);hall::require(commit_joint(state,arena,thermal,std::move(trial)),"restart pilot commit");}
    save_checkpoint(out/mode,g,state,arena,thermal,ledger);
}
ChemistryPrepared manufactured_angular_operator(const ChemistryStage&stage,std::span<const GapParticle>arena){
    ChemistryPrepared out;out.stage=stage;out.charged.assign(arena.begin(),arena.end());
    out.source_token="manufactured-angular/"+std::to_string(stage.common.accepted_revision)+"/"+std::to_string(int(stage.side));
    out.other_energy_owner="manufactured_stateless_scatterer";
    double angle=6e12*stage.local_exposure_s,c=std::cos(angle),s=std::sin(angle);auto before=charged_momentum(arena);
    for(auto&p:out.charged){auto u=p.proper_velocity_m_s;p.proper_velocity_m_s[0]=c*u[0]-s*u[1];p.proper_velocity_m_s[1]=s*u[0]+c*u[1];}
    auto after=charged_momentum(out.charged);for(int d=0;d<3;++d)out.other_owner_momentum_delta_kg_m_s[d]=before[d]-after[d];
    // A manufactured energy-conserving angular operator, not a water collision
    // rate. The stateless fixture absorbs momentum with zero authored energy.
    return out;
}
void full_time_phase_controls(){
    Mesh g;CurrentTrace trace((std::filesystem::path(__FILE__).parent_path()/"data/shot_17133_current.csv").string());
    struct Result{std::vector<GapParticle>arena;double max_collision_energy,max_charge;std::size_t steps;};
    auto initial_arena=[](){std::vector<GapParticle>arena;for(int sp=0;sp<2;++sp){GapParticle p;p.id=17+sp;
        p.position_m={.00013,-.00021,.003};p.proper_velocity_m_s={1e5,2e5,3e5};p.charge_C=sp?pic_e:-pic_e;
        p.mass_kg=sp?hall::proton_mass:pic_me;p.weight=1e7;arena.push_back(p);}return arena;};
    auto sample=[&](int refinement){Thermal thermal(g,trace.start_s());Fields fields(g,thermal.conductivity(),trace.current(trace.start_s()));auto state=fields.start(trace.start_s());auto arena=initial_arena();
        double max_energy=0,max_charge=0;int steps=0;
        for(int cycle=0;cycle<2;++cycle)for(double span:{1e-14,5e-15,7.5e-15,2.5e-15})for(int part=0;part<refinement;++part){double dt=span/refinement;
            auto trial=prepare_chemistry_em_chemistry(g,fields,state,arena,thermal,trace,dt,foil::NativePhase::full_time_work_conjugate_v1,{},manufactured_angular_operator);
            for(const auto&r:{trial.before_receipt,trial.after_receipt}){max_energy=std::max(max_energy,std::abs(r.energy_defect_J));max_charge=std::max(max_charge,r.charge_impulse_C);}
            hall::require(trial.before.stage.physical_stage_s==state.time_s&&trial.after.stage.physical_stage_s==state.time_s+dt&&
                trial.before.stage.local_exposure_s==dt/2&&trial.after.stage.local_exposure_s==dt/2,"variable chemistry phase stamps");
            // This fixture has no persistent source state or RNG. Production
            // must also preflight/commit both actual source tokens here.
            hall::require(commit_joint(state,arena,thermal,std::move(trial.candidate)),"manufactured split common commit");++steps;
        }
        hall::require(state.time_s==thermal.time_s()&&state.revision==std::uint64_t(steps),"variable split common clock");
        return Result{arena,max_energy,max_charge,std::size_t(steps)};
    };
    auto coarse=sample(1),middle=sample(2),fine=sample(4);
    auto distance=[](const auto&a,const auto&b){double result=0;for(std::size_t i=0;i<a.arena.size();++i)for(int d=0;d<3;++d)
        result=std::max(result,std::abs(a.arena[i].proper_velocity_m_s[d]-b.arena[i].proper_velocity_m_s[d]));return result;};
    double d01=distance(coarse,middle),d12=distance(middle,fine);
    check("variable-dt full-time chemistry EM chemistry converges",d12<d01,{{"coarse_middle_u_m_s",d01},{"middle_fine_u_m_s",d12},{"ratio",d01/d12},
        {"steps",{coarse.steps,middle.steps,fine.steps}},{"maximum_dt_s",1e-14},{"minimum_dt_s",6.25e-16}});
    check("frozen chemistry stages preserve nodal charge without fictitious current",std::max({coarse.max_charge,middle.max_charge,fine.max_charge})<=1e-24);
    check("collision energy remains separate from EM work",std::max({coarse.max_collision_energy,middle.max_collision_energy,fine.max_collision_energy})<=1e-20,
        {{"max_collision_energy_defect_J",std::max({coarse.max_collision_energy,middle.max_collision_energy,fine.max_collision_energy})}});
    Thermal thermal(g,trace.start_s());Fields fields(g,thermal.conductivity(),trace.current(trace.start_s()));auto state=fields.start(trace.start_s());auto arena=initial_arena();
    auto original=field_fingerprint(state),original_arena=arena_fingerprint(arena);double originalT=thermal.max_temperature_K();
    auto rejected=[&](auto fn){try{fn();return false;}catch(const std::exception&){return true;}};
    check("full-time composition rejects legacy D-K-D state",rejected([&]{prepare_chemistry_em_chemistry(g,fields,state,arena,thermal,trace,1e-14,foil::NativePhase::legacy_leapfrog_dkd,{},manufactured_angular_operator);}));
    check("reference composition refuses unhandled active layer ownership",rejected([&]{prepare_chemistry_em_chemistry(g,fields,state,arena,thermal,trace,1e-14,foil::NativePhase::full_time_work_conjugate_v1,{17},manufactured_angular_operator);}));
    ChemistryPrepare moves=[](const auto&s,auto p){auto t=manufactured_angular_operator(s,p);t.charged[0].position_m[0]+=1e-9;return t;};
    check("SSA exposure cannot advance particle positions",rejected([&]{prepare_chemistry_em_chemistry(g,fields,state,arena,thermal,trace,1e-14,foil::NativePhase::full_time_work_conjugate_v1,{},moves);}));
    ChemistryPrepare energy=[](const auto&s,auto p){auto t=manufactured_angular_operator(s,p);t.other_owner_energy_delta_J=1e-8;return t;};
    check("split collision energy mismatch rejects transaction",rejected([&]{prepare_chemistry_em_chemistry(g,fields,state,arena,thermal,trace,1e-14,foil::NativePhase::full_time_work_conjugate_v1,{},energy);}));
    auto small=arena;for(auto&p:small)p.weight=1e-5;
    ChemistryPrepare small_bad=[](const auto&s,auto p){auto t=manufactured_angular_operator(s,p);t.other_owner_energy_delta_J=1e-28;return t;};
    check("small collision energy error is not hidden by an absolute floor",rejected([&]{prepare_chemistry_em_chemistry(g,fields,state,small,thermal,trace,1e-14,foil::NativePhase::full_time_work_conjugate_v1,{},small_bad);}));
    check("rejected split leaves accepted field arena and heat untouched",field_fingerprint(state)==original&&arena_fingerprint(arena)==original_arena&&thermal.max_temperature_K()==originalT&&thermal.revision()==0);
}
int main(int argc,char**argv){try{
    hall::require(argc==2||argc==3,"output directory required");std::filesystem::path out=argv[1];
    if(argc==3){if(std::string(argv[2])=="refinement"){
        CurrentTrace trace((std::filesystem::path(__FILE__).parent_path()/"data/shot_17133_current.csv").string());
        auto result=run_refinements(trace);std::ofstream(out/"refinement.json")<<result.dump(2)<<'\n';
    }else restart_run(out,argv[2]);return 0;}
    auto begin=std::chrono::steady_clock::now();operators();particles_and_heat();surface_and_rejections();full_time_phase_controls();int failures=0;
    for(const auto&c:checks)if(!c.at("passed").get<bool>())++failures;
    Json result{{"status",failures?"failed":"passed"},{"count",checks.size()},{"failures",failures},{"checks",checks},
                {"wall_s",std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count()}};
    std::ofstream(out/"validation.json")<<result.dump(2)<<'\n';return failures?1:0;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 2;}}
