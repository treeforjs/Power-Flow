#include "surface_path_port.hpp"
#include "native_total_field_transport.hpp"
#include <iostream>
#include <chrono>
using namespace foil;
namespace fs=std::filesystem;
Json checks=Json::array();
void check(std::string name,bool ok,Json detail=Json::object()){checks.push_back({{"name",name},{"passed",ok},{"detail",detail}});}
void baseline(){
    Mesh g(45);Vec phi(g.nn),E(g.ne),B(g.nb);
    for(int n=0;n<g.nn;++n)phi[n]=std::sin(.13*n);for(int e=0;e<g.ne;++e)E[e]=std::cos(.31*e);for(int b=0;b<g.nb;++b)B[b]=std::sin(.23*b);
    auto grad=g.gradient(phi),cg=g.curl(grad);double curlgrad=maximum(cg)/(maximum(grad)/g.dx);
    auto CE=g.curl(E),dc=g.magnetic_divergence(CE);double divcurl=maximum(dc)/(maximum(CE)/g.dx);
    check("2D curl gradient incidence",curlgrad<=2e-14,{{"relative",curlgrad}});check("2D divergence curl incidence",divcurl<=2e-14,{{"relative",divcurl}});
    check("weighted curl adjoint",std::abs(dot(CE,B)-dot(E,g.curl_transpose(B)))<=1e-12*std::max(1.,std::abs(dot(CE,B))));
    Fields solver(g,1e-14,-13468.190946430152,6012.0713604700004);auto state=solver.initial(-13468.190946430152);auto original=state;Vec J(g.ne);
    auto [next,receipt]=solver.advance(state,J,state.voltage_V);
    double field_error=maximum(subtract(next.E,state.E))/maximum(state.E),magnetic_error=maximum(subtract(next.B,state.B))/maximum(state.B);
    check("supported-start static field remains stationary",field_error<=2e-13&&magnetic_error<=2e-13,{{"electric_relative",field_error},{"magnetic_relative",magnetic_error}});
    check("background magnetic field stays in plane",maximum(Vec(state.B.begin()+g.ney+g.nex,state.B.end()))==0);
    check("terminal voltage and capacitance defined",solver.capacitance_F>0&&std::abs(receipt.voltage_error_V)<=1e-8,{{"capacitance_F",solver.capacitance_F},{"voltage_error_V",receipt.voltage_error_V}});
    auto charge0=g.field_charge(state.E);
    double maxbalance=0,maxgauss=0,maxvoltage=0;double accumulated_work=0,accumulated_boundary=0,accumulated_drive=0;
    for(int k=0;k<32;++k){
        std::fill(J.begin(),J.end(),0);J[g.ex(21,20)]=1e-5;J[g.ey(22,20)]=-2e-5;J[g.ez(22,20)]=3e-5;
        auto step=solver.advance(state,J,original.voltage_V+10*std::sin(.1*(k+1)));auto q0=g.field_charge(state.E),q1=g.field_charge(step.first.E),div=g.gradient_transpose(J);double scale=0;
        for(int p=0;p<g.nn;++p)if(g.owner[p]<0&&!g.outer_node[p]){maxgauss=std::max(maxgauss,std::abs(q1[p]-q0[p]-1e-14*div[p]));scale=std::max(scale,std::abs(1e-14*div[p]));}
        const auto&r=step.second;double transfers=std::abs(r.field_change_J)+std::abs(r.particle_work_J)+std::abs(r.outward_boundary_energy_J)+std::abs(r.background_supply_work_J)+std::abs(r.circuit_supply_work_J);
        maxbalance=std::max(maxbalance,std::abs(r.balance_J)/(1e-20+1e-10*transfers));maxvoltage=std::max(maxvoltage,std::abs(r.voltage_error_V));
        check("driven current step "+std::to_string(k)+" Gauss",maxgauss<=1e-22+2e-10*scale);
        accumulated_work+=r.particle_work_J;accumulated_boundary+=r.outward_boundary_energy_J;accumulated_drive+=r.background_supply_work_J+r.circuit_supply_work_J;state=std::move(step.first);
    }
    check("finite foil field work circuit and boundary balance",maxbalance<=1,{{"fraction_of_budget",maxbalance}});
    check("finite foil circuit voltage constraint",maxvoltage<=1e-8,{{"max_error_V",maxvoltage}});
    check("all three current components feed the field",state.E[g.ex(21,20)]!=original.E[g.ex(21,20)]&&state.E[g.ey(22,20)]!=original.E[g.ey(22,20)]&&state.E[g.ez(22,20)]!=original.E[g.ez(22,20)]);
    check("finite foil discrete divB retained",maximum(subtract(g.magnetic_divergence(state.B),g.magnetic_divergence(original.B)))<=2e-10,{{"T_per_m",maximum(g.magnetic_divergence(state.B))}});
    check("PEC tangential field zero",[&]{for(int e=0;e<g.ne;++e)if(g.fixed_E[e]&&state.E[e]!=0)return false;return true;}());
    // Zero-background transverse pulse touching an outer face: its impedance
    // power must remove exactly the field energy, without an incoming drive.
    Mesh vacuum(45,.009,.025,false);Fields radiation(vacuum,1e-14,0,0);auto wave=radiation.initial(0);
    for(int j=0;j<=vacuum.ny;++j)for(int i=0;i<=vacuum.nx;++i)wave.E[vacuum.ez(i,j)]=1e4*std::exp(-std::pow(i/3.,2))*std::sin(3.141592653589793*j/vacuum.ny);
    auto w=radiation.advance(wave,Vec(vacuum.ne),0);double work=w.second.outward_boundary_energy_J;
    check("outer impedance has positive outward Poynting work",work>0&&radiation.energy(w.first)<radiation.energy(wave));
    check("outer wave loss equals field energy debit",std::abs(w.second.field_change_J+work)<=1e-20+1e-10*work,{{"outward_J",work},{"balance_J",w.second.balance_J}});
    std::cout<<Json{{"geometry",g.json()},{"dc_current_A",solver.dc_current_A},{"initial_energy_J",solver.energy(original)},{"final_energy_J",solver.energy(state)},{"max_gauss_increment_C",maxgauss},{"work_to_authored_current_J",accumulated_work},{"outward_boundary_energy_J",accumulated_boundary},{"external_drive_work_J",accumulated_drive},{"actual_particle_feedback",false}}.dump()<<'\n';
}
void particle_controls(){
    Mesh g(45);Vec E(g.ne),B(g.nb);for(int i=0;i<g.ne;++i)E[i]=std::sin(.3*i);for(int i=0;i<g.nb;++i)B[i]=std::cos(.2*i);
    double continuity=0,adjoint=0;
    for(double sign:{-1.,1.})for(V3 D:std::vector<V3>{{.3*g.dx,-.2*g.dy,.1*g.dx},{2.7*g.dx,1.4*g.dy,-.2*g.dx},{0,0,.2*g.dx}}){V3 x{.123*g.dx,-.45*g.dy,0},y=x;for(int d=0;d<3;++d)y[d]+=D[d];Path path(g,x,D);Vec J(g.ne),q0(g.nn),q1(g.nn);path.deposit(J,sign,D,1e-14);point_deposit(g,q0,x,sign);point_deposit(g,q1,y,sign);auto div=g.gradient_transpose(J);for(int n=0;n<g.nn;++n)continuity=std::max(continuity,std::abs(q1[n]-q0[n]-1e-14*div[n]));double grid=1e-14*dot(E,J),particle=sign*hall::dot(path.gather_E(E),D);adjoint=std::max(adjoint,std::abs(grid-particle)/(1+std::abs(particle)));}
    check("2D3V path charge continuity including cell crossings",continuity<=1e-12,{{"relative",continuity}});check("2D3V current gather work adjoint",adjoint<=1e-12,{{"relative",adjoint}});
    std::vector<Species>registry{{"electron",-pic_e,pic_me},{"proton",pic_e,hall::proton_mass}};System system(g,1e-14,-13468.190946430152,6012.0713604700004,registry);auto state=system.initial(-13468.190946430152);
    for(int j=0;j<16;++j)for(int i=0;i<16;++i)for(int sp=0;sp<2;++sp)system.add(state,{0,sp,{-.002+(i+.5)*.004/16,-.002+(j+.5)*.004/16,0},{0,0,sp?0.:1e5},1e18*.004*.004*.025/256});
    auto original=state;double K0=kinetic(state.particles,registry),maxwork=0,maxgauss=0,maxpath=0;Json history=Json::array();
    for(int k=0;k<48;++k){auto r=system.advance(state,state.field.voltage_V);maxwork=std::max(maxwork,std::abs(r.particle_work_defect_J));maxgauss=std::max(maxgauss,r.gauss_increment_C);maxpath=std::max(maxpath,r.path_continuity_C);history.push_back({{"time_s",state.field.time_s},{"revision",state.field.revision},{"kinetic_J",kinetic(state.particles,registry)},{"receipt",r.json()}});}
    check("finite foil particle current changes all E components",maximum(subtract(state.field.E,original.field.E))>0&&maximum(Vec(state.last_current.begin()+g.nex+g.ney,state.last_current.end()))>0);
    check("finite foil plasma current induces B change",maximum(subtract(state.field.B,original.field.B))>0);
    check("finite foil particles and field exchange actual work",std::abs(kinetic(state.particles,registry)-K0-state.cumulative_particle_work_J)<=1e-18+1e-10*std::abs(state.cumulative_particle_work_J),{{"work_J",state.cumulative_particle_work_J},{"max_step_defect_J",maxwork}});
    check("finite foil particle Gauss and path checks admitted",maxgauss<=1e-20&&maxpath<=1e-21,{{"Gauss_C",maxgauss},{"path_C",maxpath}});
    check("axial velocity retained with in-plane background B",state.particles.front().x[2]!=original.particles.front().x[2]);
    for(bool outside:{false,true}){auto boundary=system.initial(-13468.190946430152);for(int sp=0;sp<2;++sp)system.add(boundary,{0,sp,outside?V3{g.width/2-1e-9,0,0}:V3{0,-.0025+1e-9,0},outside?V3{1e6,0,0}:V3{0,-1e6,0},1000});double K=kinetic(boundary.particles,registry);auto r=system.advance(boundary,boundary.field.voltage_V);double collected=outside?boundary.escaped_energy_J:boundary.collected_energy_J;
        check(outside?"outer particle absorption retains path energy":"foil surface absorption retains path energy",boundary.particles.empty()&&r.boundary.size()==2&&std::abs(collected-K-r.field.particle_work_J)<=1e-20+1e-10*K,{{"collected_J",collected},{"flight_s",r.boundary.empty()?0:r.boundary[0].flight_s}});
        check(outside?"outer absorbing charge receipt":"foil absorbing charge receipt",outside?boundary.escaped_charge_C==0:boundary.collected_charge_C==0);
    }
    std::cout<<Json{{"family","finite_foil_particles"},{"particles",state.particles.size()},{"steps",state.field.revision},{"dt_s",1e-14},{"duration_s",48e-14},{"authored_density_cm3",1e12},{"max_particle_work_defect_J",maxwork},{"max_Gauss_increment_C",maxgauss},{"max_path_continuity_C",maxpath},{"history",history}}.dump()<<'\n';
}
void surface_controls(){
    Mesh g(45);Vec E(g.ne);for(int i=0;i<g.ne;++i)E[i]=i<g.nex?2e6:i<g.nex+g.ney?-3e6:1e5;
    CommonInterval c{4,1,2,1e-9,1e-9+1e-14,1e-9+5e-15,1e-14,"field","arena","mesh","material"};V3 origin{0,.0025,0};std::vector<SurfacePathSI>paths;
    auto add=[&](V3 from,V3 to,double t0,double t1){V3 d{};for(int k=0;k<3;++k)d[k]=to[k]-from[k];paths.push_back({from,to,t0,t1,-1e-12,-1e-12*hall::dot(V3{2e6,-3e6,1e5},d)});};
    add({0,0,0},{0,-2e-9,0},c.begin_s,c.begin_s);add({0,-2e-9,0},{1e-9,-3e-9,2e-9},c.begin_s,c.end_s);
    auto r=import_surface_paths(g,c,c,origin,17,{17},paths,E);
    check("surface zero-duration bridge uses common nonzero dt",r.zero_duration_bridges==1&&finite(r.current_load_C_m_s)&&maximum(r.current_load_C_m_s)>0);
    check("surface path current and endpoint charge are compatible",r.continuity_C<=1e-24,{{"C",r.continuity_C}});
    check("surface field work booked exactly once",std::abs(r.work_defect_J)<=1e-20,{{"J",r.work_defect_J}});
    auto rejected=[&](auto function){try{function();return false;}catch(const std::exception&){return true;}};
    check("surface path rejects duplicate ordinary motion owner",rejected([&]{import_surface_paths(g,c,c,origin,17,{},paths,E);}));
    auto wrong=c;++wrong.accepted_revision;check("surface path rejects stale field revision",rejected([&]{import_surface_paths(g,c,wrong,origin,17,{17},paths,E);}));
    auto wrong_field=c;wrong_field.field_hash="different field";check("surface path rejects different field at same clock",rejected([&]{import_surface_paths(g,c,wrong_field,origin,17,{17},paths,E);}));
    paths[0].field_work_J+=1e-10;check("surface path rejects inconsistent source field work",rejected([&]{import_surface_paths(g,c,c,origin,17,{17},paths,E);}));
}
void native_port_controls(){
    Mesh g(45);Fields f(g,1e-14,-13468.190946430152,6012.0713604700004);auto fields=f.initial(-13468.190946430152);
    NativeInterval epoch{FieldMode::maxwell_total,0,1,fields.time_s,fields.time_s+1e-14,1e-14,"arena","field","mesh","material","species",NativePhase::full_time_work_conjugate_v1};
    std::vector<GapParticle> arena(2);for(int i=0;i<2;++i){auto&p=arena[i];p.id=17+i;p.charge_C=i?pic_e:-pic_e;p.mass_kg=i?hall::proton_mass:pic_me;p.weight=1000;p.position_m={0,0,0};p.proper_velocity_m_s={1e4,2e4,3e4};}
    auto old=arena;auto tr=prepare_native_total_field_transport(g,epoch,arena,{},fields.E,fields.B);
    check("native total-field transport returns same canonical IDs",tr.updates.size()==2&&tr.updates[0].id==17&&tr.updates[1].id==18);
    check("native transport leaves accepted arena and clock unchanged",arena[0].position_m==old[0].position_m&&arena[0].proper_velocity_m_s==old[0].proper_velocity_m_s&&tr.interval==epoch);
    check("native total-field actual current work closes",std::abs(tr.work_defect_J)<=1e-20&&tr.continuity_C<=1e-24,{{"work_defect_J",tr.work_defect_J},{"continuity_C",tr.continuity_C}});
    auto excluded=prepare_native_total_field_transport(g,epoch,arena,{17},fields.E,fields.B);
    check("canonical micro-layer owner excluded from ordinary push",excluded.updates.size()==1&&excluded.updates[0].id==18);
    auto rejected=[&](auto fn){try{fn();return false;}catch(const std::exception&){return true;}};
    auto legacy=epoch;legacy.mode=FieldMode::legacy_electrostatic_plus_external;
    check("native EM port rejects electrostatic plus external mode",rejected([&]{prepare_native_total_field_transport(g,legacy,arena,{},fields.E,fields.B);}));
    auto old_phase=epoch;old_phase.particle_phase=NativePhase::legacy_leapfrog_dkd;
    check("native EM port rejects relabeled legacy D-K-D phase",rejected([&]{prepare_native_total_field_transport(g,old_phase,arena,{},fields.E,fields.B);}));
    check("native port rejects noncanonical layer ID",rejected([&]{prepare_native_total_field_transport(g,epoch,arena,{19},fields.E,fields.B);}));
}
int main(int argc,char**argv){try{hall::require(argc==2,"output directory required");fs::path out=argv[1];fs::create_directories(out);auto start=std::chrono::steady_clock::now();baseline();particle_controls();surface_controls();native_port_controls();int failures=0;for(const auto&x:checks)if(!x.at("passed").get<bool>()){++failures;std::cout<<x.dump()<<'\n';}Json result{{"status",failures?"failed":"passed"},{"count",checks.size()},{"failures",failures},{"checks",checks},{"wall_s",std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()}};std::ofstream(out/"validation.json")<<result.dump(2)<<'\n';std::cout<<result.at("status")<<' '<<checks.size()<<" checks\n";return failures?1:0;}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}}
