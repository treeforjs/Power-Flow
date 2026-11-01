#include "surface_emission.hpp"
#include "inward_layer_adapter.hpp"
#include "ss304_thermal.hpp"
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <sched.h>
#include <sys/resource.h>

namespace se=surface_emission;
namespace {
constexpr double area=1.25e-4,depth=2e-4,rho_ref=8000,mass=rho_ref*area*depth;
constexpr int cells=32;constexpr double dx=25e-6,length=cells*dx,mu0=1.25663706212e-6,pi=3.14159265358979323846;
const ss304::Model material(ss304::Qualification::evaluated_solid_1500,ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation,rho_ref);
void need(bool okay,const char* text){if(!okay)throw std::runtime_error(text);}
double dot(const se::Vector& a,const se::Vector& b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
double norm(const se::Vector& a){return std::sqrt(dot(a,a));}
se::Vector add(const se::Vector& a,const se::Vector& b,double scale=1){return {a[0]+scale*b[0],a[1]+scale*b[1],a[2]+scale*b[2]};}
se::Vector mul(se::Vector a,double s){for(double& v:a)v*=s;return a;}
se::Vector cross(const se::Vector& a,const se::Vector& b){return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
double gamma(const se::Vector& u){return std::sqrt(1+dot(u,u)/(se::c_light*se::c_light));}
void plus(se::LayerExchange& a,const se::LayerExchange& b){
    a.physical_electron_number_change+=b.physical_electron_number_change;a.particle_charge_change_C+=b.particle_charge_change_C;a.electrode_charge_change_C+=b.electrode_charge_change_C;
    a.wall_heat_J+=b.wall_heat_J;a.particle_K_change_J+=b.particle_K_change_J;a.electron_reference_change_J+=b.electron_reference_change_J;a.image_energy_change_J+=b.image_energy_change_J;a.field_work_J+=b.field_work_J;
}
struct Slot {int side=0;se::LayerState layer;};
struct State {
    std::uint64_t revision=1,next_id=1;double time=0,common_dt=1e-14;
    double drive_F=1e8,drive_B_tangent=14,drive_B_normal=27;
    std::array<double,2> specific_u{material.evaluate(1490).u_model_j_kg,material.evaluate(1490).u_model_j_kg};
    std::array<long double,2> thermal_carry{};std::array<double,2> source_reservoir{};
    std::vector<Slot> electrons;
    std::array<double,cells+1> phi{},az{};
    se::LayerExchange cumulative;long double joule=0;
    double emitted=0,returned=0,max_local_residual=0,max_residual_normalizer=0;
    std::array<double,2> emitted_side{},returned_side{};
    std::uint64_t events=0,return_events=0,paths=0,microsteps=0,reentries=0,handoffs=0;
};
double physical_y(const Slot& p){return (p.side?length:0)+p.layer.position_from_surface_m[1];}
std::array<double,cells+1> solve(double boundary,const std::array<double,cells+1>& load,double coefficient) {
    std::array<double,cells+1> answer{};std::array<double,cells-1> diag{},rhs{};
    for(int i=0;i<cells-1;++i){diag[i]=2;rhs[i]=load[i+1]*dx/coefficient;}rhs.back()+=boundary;
    for(int i=1;i<cells-1;++i){const double lower=-1/diag[i-1];diag[i]+=lower;rhs[i]-=lower*rhs[i-1];}
    answer[cells]=boundary;answer[cells-1]=rhs.back()/diag.back();
    for(int i=cells-3;i>=0;--i)answer[i+1]=(rhs[i]+answer[i+2])/diag[i];
    long double residual=0,scale=0;
    for(int i=1;i<cells;++i){const double expected=load[i]*dx/coefficient,actual=2*answer[i]-answer[i-1]-answer[i+1];residual+=std::pow(actual-expected,2);scale+=expected*expected;}
    // Include the Dirichlet RHS entry, as in an eliminated original SPD solve.
    scale+=boundary*boundary;
    need(std::sqrt(residual)<=std::max(1e-30L,2e-13L*std::sqrt(scale)),"Original finite/SPD residual bound failed in manufactured Poisson stage");return answer;
}
void deposit(std::array<double,cells+1>& load,double y,double value) {
    need(y>=0 && y<=length && std::isfinite(value),"Out-of-domain manufactured charge/current");const double s=y/dx;const int i=std::min(int(s),cells-1);const double f=s-i;load[i]+=value*(1-f);load[i+1]+=value*f;
}
double field_y(const std::array<double,cells+1>& phi,double y){const int i=std::clamp(int(y/dx),0,cells-1);return -(phi[i+1]-phi[i])/dx;}
double value(const std::array<double,cells+1>& a,double y){const double s=std::clamp(y/dx,0.,double(cells));const int i=std::min(int(s),cells-1);return a[i]*(1-(s-i))+a[i+1]*(s-i);}
struct Fields {std::array<double,cells+1> phi,az,old_az;double time,dt;};
Fields fields(const State& state) {
    std::array<double,cells+1> charge{},current{};
    for(const auto& p:state.electrons){const double q=-se::qe*p.layer.particle.weight;deposit(charge,physical_y(p),q/area);deposit(current,physical_y(p),mu0*q*p.layer.particle.proper_velocity_m_s[2]/(gamma(p.layer.particle.proper_velocity_m_s)*area));}
    const double drive=state.drive_F*(1+.02*state.time/1e-12),Bdrive=state.drive_B_tangent*(1+.01*state.time/1e-12);
    return {solve(drive*length,charge,se::eps0),solve(Bdrive*length,current,1),state.az,state.time,state.common_dt};
}
se::Stage stage(const State& state,const Fields& f,int side,double y) {
    se::Stage s;const double E=field_y(f.phi,y);s.accepted_state=s.field_accepted_state=s.thermal_accepted_state=state.revision;
    s.mesh_revision=s.material_revision=1;s.patch_id=side+1;s.t=s.stage_time=s.field_time=s.thermal_time=state.time;s.dt=state.common_dt;
    s.total_E={0,E,state.revision==1?0:-(value(f.az,y)-value(f.old_az,y))/state.common_dt};
    s.total_B_T={-field_y(f.az,y),state.drive_B_normal,0};s.outward_normal={0,side?-1.:1.,0};s.T_K=material.temperature_from_u_model_k(state.specific_u[side]);s.area_m2=area;
    s.material_temperature_phase_admitted=s.locally_planar_barrier_admitted=s.local_equilibrium_electron_supply_admitted=true;return s;
}
se::Stage self_excluded_stage(const State& state,const Fields& f,const Slot& p) {
    auto s=stage(state,f,p.side,physical_y(p));
    const bool in_field_inventory=std::any_of(state.electrons.begin(),state.electrons.end(),[&](const Slot& old){return old.layer.particle.event_id==p.layer.particle.event_id;});
    if(!in_field_inventory)return s;
    std::array<double,cells+1> one{};deposit(one,physical_y(p),-se::qe/area);
    const auto self=solve(0,one,se::eps0);s.total_E[1]-=field_y(self,physical_y(p));return s;
}
struct End {se::Vector r,u;};
End macro_step(const End& a,const se::Stage& f,double h) {
    auto u=a.u;const double g0=gamma(a.u),qm=-se::qe/se::me;
    for(int i=0;i<60;++i){const auto k=mul(f.total_B_T,qm*h/(g0+gamma(u)));const auto rhs=add(add(a.u,mul(f.total_E,qm*h)),cross(a.u,k));const auto next=mul(add(add(rhs,cross(rhs,k)),mul(k,dot(rhs,k))),1/(1+dot(k,k)));
        const double err=norm(add(next,u,-1));u=next;if(err<2e-13*std::max(1.,norm(u)))return {add(a.r,mul(add(u,a.u),h/(g0+gamma(u)))),u};}
    throw std::runtime_error("Manufactured ordinary PIC midpoint did not converge");
}
void ordinary(Slot& slot,const se::Stage& f,const se::LayerPolicy& policy,double end,se::LayerExchange& ledger,std::uint64_t& paths,std::uint64_t& reentries) {
    auto& s=slot.layer;const double duration=end-s.current_time_s;if(duration<=0)return;
    End a{s.position_from_surface_m,s.particle.proper_velocity_m_s};auto b=macro_step(a,f,duration);double h=duration;
    auto global=[&](const End& p){return (slot.side?length:0)+p.r[1];};bool hit=global(b)<=0 || global(b)>=length;
    const bool reenter=s.handoff_distance_m>0 && dot(b.r,s.outward_normal)<=s.handoff_distance_m && dot(b.u,s.outward_normal)<0;
    if(reenter){double lo=0,hi=h;
        for(int j=0;j<48;++j){const double mid=(lo+hi)/2;const auto candidate=macro_step(a,f,mid);if(dot(candidate.r,s.outward_normal)>s.handoff_distance_m)lo=mid;else hi=mid;}
        h=(lo+hi)/2;b=macro_step(a,f,h);hit=false;
    }
    if(hit){const double wall=global(b)<=0?0:length;double lo=0,hi=h;
        for(int j=0;j<48;++j){const double mid=(lo+hi)/2;const auto candidate=macro_step(a,f,mid);if((wall==0 && global(candidate)>0)||(wall==length && global(candidate)<length))lo=mid;else hi=mid;}h=(lo+hi)/2;b=macro_step(a,f,h);}
    const double K0=se::kinetic_energy_J(a.u),K1=se::kinetic_energy_J(b.u),w=s.particle.weight,work=-se::qe*dot(f.total_E,add(b.r,a.r,-1));
    need(std::abs(K1-K0-work)<1e-10*std::max({K0,K1,1e-20}),"Ordinary PIC kinetic/electric-work mismatch");plus(ledger,{0,0,0,0,(K1-K0)*w,0,0,work*w});++paths;
    s.position_from_surface_m=b.r;s.particle.proper_velocity_m_s=b.u;s.current_time_s+=h;
    if(reenter) {
        const auto trial=surface_emission_integration::reenter_image_layer(s,f,policy);
        bool duplicate_rejected=false;try{surface_emission_integration::reenter_image_layer(trial.state,f,policy);}catch(const std::invalid_argument&){duplicate_rejected=true;}
        need(duplicate_rejected,"Image tail restored twice");
        s=trial.state;plus(ledger,trial.exchange);
        ++reentries;
    }
    if(hit){const auto absorbed=se::collect_at_material_surface(s.work_function_eV,K1,w,true);plus(ledger,{absorbed.physical_electron_number,absorbed.particle_charge_C,absorbed.electrode_charge_C,absorbed.wall_heat_J,-absorbed.incoming_K_J,absorbed.electron_reference_change_J,0,0});s.status=se::LayerStatus::Recollected;}
}
struct Step {State state;std::array<double,2> emitted_J{},surface_F{},temperatures{};};
void validate_inventory(const State& state) {
    long double number=0,K=0,image=0,reference=0;
    for(const auto& p:state.electrons) {
        const auto& a=p.layer;need(a.status!=se::LayerStatus::Recollected,"Recollected packet remained resident");
        number+=a.particle.weight;K+=a.particle.weight*se::kinetic_energy_J(a.particle.proper_velocity_m_s);
        reference+=a.particle.weight*se::qe*a.work_function_eV;
        if(a.status==se::LayerStatus::Active)image-=a.particle.weight*se::qe*se::qe/(16*pi*se::eps0*dot(a.position_from_surface_m,a.outward_normal));
    }
    need(std::abs(number-state.cumulative.physical_electron_number_change)<1e-12L*std::max(1.L,number),"Actual resident inventory differs from accepted source/return transaction");
    const long double residual=state.cumulative.wall_heat_J+K+image+reference-state.cumulative.field_work_J;
    const long double scale=std::abs(state.cumulative.wall_heat_J)+K+std::abs(image)+reference+std::abs(state.cumulative.field_work_J);
    need(std::abs(residual)<1e-10L*std::max(scale,1e-30L),"Actual resident phase-space/reference energy differs from accepted ledger");
    long double thermal=0;const double u0=material.evaluate(1490).u_model_j_kg;
    for(int side=0;side<2;++side)thermal+=static_cast<long double>(mass)*(state.specific_u[side]-u0)+state.thermal_carry[side];
    const long double caloric_input=state.joule+state.cumulative.wall_heat_J;
    need(std::abs(thermal-caloric_input)<1e-12L*std::max(std::abs(caloric_input),1e-30L),"Sampled wall heat/Joule/caloric carry closure failed");
}
Step prepare_step(const State& old,const se::LayerPolicy& policy,bool bad_epoch=false) {
    Step trial;trial.state=old;auto& next=trial.state;const auto field=fields(old);next.phi=field.phi;next.az=field.az;
    std::array<double,2> heat{};const double end=old.time+old.common_dt;std::uint64_t added_paths=0;
    for(int side=0;side<2;++side) {
        const auto s=stage(old,field,side,side?length:0);const auto flux=se::prepare(s,se::Config{4.3});trial.emitted_J[side]=flux.spectrum.current_A_m2;trial.surface_F[side]=flux.spectrum.extraction_field_V_m;trial.temperatures[side]=s.T_K;
        next.source_reservoir[side]+=flux.expected_physical_electrons;
        if(next.source_reservoir[side]>=1) {
            auto birth=se::begin_layer(flux,policy,next.source_reservoir[side],.6,.35,.1,next.next_id++);
            next.source_reservoir[side]=0;next.emitted+=birth.exchange.physical_electron_number_change;next.emitted_side[side]+=birth.exchange.physical_electron_number_change;++next.events;added_paths+=birth.paths.size();plus(next.cumulative,birth.exchange);heat[side]+=birth.exchange.wall_heat_J;
            next.electrons.push_back({side,birth.state});
        }
    }
    std::vector<Slot> survivors;survivors.reserve(next.electrons.size());
    for(auto slot:next.electrons) {
        auto s=self_excluded_stage(old,field,slot);se::LayerExchange increment;
        int transitions=0;
        while(slot.layer.status!=se::LayerStatus::Recollected && slot.layer.current_time_s<end) {
            need(++transitions<=8,"Image/PIC interface transitions exceeded bounded epoch budget");
            if(slot.layer.status==se::LayerStatus::Active) {
                if(bad_epoch)s.mesh_revision=2;
                const auto advance=se::advance_layer(slot.layer,s,policy,end-slot.layer.current_time_s);slot.layer=advance.state;plus(increment,advance.exchange);next.microsteps+=advance.substeps;added_paths+=advance.paths.size();
                if(slot.layer.status==se::LayerStatus::PICHandedOff)++next.handoffs;
            }
            if(slot.layer.status==se::LayerStatus::PICHandedOff)ordinary(slot,s,policy,end,increment,added_paths,next.reentries);
        }
        heat[slot.side]+=increment.wall_heat_J;plus(next.cumulative,increment);
        if(slot.layer.status==se::LayerStatus::Recollected){next.returned+=slot.layer.particle.weight;next.returned_side[slot.side]+=slot.layer.particle.weight;++next.return_events;}else survivors.push_back(slot);
    }
    next.electrons=std::move(survivors);next.paths+=added_paths;
    for(int side=0;side<2;++side) {
        const auto m=material.evaluate(trial.temperatures[side]);const double J=(old.drive_B_tangent*(1+.01*old.time/1e-12))/(mu0*depth);
        const double Qjoule=J*J*m.resistivity_ohm_m*area*depth*old.common_dt;
        const long double input=old.thermal_carry[side]+Qjoule+heat[side];
        const auto result=material.after_deposited_energy(old.specific_u[side],mass,double(input));next.specific_u[side]=result.state.u_model_j_kg;
        next.thermal_carry[side]=input-static_cast<long double>(mass)*(next.specific_u[side]-old.specific_u[side]);next.joule+=Qjoule;
    }
    next.time=end;++next.revision;need(next.electrons.size()<=128,"Bounded source-feedback inventory exceeded128");validate_inventory(next);return trial;
}
// Versioned fieldwise binary checkpoint. Two double limbs preserve the small
// caloric carry; no struct padding, global cache, RNG or live pointer is saved.
struct Codec {
    std::string data;std::size_t at=0;
    void put(std::uint64_t n){for(int i=0;i<8;++i)data.push_back(char((n>>(8*i))&255));}
    void put(double x){put(std::bit_cast<std::uint64_t>(x));}
    std::uint64_t integer(){need(at+8<=data.size(),"Truncated layer checkpoint");std::uint64_t n=0;for(int i=0;i<8;++i)n|=std::uint64_t((unsigned char)data[at++])<<(8*i);return n;}
    double number(){return std::bit_cast<double>(integer());}
    void put(const se::Vector& a){for(double x:a)put(x);}se::Vector vector(){return {number(),number(),number()};}
    void put(const se::Ledger& x){for(double v:{x.physical_electron_number,x.particle_charge_C,x.electrode_charge_C,x.wall_heat_J,x.outgoing_K_J,x.incoming_K_J,x.electron_reference_change_J,x.image_reference_change_J,x.resolved_matching_field_work_J})put(v);}
    se::Ledger ledger(){return {number(),number(),number(),number(),number(),number(),number(),number(),number()};}
    void put(const se::LayerExchange& x){for(double v:{x.physical_electron_number_change,x.particle_charge_change_C,x.electrode_charge_change_C,x.wall_heat_J,x.particle_K_change_J,x.electron_reference_change_J,x.image_energy_change_J,x.field_work_J})put(v);}
    se::LayerExchange exchange(){return {number(),number(),number(),number(),number(),number(),number(),number()};}
    void put(long double x){const double hi=double(x);put(hi);put(double(x-hi));}long double wide(){const double hi=number(),lo=number();return static_cast<long double>(hi)+lo;}
};
std::string encode(const State& s) {
    Codec c;c.put(std::uint64_t(0x454d49544c415933ULL));c.put(s.revision);c.put(s.next_id);c.put(s.time);c.put(s.common_dt);
    c.put(s.drive_F);c.put(s.drive_B_tangent);c.put(s.drive_B_normal);
    for(double x:s.specific_u)c.put(x);
    for(long double x:s.thermal_carry)c.put(x);
    for(double x:s.source_reservoir)c.put(x);
    for(double x:s.phi)c.put(x);
    for(double x:s.az)c.put(x);
    c.put(s.cumulative);c.put(s.joule);for(double x:{s.emitted,s.returned,s.max_local_residual,s.max_residual_normalizer})c.put(x);c.put(s.events);c.put(s.return_events);c.put(s.paths);c.put(s.microsteps);c.put(s.reentries);c.put(s.handoffs);c.put(std::uint64_t(s.electrons.size()));
    for(double x:s.emitted_side)c.put(x);
    for(double x:s.returned_side)c.put(x);
    for(const auto& p:s.electrons){c.put(std::uint64_t(p.side));const auto& a=p.layer;const auto& e=a.particle;c.put(e.event_id);c.put(e.accepted_state);c.put(e.patch_id);c.put(e.weight);c.put(e.normal_metal_energy_eV);c.put(e.transverse_energy_eV);c.put(e.proper_velocity_m_s);c.put(e.ledger);c.put(a.position_from_surface_m);c.put(a.outward_normal);
        for(double x:{a.current_time_s,a.return_overlap_distance_m,a.handoff_distance_m,a.work_function_eV,a.source_energy_scale_eV,a.quantum_magnetic_energy_fraction})c.put(x);
        c.put(a.mesh_revision);c.put(a.material_revision);c.put(std::uint64_t(a.status));}
    return c.data;
}
State decode(const std::string& data) {
    Codec c{data};need(c.integer()==0x454d49544c415933ULL,"Checkpoint version");State s;s.revision=c.integer();s.next_id=c.integer();s.time=c.number();s.common_dt=c.number();
    s.drive_F=c.number();s.drive_B_tangent=c.number();s.drive_B_normal=c.number();
    for(double& x:s.specific_u)x=c.number();
    for(long double& x:s.thermal_carry)x=c.wide();
    for(double& x:s.source_reservoir)x=c.number();
    for(double& x:s.phi)x=c.number();
    for(double& x:s.az)x=c.number();
    s.cumulative=c.exchange();s.joule=c.wide();s.emitted=c.number();s.returned=c.number();s.max_local_residual=c.number();s.max_residual_normalizer=c.number();s.events=c.integer();s.return_events=c.integer();s.paths=c.integer();s.microsteps=c.integer();s.reentries=c.integer();s.handoffs=c.integer();const auto count=c.integer();need(count<=128,"Checkpoint inventory limit");
    for(double& x:s.emitted_side)x=c.number();
    for(double& x:s.returned_side)x=c.number();
    for(std::uint64_t i=0;i<count;++i){Slot p;p.side=int(c.integer());auto& a=p.layer;auto& e=a.particle;e.event_id=c.integer();e.accepted_state=c.integer();e.patch_id=c.integer();e.weight=c.number();e.normal_metal_energy_eV=c.number();e.transverse_energy_eV=c.number();e.proper_velocity_m_s=c.vector();e.ledger=c.ledger();a.position_from_surface_m=c.vector();a.outward_normal=c.vector();a.current_time_s=c.number();a.return_overlap_distance_m=c.number();a.handoff_distance_m=c.number();a.work_function_eV=c.number();a.source_energy_scale_eV=c.number();a.quantum_magnetic_energy_fraction=c.number();a.mesh_revision=c.integer();a.material_revision=c.integer();a.status=se::LayerStatus(c.integer());s.electrons.push_back(p);}
    need(c.at==data.size(),"Trailing checkpoint data");return s;
}
std::vector<std::string> split(const std::string& line){std::vector<std::string> out;std::istringstream s(line);std::string v;while(std::getline(s,v,','))out.push_back(v);return out;}
}
int main(int argc,char** argv) {
    try {
        need(argc==3,"Expected output directory and legacy field probes");const std::filesystem::path out=argv[1];std::filesystem::create_directories(out);cpu_set_t mask;CPU_ZERO(&mask);
        need(sched_getaffinity(0,sizeof(mask),&mask)==0 && CPU_COUNT(&mask)==1 && CPU_ISSET(7,&mask) && getpriority(PRIO_PROCESS,0)>=15,"CPU7/nice15 required");
        se::LayerPolicy nominal;nominal.local_linear_vacuum_extent_m=dx;
        State baseline;const auto initial=fields(baseline);baseline.az=initial.az;baseline.phi=initial.phi;
        const auto original=encode(baseline);const auto proposed=prepare_step(baseline,nominal);need(encode(baseline)==original,"Source/field/T preparation mutated accepted state");
        const auto rollback_base=encode(proposed.state);
        bool rejected=false;try{prepare_step(proposed.state,nominal,true);}catch(const std::exception& e){rejected=std::string(e.what()).find("revision mismatch")!=std::string::npos;}
        need(rejected && encode(proposed.state)==rollback_base && encode(baseline)==original,"Failed field epoch did not preserve rollback base");
        State running=baseline;std::ofstream trace(out/"common-clock-feedback.csv");trace<<std::setprecision(17)<<"revision,t_s,dt_s,T_lower_K,T_upper_K,F_lower_V_m,F_upper_V_m,J_lower_A_m2,J_upper_A_m2,inventory,emitted_physical,returned_physical,wall_heat_J,field_work_J\n";
        std::string checkpoint;
        for(int step=0;step<32;++step){const auto next=prepare_step(running,nominal);running=next.state;trace<<running.revision<<','<<running.time<<','<<running.common_dt<<','<<next.temperatures[0]<<','<<next.temperatures[1]<<','<<next.surface_F[0]<<','<<next.surface_F[1]<<','<<next.emitted_J[0]<<','<<next.emitted_J[1]<<','<<running.electrons.size()<<','<<running.emitted<<','<<running.returned<<','<<running.cumulative.wall_heat_J<<','<<running.cumulative.field_work_J<<'\n';if(step==15)checkpoint=encode(running);}
        need(running.events>0 && running.emitted>0 && !running.electrons.empty(),"Evolving solved field/T pilot produced no actual source particles");
        need(material.temperature_from_u_model_k(running.specific_u[0])>1490 && running.phi!=baseline.phi && running.az!=baseline.az,"T/E/B stages did not evolve");
        State restarted=decode(checkpoint);need(encode(restarted)==checkpoint,"Versioned fieldwise layer/thermal/source checkpoint differs");
        need(std::any_of(restarted.electrons.begin(),restarted.electrons.end(),[](const Slot& p){return p.layer.status==se::LayerStatus::Active;}),"Restart checkpoint has no actual active layer");
        for(int step=16;step<32;++step)restarted=prepare_step(restarted,nominal).state;
        need(encode(restarted)==encode(running),"Checkpoint/restart across evolving field epochs is not exact");
        std::ofstream bytes(out/"layer-checkpoint.bin",std::ios::binary);bytes.write(checkpoint.data(),checkpoint.size());bytes.close();
        std::ofstream sensitivity(out/"matching-substep-sensitivity.csv");sensitivity<<std::setprecision(17)<<"eta_energy,relative_micro_step,common_dt_s,source_events,returned_events,emitted_physical,returned_physical,resident,wall_heat_J,field_work_J,energy_residual_J,microsteps,paths,mean_net_outward_current_A_m2\n";
        for(double eta:{.1,.03,.01})for(double resolution:{.02,.01,.005}) {
            auto policy=nominal;policy.image_energy_fraction=eta;policy.maximum_relative_normal_step=resolution;
            State s=baseline;for(int i=0;i<32;++i)s=prepare_step(s,policy).state;
            const double scale=std::abs(s.cumulative.wall_heat_J)+std::abs(s.cumulative.particle_K_change_J)+std::abs(s.cumulative.electron_reference_change_J)+std::abs(s.cumulative.image_energy_change_J)+std::abs(s.cumulative.field_work_J);
            need(std::abs(s.cumulative.residual_J())<1e-10*std::max(scale,1e-30),"Common-clock accumulated source/flight/return energy mismatch");
            need(std::abs(s.emitted-s.returned-s.cumulative.physical_electron_number_change)<1e-12*std::max(1.,s.emitted),"One electron inventory/count mismatch");
            sensitivity<<eta<<','<<resolution<<','<<s.common_dt<<','<<s.events<<','<<s.return_events<<','<<s.emitted<<','<<s.returned<<','<<s.electrons.size()<<','<<s.cumulative.wall_heat_J<<','<<s.cumulative.field_work_J<<','<<s.cumulative.residual_J()<<','<<s.microsteps<<','<<s.paths<<','<<se::qe*(s.emitted-s.returned)/(area*s.time)<<'\n';
        }
        // Separate lower-field/tangential-B pilot. The nominal barrier criterion
        // remains0.2, while resolved Lorentz returns can cross many common epochs.
        std::ofstream returns(out/"magnetic-return-sensitivity.csv");
        returns<<std::setprecision(17)<<"eta_energy,relative_micro_step,common_dt_s,stages,emitted_lower,returned_lower,emitted_upper,returned_upper,resident,wall_heat_J,field_work_J,energy_residual_J,mean_net_outward_current_A_m2,image_reentries,outgoing_handoffs\n";
        for(double eta:{.1,.03,.01})for(double resolution:{.02,.01,.005}) {
            auto policy=nominal;policy.image_energy_fraction=eta;policy.maximum_relative_normal_step=resolution;
            State s;s.drive_F=1e7;s.drive_B_tangent=25;s.drive_B_normal=0;const auto first=fields(s);s.phi=first.phi;s.az=first.az;
            for(int i=0;i<200;++i)s=prepare_step(s,policy).state;
            const double scale=std::abs(s.cumulative.wall_heat_J)+std::abs(s.cumulative.particle_K_change_J)+std::abs(s.cumulative.electron_reference_change_J)+std::abs(s.cumulative.image_energy_change_J)+std::abs(s.cumulative.field_work_J);
            need(s.returned_side[0]>0,"Tangential-B common-clock pilot did not resolve an extracting-side return");
            need(std::abs(s.cumulative.residual_J())<1e-10*std::max(scale,1e-30),"Magnetic return ledger mismatch");
            returns<<eta<<','<<resolution<<','<<s.common_dt<<",200,"<<s.emitted_side[0]<<','<<s.returned_side[0]<<','<<s.emitted_side[1]<<','<<s.returned_side[1]<<','<<s.electrons.size()<<','<<s.cumulative.wall_heat_J<<','<<s.cumulative.field_work_J<<','<<s.cumulative.residual_J()<<','<<se::qe*(s.emitted-s.returned)/(area*s.time)<<','<<s.reentries<<','<<s.handoffs<<'\n';
        }
        // Nominal0.2 guard assessment on requested legacy magnitudes, never widening
        // or deleting failed samples. These remain finite-offset probes, not true
        // one-sided coupled boundary fields.
        std::ifstream probes(argv[2]);std::string line;std::getline(probes,line);std::vector<std::vector<std::string>> rows;while(std::getline(probes,line))rows.push_back(split(line));
        std::vector<std::size_t> selected;
        for(const auto& name:std::vector<std::string>{"upper_back","upper_gap","upper_left","upper_right","lower_back","lower_gap","lower_left","lower_right"}) {
            std::size_t lo=rows.size(),hi=rows.size();double low=1e300,high=0;
            for(std::size_t i=0;i<rows.size();++i)if(rows[i][0]==name){double F=0;for(int j=0;j<3;++j)F-=std::stod(rows[i][4+j])*std::stod(rows[i][9+j]);if(std::abs(F)<low){low=std::abs(F);lo=i;}if(std::abs(F)>high){high=std::abs(F);hi=i;}}
            need(lo<rows.size() && hi<rows.size(),"Missing requested probe face");selected.push_back(lo);selected.push_back(hi);
        }
        std::ofstream admission(out/"nominal-barrier-admission.csv");admission<<"face,legacy_row,T_K,signed_F_V_m,B_T,normal_quantile,admitted,reason\n";int admissions=0,failures=0;
        for(auto index:selected)for(double T:{300.,773.15,1490.}) {
            auto s=stage(baseline,initial,0,0);for(int j=0;j<3;++j){s.outward_normal[j]=std::stod(rows[index][4+j]);s.total_E[j]=std::stod(rows[index][9+j]);s.total_B_T[j]=std::stod(rows[index][12+j]);}s.T_K=T;
            const auto f=se::prepare(s,se::Config{4.3});
            for(double q:{.1,.5,.9}) {bool okay=true;std::string reason="passed";try{se::begin_layer(f,nominal,1,q,.5,.1,1);}catch(const std::exception& e){okay=false;reason=e.what();}
                admission<<std::setprecision(17)<<rows[index][0]<<','<<index<<','<<T<<','<<f.spectrum.extraction_field_V_m<<','<<norm(s.total_B_T)<<','<<q<<','<<okay<<','<<std::quoted(reason)<<'\n';if(okay)++admissions;else++failures;}
        }
        std::ofstream report(out/"integration-qualification.json");report<<std::setprecision(17)<<"{\"passed\":true,\"stable_API_modified\":false,\"physical_foil_model_admitted\":false,\"GPU_executed\":false,\"common_dt_s\":"<<baseline.common_dt<<",\"common_stages\":32,\"source_events\":"<<running.events<<",\"emitted_physical_electrons\":"<<running.emitted<<",\"returned_physical_electrons\":"<<running.returned<<",\"final_inventory_packets\":"<<running.electrons.size()<<",\"exact_field_epoch_restart\":true,\"prepare_and_failure_rollback\":true,\"matching_substep_variants\":9,\"nominal_magnetic_guard\":0.2,\"representative_probe_phase_admissions\":"<<admissions<<",\"representative_probe_phase_guard_failures\":"<<failures<<",\"probe_counts_are_not_physical_probabilities\":true,\"scope\":\"Manufactured1D Poisson/vector-potential field stages plus explicitly selected SS304 caloric scenario and boundary component. Not a native foil/Maxwell validation or measured surface model.\"}\n";
        std::cout<<"{\"passed\":true,\"emission_events\":"<<running.events<<",\"nominal_guard_failures\":"<<failures<<"}\n";
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
