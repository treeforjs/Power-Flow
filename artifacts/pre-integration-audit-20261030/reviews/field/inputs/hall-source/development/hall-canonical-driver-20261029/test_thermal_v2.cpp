#include "axial_thermal_v2.hpp"
#include "raikou/mfem_pic_arena.h"
#include <iostream>
#include <type_traits>

namespace hd=hall_driver;
using axial::Mesh;using axial::Fields;using axial::State;using axial::Vec;
Json checks=Json::array();
void check(bool ok,const std::string&name,Json detail=Json::object()){
    checks.push_back({{"name",name},{"passed",ok},{"detail",detail}});
    if(!ok)throw std::runtime_error(name);
}
template<class F>void rejects(F f,const std::string&name){bool caught=false;try{f();}catch(const std::exception&){caught=true;}check(caught,name);}
Json states(const hd::Thermal&t){Json out=Json::array();for(const auto&s:t.states())out.push_back({
    {"T",s.T},{"time",s.time},{"revision",s.revision},{"cumulative",s.cumulative},{"watermarks",s.watermarks}});return out;}
Json field_json(const State&s){return {{"schema","v2-thermal-field-checkpoint-fixture-v1"},
    {"E",s.E},{"B",s.B},{"E_low",s.E_low},{"B_low",s.B_low},{"time_s",s.time_s},{"revision",s.revision},
    {"material",hd::Thermal::material_binding()},{"canonical_arena_checkpoint",false}};}
void save(const std::filesystem::path&path,const State&s,const hd::Thermal&t){
    std::filesystem::create_directories(path);std::ofstream(path/"field.json")<<field_json(s).dump()<<'\n';t.save(path/"thermal");
}
State load(const std::filesystem::path&path,const Mesh&g,hd::Thermal&t){
    Json j;std::ifstream(path/"field.json")>>j;
    hall::require(j.at("schema")=="v2-thermal-field-checkpoint-fixture-v1"&&j.at("material")==hd::Thermal::material_binding(),"v2 fixture schema/material");
    State s{j.at("E").get<Vec>(),j.at("B").get<Vec>(),j.at("time_s"),j.at("revision"),j.at("E_low").get<Vec>(),j.at("B_low").get<Vec>()};
    hall::require(s.E.size()==std::size_t(g.ne)&&s.E_low.size()==s.E.size()&&s.B.size()==std::size_t(g.nb)&&s.B_low.size()==s.B.size(),"v2 fixture dimensions");
    t.load(path/"thermal");hall::require(t.time_s()==s.time_s&&t.revision()==s.revision,"v2 fixture clock");return s;
}
double step_dt(unsigned i){return std::array<double,4>{1e-14,2.5e-15,7.5e-15,5e-15}.at(i%4);}
struct Step {State field;axial::Receipt receipt;hd::ThermalTrial heat;};
Step prepare(const Mesh&g,const Fields&f,const State&s,hd::Thermal&t,const hd::CurrentTrace&trace,double dt){
    hall::require(s.time_s==t.time_s()&&s.revision==t.revision(),"v2 fixture common accepted clock");
    auto proposal=f.prepare(s,Vec(g.ne),t.conductivity(),dt,trace.current(s.time_s+dt/2));
    auto heat=t.prepare(proposal.second,dt,field_json(proposal.first).dump());
    double balance=proposal.second.balance_J+heat.thermal_delta_J-proposal.second.joule_J;
    double scale=std::abs(proposal.second.field_delta_J)+std::abs(heat.thermal_delta_J)+
        std::abs(proposal.second.outer_loss_J)+std::abs(proposal.second.outer_supply_J)+std::abs(proposal.second.circuit_work_J);
    hall::require(std::abs(balance)<=2e-20+1e-10*scale+heat.roundoff_bound_J,"v2 field thermal ledger");
    return {std::move(proposal.first),std::move(proposal.second),std::move(heat)};
}
void commit(State&s,hd::Thermal&t,Step&&p){hall::require(t.preflight(p.heat),"v2 fixture heat preflight");
    if(!t.commit(std::move(p.heat)))std::terminate();using std::swap;swap(s,p.field);}
void benchmark(const std::filesystem::path&out,const hd::CurrentTrace&trace){
    static_assert(!std::is_same_v<electrode_heat_v2::Count,ArenaCount>);
    check(true,"v2 thermal and canonical v6.2 Count namespaces coexist without conversion");
    Mesh g;hd::Thermal t(g,trace.start_s());const auto sigma0=t.conductivity();Fields fields(g,sigma0,trace.current(trace.start_s()));auto s=fields.start(trace.start_s());
    check(t.material_binding().at("minimum_temperature_K")==298.2&&t.material_binding().at("maximum_temperature_K")==1500.,"explicit v2 evaluated-solid domain");
    const auto original=states(t),original_field=field_json(s);
    auto p=prepare(g,fields,s,t,trace,step_dt(0));
    check(states(t)==original&&field_json(s)==original_field,"field and all thermal slices unchanged by prepare");
    check(std::abs(p.heat.admitted_joule_J-p.receipt.joule_J)<=1e-22+2e-12*std::abs(p.receipt.joule_J),"same actual midpoint Joule operand admitted once");
    auto peer=prepare(g,fields,s,t,trace,step_dt(0));
    hd::Thermal foreign(g,trace.start_s());check(!foreign.preflight(p.heat),"foreign heat aggregate rejects token");
    auto bad=p.receipt;bad.cell_joule_W.front()=1;
    rejects([&]{t.prepare(bad,step_dt(0),"vacuum-power-fixture");},"power in vacuum rejected");
    bad=p.receipt;bad.cell_joule_W.front()=-1e-100;
    rejects([&]{t.prepare(bad,step_dt(0),"negative-power-fixture");},"negative power has no admission floor");
    bad=p.receipt;bad.cell_joule_W.pop_back();
    rejects([&]{t.prepare(bad,step_dt(0),"wrong-dimensions-fixture");},"wrong field cell dimension rejected");
    bad=p.receipt;bad.joule_J*=2;
    rejects([&]{t.prepare(bad,step_dt(0),"different-joule-fixture");},"changed Joule total rejected");
    bad=p.receipt;for(auto&v:bad.cell_joule_W)v*=1e20;bad.joule_J*=1e20;
    rejects([&]{t.prepare(bad,step_dt(0),"outside-domain-fixture");},"out-of-domain trial rejected");
    rejects([&]{t.prepare(p.receipt,2e-14,"dt-fixture");},"10 fs guard retained");
    check(states(t)==original&&field_json(s)==original_field,"failed trials preserve all accepted heat and fields");
    double joule=p.receipt.joule_J,thermal=p.heat.thermal_delta_J,residual=p.heat.numerical_residual_J;
    commit(s,t,std::move(p));
    check(!t.preflight(peer.heat)&&!t.commit(std::move(peer.heat)),"stale thermal proposal cannot publish");
    check(!t.preflight(p.heat)&&!t.commit(std::move(p.heat)),"consumed thermal proposal cannot publish twice");
    Json history=Json::array();
    for(unsigned i=1;i<16;++i){auto q=prepare(g,fields,s,t,trace,step_dt(i));joule+=q.receipt.joule_J;thermal+=q.heat.thermal_delta_J;residual+=q.heat.numerical_residual_J;
        history.push_back({{"dt_s",step_dt(i)},{"Joule_J",q.receipt.joule_J},{"thermal_delta_J",q.heat.thermal_delta_J},{"thermal_residual_J",q.heat.numerical_residual_J}});
        commit(s,t,std::move(q));}
    check(t.revision()==16&&s.revision==16&&t.time_s()==s.time_s,"all eight heat owners and fields share variable-step endpoint clock");
    check(t.max_temperature_K()>300&&joule>0&&t.energy_J()>0,"v2 enthalpy owner heats from actual Maxwell Joule",
        {{"joule_J",joule},{"thermal_delta_J",thermal},{"thermal_residual_J",residual},{"max_T_K",t.max_temperature_K()}});
    auto sigma=t.conductivity();check(axial::maximum(axial::subtract(sigma,sigma0))>0,"v2 accepted temperature changes conductivity");
    auto warm=fields.prepare(s,Vec(g.ne),sigma,step_dt(0),trace.current(s.time_s+step_dt(0)/2));
    auto cold=fields.prepare(s,Vec(g.ne),sigma0,step_dt(0),trace.current(s.time_s+step_dt(0)/2));
    check(axial::maximum(axial::subtract(warm.first.E,cold.first.E))>0,"v2 conductivity affects next actual Maxwell solution");
    save(out/"saved",s,t);hd::Thermal restored(g,trace.start_s());auto loaded=load(out/"saved",g,restored);
    check(states(restored)==states(t)&&field_json(loaded)==field_json(s),"same-process exact v2 aggregate restore");
    auto damaged=out/"damaged";std::filesystem::copy(out/"saved"/"thermal",damaged,std::filesystem::copy_options::recursive);
    {std::ofstream f(damaged/"7.bin",std::ios::binary|std::ios::trunc);f<<"corrupt final slice";}
    hd::Thermal rejectee(g,trace.start_s());auto base=states(rejectee);
    rejects([&]{rejectee.load(damaged);},"corrupt final heat slice rejected");
    check(states(rejectee)==base,"corrupt final slice leaves earlier accepted slices untouched");
    rejectee.load(out/"saved"/"thermal");check(states(rejectee)==states(t),"valid restore after rejected aggregate remains exact");
    std::ofstream(out/"history.json")<<Json{{"material",hd::Thermal::material_binding()},{"history",history}}.dump(2)<<'\n';
}
int main(int argc,char**argv){MPI_Init(&argc,&argv);int rc=0;try{
    hall::require(argc==3,"output and mode required");std::filesystem::path out=argv[1];std::string mode=argv[2];
    hd::CurrentTrace trace((std::filesystem::path(__FILE__).parent_path()/"../hall-finite-foil-em-20261029/data/shot_17133_current.csv").string());
    if(mode=="benchmark"){benchmark(out,trace);std::ofstream(out/"validation.json")<<Json{{"status","passed"},{"count",checks.size()},{"checks",checks}}.dump(2)<<'\n';}
    else{Mesh g;hd::Thermal t(g,trace.start_s());Fields f(g,t.conductivity(),trace.current(trace.start_s()));auto s=f.start(trace.start_s());
        unsigned begin=0,end=16;if(mode=="split")end=7;else if(mode=="resume"){s=load(out/"split",g,t);begin=7;}else hall::require(mode=="full","unknown mode");
        for(unsigned i=begin;i<end;++i){auto p=prepare(g,f,s,t,trace,step_dt(i));commit(s,t,std::move(p));}save(out/mode,s,t);}
    std::cout<<"passed "<<mode<<'\n';
    }catch(const std::exception&e){std::cerr<<e.what()<<'\n';rc=1;}MPI_Finalize();return rc;}
