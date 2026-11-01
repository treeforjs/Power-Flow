#include "native_gas_motion.hpp"
#include "chem_v3_native.hpp"
#include <iostream>
#include <limits>
namespace a=native_gas_motion; namespace n=a::n; namespace g=a::g;
unsigned checks=0;
void check(bool x,const char*s){++checks;if(!x)throw std::runtime_error(s);}
template<class F> void refuses(F f){bool yes=false;try{f();}catch(const std::exception&){yes=true;}check(yes,"missing refusal");}
struct Fixture {
  g::Configuration cfg;
  std::vector<n::ArenaRateEntry> before,start,after;
  std::vector<n::FullTimeStoredPacketEnergy> e0,e1,e2;
  n::FullTimeBirthReceipt births;
  n::FullTimeSinkLedger sinks;
  bool merge=false;
  n::ArenaRateEntry row(std::uint64_t id,std::uint32_t sp,double w,double x,double u) {
    n::ArenaRateEntry r;r.packet_id="packet/"+std::to_string(id);r.state=cfg.species.at(sp).state;
    r.native.id=id;r.native.weight=w;r.native.mass_kg=cfg.species[sp].mass_kg;
    r.native.charge_C=cfg.species[sp].charge*n::pic_e;
    r.native.position_m={x,.5,.5};r.native.proper_velocity_m_s={u,0,0};return r;
  }
  void sink(const n::ArenaRateEntry&r,n::FullTimeSinkKind kind,n::FullTimeSinkSurface surface,double incident,double U) {
    n::FullTimeSinkExposure x;x.packet_id=r.packet_id;x.state=r.state;x.weight=r.native.weight;
    x.charge_C=r.native.charge_C;x.mass_kg=r.native.mass_kg;x.start_position_m=r.native.position_m;x.internal_energy_J=U;
    auto&k=x.record;k.native_id=r.native.id;k.kind=kind;k.surface=surface;k.axis=0;k.plane_m=1;
    k.fraction=.5;k.segment_m={1,0,0};k.hit_point_m={1,.5,.5};k.incident_proper_velocity_m_s={incident,0,0};
    k.flight_s=1e-9;k.discarded_flight_s=1e-9;sinks.records.push_back(k);sinks.exposed.push_back(x);++sinks.carriers;
  }
  Fixture(){
    cfg.preset="manufactured-motion-records";cfg.transport=g::Transport::carrier_shares;
    cfg.grid={{0,0,0},{1,1,1},{2,1,1},.5};
    cfg.species={{"water",2.9915e-26,0,true},{"electron",9.109383713868749e-31,-1,false}};
    cfg.chemistry_blob_sha256=std::string(64,'a');cfg.chemistry_registry_sha256=std::string(64,'b');
    before={row(1,0,2,.25,1000),row(2,0,6,.25,2000),row(3,1,1,.25,3000)};
    start=before;start.push_back(row(4,0,4,.25,4000));start.push_back(row(5,1,3,.25,5000));
    after={start[0],start[3],start[4]};after[0].native.position_m[0]=.75;after[1].native.position_m[0]=.75;
    e0={{"packet/2",6,1e-19}};e1=e0;e1.push_back({"packet/4",4,8e-20});e2={e1.back()};
    births.accepted_step=0;births.entries={{"packet/4","water",4,4,false},{"packet/5","electron",5,3,false}};
    n::chem8::ESum sum;sum.add(8e-20);const auto native=n::chem8::exact_value(sum);
    births.internal_energy_J=native.value;births.internal_energy_J_digits=native.digits;
    sink(start[1],n::FullTimeSinkKind::escaped_exterior,n::FullTimeSinkSurface::axial,2500,1e-19);
    sink(start[2],n::FullTimeSinkKind::collected_at_conductor,n::FullTimeSinkSurface::upper_foil,9000,0);
  }
  a::View view()const{return {0,before,start,after,e0,e1,e2,&births,&sinks,merge};}
  a::Records prepare()const{return a::prepare(cfg,view(),100);}
};
int main(){try{
  Fixture f;const auto r=f.prepare();
  check(r.births.size()==2 && r.sinks.size()==2 && r.movers.size()==3,"source/sink/mover counts");
  check(r.births[0].native_id==4 && r.births[0].internal_J==8e-20 && r.births[1].species==1,"source order/energy");
  check(r.sinks[0].negative==1 && r.sinks[0].u[0]==2500 && r.sinks[1].u[0]==9000,"sink incident velocity");
  check(r.sink_destinations[0].kind==n::FullTimeSinkKind::escaped_exterior &&
        r.sink_destinations[1].kind==n::FullTimeSinkKind::collected_at_conductor,"escape versus conductor retained");
  check(r.movers[0].native_id==1 && r.movers[1].native_id==2 && r.movers[2].native_id==4,"mover identity/order");
  for(const auto&m:r.movers)check(m.cell_count==12 && m.from_cell==0,"postbirth source-cell census");
  check(r.movers[0].to_cell==1 && r.movers[1].to_cell==pfg::kNone && r.movers[2].to_cell==1,"actual cell/sink destinations");
  for(double x:{0.,1e-300,8e-20,1.,1e100}){
    n::chem8::ESum s;s.add(x);g::Exact e;e.add_product(x,1);
    check(a::detail::native_digits(e)==n::chem8::exact_value(s).digits,"independent native/forms internal codec");
  }
  // Real gas owner: transport exactly preserves all five deposit quantities
  // once population exits are included, independently of flow ledger K/P.
  g::Owner gas(f.cfg);gas.begin_step(0);g::FlowRecord seed{};seed.kind=pfg::kFlowDebit;seed.negative=1;
  seed.species=0;seed.grid_cell=0;seed.native_id=99;seed.weight=1;seed.u[0]=8000;seed.internal_J=3e-19;
  gas.take(0,1,{},std::span(&seed,1));
  std::array<g::Exact,5> initial;
  for(unsigned q=0;q<5;++q)initial[q]=gas.entry_total(q,true);
  gas.apply_flows(0,r.births);gas.transport(0,r.movers);gas.apply_flows(0,r.sinks);
  for(unsigned q=0;q<5;++q){auto total=gas.entry_total(q,true);const auto d=gas.global_digits(gas.layout().wall(0,q),true);total.add_digits(d.data());
    check(total==initial[q],"exact transported deposit plus exits");}
  unsigned empty_cell_deposits=0;
  for(const auto&e:gas.entries(true))if(e.cell==0)++empty_cell_deposits;
  auto born_count=gas.global_digits(gas.layout().flow(pfg::kFlowBirth,0,0),true);g::Exact b;b.add_digits(born_count.data());check(b.value()==4,"gas birth count once");
  auto sunk_count=gas.global_digits(gas.layout().flow(pfg::kFlowSink,0,0),true);g::Exact s;s.add_digits(sunk_count.data());check(s.value()==-6,"gas sink count once");
  gas.take(0,2,{},{});auto p=gas.prepare_commit(gas.candidate_sha256());check(gas.try_publish(p),"gas exact accepted publication");
  const auto encoded=gas.encode();auto restored=g::Owner::restore(encoded,f.cfg,gas.accepted_sha256());check(restored->encode()==encoded,"nonempty gas restore");
  // Input-order independence, including charged rows which create no movers.
  auto perm=f;std::reverse(perm.before.begin(),perm.before.end());std::reverse(perm.start.begin(),perm.start.end());
  std::reverse(perm.after.begin(),perm.after.end());std::reverse(perm.births.entries.begin(),perm.births.entries.end());
  const auto pr=perm.prepare();check(pr.movers.size()==r.movers.size() && pr.births.size()==r.births.size(),"permuted cardinality");
  check(std::memcmp(pr.movers.data(),r.movers.data(),r.movers.size()*sizeof(g::MoverRecord))==0,"permuted exact movers");
  check(std::memcmp(pr.births.data(),r.births.data(),r.births.size()*sizeof(g::FlowRecord))==0,"permuted exact births");
  auto bad=[&](auto modify){auto x=f;modify(x);refuses([&]{(void)x.prepare();});};
  bad([](auto&x){x.start[0].native.weight=3;});
  bad([](auto&x){x.after[0].native.weight=3;});
  bad([](auto&x){x.start[0].native.position_m[0]=.3;});
  bad([](auto&x){x.after[0].state="unknown";});
  bad([](auto&x){x.after[0].native.mass_kg*=2;});
  bad([](auto&x){x.after[0].native.charge_C=n::pic_e;});
  bad([](auto&x){x.after[0].native.id=44;});
  bad([](auto&x){x.after.push_back(x.after[0]);});
  bad([](auto&x){x.after[0].native.position_m[0]=2;});
  bad([](auto&x){x.start[0].native.proper_velocity_m_s[0]=std::numeric_limits<double>::infinity();});
  bad([](auto&x){x.births.accepted_step=1;});
  bad([](auto&x){x.births.entries[0].weight=3;});
  bad([](auto&x){x.births.entries[0].returned=true;});
  bad([](auto&x){x.births.entries[0].packet_id="wrong";});
  bad([](auto&x){x.births.entries.pop_back();});
  bad([](auto&x){x.births.wall_entries.push_back(x.births.entries[0]);});
  bad([](auto&x){x.births.internal_energy_J*=2;});
  bad([](auto&x){x.births.internal_energy_J_digits="0";});
  bad([](auto&x){x.births.internal_energy_J_digits.clear();});
  bad([](auto&x){x.e1[1].original_birth_weight=5;});
  bad([](auto&x){x.e1[1].packet_id="absent";});
  bad([](auto&x){x.e1[0].energy_J*=2;});
  bad([](auto&x){x.e2[0].energy_J*=2;});
  bad([](auto&x){x.sinks.accepted_step=1;});
  bad([](auto&x){x.sinks.exposed[0].weight=3;});
  bad([](auto&x){x.sinks.exposed[0].start_position_m[0]=.2;});
  bad([](auto&x){x.sinks.exposed[0].internal_energy_J*=2;});
  bad([](auto&x){x.sinks.records[0].incident_proper_velocity_m_s[0]*=2;});
  bad([](auto&x){x.sinks.records[0].fraction=x.sinks.exposed[0].record.fraction=2;});
  bad([](auto&x){x.sinks.records[0].surface=x.sinks.exposed[0].record.surface=n::FullTimeSinkSurface::lower_foil;});
  bad([](auto&x){x.sinks.records[0].hit_point_m[1]=x.sinks.exposed[0].record.hit_point_m[1]=std::numeric_limits<double>::infinity();});
  bad([](auto&x){x.sinks.records.pop_back();});
  bad([](auto&x){x.sinks.exposed.pop_back();});
  bad([](auto&x){x.sinks.carriers=1;});
  bad([](auto&x){std::reverse(x.sinks.exposed.begin(),x.sinks.exposed.end());std::reverse(x.sinks.records.begin(),x.sinks.records.end());});
  bad([](auto&x){x.merge=true;});
  bad([](auto&x){x.cfg.transport=g::Transport::none;});
  refuses([&]{(void)a::prepare(f.cfg,f.view(),1);});
  std::cout<<"{\"passed\":true,\"checks\":"<<checks<<",\"empty_origin_cells_with_deposit\":"<<empty_cell_deposits<<",\"scope\":\"manufactured checked mapping and actual gas owner; no native motion authority\"}\n";
  return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}
