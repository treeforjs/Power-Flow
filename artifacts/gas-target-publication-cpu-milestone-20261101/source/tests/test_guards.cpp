#include "pf_gas_v1/gas_owner.hpp"
#include "allocation_probe.hpp"
#include <array>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace pf_gas_v1;
int checks=0;
void check(bool v,const char *label) { if(!v) throw std::runtime_error(label); ++checks; }
template<class F> void refuses(F f,const char *label) { bool ok=false;try{f();}catch(const Refusal&){ok=true;}check(ok,label); }
Configuration config(Execution e) {
  Configuration c;c.feedback=true;c.species={{"test-neutral",1.0,0,true}};
  c.grid.low={0,0,0};c.grid.high={2,1,1};c.grid.n={2,1,1};c.grid.volume=1;c.execution=e;
  c.device_uuid="GPU-a7b9b458-167f-cd57-fec9-3dea0995b6e9";
  c.nvrtc_library="/home/treefor/projects/runtimes/pic-cuda-gather-investigation-20261029/tools/nvrtc/nvidia/cuda_nvrtc/lib/libnvrtc.so.12";
  return c;
}
FlowRecord debit(unsigned cell,double u,bool negative) {
  FlowRecord f{};f.kind=pfg::kFlowDebit;f.negative=negative;f.grid_cell=cell;f.species=0;f.native_id=cell+1;f.weight=1;f.u[0]=u;return f;
}
CellTargetsAbi cell(double s=10) { CellTargetsAbi c{};c.n_targets=1;c.n[0]=1;c.species[0]=0;c.s[0]=s;return c; }
void slot_tests() {
  auto trial=[](double n,double v,double m,std::array<double,5> val,std::array<double,3> d,double s,bool reject) {
    const auto before=d;const double bs=s;int clamped=7;
    int r=pfg::feedback_slot(n,v,m,val.data(),d.data(),&s,&clamped);
    check((r<0)==reject,"slot outcome");check(clamped==0,"never clamp");
    if(reject)check(std::memcmp(before.data(),d.data(),sizeof d)==0 && std::memcmp(&bs,&s,8)==0,"refusal atomic");
  };
  trial(1,1,1,{-2,0,0,0,0},{0,0,0},1,true);
  trial(1,1,1,{-1.5,0,0,0,0},{0,0,0},1,false);
  trial(0,1,1,{1,0,0,0,0},{0,0,0},1,true);
  trial(0,1,1,{0,0,0,0,1},{0,0,0},1,true);
  trial(0,1,1,{0,0,0,0,0},{0,0,0},1,false);
  trial(1,1,1,{1,0,0,0,0},{0,0,0},1,false);
  trial(1,1,1,{0,3,0,0,0},{0,0,0},1,true);
  trial(1,1,1,{0,0,0,0,0},{0,0,0},-1,true);
  trial(1,1,1,{0,0,0,0,0},{0,0,0},1e308,true);
  const double nan=std::numeric_limits<double>::quiet_NaN(), inf=std::numeric_limits<double>::infinity();
  trial(nan,1,1,{0,0,0,0,0},{0,0,0},1,true);
  trial(1,1,1,{nan,0,0,0,0},{0,0,0},1,true);
  trial(1,1,1,{0,0,0,0,inf},{0,0,0},1,true);
  trial(1,1,1,{0,0,0,0,0},{0,nan,0},1,true);
  trial(1e308,1e308,1,{1,0,0,0,0},{0,0,0},1,true);
  trial(1,1,1,{1,1e308,0,0,0},{0,0,0},1,true);
}
void owner_tests(Execution execution) {
  const auto cfg=config(execution);
  auto bad_cfg=cfg;bad_cfg.grid.volume=2;
  refuses([&]{Owner bad(bad_cfg);},"inconsistent cell volume");
  bad_cfg=cfg;bad_cfg.grid.n={0xffffffffu,2,1};
  refuses([&]{Owner bad(bad_cfg);},"grid cell ID overflow");
  Owner o(cfg);const auto genesis=o.accepted_sha256();
  o.begin_step(1);
  auto outside=debit(2,2,true);
  refuses([&]{o.take(1,1,{},{&outside,1});},"out-of-grid debit");
  refuses([&]{o.commit(o.candidate_sha256());},"empty early commit");
  refuses([&]{o.take(1,2,{},{});},"stage two first");
  std::array<CellTargetsAbi,1> one={cell()};std::array<unsigned,1> a={0};
  refuses([&]{o.feedback(1,2,a,one);},"feedback stage order");
  o.take(1,1,{},{});
  refuses([&]{o.take(1,1,{},{});},"duplicate stage one");
  refuses([&]{o.commit(o.candidate_sha256());},"one-stage commit");
  auto f=debit(0,2,true);
  refuses([&]{o.apply_flows(1,{&f,1});},"debit bypass");
  const auto d=o.take(1,2,{},{});
  refuses([&]{o.take(1,2,{},{});},"duplicate stage two");o.commit(d);
  check(o.accepted_step()==1,"accepted step");
  refuses([&]{o.begin_step(0);},"backward step");
  refuses([&]{o.begin_step(1);},"replayed step");
  refuses([&]{o.begin_step(3);},"skipped step");
  auto restored=Owner::restore(o.encode(),cfg,o.accepted_sha256());
  check(restored->accepted_sha256()==o.accepted_sha256(),"restore identity");
  o.begin_step(2);o.take(2,1,{},{});o.abandon();
  check(o.accepted_sha256()==restored->accepted_sha256(),"abandon");
  // Actual state deposits: first cell changes, second fails; caller receives neither.
  o.begin_step(2);std::array<FlowRecord,2> flows={debit(0,2,true),debit(1,100,false)};
  o.take(2,1,{},flows);
  std::array<CellTargetsAbi,2> cells={cell(),cell()};const auto before=cells;std::array<unsigned,2> atlas={0,1};
  const auto accepted=o.accepted_sha256(), candidate=o.candidate_sha256();
  refuses([&]{o.feedback(2,2,atlas,cells);},"nonrealizable second slot");
  check(std::memcmp(cells.data(),before.data(),sizeof cells)==0,"whole-span refusal atomic");
  check(o.accepted_sha256()==accepted && o.candidate_sha256()==candidate,"feedback no state mutation");
  o.abandon();
  o.begin_step(2);o.take(2,1,{},std::span<const FlowRecord>(flows.data(),1));
  auto before_one=one;o.feedback(2,2,a,one);check(one[0].drift[0][0]==2 && one[0].s[0]>0,"valid feedback changes drift");
  check(std::memcmp(one.data(),before_one.data(),sizeof one)!=0,"valid feedback publishes");
  one={cell()};one[0].n[0]=0;before_one=one;
  refuses([&]{o.feedback(2,2,a,one);},"zero-population residual");
  check(std::memcmp(one.data(),before_one.data(),sizeof one)==0,"zero-population atomic");
  cells={cell(),cell()};atlas={0,0};refuses([&]{o.feedback(2,2,atlas,cells);},"duplicate atlas");
  a={2};one={cell()};refuses([&]{o.feedback(2,2,a,one);},"out of grid");
  a={0};one[0].species[0]=1;refuses([&]{o.feedback(2,2,a,one);},"unknown target");
  one={cell()};one[0].n_targets=2;one[0].species[1]=0;refuses([&]{o.feedback(2,2,a,one);},"duplicate target");
  o.take(2,2,{},{});o.commit(o.candidate_sha256());
  check(o.entries().size()==1,"nonempty committed deposit");
  auto nonempty=Owner::restore(o.encode(),cfg,o.accepted_sha256());
  check(nonempty->accepted_sha256()==o.accepted_sha256(),"nonempty restore identity");
  o.begin_step(3);RecoilRecord bad{};bad.accepted_step=3;bad.stage=1;bad.ev.process=0;bad.ev.n_out=1;bad.ev.weight=1;
  refuses([&]{o.take(3,1,{&bad,1},{});},"empty process table rejects dummy process");o.abandon();
  for(Owner *continued : {&o,nonempty.get()}) {
    continued->begin_step(3);continued->take(3,1,{},{});continued->take(3,2,{},{});
    continued->commit(continued->candidate_sha256());
  }
  check(nonempty->accepted_sha256()==o.accepted_sha256(),"nonempty continued restore identity");
  check(genesis!=o.accepted_sha256(),"committed state not genesis");
  if(execution==Execution::host) {
    Owner x(cfg),y(cfg);x.begin_step(1);y.begin_step(1);
    auto f0=debit(0,2,true),f1=debit(1,2,true);x.prepare_keys({},{&f0,1},{});y.prepare_keys({},{&f1,1},{});
    auto p=x.accumulate({},{&f0,1});refuses([&]{y.merge(1,1,{&p,1});},"same-size wrong-key partial");
    auto altered=p;altered.global[0]=1;refuses([&]{x.merge(1,1,{&altered,1});},"altered partial");
    std::array<Partial,2> duplicate={p,p};
    refuses([&]{x.merge(1,1,duplicate);},"duplicated partial payload");
    refuses([&]{x.merge(1,1,{});},"missing partial payload");
    x.merge(1,1,{&p,1});refuses([&]{x.merge(1,1,{&p,1});},"duplicate partial stage");
    refuses([&]{x.merge(1,2,{&p,1});},"stale-stage partial");
    x.prepare_keys({}, {}, {});
    auto p2=x.accumulate({},{});x.merge(1,2,{&p2,1});x.commit(x.candidate_sha256());
    x.begin_step(2);refuses([&]{x.merge(2,1,{&p,1});},"stale-step partial");
    FlowRecord birth=f0;birth.kind=pfg::kFlowBirth;x.prepare_keys({},{&birth,1},{});auto pb=x.accumulate({},{&birth,1});
    refuses([&]{x.merge(2,1,{&pb,1});},"external flow as chemistry partial");x.abandon();
  }
  if(execution==Execution::host) {
    Owner direct(cfg),split(cfg);direct.begin_step(0);split.begin_step(0);
    std::array<FlowRecord,2> batch={debit(0,2,true),debit(1,3,true)};
    const auto expected=direct.take(0,1,{},batch);split.prepare_keys({},batch,{});
    std::array<Partial,3> parts={split.accumulate({},std::span<const FlowRecord>(batch).subspan(1)),
                                split.accumulate({},{}),split.accumulate({},std::span<const FlowRecord>(batch).first(1))};
    check(split.merge(0,1,parts)==expected,"partition and empty-shard exact identity");
    direct.take(0,2,{},{});split.take(0,2,{},{});
    direct.commit(direct.candidate_sha256());split.commit(split.candidate_sha256());
    check(direct.accepted_sha256()==split.accepted_sha256(),"partition committed identity");
    Owner other(cfg);other.begin_step(0);auto different=batch;different[0].weight=2;
    other.take(0,1,{},different);other.take(0,2,{},{});other.commit(other.candidate_sha256());
    direct.begin_step(1);other.begin_step(1);direct.prepare_keys({},batch,{});other.prepare_keys({},batch,{});
    auto old_history=direct.accumulate({},batch);
    refuses([&]{other.merge(1,1,{&old_history,1});},"same-key foreign accepted history");
  }
  std::cout<<"state "<<o.accepted_sha256()<<"\n";
}
#include "motion_inputs.inc"
#include "publication.inc"
int main(int argc,char**argv) {
  try {slot_tests();Execution e=Execution::host;
    if(argc>1)e=std::string(argv[1])=="device"?Execution::device:Execution::device_verified;
    owner_tests(e);motion_inputs(e);publication_tests(e);std::cout<<"PASS "<<checks<<" checks\n";return 0;
  }catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<" after "<<checks<<" checks\n";return 1;}
}
