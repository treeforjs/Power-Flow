CensusRow row(std::uint64_t id, double x=0.25, double weight=1) {
  CensusRow r{};r.native_id=id;r.species=0;r.position_m[0]=x;r.position_m[1]=0.5;r.position_m[2]=0.5;r.weight=weight;
  return r;
}
void motion_inputs(Execution execution) {
  auto cfg=config(execution);cfg.transport=Transport::carrier_shares;
  auto a=row(1),b=row(2);std::array<CensusRow,2> before={a,b},after={a,b};
  after[0].position_m[0]=1.25;
  auto valid=movers(cfg,before,after,{});
  check(valid.size()==1 && valid[0].native_id==1 && valid[0].from_cell==0 &&
        valid[0].to_cell==1 && valid[0].weight==1 && valid[0].cell_count==2,"valid mover");
  std::reverse(before.begin(),before.end());std::reverse(after.begin(),after.end());
  auto permuted=movers(cfg,before,after,{});
  check(std::memcmp(valid.data(),permuted.data(),sizeof(MoverRecord))==0,"motion row order independent");
  before={a,b};after={a,b};
  const std::uint64_t one=1,two=2,unknown=3;
  check(movers(cfg,before,{&b,1},{&one,1}).front().to_cell==pfg::kNone,"valid sink");
  check(movers(cfg,before,after,{}).empty(),"stationary census");
  auto duplicate=before;duplicate[1]=duplicate[0];
  refuses([&]{census(cfg,duplicate);},"duplicate before census");
  refuses([&]{movers(cfg,before,duplicate,{});},"duplicate after census");
  refuses([&]{movers(cfg,before,after,{&one,1});},"after plus sink overlap");
  refuses([&]{movers(cfg,before,after,{&unknown,1});},"unknown sink");
  std::array<std::uint64_t,2> dup_sink={two,two};
  refuses([&]{movers(cfg,before,{&a,1},dup_sink);},"duplicate sink");
  refuses([&]{movers(cfg,before,{&a,1},{});},"missing survivor");
  auto changed=after;changed[1].native_id=3;
  refuses([&]{movers(cfg,before,changed,{});},"unknown survivor");
  changed=after;changed[0].weight=2;
  refuses([&]{movers(cfg,before,changed,{});},"motion changes weight");
  auto cfg2=cfg;cfg2.species.push_back({"charged",1,1,false});changed=after;changed[0].species=1;
  refuses([&]{movers(cfg2,before,changed,{});},"motion changes species");
  changed=after;changed[0].position_m[0]=2.25;
  refuses([&]{movers(cfg,before,changed,{});},"uncollected survivor outside grid");
  for(double bad : {0.0,-1.0,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}) {
    auto r=a;r.weight=bad;refuses([&]{census(cfg,{&r,1});},"invalid census weight");
  }
  auto r=a;r.proper_velocity_m_s[0]=std::numeric_limits<double>::quiet_NaN();
  refuses([&]{census(cfg,{&r,1});},"nonfinite census velocity");
  r=a;r.position_m[0]=std::numeric_limits<double>::quiet_NaN();
  refuses([&]{census(cfg,{&r,1});},"nonfinite census position");
  auto bad_cfg=cfg;bad_cfg.grid.n[0]=0;refuses([&]{census(bad_cfg,before);},"utility validates grid");
  bad_cfg=cfg;bad_cfg.grid.volume=2;refuses([&]{census(bad_cfg,before);},"utility validates volume");
  bad_cfg=cfg;bad_cfg.species[0].mass_kg=0;refuses([&]{census(bad_cfg,before);},"utility validates mass");
  Exact exact;exact.add_product(2,3);const auto digits=exact.normalised();
  refuses([&]{exact.add_product(std::numeric_limits<double>::infinity(),1);},"exact rejects nonfinite term");
  check(exact.normalised()==digits,"failed exact add atomic");
  refuses([&]{exact.add_triple(0x1p-1000,0x1p-1000,0x1p-1000);},"exact rejects below-grid term");
  check(exact.normalised()==digits,"failed triple atomic");
  Exact huge;huge.add_product(1e308,2);refuses([&]{huge.value();},"exact rounded overflow refuses");
  Owner o(cfg);o.begin_step(1);auto debit_record=debit(0,2,true);o.take(1,1,{},{&debit_record,1});
  const auto candidate=o.candidate_sha256(),accepted=o.accepted_sha256();
  MoverRecord m{};m.native_id=1;m.species=0;m.from_cell=0;m.to_cell=1;m.weight=1;m.cell_count=2;
  auto reject_transport=[&](std::span<const MoverRecord> v,const char *label){
    refuses([&]{o.transport(1,v);},label);
    check(candidate==o.candidate_sha256() && accepted==o.accepted_sha256(),"invalid transport atomic");
  };
  std::array<MoverRecord,2> mm={m,m};reject_transport(mm,"duplicate mover ID");
  mm[1].native_id=2;mm[1].cell_count=3;reject_transport(mm,"inconsistent cell count");
  mm={m,m};mm[1].native_id=2;mm[0].cell_count=mm[1].cell_count=1;reject_transport(mm,"overdraw count");
  auto invalid=m;invalid.weight=std::numeric_limits<double>::infinity();reject_transport({&invalid,1},"nonfinite mover weight");
  invalid=m;invalid.from_cell=2;reject_transport({&invalid,1},"invalid mover origin");
  mm={m,m};mm[1].native_id=2;mm[1].to_cell=pfg::kNone;
  o.transport(1,mm);o.take(1,2,{},{});o.commit(o.candidate_sha256());
  auto restored=Owner::restore(o.encode(),cfg,o.accepted_sha256());
  check(restored->accepted_sha256()==o.accepted_sha256(),"transport nonempty restore");
  check(o.entries().size()==1 && o.entries()[0].cell==1,"transport source cleared and destination retained");
  std::cout<<"motion_state "<<o.accepted_sha256()<<"\n";
}
