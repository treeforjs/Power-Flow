// Synthetic admitted event with unchanged input/output momenta; only the
// target's numerical destination is reservoir. This isolates stock bookkeeping.
void reservoir_release_tests(Execution execution) {
  auto cfg=config(execution);cfg.species[0].mass_kg=2.99e-26;
  cfg.processes={{0,2,pfg::kNone,1}};
  auto rec=[](std::uint64_t step){
    RecoilRecord r{};r.accepted_step=step;r.stage=1;r.grid_cell=0;
    auto& e=r.ev;e.process=0;e.weight=2;e.target_species=0;e.n_out=2;
    e.u_in[0]=3000;e.u_in[1]=-2000;e.u_in[2]=1000;
    e.ut[0]=-4000;e.ut[1]=2500;e.ut[2]=1500;
    e.out_species[0]=e.out_species[1]=0;
    e.out_dest[0]=pfg::kDestSurvivor;e.out_dest[1]=pfg::kDestReservoir;
    for(unsigned a=0;a<3;++a){e.out_u[0][a]=e.u_in[a];e.out_u[1][a]=e.ut[a];}
    return r;
  };
  const auto event=rec(0);
  auto release=[&](double weight){FlowRecord f{};f.kind=pfg::kFlowRelease;f.species=0;f.grid_cell=0;f.native_id=100;f.weight=weight;for(unsigned a=0;a<3;++a)f.u[a]=event.ev.ut[a];return f;};
  auto stock=[&](const Owner& o,double weight,bool candidate){
    Exact energy;energy.add_product(weight,pfg::kinetic(cfg.species[0].mass_kg,event.ev.ut));
    check(o.global_digits(pfg::G_RES_E,candidate)==energy.normalised(),"actual pending kinetic stock equals mirror");
    for(unsigned a=0;a<3;++a){Exact p;p.add_triple(weight,cfg.species[0].mass_kg,event.ev.ut[a]);check(o.global_digits(pfg::G_RES_P+a,candidate)==p.normalised(),"actual pending momentum stock equals mirror");}
  };
  Owner o(cfg);const auto before=o.encode();o.begin_step(0);
  o.take(0,1,{&event,1},{});stock(o,2,true);
  auto partial=release(.5);o.take(0,2,{},{&partial,1});stock(o,1.5,true);
  auto prepared=o.prepare_commit(o.candidate_sha256());check(o.try_publish(prepared),"partial reservoir release published");stock(o,1.5,false);
  auto restored=Owner::restore(o.encode(),cfg,o.accepted_sha256());
  for(Owner* owner:{&o,restored.get()}){
    owner->begin_step(1);owner->take(1,1,{},{});
    std::array<FlowRecord,2> remainder={release(.25),release(1.25)};remainder[1].native_id=101;
    owner->take(1,2,{},remainder);stock(*owner,0,true);
    owner->commit(owner->candidate_sha256());stock(*owner,0,false);
  }
  check(o.encode()==restored->encode(),"pending-reservoir restore and exact split release continuation agree");
  const auto accepted=o.encode();o.begin_step(2);const auto another=rec(2);o.take(2,1,{&another,1},{});
  auto all=release(2);o.take(2,2,{},{&all,1});stock(o,0,true);o.abandon();
  check(o.encode()==accepted && before!=accepted,"abandoned populated release leaves accepted inventory unchanged");
  // Birth/sink/transaction flows are distinct and must never drain this mirror.
  Owner unrelated(cfg);unrelated.begin_step(0);unrelated.take(0,1,{&event,1},{});
  for(auto kind:{pfg::kFlowBirth,pfg::kFlowSink,pfg::kFlowTransaction}){
    auto f=release(.5);f.kind=kind;f.negative=kind==pfg::kFlowSink;unrelated.apply_flows(0,{&f,1});stock(unrelated,2,true);
  }
  unrelated.take(0,2,{},{&all,1});stock(unrelated,0,true);unrelated.abandon();
  // Algebraic signed population flow convention is retained (no new pathway).
  Owner inverse(cfg);inverse.begin_step(0);inverse.take(0,1,{},{});auto back=release(.5);back.negative=1;
  inverse.take(0,2,{},{&back,1});stock(inverse,.5,true);inverse.abandon();
  std::cout<<"reservoir_release_state "<<o.accepted_sha256()<<"\n";
}
