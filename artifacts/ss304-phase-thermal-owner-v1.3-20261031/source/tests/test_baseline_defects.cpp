// Read-only reproduction against historical v1.2 headers. Passing means the old defect was observed.
#include "check.hpp"
#include "/home/treefor/projects/runtimes/ss304-phase-thermal-owner-v1.2-20261031/source/phase_owner.hpp"
using namespace ss304_phase_v1_2;
static Config cfg(){Config c;c.slices={{0,0,-0.001,0}};c.width=.002;c.thickness=.0002;c.height=.0002;c.nx_em=2;c.ny_em=1;c.y_edges={0,.00005,.0001,.00015,.0002};c.evap_mode=0;return c;}
static StepInput in(double p){StepInput x;x.dt=1e-12;x.em_power_W={p,p};return x;}
int main(int argc,char**argv){
    auto c=cfg();
    {Owner o(c),ref(c);auto a=o.prepare(in(1000));auto b=o.prepare(in(2000));ref.advance(in(2000));
     bool accepted=o.commit(std::move(a));auto x=o.state(),y=ref.state();bool cells=true;
     for(size_t k=0;k<x.E.size();++k)cells=cells&&px_eq(x.E[k],y.E[k]);
     chk::check(accepted&&cells&&!px_eq(o.ledger().joule,ref.ledger().joule),"v1.2 reproduced: old trial commits newer cells with older Joule ledger");}
    {Owner o(c);auto q=o.quote(1e-12).sha256, h=o.state_sha256();auto x=in(0);x.back_pressure_Pa.assign(o.n_patches(),1000);x.loads.push_back({0,0,0,0,"","wall",1});
     bool refused=false;try{o.prepare(x);}catch(const PhaseRefusal&){refused=true;}
     chk::check(refused&&o.quote(1e-12).sha256!=q&&o.state_sha256()==h,"v1.2 reproduced: refusal changes pressure quote without changing state digest");}
    {Owner o(c),ref(c);o.accumulate({1000,1000},{},1e-12);StepInput w;w.use_accumulated=true;auto t=o.prepare(w);o.accumulate({2000,2000},{},1e-12);
     ref.accumulate({1000,1000},{},1e-12);ref.advance(w);
     chk::check(o.commit(std::move(t))&&o.window_steps()==0&&px_eq(o.ledger().joule,ref.ledger().joule),"v1.2 reproduced: stale commit clears newer accumulated Joule contribution");}
    return chk::finish("test_baseline_defects",argc>1?argv[1]:".");
}
