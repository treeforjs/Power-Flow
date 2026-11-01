#pragma once
#include "check.hpp"
#include "phase_owner.hpp"
#include <filesystem>
#include <functional>

namespace transaction_tests {
using namespace ss304_phase_v1_3;
using Factory = std::function<std::unique_ptr<Engine>()>;
inline Config config() {
    Config c; c.slices = {{0,0,-0.001,0.0}};
    c.width=0.002; c.thickness=0.0002; c.height=0.0002;
    c.nx_em=2; c.ny_em=1; c.rx=1; c.y_edges={0,0.00005,0.0001,0.00015,0.0002};
    c.evap_mode=0; return c;
}
inline StepInput input(double p=1000) { StepInput x; x.dt=1e-12; x.em_power_W={p,p}; return x; }
inline bool refused(const std::function<void()>& f) { try {f();} catch(const PhaseRefusal&){return true;} return false; }
inline std::string run(const Factory& make, const std::filesystem::path& out) {
    std::filesystem::create_directories(out);
    const auto c=config();
    {
        Owner o(c,make()), ref(c,make()); auto a=input(), b=input(2000);
        auto ta=o.prepare(a); auto tb=o.prepare(b);
        chk::check(!o.preflight(ta) && !o.commit(std::move(ta)), "older trial refused after a different prepare reused its slot");
        chk::check(o.preflight(tb) && o.commit(std::move(tb)), "latest trial commits after replacement");
        ref.advance(b);
        chk::check(o.state_sha256()==ref.state_sha256() && o.ledger_closes(), "replacement publishes matching cells, conductivity and ledger");
        const auto before=o.state_sha256(); auto good=o.prepare(a); auto bad=input(1e20);
        chk::check(refused([&]{o.prepare(bad);}), "post-engine temperature refusal reached");
        chk::check(!o.preflight(good) && !o.commit(std::move(good)), "failed preparation also invalidates prior scratch-backed trial");
        chk::check(o.state_sha256()==before, "failed replacement leaves accepted cells and ledger unchanged");
    }
    {
        Owner o(c,make()); const auto before=o.state_sha256(), quote=o.quote(1e-12).sha256;
        auto bad=input(); bad.back_pressure_Pa.assign(o.n_patches(),1000);
        bad.loads.push_back({0,0,0,0,"","wall",1});
        chk::check(refused([&]{o.prepare(bad);}), "pressure-changing proposal refused on invalid producer");
        chk::check(o.state_sha256()==before && o.quote(1e-12).sha256==quote, "refusal preserves pressure, accepted digest and quote");
        auto good=input(); good.back_pressure_Pa=bad.back_pressure_Pa;
        auto t=o.prepare(good);
        chk::check(o.state_sha256()==before && o.quote(1e-12).sha256==quote, "uncommitted pressure proposal has no accepted-state effect");
        chk::check(o.commit(std::move(t)), "pressure publishes with successful commit");
        Owner ref(c,make()); ref.set_back_pressure(good.back_pressure_Pa); ref.advance(input());
        chk::check(o.state_sha256()==ref.state_sha256() && o.quote(1e-12).sha256==ref.quote(1e-12).sha256, "committed pressure equals explicitly set pressure evolution");
        auto pending=o.prepare(input()); o.set_back_pressure(std::vector<double>(o.n_patches(),2000));
        chk::check(!o.preflight(pending), "explicit pressure change invalidates pending trial");
        auto cp=out/"pressure.chk"; o.save(cp); Owner restored(c,make()); restored.load(cp);
        chk::check(o.state_sha256()==restored.state_sha256() && o.quote(1e-12).sha256==restored.quote(1e-12).sha256, "pressure is included in checkpoint and state identity");
        restored.advance(input()); o.advance(input());
        chk::check(o.state_sha256()==restored.state_sha256(), "pressure checkpoint continues identically");
        const auto accepted=o.state_sha256(); const auto aq=o.quote(1e-12).sha256;
        auto pq=o.quote(1e-12,std::vector<double>(o.n_patches(),3000));
        chk::check(pq.sha256!=aq && o.state_sha256()==accepted && o.quote(1e-12).sha256==aq, "proposed-pressure quote does not publish pressure");
    }
    {
        Owner o(c,make()), ref(c,make()); const std::vector<double> a={1000,2000}, b={3000,4000};
        o.accumulate(a,{},1e-12); StepInput w; w.use_accumulated=true;
        auto t=o.prepare(w); o.accumulate(b,{},1e-12);
        const auto before=o.state_sha256();
        chk::check(!o.preflight(t) && !o.commit(std::move(t)), "extra accumulation invalidates prepared window");
        chk::check(o.window_steps()==2 && o.state_sha256()==before, "stale commit preserves both pending contributions");
        chk::check(refused([&]{o.prepare(input());}) && o.state_sha256()==before, "direct step cannot bypass pending window");
        chk::check(refused([&]{o.save(out/"pending.chk");}) && !std::filesystem::exists(out/"pending.chk"), "pending window checkpoint is refused");
        o.advance(w); ref.accumulate(a,{},1e-12); ref.accumulate(b,{},1e-12); ref.advance(w);
        chk::check(o.state_sha256()==ref.state_sha256() && o.ledger_closes() && o.window_steps()==0, "new trial consumes all Joule contributions exactly once");
        auto cp=out/"consumed.chk"; o.save(cp); Owner restored(c,make()); restored.load(cp);
        chk::check(restored.state_sha256()==o.state_sha256(), "consumed window checkpoint restores identity");
        Owner pending(c,make()); pending.accumulate(a,{},1e-12); const auto ph=pending.state_sha256();
        chk::check(refused([&]{pending.load(cp);}) && pending.state_sha256()==ph, "restart cannot discard an outstanding accumulation window");
    }
    {
        Owner a(c,make()), b(c,make());
        a.accumulate({1000,2000},{},1e-12); b.accumulate({2000,1000},{},1e-12);
        chk::check(a.state_sha256()!=b.state_sha256(), "pending digest binds spatial Joule distribution, not just total or count");
        StepInput w; w.use_accumulated=true;
        auto ta=a.prepare(w),tb=b.prepare(w);
        chk::check(ta.record().inputs_sha256!=tb.record().inputs_sha256, "step input digest binds pending Joule distribution");
        Owner x(c,make()), y(c,make()); std::vector<double> f(12),g(12); f[2]=g[2]=1e8; f[3]=1; g[3]=2;
        x.accumulate({0,0},f,1e-12); y.accumulate({0,0},g,1e-12);
        chk::check(x.state_sha256()!=y.state_sha256(), "pending digest binds Lorentz impulse");
        auto tx=x.prepare(w),ty=y.prepare(w);
        chk::check(tx.record().inputs_sha256!=ty.record().inputs_sha256, "step input digest binds Lorentz impulse even with flow disabled");
        Owner z(c,make()), v(c,make()); z.accumulate({NAN,0},{},1e-12); v.accumulate({0,0},{},1e-12);
        chk::check(z.state_sha256()!=v.state_sha256(), "pending digest binds invalid-input flags");
    }
    Owner result(c,make()); result.set_back_pressure(std::vector<double>(result.n_patches(),500));
    result.accumulate({1200,3400},{},1e-12); StepInput w; w.use_accumulated=true; result.advance(w);
    return result.state_sha256();
}
}
