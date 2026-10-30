#include "water_parent_bridge.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>
using namespace raikou::pic_water;
namespace pc=raikou::pic_collisions;
int checks=0;
void check(bool ok,const char* message){++checks;if(!ok)throw pc::Error(message);}
template<class F>void rejects(F f,const char* message){bool caught=false;try{f();}catch(const std::exception&){caught=true;}check(caught,message);}
bool close(double a,double b,double tolerance=2e-13){return std::abs(a-b)<=tolerance*std::max({std::abs(a),std::abs(b),1e-300});}
int main(int argc,char**argv){try{
    pc::require(argc==5,"dataset, kernel hash, data hash, result path required");
    const auto data=std::make_shared<const Data>(Data::read(argv[1]));
    const std::string hash(64,'a');constexpr double mn=2.9915e-26;
    ParentBridge bridge(data,parent_model,argv[3],argv[2],mn);
    const double volume=25e-6*25e-6*.025;
    const std::vector<pc::Cell> cells{{7,volume}};
    auto electron=[](double E){const double q=E*pc::qe/(me*pc::c*pc::c);return pc::Packet{"e0","electron",7,.25,me,{pc::c*std::sqrt(q*(q+2)),0,0}};};
    pc::Packet neutral{"water0","H2O",7,64,mn,{0,0,0}};
    rejects([&]{ParentBridge wrong(data,"",argv[3],argv[2],mn);},"conditional model must be selected");
    check(pc::candidates({electron(12)},{neutral},cells,{bridge.channel()}).empty(),"below-threshold zero does not invent event");
    rejects([&]{pc::candidates({electron(13)},{neutral},cells,{bridge.channel()});},"threshold gap not interpolated away");
    rejects([&]{pc::candidates({electron(8000)},{neutral},cells,{bridge.channel()});},"8keV unsupported parent tail refused");
    std::mt19937_64 rng(8131);std::uniform_real_distribution<double> uniform(0,1);
    const std::vector<Vec> boosts{{0,0,0},{312,-840,11},{.00005*pc::c,.00002*pc::c,-.00003*pc::c},{-.0001*pc::c,.00005*pc::c,0},{.001*pc::c,0,0}};
    double max_energy_fraction=0,max_momentum_fraction=0;int events=0;
    for(double E:{14.,20.,100.,500.,999.})for(const auto& boost:boosts)for(int sample=0;sample<20;++sample){
        auto e=electron(E),n=neutral;e.proper_velocity_m_s=to_lab(e.proper_velocity_m_s,boost);n.proper_velocity_m_s=boost;
        const auto before_e=e,before_n=n;
        const auto rows=pc::candidates({e},{n},cells,{bridge.channel()},1,1,1);check(rows.size()==1,"positive actual parent cross section");
        check(close(rows[0].energy_eV,E),"invariant energy under common boost");
        const double z=2*uniform(rng)-1,phi=2*std::numbers::pi*uniform(rng),s=std::sqrt(1-z*z);
        const auto out=bridge.sample(rows[0],e,n,{s*std::cos(phi),s*std::sin(phi),z});++events;
        check(out.particles[0].charge_e+out.particles[1].charge_e+out.particles[2].charge_e==-1,"charge conserved");
        check(out.particles[2].state=="H2O+"&&out.particles[2].reference_J==threshold_eV*pc::qe,"parent retains2H1O and declared reference");
        check(out.projectile_remainder.require_single_weight()==0&&out.target_remainder.require_single_weight()==63.75,"finite one-quarter parent debit");
        check(out.particles[0].weight==.25&&out.particles[1].weight==.25&&out.particles[2].weight==.25,"same reacting cohort for all products");
        check(e.proper_velocity_m_s==before_e.proper_velocity_m_s&&n.proper_velocity_m_s==before_n.proper_velocity_m_s&&e.weight==before_e.weight&&n.weight==before_n.weight,"candidate leaves native owners unchanged");
        check(std::abs(out.energy_defect_J)<=out.energy_roundoff_limit_J,"lab energy accounted once");
        max_energy_fraction=std::max(max_energy_fraction,std::abs(out.energy_defect_J)/out.energy_roundoff_limit_J);
        for(double d:out.momentum_defect_kg_m_s){check(std::abs(d)<=out.momentum_roundoff_limit_kg_m_s,"lab momentum closes");max_momentum_fraction=std::max(max_momentum_fraction,std::abs(d)/out.momentum_roundoff_limit_kg_m_s);}
    }
    auto e=electron(100);auto rows=pc::candidates({e},{neutral},cells,{bridge.channel()});
    const double represented_rate=e.weight*neutral.weight*rows[0].sigma_m2*rows[0].flux_m_s/volume;
    check(close(rows[0].rate_s*rows[0].event_weight,represented_rate),"finite physical represented rate with depth exactly once");
    auto huge=neutral;huge.weight=1e20;const auto bigrows=pc::candidates({e},{huge},cells,{bridge.channel()});
    const auto big=bridge.sample(bigrows[0],e,huge,{0,1,0});
    check(big.target_remainder.high==1e20&&big.target_remainder.low==-.25,"low neutral debit survives lab-product bridge");
    rejects([&]{big.target_remainder.require_single_weight();},"bridge never silently collapses neutral reservoir");
    auto bad=rows[0];bad.target_before+=1;rejects([&]{bridge.sample(bad,e,neutral,{0,1,0});},"stale weight refused");
    bad=rows[0];bad.energy_eV*=2;rejects([&]{bridge.sample(bad,e,neutral,{0,1,0});},"stale energy refused");
    bad=rows[0];bad.data_sha256=hash;rejects([&]{bridge.sample(bad,e,neutral,{0,1,0});},"foreign dataset refused");
    rejects([&]{bridge.sample(rows[0],e,neutral,{0,2,0});},"nonunit product axis refused");
    auto wrong=e;wrong.mass_kg*=2;rejects([&]{bridge.sample(rows[0],wrong,neutral,{0,1,0});},"wrong projectile mass refused");
    auto wrong_target=neutral;wrong_target.mass_kg*=2;rejects([&]{bridge.sample(rows[0],e,wrong_target,{0,1,0});},"changed catalog target mass refused");
    for(const Vec boost:{Vec{0,1e5*pc::c,0},Vec{.8*pc::c,0,0}}){auto fast_e=electron(14),fast_n=neutral;fast_e.proper_velocity_m_s=to_lab(fast_e.proper_velocity_m_s,boost);fast_n.proper_velocity_m_s=boost;auto fast=pc::candidates({fast_e},{fast_n},cells,{bridge.channel()});rejects([&]{bridge.sample(fast[0],fast_e,fast_n,{0,1,0});},"numerically unresolved bulk boost refused");}
    auto tiny_e=e,tiny_n=neutral;tiny_e.weight=1e-300;tiny_n.weight=64e-300;const auto tiny_rows=pc::candidates({tiny_e},{tiny_n},cells,{bridge.channel()});
    rejects([&]{bridge.sample(tiny_rows[0],tiny_e,tiny_n,{0,0,1+4e-13});},"tiny weight cannot hide unweighted conservation defect");
    rejects([&]{bridge.sample(tiny_rows[0],tiny_e,tiny_n,{0,1,0});},"tiny weighted receipt underflow refused");
    // A constant-population optical-depth test only, not physical evolution.
    // Actual sigma/flux/volume are used; no fabricated rate replaces the data.
    pc::Scheduler clock(51);const double rate=rows[0].rate_s;
    pc::Clock stage{1,0,0,20/rate,hash,hash,hash,10/rate,0};
    long double sum=0;constexpr int trials=20000;
    for(int i=0;i<trials;++i){const auto proposal=clock.prepare(stage,rows);pc::require(proposal.event().has_value(),"unexpected optical-depth truncation");sum+=proposal.event_subflow_elapsed_s()*rate;pc::require(clock.try_publish(proposal),"scheduler publish");}
    const double mean=static_cast<double>(sum/trials);check(std::abs(mean-1)<6/std::sqrt(trials),"actual-rate exponential optical depth mean");
    pc::Clock tiny{1,0,0,1e-9,hash,hash,hash,.5e-9,0};
    const auto noevent=clock.prepare(tiny,rows);check(!noevent.event(),"actual dilute finite interval no artificial forced reaction");
    std::ofstream out(argv[4]);out<<std::setprecision(17)<<"{\"passed\":true,\"checks\":"<<checks<<",\"lab_product_events\":"<<events<<",\"constant_rate_trials\":"<<trials<<",\"mean_optical_depth\":"<<mean<<",\"max_energy_bound_fraction\":"<<max_energy_fraction<<",\"max_momentum_bound_fraction\":"<<max_momentum_fraction<<",\"native_inventory_mutation\":false,\"native_field_evolution\":false,\"full_fragmentation\":false,\"parent_model\":\""<<parent_model<<"\",\"tested_incident_eV\":[14,20,100,500,999],\"number_density_display_unit\":\"cm^-3\"}\n";
    pc::require(bool(out),"result write");std::cout<<checks<<" water bridge checks passed\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
