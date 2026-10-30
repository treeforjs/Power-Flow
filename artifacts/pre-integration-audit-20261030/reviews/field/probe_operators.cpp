#include "axial_fields.hpp"
#include "axial_particles.hpp"
#include <iostream>

// Narrow operator observations only: no canonical arena, SSA, heat owner,
// source count/ID allocation, joint proposal, accepted commit or restart.
int main(){try{
    axial::Mesh g;axial::Vec sigma(g.cells.size());
    for(std::size_t c=0;c<sigma.size();++c)if(g.cells[c].body>=0)sigma[c]=1/(71.9e-8);
    const double t=7.4978805549199997e-10,dt=1e-14,I=6012.0713604700004;
    axial::Fields fields(g,sigma,I);const auto initial=fields.start(t);
    auto divergence=[&](const axial::State&s){auto d=g.magnetic_divergence(s.B),low=g.magnetic_divergence(s.B_low);
        for(std::size_t i=0;i<d.size();++i)d[i]+=low[i];return axial::maximum(d);};
    const auto qfield=g.field_charge(initial.E);axial::Vec qauthored(g.nn);
    axial::point_deposit(g,qauthored,{.00021,-.00017,.001},-pic_e*1e9);
    double initial_vacuum=0,authored_mismatch=0,authored_max=0;
    for(int n=0;n<g.nn;++n)if(g.node_body[n]<0&&!g.short_node[n]&&!g.outer_node[n]){
        initial_vacuum=std::max(initial_vacuum,std::abs(qfield[n]));
        authored_mismatch=std::max(authored_mismatch,std::abs(qfield[n]-qauthored[n]));
        authored_max=std::max(authored_max,std::abs(qauthored[n]));}
    const auto step=fields.prepare(initial,axial::Vec(g.ne),sigma,dt,I+.008);
    std::size_t differing=0;double max_difference=0;
    for(int e=0;e<g.ne;++e){const double midpoint=.5*(initial.E[e]+step.first.E[e]);
        if(midpoint!=step.second.mid_E[e])++differing;
        max_difference=std::max(max_difference,std::abs(midpoint-step.second.mid_E[e]));}
    auto perturbed=initial;perturbed.B[g.bx(22,22,g.kshort+1)]+=.01;
    Json bad{{"divB_before_T_per_m",divergence(perturbed)}};
    try{auto trial=fields.prepare(perturbed,axial::Vec(g.ne),sigma,dt,I);
        bad["accepted_by_Fields_prepare"]=true;bad["divB_after_T_per_m"]=divergence(trial.first);
        bad["gauss_increment_C"]=trial.second.gauss_increment_C;bad["energy_balance_J"]=trial.second.balance_J;
    }catch(const std::exception&e){bad["accepted_by_Fields_prepare"]=false;bad["rejection"]=e.what();}
    Json work=Json::array();const double qm=-pic_e/pic_me;
    for(double factor:{0.,.1,1.,10.})for(foil::V3 E:std::vector<foil::V3>{{0,0,0},{2e6,-3e6,1e6}}){
        foil::V3 u{factor*pic_c,.2*factor*pic_c,-.1*factor*pic_c},B{.4,-.3,.2};
        auto [next,D]=foil::work_kick(u,qm,E,B,dt);
        raikou::mfem_backend::GapParticle p;p.mass_kg=pic_me;p.proper_velocity_m_s=u;
        auto q=p;q.proper_velocity_m_s=next;
        const double K0=pic_kinetic_energy_J(p),K1=pic_kinetic_energy_J(q),W=-pic_e*hall::dot(E,D);
        work.push_back({{"proper_velocity_scale_over_c",factor},{"E",E},
            {"momentum_identical_to_native_Boris",next==pic_boris(u,qm,E,B,dt)},
            {"delta_K_J",K1-K0},{"work_J",W},{"work_defect_J",K1-K0-W},
            {"floating_operand_bound_J",128*std::numeric_limits<double>::epsilon()*(K0+K1+std::abs(W))},
            {"displacement_speed_over_c",std::sqrt(hall::dot(D,D))/(dt*pic_c)}});
    }
    std::cout<<Json{{"kind","narrow_operator_audit_not_full_coupled_qualification"},
        {"initial_divB_T_per_m",divergence(initial)},
        {"initial_vacuum_field_charge_max_C",initial_vacuum},
        {"manufactured_authored_charge_max_C",authored_max},
        {"initial_field_vs_manufactured_authored_charge_max_C",authored_mismatch},
        {"authored_particle_or_ID_created",false},
        {"midpoint",{{"different_edges",differing},{"edges",g.ne},{"max_abs_difference_V_per_m",max_difference},
            {"comparison","arithmetic high-word midpoint used by Hall joint vs actual Maxwell receipt.mid_E"}}},
        {"invalid_divB_input",bad},{"relativistic_single_kicks",work},
        {"actual_SSA",false},{"canonical_arena",false},{"source_sampling",false},{"heat_owner",false},
        {"common_joint_execution",false},{"accepted_commit",false}}.dump(2)<<'\n';
    return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
