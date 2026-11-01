// Isolated diagnostic channel, not a complete H2 encounter group or native binding.
#pragma once
#include <cmath>
#include <stdexcept>
#include "chem_event.h"
#include "generated_channel.hpp"

namespace h2diag {
enum class Mode { refused, diagnostic_primary_cdf_iso };
enum class EnergyFrame { target_rest, laboratory };
struct Result {
    pfc::V3 incoming[2], outgoing[3], target_u;
    double kinetic_in[2], kinetic_out[3], sigma_m2, internal_eV, target_rest_eV;
    double quantile, quantile_hi, fragment_KER_eV, kinematic_loss_eV;
};
inline void require(bool ok, const char* why) { if (!ok) throw std::runtime_error(why); }
inline double threshold() {
    const double m=species[0].mass_eV, M=species[1].mass_eV;
    const double loss=reference_cost_eV+2*fragment_eV[0];
    return loss*(2*(m+M)+loss)/(2*M);
}
inline double cross_section(double eV, EnergyFrame frame, Mode mode=Mode::refused) {
    require(mode==Mode::diagnostic_primary_cdf_iso,"explicit diagnostic mode required");
    require(frame==EnergyFrame::target_rest,"target-rest energy required");
    require(std::isfinite(eV) && eV>=energies[0] && eV<=energies[point_count-1],"unsupported source energy");
    const int i=pfc::grid_locate(energies,point_count,eV);
    const double s=pfc::row_eval(cross_sections_m2,energies,point_count,i,eV);
    require(std::isfinite(s) && s>0,"nonpositive supported source rate");
    return s;
}
inline Result sample(double eV, EnergyFrame frame, pfc::V3 incident_direction,
                     pfc::V3 target_u, const pfc::ChemKey& key, pfc_u32 candidate,
                     Mode mode=Mode::refused) {
    Result out{};
    out.sigma_m2=cross_section(eV,frame,mode); // All admission checks precede the stateless RNG.
    require(pfc::chem_key_valid(key,candidate),"RNG domain");
    require(key.purpose==PFC_PURPOSE_FIXTURE_TEST,"diagnostic RNG purpose required");
    require(std::isfinite(incident_direction.x) && std::isfinite(incident_direction.y) &&
            std::isfinite(incident_direction.z) && std::fabs(pfc::norm2(incident_direction)-1)<=0x1p-45,
            "unit incident direction required");
    require(std::isfinite(target_u.x) && std::isfinite(target_u.y) && std::isfinite(target_u.z),"target velocity");
    const double m=species[0].mass_eV, M=species[1].mass_eV;
    const pfc::V3 wt=pfc::w_from_u(target_u);
    const double gt=std::sqrt(1+pfc::norm2(wt));
    require(std::isfinite(gt),"target boost overflow");
    const pfc::V3 pp=pfc::scale(incident_direction,std::sqrt(eV*(eV+2*m)));
    out.incoming[0]=pfc::boost_out(wt,gt,pp,m+eV);
    out.incoming[1]=pfc::p_from_u(M,target_u);
    out.kinetic_in[0]=pfc::kinetic_from_p2(m,pfc::norm2(out.incoming[0]));
    out.kinetic_in[1]=pfc::ke_from_u(M,target_u);
    // Stable total collision-CM kinetic energy, retaining finite target recoil.
    const double total=m+M;
    const double Kcm=2*M*eV/(std::sqrt(total*total+2*M*eV)+total);
    const double cap=(Kcm-reference_cost_eV)/2;
    require(std::isfinite(cap) && cap>fragment_eV[0],"no supported source fragment phase space");
    double upper=quantiles[cdf_count-1];
    if(cap<fragment_eV[cdf_count-1]) {
        const int j=pfc::grid_locate(fragment_eV,cdf_count,cap);
        upper=pfc::row_eval(quantiles,fragment_eV,cdf_count,j,cap);
    }
    require(upper>quantiles[0] && upper<=quantiles[cdf_count-1],"empty conditional probability support");
    static_assert(PFC_RNG_LAYOUT_VERSION==2 && PFC_W_H2_PARTITION>=PFC_W_KIN_END,"disjoint qualified RNG layout required");
    const double u=pfc::chem_u01(key,candidate,PFC_W_H2_PARTITION);
    out.quantile=quantiles[0]+u*(upper-quantiles[0]);out.quantile_hi=upper;
    require(out.quantile>quantiles[0] && out.quantile<upper,"conditional draw lost open support");
    const int j=pfc::grid_locate(quantiles,cdf_count,out.quantile);
    const double fragment=pfc::row_eval(fragment_eV,quantiles,cdf_count,j,out.quantile);
    require(fragment>=fragment_eV[0] && fragment<cap && fragment<=fragment_eV[cdf_count-1],"fragment outside retained support");
    const double loss_eV=reference_cost_eV+2*fragment;
    require(loss_eV<Kcm,"rounded event has no electron phase space");
    out.kinematic_loss_eV=loss_eV;out.fragment_KER_eV=loss_eV-reference_cost_eV;
    pfc::TableView table{};table.n_species=3;table.species=species;
    pfc::GroupRec group{};
    group.projectile=0;group.target=1;group.xi_off=0xffffffffu;
    pfc::ProcessRec process{};
    process.kind=PFC_P_DISSOC;process.angular=PFC_ANG_ISO;process.tag=4;
    process.n_out=3;process.n_chain=2;
    process.out_species[0]=0;process.out_species[1]=2;process.out_species[2]=2;
    process.out_dest[0]=PFC_DEST_SURVIVOR;
    process.out_dest[1]=PFC_DEST_RESERVOIR;process.out_dest[2]=PFC_DEST_RESERVOIR;
    process.loss_eV=loss_eV;process.q_eV=reference_cost_eV;process.chain_q_eV=reference_cost_eV;
    process.kinetic_fraction=1;process.internal_sink=PFC_SINK_NONE;
    double mk[PFC_MAX_OUT];pfc::V3 momenta[PFC_MAX_OUT];
    require(pfc::event_kinematics(table,process,group,0,key,candidate,m,M,pp,eV,wt,gt,
                                momenta,mk,&out.internal_eV),"kinematically closed");
    out.target_rest_eV=eV;out.target_u=target_u;
    for (unsigned i=0;i<3;++i) {
        out.outgoing[i]=momenta[i];
        const auto u=pfc::u_from_p(mk[i],momenta[i]);
        out.kinetic_out[i]=pfc::ke_from_u(species[process.out_species[i]].mass_eV,u);
        require(std::isfinite(out.kinetic_out[i]) && out.kinetic_out[i]>=0,"invalid product energy");
    }
    // The existing per-event chemical-energy tolerance is retained, not enlarged.
    const double defect=out.kinetic_in[0]+out.kinetic_in[1]-out.kinetic_out[0]-out.kinetic_out[1]-out.kinetic_out[2]-reference_cost_eV-out.internal_eV;
    const double scale=out.kinetic_in[0]+out.kinetic_in[1]+out.kinetic_out[0]+out.kinetic_out[1]+out.kinetic_out[2]+std::fabs(reference_cost_eV)+std::fabs(out.internal_eV);
    require(std::fabs(defect)<=0x1p-40*scale+1e-300,"event energy closure");
    require(out.internal_eV==0,"ground atomic products cannot absorb unspecified internal energy");
    return out;
}
} // namespace h2diag
