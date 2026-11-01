#pragma once
#include "root_api.hpp"
#include "pf_gas_v1/gas_owner.hpp"
#include <bit>
#include <map>

// Independent stock operands. No component energy-residual value is used.
// The caller retains the actual Hall candidate and publication capability.
namespace powerflow::common_conservation_v1 {
namespace h=::hall_corrective_v3;
namespace g=::pf_gas_v1;
namespace eh=::electrode_heat_v2;
inline void need(bool x,const char*s){if(!x)throw std::invalid_argument(std::string("common conservation: ")+s);}
inline bool bits(double a,double b){return std::bit_cast<std::uint64_t>(a)==std::bit_cast<std::uint64_t>(b);}
inline g::Exact difference(const g::Exact&a,const g::Exact&b){auto r=a;r.merge(b,1);return r;}
struct FieldStock {g::Exact electric_J,magnetic_J;g::Exact total()const{auto r=electric_J;r.merge(magnetic_J);return r;}};
inline FieldStock field_stock(const axial::Mesh&mesh,const axial::State&state){
  need(state.E.size()==std::size_t(mesh.ne)&&state.B.size()==std::size_t(mesh.nb)&&
       state.E_low.size()==state.E.size()&&state.B_low.size()==state.B.size(),"complete compensated field state required");
  FieldStock out;
  auto add=[](g::Exact&sum,const std::vector<double>&hi,const std::vector<double>&lo,const std::vector<double>&metric){
    need(hi.size()==metric.size(),"field metric shape");
    for(std::size_t i=0;i<hi.size();++i){need(std::isfinite(hi[i])&&std::isfinite(lo[i])&&std::isfinite(metric[i])&&metric[i]>=0,"field stock operands");
      const double half=.5*metric[i];need(metric[i]==0||(half>0&&2*half==metric[i]),"field metric underflow");
      sum.add_triple(half,hi[i],hi[i]);sum.add_triple(metric[i],hi[i],lo[i]);sum.add_triple(half,lo[i],lo[i]);}
  };
  add(out.electric_J,state.E,state.E_low,mesh.M);add(out.magnetic_J,state.B,state.B_low,mesh.H);return out;
}
struct ThermalSliceStock {
  g::Exact before_J,after_J,stable_delta_J;
  double accepted_ordered_J=0,candidate_ordered_J=0;
  std::size_t nodes=0;
};
inline std::vector<double> lumped_control_areas(const eh::Geometry&geometry){
  need(geometry.nx>0&&geometry.ny>0&&geometry.nx<=1024&&geometry.ny<=1024,"bounded Q1 thermal geometry");
  mfem::Mesh mesh=mfem::Mesh::MakeCartesian2D(geometry.nx,geometry.ny,mfem::Element::QUADRILATERAL,true,geometry.width,geometry.thickness);
  mfem::H1_FECollection collection(1,2);mfem::FiniteElementSpace space(&mesh,&collection);
  std::vector<double>x,y;
  for(int v=0;v<mesh.GetNV();++v){auto p=mesh.GetVertex(v);p[0]+=geometry.x0;p[1]+=geometry.y0;x.push_back(p[0]);y.push_back(p[1]);}
  for(auto*axis:{&x,&y}){std::sort(axis->begin(),axis->end());axis->erase(std::unique(axis->begin(),axis->end()),axis->end());}
  need(x.size()==std::size_t(geometry.nx+1)&&y.size()==std::size_t(geometry.ny+1),"thermal axes");
  std::vector<double>dx(x.size()),dy(y.size()),area(space.GetVSize());std::vector<bool>seen(area.size());
  for(std::size_t i=0;i<x.size();++i)dx[i]=.5*((i?x[i]-x[i-1]:0)+(i+1<x.size()?x[i+1]-x[i]:0));
  for(std::size_t i=0;i<y.size();++i)dy[i]=.5*((i?y[i]-y[i-1]:0)+(i+1<y.size()?y[i+1]-y[i]:0));
  for(int v=0;v<mesh.GetNV();++v){auto p=mesh.GetVertex(v);auto i=std::lower_bound(x.begin(),x.end(),p[0])-x.begin();auto j=std::lower_bound(y.begin(),y.end(),p[1])-y.begin();
    mfem::Array<int>d;space.GetVertexDofs(v,d);need(d.Size()==1&&d[0]>=0&&std::size_t(d[0])<area.size()&&!seen[d[0]],"Q1 scalar nodal ownership");seen[d[0]]=true;area[d[0]]=dx[i]*dy[j];}
  need(std::all_of(seen.begin(),seen.end(),[](bool v){return v;}),"thermal node coverage");return area;
}
inline ThermalSliceStock thermal_slice_stock(const eh::Coupler&accepted,const eh::State&candidate){
  const auto&geometry=accepted.geometry();const auto&before=accepted.state();const auto area=lumped_control_areas(geometry);
  need(before.T.size()==area.size()&&candidate.T.size()==area.size(),"thermal state/node mismatch");
  // This is exactly the pinned material admitted by Hall ThermalV3, not a new law.
  const auto law=eh::conditional_ss304_v2(ss304v2::Qualification::evaluated_solid_1500,
    ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation,8000.,eh::material_headers_sha256);
  ThermalSliceStock out;out.nodes=area.size();
  for(std::size_t i=0;i<area.size();++i){const double h0=law.energy_reference(before.T[i]),h1=law.energy_reference(candidate.T[i]),dh=law.delta_h(candidate.T[i],before.T[i]);
    out.before_J.add_triple(area[i],h0,geometry.height);out.after_J.add_triple(area[i],h1,geometry.height);out.stable_delta_J.add_triple(area[i],dh,geometry.height);
    // Independent construction must reproduce the accepted public evaluator's
    // exact operation order before candidate evaluation is usable.
    out.accepted_ordered_J+=area[i]*h0*geometry.height;out.candidate_ordered_J+=area[i]*h1*geometry.height;
  }
  need(bits(out.accepted_ordered_J,accepted.stored_energy_J()),"independent thermal operands differ from accepted Coupler evaluator");return out;
}
struct FieldThermalAssessment {
  FieldStock field_before,field_after;
  g::Exact field_delta_J,thermal_before_J,thermal_after_J,thermal_delta_J,thermal_stable_delta_J,thermal_reference_rounding_J;
  g::Exact circuit_input_J,outer_input_J,outer_loss_J;
  std::size_t thermal_slices=0,thermal_nodes=0;
};
inline FieldThermalAssessment field_thermal_assessment(const h::Owner&before,const h::FieldHeatProposal&candidate,double dt){
  need(candidate.configuration_sha256()==before.configuration_sha256()&&candidate.endpoint_revision()==before.revision()+1,"actual Hall candidate configuration/revision");
  const auto&mesh=before.mesh();const auto&state=candidate.candidate_field();const auto&heat=candidate.heat();
  need(heat.bodies.size()==before.thermal().slices(),"all thermal candidate slices required");
  FieldThermalAssessment out;out.field_before=field_stock(mesh,before.field());out.field_after=field_stock(mesh,state);out.field_delta_J=difference(out.field_after.total(),out.field_before.total());
  for(int k=mesh.kshort;k<mesh.kfeed;++k)for(int body=0;body<2;++body){const auto i=before.thermal().index(body,k);const auto&after=heat.bodies.at(i).state();
    need(after.time==state.time_s&&after.revision==state.revision,"thermal candidate clock");const auto slice=thermal_slice_stock(before.thermal().slice(body,k),after);
    out.thermal_before_J.merge(slice.before_J);out.thermal_after_J.merge(slice.after_J);out.thermal_stable_delta_J.merge(slice.stable_delta_J);++out.thermal_slices;out.thermal_nodes+=slice.nodes;
  }
  out.thermal_delta_J=difference(out.thermal_after_J,out.thermal_before_J);
  out.thermal_reference_rounding_J=difference(out.thermal_delta_J,out.thermal_stable_delta_J);
  need(std::isfinite(dt)&&dt>0&&before.time_s()+dt==candidate.end_s(),"field candidate duration");
  const auto&initial=before.model().initial_fields();const auto&E=candidate.midpoint().E();
  const double I=before.model().drive().current(before.time_s()+.5*dt);
  need(initial.feed.size()==E.size()&&initial.background_force.size()==E.size()&&mesh.Y.size()==E.size(),"external field operand shape");
  for(std::size_t i=0;i<E.size();++i){
    // These are independent terminal/exterior operand sums, not field residuals.
    // Products use the exact product of three binary64 operands; I*feed and
    // Y*E are explicit model operator operands retained as binary64.
    out.circuit_input_J.add_triple(-dt,I*initial.feed[i],E[i]);
    out.outer_input_J.add_triple(dt,initial.background_force[i],E[i]);
    out.outer_loss_J.add_triple(dt,mesh.Y[i]*E[i],E[i]);
  }
  return out;
}
inline FieldThermalAssessment field_thermal_assessment(const h::Owner&before,const h::JointProposal&candidate){
  return field_thermal_assessment(before,candidate.field_heat(),candidate.motion().binding().dt_s);
}
} // namespace powerflow::common_conservation_v1
