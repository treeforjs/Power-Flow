// Minimal reproducer: two residual A08 gaps in heat/source owner v4 (c65ac51),
// on the audit's own 300 K, F=0 source fixture.
//  (1) An admitted tiny-weight birth (w=1e-301, all LayerExchange operands
//      nonzero) is handed off; the ordinary-handoff current arrays contain
//      entries that are nonzero at unit weight but exactly zero at weight w.
//      No Refusal is raised, and the component trial preflights.
//  (2) prepare_patch admits area fraction 3.1622776601683792e-245 with a
//      nonzero scalar expected count about twice the exact finite-log count.
// Prints REPRODUCED/NOT-REPRODUCED per gap; exit 0 when the program ran.
// Same compiler flags, include paths and libraries as the v4 test target.
#include "hall-owner/development/hall-finite-foil-em-20261029/axial_surface_path_port.hpp"
#include "source_owner.hpp"
#include "ss304_heat_material_v2.hpp"
#include <iostream>
using namespace electrode_heat_v3;
namespace s4=electrode_heat_v4;
int main(){
 try{
  std::cout<<std::setprecision(17);
  eh::Geometry geo;geo.nx=4;geo.ny=2;geo.revision="thermal-mesh-v3";geo.electrode="upper";
  auto mat=eh::conditional_ss304_v2(ss304v2::Qualification::source_conditioned_liquid_3000,ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation,8000,eh::material_headers_sha256);
  eh::Patch patch;patch.id="lower-face";patch.face=eh::Face::LowY;patch.source_hash="v3-source-scenario";
  se::LayerPolicy pol;pol.local_linear_vacuum_extent_m=25e-6;
  s4::Owner o(geo,mat,{patch},{{"upper-foil-lower-surface","conditional-single-layer",.005*.025}},{{se::Config{4.3},pol,"minimal-reproducer-marks"}},300,"v3-source-scenario",5e-9);
  axial::Mesh g(45,4,.009,1);foil::Vec E(g.ne),B(g.nb); // F=0, B=0
  TrialFrame f;f.role=FieldRole::TotalMidpoint;f.accepted_revision=o.thermal_state().revision;f.trial_stage=2;f.mesh_revision=1;f.material_revision=2;f.begin_s=o.thermal_state().time;f.dt_s=1e-14;f.end_s=f.begin_s+f.dt_s;f.field_stage_s=f.begin_s+.5e-14;f.arena_hash="accepted-arena-revision-0";f.field_hash="immutable-midpoint-trial-2";f.mesh_hash="Hall-axial-sealed-mesh";f.material_hash="global-material-v2-map";f.species_hash="actual-source-v1-electron";
  auto field=std::make_shared<const s4::Field>(g,f,E,B);s4::FinalFieldView view{field->frame(),g,E,B};
  PatchBinding b;b.thermal_revision=o.thermal_state().revision;b.source_accepted_state=b.thermal_revision+1;b.patch_id=7;b.thermal_mesh_identity="thermal-mesh-v3";b.thermal_material_hash=eh::material_headers_sha256;b.material_phase_admitted=b.planar_barrier_admitted=b.electron_supply_admitted=true;b.temperature=HeldTemperature::AcceptedPatchMean;b.extrusion_center_z_m=0;b.extrusion_height_m=geo.height;
  auto prep=[&](double area){return o.prepare_patch(o.make_epoch(1e-14,"trial-midpoint-field"),0,b,field,{.0001,.0025,.0002},Exposure{f.begin_s,1e-14,area,"tiny-positive-area"});};
  auto birth=[&](const s4::PreparedPatch&p,double w){BirthSample s;s.physical_weight=w;s.normal_uniform=.6;s.transverse_uniform=.35;s.azimuth_uniform=.1;s.time_uniform=.23;s.event_id=123;s.canonical_particle_id=42;s.count_owner_revision=0;s.count_owner_base_hash="unchanged-accepted-source-count-owner";s.count_owner_token="staged-canonical-fractional-count-receipt";s.exposure_key=p.exposure_key();s.time_policy=BirthTime::ExposureBeginLumped;s.weight_provenance=WeightProvenance::StagedCanonicalOwner;return o.begin_source(p,s,{"source-v3","e-","electrode-electron-reservoir","native-layer","canonical-electron-42",1,1});};
  auto electron=[](double w){s4::Particle c;c.id=42;c.mass_kg=se::me;c.charge_C=-se::qe;c.weight=w;return c;};
  // (1)
  auto p1=prep(1);const double w=1e-301;
  auto unit=birth(p1,1);auto tiny=birth(p1,w);
  bool all_nonzero=true;for(double x:{tiny.exchange().particle_charge_change_C,tiny.exchange().wall_heat_J,tiny.exchange().particle_K_change_J})all_nonzero=all_nonzero&&x!=0;
  auto hu=o.handoff(unit,view,electron(1));auto ht=o.handoff(tiny,view,electron(w));
  int zeroed=0,unit_nonzero=0;long double lost=0;
  for(std::size_t i=0;i<hu.current().current_load_C_m_s.size();++i){double u=hu.current().current_load_C_m_s[i];if(u!=0){++unit_nonzero;if(ht.current().current_load_C_m_s[i]==0){++zeroed;lost+=std::fabs((long double)u*w);}}}
  mfem::ConstantCoefficient zero(0);std::vector<s4::Candidate>sources{tiny};std::vector<s4::Handoff>hands{ht};
  auto trial=o.prepare(o.make_epoch(1e-14,"trial-midpoint-field"),zero,field,view,sources,hands,{{}});
  std::cout<<"gap1 weight="<<w<<" birth_operands_nonzero="<<all_nonzero<<" handoff_current_entries_nonzero_at_unit_weight="<<unit_nonzero<<" zeroed_at_weight_w="<<zeroed<<" exact_magnitude_lost="<<double(lost)<<" trial_preflight="<<o.preflight(trial,view)<<"\n";
  std::cout<<(zeroed>0&&o.preflight(trial,view)?"REPRODUCED":"NOT-REPRODUCED")<<" gap1 admitted handoff current silently zeroed\n";
  // (2)
  const double a=3.1622776601683792e-245;auto p2=prep(a);const long double reference=(long double)p1.expected_number()*a;
  std::cout<<"gap2 area_fraction="<<a<<" admitted_scalar_count="<<p2.expected_number()<<" exp_log_count="<<double(std::exp((long double)p2.log_expected_number()))<<" area1_count_times_fraction="<<double(reference)<<" ratio="<<double(p2.expected_number()/reference)<<"\n";
  std::cout<<(std::fabs(p2.expected_number()/reference-1)>1e-9?"REPRODUCED":"NOT-REPRODUCED")<<" gap2 imprecise nonzero exposure scalar admitted\n";
  return 0;
 }catch(const std::exception&e){std::cerr<<"error: "<<e.what()<<"\n";return 1;}
}
