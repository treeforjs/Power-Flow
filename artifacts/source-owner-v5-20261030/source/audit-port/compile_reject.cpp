// Compile-reject evidence for the audit source probes that cannot be written
// against owner v5 (cases unchanged from the source-v4 port). Built with
// -fsyntax-only and -DREJECT_CASE=N using the v5 component flags and include paths.
//   case 0      control: the same objects used read-only; must compile.
//   cases 1..11 the direct v3-style mutation/forgery; each must be rejected,
//               with the first diagnostic located in this file.
#include "hall-owner/development/hall-finite-foil-em-20261029/axial_surface_path_port.hpp"
#include "source_owner.hpp"
#include "ss304_heat_material_v2.hpp"
#ifndef REJECT_CASE
#error "define REJECT_CASE"
#endif
namespace s5=electrode_heat_v5;
namespace v3=electrode_heat_v3;
namespace se=surface_emission;
void probe(s5::Owner&o,s5::Candidate&cand,s5::PreparedPatch&p,v3::Candidate&v3cand,v3::PreparedPatch&v3p,const std::shared_ptr<const s5::Field>&f,const s5::FinalFieldView&v,s5::Trial&t){
#if REJECT_CASE==0
 double w=cand.state().native.particle.weight;(void)w;auto copy=cand;(void)copy;const auto&key=p.exposure_key();(void)key;
 auto receipt=o.import(f,v,std::span<const s5::Candidate>(&cand,1),{42});(void)receipt;
 v3cand.state.native.particle.weight*=2;v3p.native.spectrum=se::evaluate(1600,1e8,se::Config{4.3});(void)t;
 double charge=cand.exchange().particle_charge_change_C;(void)charge;
#elif REJECT_CASE==1  // A02/C3 weight: mutated.state.native.particle.weight*=2
 auto mutated=cand;mutated.state().native.particle.weight*=2;
#elif REJECT_CASE==2  // A02/C3 weight through the private payload
 auto mutated=cand;mutated.d_->native.state.native.particle.weight*=2;
#elif REJECT_CASE==3  // A02/C3 weight: forge a v5 candidate from an edited v3 candidate
 v3cand.state.native.particle.weight*=2;s5::Candidate forged(v3cand);(void)forged;
#elif REJECT_CASE==4  // A02/C3 weight: import an edited v3 candidate through the v5 owner
 v3cand.state.native.particle.weight*=2;auto receipt=o.import(f,v,std::span<const v3::Candidate>(&v3cand,1),{42});(void)receipt;
#elif REJECT_CASE==5  // A02/C3 spectrum: altered.native.spectrum=se::evaluate(1600,...)
 auto altered=p;altered.native.spectrum=se::evaluate(1600,1e8,se::Config{4.3});
#elif REJECT_CASE==6  // A02/C3 spectrum through the private payload
 p.d_->native.native.spectrum=se::evaluate(1600,1e8,se::Config{4.3});
#elif REJECT_CASE==7  // A02/C3 spectrum: forge a v5 prepared patch from an edited v3 patch
 v3p.native.spectrum=se::evaluate(1600,1e8,se::Config{4.3});s5::PreparedPatch forged(v3p);(void)forged;
#elif REJECT_CASE==8  // A02/C3 mesh: mutate the sealed field's private mesh copy
 f->mesh_->dx*=1.001;
#elif REJECT_CASE==9  // A02/C3 owner address: relocate an owner (nonmovable)
 s5::Owner moved(std::move(o));(void)moved;
#elif REJECT_CASE==10 // duplicate an opaque trial proposal
 s5::Trial duplicate(t);(void)duplicate;
#elif REJECT_CASE==11 // A08/C4: zero a sealed weighted charge after admission
 cand.exchange().particle_charge_change_C=0;
#else
#error "unknown REJECT_CASE"
#endif
}
int main(){return 0;}
