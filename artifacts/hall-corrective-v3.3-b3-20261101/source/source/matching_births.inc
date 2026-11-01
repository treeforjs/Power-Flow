// Genuine native B3 evidence completes an ordinary private host motion. It
// neither changes the particle push nor fabricates an electrode assignment.
namespace hall_corrective_v3 {
void MotionAssembly::apply_matching_births(hall_driver::MotionReceipt&r,const MotionContext&c,
 std::span<const native::ArenaRateEntry>census,const native::PreparedBirths&prepared,const std::vector<double>&E){
 const auto&receipt=prepared.receipt();const auto&flights=prepared.flights();const auto&g=*c.mesh;const auto&o=c.options;
 hall::require(r.forms_v7_&&r.bridges_.empty()&&r.extra_current_A_m_.empty()&&r.extra_work_J_==0,"B3 requires fresh host native-v7 motion; cannot apply twice");
 hall::require(receipt.accepted_step==r.binding_.accepted_step&&flights.size()==r.births_.size(),"B3 authentic birth capability/common step mismatch");
 for(std::size_t k=0;k<flights.size();++k){const auto&a=flights[k];const auto&b=r.births_[k];hall::require(a.native_id==b.native_id&&a.birth_time_s==b.birth_time_s&&a.flight_s==b.flight_s,"B3 native versus Hall newborn flight mismatch");}
 if(receipt.bridges.empty())return;
 hall::require(E.size()==std::size_t(g.ne)&&axial::finite(E)&&receipt.wall_entries.empty(),"B3 field dimensions or unsupported wall births");
 for(const auto&b:receipt.entries)hall::require(!b.returned,"B3 returned-subcell births are not qualified");
 std::map<std::uint64_t,const native::ArenaRateEntry*>rows;for(const auto&x:census)hall::require(rows.emplace(x.native.id,&x).second,"B3 duplicate census ID");
 std::map<std::uint64_t,NewbornRecord*>born;for(auto&x:r.births_)hall::require(born.emplace(x.native_id,&x).second,"B3 duplicate newborn ID");
 const native::PicVector shift{-o.width_m/2,-o.gap_m/2,0};hall::require(r.frame_.native_to_field_m==shift,"B3 field translation differs from actual kinetic geometry");
 r.extra_current_A_m_.assign(g.ne,0.);std::uint64_t last=0;
 for(const auto&b:receipt.bridges){
  hall::require(b.native_id>last&&rows.contains(b.native_id)&&born.contains(b.native_id),"B3 bridge order or missing actual newborn census");last=b.native_id;
  const auto&row=*rows.at(b.native_id);const auto&p=row.native;auto&birth=*born.at(b.native_id);
  const auto entry=std::find_if(receipt.entries.begin(),receipt.entries.end(),[&](const auto&x){return x.native_id==b.native_id;});
  hall::require(entry!=receipt.entries.end()&&entry->packet_id==row.packet_id&&entry->state==row.state&&entry->weight==p.weight&&birth.packet_id==row.packet_id&&birth.state==row.state&&birth.start_m==p.position_m,"B3 native source receipt differs from actual Hall census");
  hall::require(birth.origin<0&&birth.bridge_origin<0&&p.charge_C==-native::pic_e&&p.weight>0&&std::isfinite(p.weight)&&b.charge_C==p.charge_C*p.weight,"B3 requires an interior matching-plane electron of the actual macro-charge");
  native::PicVector end=p.position_m;for(int k=0;k<3;++k)end[k]+=shift[k];
  int axis=-1;for(int k=0;k<3;++k){hall::require(std::isfinite(b.field_start_m[k])&&std::isfinite(b.field_displacement_m[k])&&b.field_displacement_m[k]==end[k]-b.field_start_m[k],"B3 bridge endpoint differs from actual newborn");if(b.field_displacement_m[k]!=0){hall::require(axis<0,"B3 bridge is not a normal surface displacement");axis=k;}}
  hall::require(axis==0||axis==1,"B3 bridge must leave a foil side, not a target end");
  const double sign=b.field_displacement_m[axis]>0?1.:-1.;int body=-1;
  for(int k=0;k<2;++k){const auto&foil=k==0?o.lower_foil:o.upper_foil;const double plane=axis==0?(sign>0?foil.right_m:foil.left_m):(sign>0?foil.top_m:foil.bottom_m);const int tangent=1-axis;
   const double lo=axis==0?foil.bottom_m:foil.left_m,hi=axis==0?foil.top_m:foil.right_m;
   if(b.field_start_m[axis]==plane+shift[axis]&&b.field_start_m[tangent]>=lo+shift[tangent]&&b.field_start_m[tangent]<=hi+shift[tangent]&&b.field_start_m[2]>o.short_z_m&&b.field_start_m[2]<o.feed_z_m){hall::require(body<0,"B3 ambiguous emitting foil");body=k;}}
  hall::require(body>=0&&p.proper_velocity_m_s[axis]*sign>0,"B3 start is not an actual emitting foil face with outward phase space");
  birth.bridge_origin=body;birth.bridge_start_field_m=b.field_start_m;
  axial::Path path(g,b.field_start_m,b.field_displacement_m);path.deposit(r.extra_current_A_m_,b.charge_C,b.field_displacement_m,r.binding_.dt_s);
  const double work=b.charge_C*hall::dot(path.gather_E(E),b.field_displacement_m);hall::require(std::isfinite(work),"B3 bridge work overflow");r.bridge_work_J_.push_back(work);
  // Initial charge of the combined source+ordinary path is at the conductor
  // face; ordinary trajectory and its endpoint remain at the matching plane.
  axial::point_deposit(g,r.initial_charge_C_,end,-b.charge_C);axial::point_deposit(g,r.initial_charge_C_,b.field_start_m,b.charge_C);
  r.bridges_.push_back(b);
 }
 for(std::size_t k=0;k<r.current_A_m_.size();++k)r.current_A_m_[k]+=r.extra_current_A_m_[k];
 const auto extra=pf_v7::field_work_J(r.binding_.dt_s,E,r.extra_current_A_m_),total=pf_v7::field_work_J(r.binding_.dt_s,E,r.current_A_m_);
 hall::require(extra.status==0&&total.status==0,"B3 exact current work operands");r.extra_work_J_=extra.value;r.budget_.field_work_J=total.value;r.budget_.work_defect_J=r.budget_.kinetic_change_J-(total.value-extra.value);
 std::vector<int>ea(g.ne),eb(g.ne);for(int k=0;k<g.ne;++k){ea[k]=g.edges[k].a;eb[k]=g.edges[k].b;}
 motion_work_v1::field_scales(r.work_,r.binding_.dt_s,E,r.current_A_m_,ea,eb,g.nn);r.work_.tolerance_J=motion_work_v1::tolerance(r.work_);r.work_.v31_tolerance_J=1e-20+1e-10*(std::abs(r.budget_.kinetic_change_J)+std::abs(total.value));
 hall::require(motion_work_v1::admits(r.work_,r.budget_.work_defect_J),"B3 ordinary-motion work bound");
 foil::Vec gross(g.nn);for(const auto&row:census){auto x=row.native.position_m;for(int k=0;k<3;++k)x[k]+=shift[k];auto i=born.find(row.native.id);if(i!=born.end()&&i->second->bridge_origin>=0)x=i->second->bridge_start_field_m;axial::point_deposit(g,gross,x,std::abs(row.native.charge_C*row.native.weight));}
 auto div=g.gradient_transpose(r.current_A_m_);r.budget_.continuity_C=0;for(int k=0;k<g.nn;++k)r.budget_.continuity_C=std::max(r.budget_.continuity_C,std::abs(r.endpoint_charge_C_[k]-r.initial_charge_C_[k]-r.binding_.dt_s*div[k]));
 hall::require(axial::finite(r.current_A_m_)&&std::isfinite(r.budget_.continuity_C)&&r.budget_.continuity_C<=1e-12*std::max(axial::maximum(gross),1e-30),"B3 combined source/ordinary charge continuity mismatch");
 r.record_v7_sha256_.clear();r.record_v3_sha256_.clear();r.seal_sha256_.clear();seal_v7(r);
}
}
