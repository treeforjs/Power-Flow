// Exact field serialization, no C++ object padding. Full source state and gas
// records are part of the same digest/checkpoint as the native thermal owner.
struct Writer {
 std::ostringstream out{std::ios::out|std::ios::binary};
 template<class T>void put(const T& x){out.write(reinterpret_cast<const char*>(&x),sizeof(x));}
 void str(const std::string& s){put(std::uint64_t(s.size()));out.write(s.data(),s.size());}
 template<class T>void array(const T& v){for(const auto& x:v)put(x);}
};
struct Reader {
 std::istringstream in;
 explicit Reader(const std::string& s):in(s,std::ios::in|std::ios::binary){}
 template<class T>void get(T& x){in.read(reinterpret_cast<char*>(&x),sizeof(x));need(bool(in),"joint checkpoint truncated");}
 std::uint64_t count(std::uint64_t bound=10000000){std::uint64_t n=0;get(n);need(n<=bound,"joint checkpoint count exceeds bound");return n;}
 std::string str(){auto n=count(1000000000);std::string s(n,'\0');in.read(s.data(),n);need(bool(in),"joint checkpoint string truncated");return s;}
 template<class T>void array(T& v){for(auto& x:v)get(x);}
};
std::string payload(Holder& h,int which){
 const auto& a=auxiliary(h,which);const auto& s=state(h,which);Writer f;
 f.put(std::uint64_t(0x4a4f494e54563034ULL));f.str(h.electron_state);f.str(h.electron_parent);f.str(eh::source_kernel_sha256);f.str(h.owner->serialize_state(&s));
 f.put(a.source_sequence);f.array(a.cumulative_charge_number);f.put(std::uint64_t(a.gas_inventory.size()));for(auto n:a.gas_inventory){f.put(n.hi);f.put(n.lo);}
 f.put(std::uint64_t(a.layers.size()));for(const auto& [id,layer]:a.layers){f.put(id);eh::write_bound_layer_state(f.out,layer);}
 f.put(std::uint64_t(a.paths.size()));for(const auto& path:a.paths){f.put(path.layer_id);f.array(path.native.from_m);f.array(path.native.to_m);for(double x:{path.native.start_time_s,path.native.end_time_s,path.native.charge_C,path.native.field_work_J})f.put(x);}
 f.put(std::uint64_t(a.gas.size()));for(const auto& p:a.gas){for(const auto* s:{&p.patch,&p.state_id,&p.parent_id,&p.event_id,&p.channel_id})f.str(*s);f.put(p.source_revision);for(double x:{p.weight,p.T,p.time})f.put(x);f.array(p.x);f.array(p.v);f.array(p.proper_velocity_m_s);f.array(p.wall_support_impulse_N_s);for(double x:{p.mass_kg,p.charge_C,p.kinetic_J,p.classical_sample_K_J,p.relativistic_conversion_J,p.internal_J,p.reference_J,p.adsorbed_reference_change_J})f.put(x);}
 need(bool(f.out),"joint state serialization failed");return f.out.str();
}
void restore(Holder& h,const std::string& bytes){
 Reader f(bytes);std::uint64_t magic=0;f.get(magic);need(magic==0x4a4f494e54563034ULL,"joint v4 checkpoint magic");need(f.str()==h.electron_state&&f.str()==h.electron_parent&&f.str()==eh::source_kernel_sha256,"joint source species/kernel identity changed");auto thermal=f.str();Aux candidate;f.get(candidate.source_sequence);f.array(candidate.cumulative_charge_number);for(double x:candidate.cumulative_charge_number)need(eh::finite(x),"nonfinite source cumulative");auto n=f.count(100000);need(n==h.patch_count,"gas inventory patch dimension");candidate.gas_inventory.resize(n);for(auto& x:candidate.gas_inventory){f.get(x.hi);f.get(x.lo);need(eh::canonical(x),"invalid canonical gas inventory");}
 n=f.count();for(std::uint64_t i=0;i<n;++i){std::uint64_t id=0;f.get(id);auto layer=eh::read_bound_layer_state(f.in);need(id==layer.native.particle.event_id&&!candidate.layers.count(id)&&layer.patch_index<h.patch_count&&layer.mesh_identity==h.mesh&&layer.material_hash==h.material,"checkpoint layer identity");eh::admit_origin(*h.owner,layer.patch_index,layer.surface_origin_m);need(layer.native.outward_normal==h.owner->normal(layer.patch_index),"checkpoint layer normal differs from physical patch");candidate.layers.emplace(id,std::move(layer));}
 n=f.count();candidate.paths.resize(n);for(auto& path:candidate.paths){f.get(path.layer_id);f.array(path.native.from_m);f.array(path.native.to_m);for(double* x:{&path.native.start_time_s,&path.native.end_time_s,&path.native.charge_C,&path.native.field_work_J})f.get(*x);need(candidate.layers.count(path.layer_id)&&path.native.start_time_s>=0&&path.native.end_time_s>=path.native.start_time_s,"checkpoint current path identity/clock");for(double x:path.native.from_m)need(eh::finite(x),"nonfinite path origin");for(double x:path.native.to_m)need(eh::finite(x),"nonfinite path endpoint");for(double x:{path.native.start_time_s,path.native.end_time_s,path.native.charge_C,path.native.field_work_J})need(eh::finite(x),"nonfinite current/work path");}
 n=f.count();candidate.gas.resize(n);std::vector<eh::Count> reconstructed(h.patch_count);for(auto& p:candidate.gas){for(auto* s:{&p.patch,&p.state_id,&p.parent_id,&p.event_id,&p.channel_id})*s=f.str();f.get(p.source_revision);for(double* x:{&p.weight,&p.T,&p.time})f.get(*x);f.array(p.x);f.array(p.v);f.array(p.proper_velocity_m_s);f.array(p.wall_support_impulse_N_s);for(double* x:{&p.mass_kg,&p.charge_C,&p.kinetic_J,&p.classical_sample_K_J,&p.relativistic_conversion_J,&p.internal_J,&p.reference_J,&p.adsorbed_reference_change_J})f.get(*x);auto index=patch_index(h,p.patch);const auto& patch=h.owner->patch(index);need(patch.desorption_enabled&&p.state_id==patch.state_id&&p.parent_id==patch.parent_id&&p.mass_kg==patch.mass_kg&&p.charge_C==0&&p.reference_J==0&&p.weight>0&&p.kinetic_J>=0&&p.internal_J>=0&&!p.event_id.empty()&&p.channel_id==(patch.coverage_temkin?"conditional-oxide-water-Temkin-frozen-coverage":"authored-first-order-desorption"),"gas checkpoint state/reference mismatch");for(double x:{p.weight,p.T,p.time,p.mass_kg,p.kinetic_J,p.classical_sample_K_J,p.relativistic_conversion_J,p.internal_J,p.adsorbed_reference_change_J})need(eh::finite(x),"nonfinite gas checkpoint");for(const auto* a:{&p.x,&p.v,&p.proper_velocity_m_s,&p.wall_support_impulse_N_s})for(double x:*a)need(eh::finite(x),"nonfinite gas kinematics");double u2=0,v2=0;for(int axis=0;axis<3;++axis){u2+=p.proper_velocity_m_s[axis]*p.proper_velocity_m_s[axis];v2+=p.v[axis]*p.v[axis];const double expected=-p.weight*p.mass_kg*p.proper_velocity_m_s[axis];need(std::abs(p.wall_support_impulse_N_s[axis]-expected)<=256*eh::eps*(std::abs(expected)+std::abs(p.wall_support_impulse_N_s[axis])),"gas support impulse not conjugate to proper momentum");}
 const double expected_K=p.mass_kg*u2/(std::sqrt(1+u2/(299792458.*299792458.))+1)*p.weight;
 need(eh::finite(expected_K)&&std::abs(p.kinetic_J-expected_K)<=256*eh::eps*(std::abs(p.kinetic_J)+std::abs(expected_K)),"gas kinetic energy differs from actual proper velocity");
 need(v2<299792458.*299792458.,"checkpoint gas physical velocity not subluminal");const double gamma=std::sqrt(1+u2/(299792458.*299792458.));for(int axis=0;axis<3;++axis){const double velocity=p.proper_velocity_m_s[axis]/gamma;need(std::abs(velocity-p.v[axis])<=256*eh::eps*(std::abs(velocity)+std::abs(p.v[axis])),"gas velocity/proper velocity mismatch");}
 const double classical=.5*p.mass_kg*v2*p.weight;need(std::abs(p.classical_sample_K_J-classical)<=256*eh::eps*(std::abs(classical)+std::abs(p.classical_sample_K_J))&&std::abs(p.relativistic_conversion_J-(p.kinetic_J-p.classical_sample_K_J))<=256*eh::eps*(std::abs(p.kinetic_J)+std::abs(p.classical_sample_K_J)),"gas conversion diagnostic mismatch");
 reconstructed[index]=eh::add(reconstructed[index],p.weight);}
 for(size_t i=0;i<h.patch_count;++i)need(reconstructed[i].hi==candidate.gas_inventory[i].hi&&reconstructed[i].lo==candidate.gas_inventory[i].lo,"gas packet/canonical inventory mismatch");need(f.in.peek()==std::char_traits<char>::eof(),"joint checkpoint trailing bytes");
 Reader clocks(thermal);std::uint64_t thermal_magic=0,revision=0;double time=0;clocks.get(thermal_magic);clocks.str();clocks.get(time);clocks.get(revision);need(eh::finite(time),"checkpoint thermal clock nonfinite");
 for(const auto& [id,l]:candidate.layers)need(l.native.current_time_s<=time&&(l.native.status!=se::LayerStatus::Active||l.native.current_time_s==time),"checkpoint active-layer/common clock mismatch");
 for(const auto& p:candidate.paths)need(p.native.end_time_s<=time,"checkpoint path after accepted clock");
 for(const auto& p:candidate.gas)need(p.time<=time&&p.source_revision<=revision,"checkpoint gas after accepted state");
 double active_number=0,number_activity=0;for(const auto& [id,l]:candidate.layers){number_activity+=l.native.particle.weight;if(l.native.status!=se::LayerStatus::Recollected)active_number+=l.native.particle.weight;}
 const std::array<double,3> expected_charge{-se::qe*active_number,se::qe*active_number,active_number};for(int i=0;i<3;++i){const double activity=i<2?se::qe*number_activity:number_activity;need(std::abs(candidate.cumulative_charge_number[i]-expected_charge[i])<=256*eh::eps*(std::abs(candidate.cumulative_charge_number[i])+std::abs(expected_charge[i])+activity),"source charge/number cumulative differs from retained electrons");}
 // Read cumulative Wfield from the serialized thermal state without accepting it.
 auto nt=clocks.count();for(std::uint64_t i=0;i<nt;++i){double value;clocks.get(value);}auto ni=clocks.count();for(std::uint64_t i=0;i<2*ni;++i){double value;clocks.get(value);}std::array<double,19> cumulative{};clocks.array(cumulative);
 double path_work=0,path_activity=0;for(const auto& p:candidate.paths){path_work+=p.native.field_work_J;path_activity+=std::abs(p.native.field_work_J);}need(std::abs(path_work-cumulative[18])<=(1e-10+256*eh::eps)*(path_activity+std::abs(cumulative[18])),"checkpoint source field work differs from weighted paths");
 // restore_state validates all thermal bytes before its noexcept state swap;
 // candidate source/gas validation and allocations above precede that swap.
 h.owner->restore_state(thermal);using std::swap;swap(h.aux,candidate);
}
