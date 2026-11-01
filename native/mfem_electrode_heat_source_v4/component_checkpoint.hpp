#pragma once
#include "source_owner.hpp"
namespace electrode_heat_v4 {
namespace detail {
struct Reader {
 const std::string&b;std::size_t pos=0;
 std::uint64_t u(){eh::need(pos<=b.size()&&b.size()-pos>=8,"v4 checkpoint truncated integer");std::uint64_t n=0;for(int i=0;i<8;++i)n|=std::uint64_t(static_cast<unsigned char>(b[pos++]))<<(8*i);return n;}
 double d(){return std::bit_cast<double>(u());}
 bool flag(){auto n=u();eh::need(n<=1,"v4 checkpoint bool encoding");return n;}
 std::string s(std::size_t cap=1000000){auto n=u();eh::need(n<=cap&&pos<=b.size()&&n<=b.size()-pos,"v4 checkpoint string cap/truncation");auto out=b.substr(pos,n);pos+=n;return out;}
};
struct TempFile {std::string path;~TempFile(){std::error_code e;std::filesystem::remove(path,e);}};
inline std::string read_file(const std::string&path,std::size_t cap){eh::need(std::filesystem::file_size(path)<=cap,"v4 component checkpoint file cap");std::ifstream f(path,std::ios::binary);eh::need(bool(f),"v4 checkpoint open");return {std::istreambuf_iterator<char>(f),{}};}
}
inline void Owner::save(const std::string&path)const {
 detail::TempFile temporary{path+".heat-stage"};heat_->save(temporary.path);auto h=detail::read_file(temporary.path,64*1024*1024);
 detail::Encoding e;e.s("heat-source-v4-component-1");e.s(signature_);e.s(h);for(double x:neutral_cumulative_)e.d(x);e.u(inventory_.size());for(const auto&n:inventory_)detail::count(e,n);e.u(layers_.size());
 for(const auto&[id,r]:layers_){e.u(id);detail::layer(e,r.state);detail::identity(e,r.identity);detail::policy(e,r.policy);e.d(r.slice_center_z_m);e.d(r.extrusion_height_m);e.u(r.qualification);e.s(r.energy_owner);}
 e.s(e.sha());eh::need(e.bytes.size()<=64*1024*1024,"v4 component checkpoint byte cap");std::ofstream f(path+".tmp",std::ios::binary);f.write(e.bytes.data(),e.bytes.size());f.close();eh::need(bool(f),"v4 component checkpoint write failed");std::filesystem::rename(path+".tmp",path);
}
inline void Owner::load(const std::string&path){
 eh::need(generation_==0&&heat_->state().revision==0,"v4 restart requires fresh component owner");auto bytes=detail::read_file(path,64*1024*1024);detail::Reader r{bytes};eh::need(r.s()=="heat-source-v4-component-1"&&r.s()==signature_,"v4 checkpoint source/material/area/layer/policy identity mismatch");auto heat_bytes=r.s(64*1024*1024);std::array<double,3>neutral_cumulative;for(double&x:neutral_cumulative){x=r.d();eh::need(std::isfinite(x)&&x>=0,"v4 neutral cumulative energy encoding");}
 auto n=r.u();eh::need(n==inventory_.size(),"v4 exact source inventory dimensions");std::vector<Count>inventory;inventory.reserve(n);
 for(std::size_t j=0;j<n;++j){auto count=r.u();eh::need(count<=Count::maximum_carriers,"v4 exact count carrier cap");std::vector<double>w;w.reserve(count);for(std::size_t k=0;k<count;++k)w.push_back(r.d());Count value=Count::from_positive_carriers(w);eh::need(value.positive_carriers()==w,"v4 noncanonical count checkpoint");inventory.push_back(std::move(value));}
 n=r.u();eh::need(n<=4096,"v4 Active source checkpoint cap");std::map<std::uint64_t,LayerRecord>layers;
 for(std::size_t j=0;j<n;++j){auto id=r.u();eh::need(id>0,"v4 canonical layer ID");auto layer_bytes=r.s();std::istringstream f(layer_bytes,std::ios::binary);LayerRecord rec;rec.state=eh::read_bound_layer_state(f);eh::need(f.peek()==std::char_traits<char>::eof(),"v4 trailing layer payload");rec.identity={r.s(),r.s(),r.d()};auto&p=rec.policy;p.local_linear_vacuum_extent_m=r.d();p.image_force_fraction=r.d();p.image_energy_fraction=r.d();p.maximum_quantum_magnetic_energy_fraction=r.d();p.maximum_gyroangle=r.d();p.maximum_relative_normal_step=r.d();auto steps=r.u();eh::need(steps>=1&&steps<=16384,"v4 layer substep bound");p.maximum_substeps=int(steps);p.PIC_already_applies_microscopic_single_electron_image_force=r.flag();p.own_return_unit_accommodation=r.flag();rec.slice_center_z_m=r.d();rec.extrusion_height_m=r.d();rec.qualification=r.flag();rec.energy_owner=r.s();eh::need(!rec.energy_owner.empty(),"v4 layer energy owner missing");v3::admit_policy(p);
  auto index=rec.state.patch_index;eh::need(index<patches_.size()&&rec.identity==identities_[index]&&detail::policy_hash(p)==detail::policy_hash(settings_[index].layer)&&std::isfinite(rec.slice_center_z_m)&&rec.extrusion_height_m==geometry_.height&&rec.state.native.status==se::LayerStatus::Active,"v4 layer patch/area/material/policy mismatch");
  double w=rec.state.native.particle.weight;detail::require_product(-se::qe,w,-se::qe*w,"restored source charge C");double K=se::kinetic_energy_J(rec.state.native.particle.proper_velocity_m_s);detail::require_product(K,w,K*w,"restored source kinetic J");eh::need(layers.emplace(id,std::move(rec)).second,"v4 duplicate canonical layer ID");
 }
 auto checksum_position=r.pos;auto checksum=r.s(64);eh::need(r.pos==bytes.size(),"v4 checkpoint trailing bytes");detail::Encoding checked;checked.bytes=bytes.substr(0,checksum_position);eh::need(checksum==checked.sha(),"v4 checkpoint checksum");
 detail::TempFile temporary{path+".private-heat-load"};{std::ofstream f(temporary.path,std::ios::binary);f.write(heat_bytes.data(),heat_bytes.size());f.close();eh::need(bool(f),"v4 private checkpoint staging failed");}
 auto heat=new_heat();heat->load(temporary.path);for(const auto&n:heat->state().inventory)eh::need(n.hi==0&&n.lo==0,"v4 heat engine cannot carry a shadow source inventory");
 std::set<std::string>energy_owners;for(const auto&[id,rec]:layers){eh::need(energy_owners.insert(rec.energy_owner).second,"v4 duplicate layer energy owner");eh::need(rec.state.native.current_time_s==heat->state().time&&rec.state.mesh_identity==geometry_.revision&&rec.state.material_hash==material_.data_hash,"v4 layer/heat accepted clock/identity mismatch");auto x=rec.state.surface_origin_m;x[2]-=rec.slice_center_z_m;eh::admit_origin(*heat,rec.state.patch_index,x);}
 heat_.swap(heat);inventory_.swap(inventory);layers_.swap(layers);neutral_cumulative_.swap(neutral_cumulative);++generation_;
}
}
