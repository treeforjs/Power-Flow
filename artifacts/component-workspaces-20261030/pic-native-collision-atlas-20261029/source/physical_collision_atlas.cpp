#include "physical_collision_atlas.hpp"
#include "raikou/json.hpp"
#include <openssl/evp.h>
#include <cstring>
#include <iomanip>
#include <set>

namespace raikou::collision_geometry {
namespace {
using Json=nlohmann::json;
class Hash {
 std::unique_ptr<EVP_MD_CTX,decltype(&EVP_MD_CTX_free)> ctx_{EVP_MD_CTX_new(),EVP_MD_CTX_free};
public:
 Hash(){pc::require(bool(ctx_)&&EVP_DigestInit_ex(ctx_.get(),EVP_sha256(),nullptr)==1,"collision SHA initialization");}
 void bytes(const void*p,std::size_t n){pc::require(EVP_DigestUpdate(ctx_.get(),p,n)==1,"collision SHA update");}
 void integer(std::uint64_t x){unsigned char b[8];for(int i=0;i<8;++i)b[7-i]=static_cast<unsigned char>(x>>(8*i));bytes(b,8);}
 void number(double x){std::uint64_t b;std::memcpy(&b,&x,8);integer(b);}
 void text(std::string_view x){integer(x.size());bytes(x.data(),x.size());}
 std::string finish(){unsigned char out[EVP_MAX_MD_SIZE];unsigned n=0;pc::require(EVP_DigestFinal_ex(ctx_.get(),out,&n)==1&&n==32,"collision SHA finalization");std::ostringstream s;for(unsigned i=0;i<n;++i)s<<std::hex<<std::setw(2)<<std::setfill('0')<<unsigned(out[i]);return s.str();}
};
std::string sha(const std::string&s){Hash h;h.bytes(s.data(),s.size());return h.finish();}
double midpoint(double lo,double hi){const double m=lo+(hi-lo)*.5;pc::require(std::isfinite(m)&&lo<m&&m<hi,"unresolved geometric midpoint");return m;}
bool overlap(double a0,double a1,double b0,double b1){return std::min(a1,b1)>std::max(a0,b0);}
double positive_volume(long double v){const double out=static_cast<double>(v);pc::require(std::isfinite(v)&&std::isfinite(out)&&out>0,"unresolved physical collision volume");return out;}
}
FrameBinding::FrameBinding(std::string canonical):canonical_(std::move(canonical)){
 pc::require(pc::digest(canonical_),"unbound canonical frame");Hash h;h.text("canonical-only-collision-frame/v1");h.text(canonical_);identity_=h.finish();
}
FrameBinding::FrameBinding(const Geometry&g,std::string canonical,std::string native):canonical_(std::move(canonical)),native_(std::move(native)){
 pc::require(pc::digest(canonical_)&&pc::digest(native_),"unbound canonical/native frame");
 for(int d=0;d<2;++d){pc::require(std::isfinite(g.low_m[d])&&std::isfinite(g.high_m[d])&&g.low_m[d]<g.high_m[d],"frame geometry");origin_[d]=g.low_m[d]+(g.high_m[d]-g.low_m[d])*.5;pc::require(std::isfinite(origin_[d]),"unrepresentable centered XY frame");}
 Hash h;h.text("canonical-to-native-centered-xy/v1");h.text(canonical_);h.text(native_);for(double x:origin_)h.number(x);h.text("subtract x/y once; z unchanged; distinct point types");identity_=h.finish();
}
NativeCenteredXYPointSI FrameBinding::native_field_point(CanonicalPointSI p)const{
 pc::require(has_native_centered_xy_map(),"canonical-only frame has no raw native-field transform");
 for(double x:p.value)pc::require(std::isfinite(x),"nonfinite canonical point");
 return {{p.value[0]-origin_[0],p.value[1]-origin_[1],p.value[2]}};
}
CollisionAtlas::CollisionAtlas(Geometry g,FrameBinding f,AtlasLimits limits):geometry_(std::move(g)),frame_(std::move(f)),limits_(limits){
 pc::require(limits_.geometric_bins>0&&limits_.geometric_bins<=1000000&&limits_.xy_nodes>0&&limits_.xy_nodes<=200000&&limits_.xy_leaves>0&&limits_.xy_leaves<=100000,"collision geometry resource caps");
 for(int d=0;d<3;++d)pc::require(std::isfinite(geometry_.low_m[d])&&std::isfinite(geometry_.high_m[d])&&geometry_.low_m[d]<geometry_.high_m[d],"physical collision extent");
 pc::require(std::isfinite(geometry_.physical_depth_m)&&geometry_.physical_depth_m==geometry_.high_m[2]-geometry_.low_m[2],"physical collision depth differs from axial height");
 if(frame_.has_native_centered_xy_map())for(int d=0;d<2;++d)pc::require(frame_.native_xy_origin_m()[d]==geometry_.low_m[d]+(geometry_.high_m[d]-geometry_.low_m[d])*.5,"foreign frame binding geometry");
 pc::require(geometry_.foils.size()<=64,"solid rectangle cap");
 for(const auto&s:geometry_.foils)pc::require(std::isfinite(s.left_m)&&std::isfinite(s.right_m)&&std::isfinite(s.bottom_m)&&std::isfinite(s.top_m)&&geometry_.low_m[0]<s.left_m&&s.left_m<s.right_m&&s.right_m<geometry_.high_m[0]&&geometry_.low_m[1]<s.bottom_m&&s.bottom_m<s.top_m&&s.top_m<geometry_.high_m[1],"invalid finite foil geometry");
 for(std::size_t i=0;i<geometry_.foils.size();++i)for(std::size_t j=0;j<i;++j){const auto&a=geometry_.foils[i];const auto&b=geometry_.foils[j];pc::require(!overlap(a.left_m,a.right_m,b.left_m,b.right_m)||!overlap(a.bottom_m,a.top_m,b.bottom_m,b.top_m),"overlapping solid rectangles");}
}
CollisionAtlas CollisionAtlas::legacy_uniform(Geometry g,FrameBinding f,std::array<std::uint32_t,3>n,AtlasLimits limits){
 CollisionAtlas a(std::move(g),std::move(f),limits);a.uniform_=true;a.counts_=n;std::size_t total=1;long double volume=1;
 for(int d=0;d<3;++d){pc::require(n[d]>0&&total<=limits.geometric_bins/n[d],"collision grid cell cap");total*=n[d];volume*=(static_cast<long double>(a.geometry_.high_m[d])-a.geometry_.low_m[d])/n[d];}
 a.geometric_bins_=total;a.vacuum_bins_=total;a.uniform_volume_=positive_volume(volume);
 auto aligned=[&](double x,int d){const long double q=(static_cast<long double>(x)-a.geometry_.low_m[d])*n[d]/(static_cast<long double>(a.geometry_.high_m[d])-a.geometry_.low_m[d]);const long double nearest=std::round(q);pc::require(nearest>=0&&nearest<=n[d]&&std::abs(q-nearest)<=64*std::numeric_limits<double>::epsilon()*std::max(1.L,std::abs(q)),"collision grid cuts a foil; exact vacuum clipping is not implemented");return static_cast<std::uint64_t>(nearest);};
 for(const auto&s:a.geometry_.foils){a.solid_indices_.push_back({aligned(s.left_m,0),aligned(s.right_m,0),aligned(s.bottom_m,1),aligned(s.top_m,1)});const auto&r=a.solid_indices_.back();pc::require(r[0]<r[1]&&r[2]<r[3],"unresolved aligned solid width");a.vacuum_bins_-=(r[1]-r[0])*(r[3]-r[2])*n[2];}
 a.finish_identity(nullptr);return a;
}
std::uint32_t CollisionAtlas::split(Box b,const GradedOptions&o,unsigned depth){
 pc::require(nodes_.size()<limits_.xy_nodes&&depth<100,"collision transverse tree cap");const auto node=static_cast<std::uint32_t>(nodes_.size());nodes_.push_back({});double tx=o.bulk_max_m,ty=o.bulk_max_m;
 for(const auto&s:geometry_.foils){const bool xo=overlap(b.x0,b.x1,s.left_m,s.right_m),yo=overlap(b.y0,b.y1,s.bottom_m,s.top_m);const double dx=std::max({s.left_m-b.x1,b.x0-s.right_m,0.}),dy=std::max({s.bottom_m-b.y1,b.y0-s.top_m,0.});
  if(xo&&(b.y1<=s.bottom_m||b.y0>=s.top_m)){ty=std::min(ty,o.normal_target_m+o.distance_growth*dy);tx=std::min(tx,o.broad_tangent_max_m);}
  else if(yo&&(b.x1<=s.left_m||b.x0>=s.right_m)){tx=std::min(tx,o.normal_target_m+o.distance_growth*dx);ty=std::min(ty,o.side_tangent_max_m);}
  else{const double t=o.normal_target_m+o.distance_growth*std::hypot(dx,dy);tx=std::min(tx,t);ty=std::min(ty,t);}
 }
 const double rx=(b.x1-b.x0)/tx,ry=(b.y1-b.y0)/ty;
 if(std::max(rx,ry)<=1){pc::require(leaves_.size()<limits_.xy_leaves,"collision transverse leaf cap");nodes_[node].leaf=static_cast<std::uint32_t>(leaves_.size());leaves_.push_back({b,0,0});return node;}
 const int axis=rx>=ry?0:1;const double cut=axis==0?midpoint(b.x0,b.x1):midpoint(b.y0,b.y1);Box left=b,right=b;if(axis==0){left.x1=cut;right.x0=cut;}else{left.y1=cut;right.y0=cut;}const auto l=split(left,o,depth+1),r=split(right,o,depth+1);nodes_[node].axis=axis;nodes_[node].cut=cut;nodes_[node].left=l;nodes_[node].right=r;return node;
}
CollisionAtlas CollisionAtlas::graded(Geometry g,FrameBinding f,GradedOptions o,AtlasLimits limits){
 CollisionAtlas a(std::move(g),std::move(f),limits);for(double x:{o.normal_target_m,o.distance_growth,o.broad_tangent_max_m,o.side_tangent_max_m,o.bulk_max_m})pc::positive(x,"positive finite graded geometry option");pc::require(o.axial_bins>0&&o.axial_bins<=limits.geometric_bins,"axial bin cap");pc::require((o.near_axial_bins==0&&o.axial_refinement_band_m==0)||(o.near_axial_bins>0&&o.near_axial_bins<=limits.geometric_bins&&std::isfinite(o.axial_refinement_band_m)&&o.axial_refinement_band_m>0),"complete local axial family required");
 std::vector<double>xs{a.geometry_.low_m[0],a.geometry_.high_m[0]},ys{a.geometry_.low_m[1],a.geometry_.high_m[1]};for(const auto&s:a.geometry_.foils){xs.push_back(s.left_m);xs.push_back(s.right_m);ys.push_back(s.bottom_m);ys.push_back(s.top_m);}for(auto*p:{&xs,&ys}){std::sort(p->begin(),p->end());p->erase(std::unique(p->begin(),p->end()),p->end());}
 for(std::size_t j=0;j+1<ys.size();++j)for(std::size_t i=0;i+1<xs.size();++i){Box b{xs[i],xs[i+1],ys[j],ys[j+1]};bool solid=false;for(const auto&s:a.geometry_.foils)solid|=b.x0>=s.left_m&&b.x1<=s.right_m&&b.y0>=s.bottom_m&&b.y1<=s.top_m;if(!solid)a.roots_.push_back({b,a.split(b,o,0)});}
 std::vector<std::uint32_t>zcounts;zcounts.reserve(a.leaves_.size());
 for(auto&leaf:a.leaves_){double distance=std::numeric_limits<double>::infinity();for(const auto&s:a.geometry_.foils)distance=std::min(distance,std::hypot(std::max({s.left_m-leaf.box.x1,leaf.box.x0-s.right_m,0.}),std::max({s.bottom_m-leaf.box.y1,leaf.box.y0-s.top_m,0.})));const auto n=o.near_axial_bins&&distance<=o.axial_refinement_band_m?o.near_axial_bins:o.axial_bins;pc::require(n<=limits.geometric_bins-a.geometric_bins_,"3D collision bin cap; no partial grid or automatic coarsening");leaf.offset=a.geometric_bins_;a.geometric_bins_+=n;zcounts.push_back(n);}
 a.vacuum_bins_=a.geometric_bins_;std::vector<std::uint32_t>families=zcounts;std::sort(families.begin(),families.end());families.erase(std::unique(families.begin(),families.end()),families.end());
 for(auto n:families){std::vector<double>axis;axis.reserve(static_cast<std::size_t>(n)+1);for(std::uint32_t k=0;k<=n;++k){const double x=k==0?a.geometry_.low_m[2]:k==n?a.geometry_.high_m[2]:static_cast<double>(static_cast<long double>(a.geometry_.low_m[2])+(static_cast<long double>(a.geometry_.high_m[2])-a.geometry_.low_m[2])*k/n);pc::require(std::isfinite(x)&&(axis.empty()||axis.back()<x),"unresolved axial bin");axis.push_back(x);}a.z_families_.push_back(std::move(axis));}
 for(std::size_t i=0;i<a.leaves_.size();++i)a.leaves_[i].z_family=static_cast<std::uint32_t>(std::lower_bound(families.begin(),families.end(),zcounts[i])-families.begin());
 a.finish_identity(&o);return a;
}
void CollisionAtlas::require_point(CanonicalPointSI p)const{
 for(int d=0;d<3;++d)pc::require(std::isfinite(p.value[d])&&p.value[d]>=geometry_.low_m[d]&&p.value[d]<=geometry_.high_m[d],"carrier outside physical collision domain");
 for(const auto&s:geometry_.foils){const bool vertical=(p.value[0]==s.left_m||p.value[0]==s.right_m)&&p.value[1]>=s.bottom_m&&p.value[1]<=s.top_m;const bool horizontal=(p.value[1]==s.bottom_m||p.value[1]==s.top_m)&&p.value[0]>=s.left_m&&p.value[0]<=s.right_m;pc::require(!vertical&&!horizontal,"exact material surface requires explicit common-stage side ownership");pc::require(!(p.value[0]>s.left_m&&p.value[0]<s.right_m&&p.value[1]>s.bottom_m&&p.value[1]<s.top_m),"carrier inside solid");}
}
std::uint64_t CollisionAtlas::locate(CanonicalPointSI p)const{
 require_point(p);
 if(uniform_){std::array<std::uint64_t,3>i{};for(int d=0;d<3;++d){const long double f=(static_cast<long double>(p.value[d])-geometry_.low_m[d])/(static_cast<long double>(geometry_.high_m[d])-geometry_.low_m[d]);i[d]=std::min<std::uint64_t>(counts_[d]-1,static_cast<std::uint64_t>(f*counts_[d]));}const auto id=i[0]+counts_[0]*(i[1]+counts_[1]*i[2]);physical_volume(id);return id;}
 const Root*root=nullptr;for(const auto&r:roots_){const auto&b=r.box;if(p.value[0]>=b.x0&&(p.value[0]<b.x1||b.x1==geometry_.high_m[0])&&p.value[1]>=b.y0&&(p.value[1]<b.y1||b.y1==geometry_.high_m[1])){root=&r;break;}}pc::require(root,"collision vacuum partition hole");std::uint32_t idx=root->node;while(nodes_[idx].axis>=0){const auto&n=nodes_[idx];idx=p.value[n.axis]<n.cut?n.left:n.right;}const auto&leaf=leaves_[nodes_[idx].leaf];const auto&z=z_families_[leaf.z_family];const auto k=std::min<std::size_t>(z.size()-2,static_cast<std::size_t>(std::upper_bound(z.begin(),z.end(),p.value[2])-z.begin()-1));return leaf.offset+k;
}
double CollisionAtlas::physical_volume(std::uint64_t id)const{
 pc::require(id<geometric_bins_,"unknown collision bin identity");
 if(uniform_){const auto i=id%counts_[0],j=(id/counts_[0])%counts_[1];for(const auto&r:solid_indices_)pc::require(!(i>=r[0]&&i<r[1]&&j>=r[2]&&j<r[3]),"solid collision bin identity");return uniform_volume_;}
 const auto next=std::upper_bound(leaves_.begin(),leaves_.end(),id,[](std::uint64_t value,const Leaf&leaf){return value<leaf.offset;});pc::require(next!=leaves_.begin(),"collision leaf index");const auto&leaf=*std::prev(next);const auto&z=z_families_[leaf.z_family];const auto k=id-leaf.offset;pc::require(k+1<z.size(),"collision axial index");return positive_volume((static_cast<long double>(leaf.box.x1)-leaf.box.x0)*(static_cast<long double>(leaf.box.y1)-leaf.box.y0)*(static_cast<long double>(z[k+1])-z[k]));
}
void CollisionAtlas::finish_identity(const GradedOptions*o){
 if(uniform_){legacy_identity_=sha(Json{{"schema","physical-cartesian-collision-grid-v1"},{"low",geometry_.low_m},{"high",geometry_.high_m},{"counts",counts_},{"solid_xy_index_ranges",solid_indices_},{"volume_m3",uniform_volume_}}.dump());}
 Hash h;h.text("native-physical-collision-atlas/v1");h.text(frame_.identity());h.text(uniform_?"legacy-uniform-exact-root-arithmetic":"graded-rectangle-binary-split");for(double x:geometry_.low_m)h.number(x);for(double x:geometry_.high_m)h.number(x);h.number(geometry_.physical_depth_m);h.integer(geometry_.foils.size());for(const auto&s:geometry_.foils)for(double x:{s.left_m,s.right_m,s.bottom_m,s.top_m})h.number(x);h.integer(geometric_bins_);h.integer(vacuum_bins_);h.text("exact material faces refused; upper external faces included in last geometric bin; common-stage operator admission separate");
 if(uniform_)h.text(legacy_identity_);else{for(double x:{o->normal_target_m,o->distance_growth,o->broad_tangent_max_m,o->side_tangent_max_m,o->bulk_max_m,o->axial_refinement_band_m})h.number(x);h.integer(o->axial_bins);h.integer(o->near_axial_bins);h.integer(roots_.size());for(const auto&r:roots_){for(double x:{r.box.x0,r.box.x1,r.box.y0,r.box.y1})h.number(x);h.integer(r.node);}h.integer(nodes_.size());for(const auto&n:nodes_){h.integer(n.left);h.integer(n.right);h.integer(n.leaf);h.integer(static_cast<std::uint64_t>(n.axis+1));h.number(n.cut);}h.integer(leaves_.size());for(const auto&leaf:leaves_){for(double x:{leaf.box.x0,leaf.box.x1,leaf.box.y0,leaf.box.y1})h.number(x);h.integer(leaf.z_family);h.integer(leaf.offset);}h.integer(z_families_.size());for(const auto&z:z_families_){h.integer(z.size());for(double x:z)h.number(x);}h.text("volume: long-double product of actual x/y/z bound differences then checked double conversion");}
 identity_=h.finish();
}
OccupiedCellView CollisionAtlas::occupied_impl(const std::vector<CarrierView>&entries,const CanonicalCensusBinding&binding,CensusLimits limits)const{
 pc::require(limits.carriers<=1000000&&limits.occupied_cells<=1000000&&limits.string_bytes<=256*1024*1024&&entries.size()<=limits.carriers&&binding.expected_native_ids.size()==entries.size()&&pc::digest(binding.accepted_state_prefix),"unbound/incomplete collision census or resource cap");
 auto expected=binding.expected_native_ids;std::sort(expected.begin(),expected.end());pc::require(std::adjacent_find(expected.begin(),expected.end())==expected.end()&&(expected.empty()||expected.front()>0),"invalid canonical expected native identities");
 OccupiedCellView out;out.accepted_state_prefix=binding.accepted_state_prefix;out.geometry_identity=identity_;out.entry_cell_ids.reserve(entries.size());std::vector<std::uint64_t>ids;ids.reserve(entries.size());Hash hash;hash.text("complete-positive-carrier-collision-census/v1");hash.text(identity_);hash.text(binding.accepted_state_prefix);hash.integer(entries.size());
 for(const auto&e:entries){pc::require(e.id>0&&!e.packet.empty()&&!e.state.empty(),"missing collision carrier identity");pc::positive(e.weight,"nonpositive native rate-carrier weight");pc::positive(e.mass,"native carrier mass");pc::require(std::isfinite(e.charge)&&std::isfinite(e.charge*e.weight),"nonfinite native carrier charge");for(double x:*e.proper_velocity)pc::require(std::isfinite(x),"nonfinite native proper velocity");const auto cell=locate({*e.position});ids.push_back(e.id);out.entry_cell_ids.push_back(cell);hash.integer(e.id);hash.text(e.packet);hash.text(e.state);for(double x:*e.position)hash.number(x);for(double x:*e.proper_velocity)hash.number(x);for(double x:{e.weight,e.mass,e.charge})hash.number(x);hash.integer(cell);}
 std::sort(ids.begin(),ids.end());pc::require(ids==expected,"collision census omits/duplicates/changes a canonical positive native carrier");auto occupied=out.entry_cell_ids;std::sort(occupied.begin(),occupied.end());occupied.erase(std::unique(occupied.begin(),occupied.end()),occupied.end());pc::require(occupied.size()<=limits.occupied_cells,"occupied collision cell cap");out.cells.reserve(occupied.size());for(auto id:occupied)out.cells.push_back({id,physical_volume(id)});out.census_identity=hash.finish();return out;
}
} // namespace raikou::collision_geometry
