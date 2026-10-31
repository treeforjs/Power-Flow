#pragma once
// Source owner v6.1: immutable mesh identity.
//
// v6 recomputed the whole-mesh digest (detail::mesh_hash, about 9 MB of encoding on the
// 200 um Hall mesh, about 10 ms) in every Field constructor and in every Field::matches
// call. v6.1 takes one immutable snapshot per mesh object and hashes it once:
//
//   MeshIdentity::snapshot(g)  private const copy of g (never mutated; mutation is not
//                              expressible), its v6 digest computed once, a process-unique
//                              generation number.
//   MeshIdentity::of(g)        the identity last taken of the object at &g, reused only
//                              when g is still exactly equal to that snapshot (every
//                              hashed member, doubles compared bit for bit); otherwise a
//                              new snapshot (new digest, new generation).
//   same_as(g)                 true when g IS the snapshot object (pointer identity of an
//                              immutable object), else the exact member-wise comparison of g
//                              with the snapshot. This is equivalent to "mesh_hash(g) equals
//                              the cached digest" (the digest is a function of exactly these
//                              members) without the collision caveat, and it detects every
//                              in-place mutation of a caller's nonconst mesh (A02 probes:
//                              dx, zmin, VE[0], edges[0].length).
// The digest value is the v6 mesh_hash value (same encoding, same bytes).
#include "admission.hpp"
#include <atomic>
#include <memory>
#include <mutex>
namespace electrode_heat_v6_1 {
namespace detail {
inline std::string mesh_hash(const axial::Mesh&g){Encoding e;e.bytes.reserve(std::size_t(16)<<20);e.s("source-v6-actual-axial-mesh");
 for(int x:{g.nx,g.ny,g.nz,g.nz_foil,g.kshort,g.kfeed,g.nn,g.ne,g.nb,g.nex,g.ney,g.nbx,g.nby})e.u(x);
 for(double x:{g.side,g.depth,g.dx,g.dy,g.dz,g.zmin})e.d(x);
 for(const auto*v:{&g.M,&g.H,&g.Y,&g.VE,&g.VB,&g.VN})e.doubles(*v);
 e.u(g.edges.size());for(const auto&x:g.edges){e.u(x.a);e.u(x.b);e.u(x.axis);e.d(x.length);}
 e.u(g.C.size());for(const auto&x:g.C){for(auto i:x.index)e.u(i);e.doubles(x.coefficient);e.u(x.count);}
 e.u(g.cells.size());for(const auto&x:g.cells){for(auto i:x.edges)e.u(i);for(int i:{x.i,x.j,x.k,x.body})e.u(i);e.d(x.volume);}
 e.u(g.node_body.size());for(auto i:g.node_body)e.u(i);
 for(const auto*v:{&g.short_node,&g.fixed_E,&g.fixed_A,&g.outer_node}){e.u(v->size());for(bool x:*v)e.u(x);}
 return e.sha();}
inline bool same_bits(double a,double b){return std::bit_cast<std::uint64_t>(a)==std::bit_cast<std::uint64_t>(b);}
// Bitwise equality of two double arrays (the v6 encoding compares exactly these bytes).
template<class V>inline bool same_doubles(const V&a,const V&b){return a.size()==b.size()&&(a.size()==0||std::memcmp(a.data(),b.data(),a.size()*sizeof(double))==0);}
// Exact equality of every member that mesh_hash encodes.
inline bool same_mesh(const axial::Mesh&a,const axial::Mesh&b){
 if(&a==&b)return true;
 const int ia[]={a.nx,a.ny,a.nz,a.nz_foil,a.kshort,a.kfeed,a.nn,a.ne,a.nb,a.nex,a.ney,a.nbx,a.nby},ib[]={b.nx,b.ny,b.nz,b.nz_foil,b.kshort,b.kfeed,b.nn,b.ne,b.nb,b.nex,b.ney,b.nbx,b.nby};
 for(std::size_t i=0;i<std::size(ia);++i)if(ia[i]!=ib[i])return false;
 const double da[]={a.side,a.depth,a.dx,a.dy,a.dz,a.zmin},db[]={b.side,b.depth,b.dx,b.dy,b.dz,b.zmin};
 for(std::size_t i=0;i<std::size(da);++i)if(!same_bits(da[i],db[i]))return false;
 if(!same_doubles(a.M,b.M)||!same_doubles(a.H,b.H)||!same_doubles(a.Y,b.Y)||!same_doubles(a.VE,b.VE)||!same_doubles(a.VB,b.VB)||!same_doubles(a.VN,b.VN))return false;
 if(a.edges.size()!=b.edges.size())return false;
 for(std::size_t i=0;i<a.edges.size();++i){const auto&x=a.edges[i],&y=b.edges[i];if(x.a!=y.a||x.b!=y.b||x.axis!=y.axis||!same_bits(x.length,y.length))return false;}
 if(a.C.size()!=b.C.size())return false;
 for(std::size_t i=0;i<a.C.size();++i){const auto&x=a.C[i],&y=b.C[i];if(x.index!=y.index||x.count!=y.count||!same_doubles(x.coefficient,y.coefficient))return false;}
 if(a.cells.size()!=b.cells.size())return false;
 for(std::size_t i=0;i<a.cells.size();++i){const auto&x=a.cells[i],&y=b.cells[i];if(x.edges!=y.edges||x.i!=y.i||x.j!=y.j||x.k!=y.k||x.body!=y.body||!same_bits(x.volume,y.volume))return false;}
 return a.node_body==b.node_body&&a.short_node==b.short_node&&a.fixed_E==b.fixed_E&&a.fixed_A==b.fixed_A&&a.outer_node==b.outer_node;}
}
class MeshIdentity final {
 struct Key {};
 std::shared_ptr<const axial::Mesh>mesh_;std::string sha_;std::uint64_t generation_=0;
 static std::uint64_t next_generation(){static std::atomic<std::uint64_t>n{0};return ++n;}
 struct Registry {std::mutex m;std::vector<std::pair<const axial::Mesh*,std::weak_ptr<const MeshIdentity>>>entries;};
 static Registry&registry(){static Registry r;return r;}
public:
 MeshIdentity(Key,const axial::Mesh&g):mesh_(std::make_shared<const axial::Mesh>(g)),sha_(detail::mesh_hash(*mesh_)),generation_(next_generation()){}
 MeshIdentity(const MeshIdentity&)=delete;MeshIdentity&operator=(const MeshIdentity&)=delete;
 // A new immutable snapshot of g, hashed once.
 static std::shared_ptr<const MeshIdentity>snapshot(const axial::Mesh&g){return std::make_shared<const MeshIdentity>(Key{},g);}
 // The identity of the object at &g; reused only when g is exactly equal to it.
 static std::shared_ptr<const MeshIdentity>of(const axial::Mesh&g){
  auto&r=registry();
  {std::lock_guard<std::mutex>lock(r.m);
   for(auto&[address,weak]:r.entries)if(address==&g)if(auto id=weak.lock();id&&id->same_as(g))return id;}
  auto id=snapshot(g);
  std::lock_guard<std::mutex>lock(r.m);
  std::erase_if(r.entries,[&](const auto&e){return e.first==&g||e.second.expired();});
  if(r.entries.size()>=64)r.entries.erase(r.entries.begin());
  r.entries.emplace_back(&g,id);return id;}
 const axial::Mesh&mesh()const noexcept{return *mesh_;}
 const std::shared_ptr<const axial::Mesh>&shared()const noexcept{return mesh_;}
 const std::string&sha256()const noexcept{return sha_;}
 std::uint64_t generation()const noexcept{return generation_;}
 // Identity of the immutable snapshot, or exact equality of every hashed member.
 bool same_as(const axial::Mesh&g)const{return &g==mesh_.get()||detail::same_mesh(*mesh_,g);}
};
}
