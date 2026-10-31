#pragma once
// Canonical carrier charge (unchanged arithmetic from 8dddec0 census_charge.hpp)
// and the v2 deposit-rounding support of the endpoint continuity guard.
#include "field_model.hpp"
#include "motion_receipt.hpp"

namespace hall_corrective_v2 {
namespace native=raikou::mfem_backend;
inline std::string Model::motion_mesh_sha256(const native::GapPicOptions&o)const{
    return ::hall_corrective_v2::motion_mesh_sha256(mesh_,o,mesh_json_);
}
inline void deposit_charge(const Mesh&g,FreeCharge&q,native::PicVector x,Words charge){
    need(canonical(charge),"invalid deposited charge words");
    native::PicVector spacing{g.dx,g.dy,g.dz},origin{-g.side/2,-g.side/2,g.zmin},a{};
    std::array<int,3>counts{g.nx,g.ny,g.nz},cell{};
    for(int d=0;d<3;++d){need(std::isfinite(x[d]),"nonfinite charge position");double v=(x[d]-origin[d])/spacing[d];
        need(v>=-1e-10&&v<=counts[d]+1e-10,"canonical charge outside mesh");v=std::clamp(v,0.,double(counts[d]));
        cell[d]=std::min(int(std::floor(v)),counts[d]-1);a[d]=v-cell[d];}
    for(int k:{0,1})for(int j:{0,1})for(int i:{0,1}){
        const int n=g.node(cell[0]+i,cell[1]+j,cell[2]+k);
        const auto v=charge*(i?a[0]:1-a[0])*(j?a[1]:1-a[1])*(k?a[2]:1-a[2]);q.net_C[n]=q.net_C[n]+v;
        const long double gross=(long double)q.gross_C[n]+magnitude(v);
        q.gross_C[n]=double(gross);if((long double)q.gross_C[n]<gross)q.gross_C[n]=std::nextafter(q.gross_C[n],INFINITY);
    }
}
// Read actual native positive carriers individually. No Count2 conversion,
// packet allocation, ID assignment, population repair or rounded count view.
// motion_mesh is the owner's precomputed motion mesh identity (geometry guard).
inline FreeCharge charge_of_particles(const Mesh&g,const native::GapPicOptions&o,
    std::span<const native::GapParticle*const>particles,std::size_t cap=100000){
    need(particles.size()<=cap,"charge census exceeds motion carrier cap");
    auto policy=hall_driver::declared_kinetic_policy(o);
    std::vector<const native::GapParticle*>ordered(particles.begin(),particles.end());
    for(auto p:ordered)need(p!=nullptr,"null canonical carrier");
    std::sort(ordered.begin(),ordered.end(),[](auto a,auto b){return a->id<b->id;});
    FreeCharge out{WordVec(g.nn),Vec(g.nn)};std::uint64_t previous=0;
    for(auto p:ordered){need(p&&p->id>previous,"duplicate/missing canonical carrier ID");previous=p->id;
        need(std::isfinite(p->weight)&&p->weight>0&&std::isfinite(p->charge_C),"invalid canonical charge operand");
        hall_driver::require_kinetic_interior(policy,p->position_m);auto x=p->position_m;x[0]-=o.width_m/2;x[1]-=o.gap_m/2;
        const auto q=Words(p->weight)*p->charge_C;need(canonical(q)&&(p->charge_C==0||q.hi!=0||q.lo!=0),"unresolved canonical charge");
        deposit_charge(g,out,x,q);}
    free_charge_guard(g,out);return out;
}
inline FreeCharge charge_of_census(const Mesh&g,const native::GapPicOptions&o,std::span<const native::ArenaRateEntry>census,std::size_t cap=100000){
    require_census_identity(census,cap);std::vector<const native::GapParticle*>particles;particles.reserve(census.size());
    for(const auto&row:census)particles.push_back(&row.native);return charge_of_particles(g,o,particles,cap);
}
inline FreeCharge charge_of_endpoints(const Mesh&g,const native::GapPicOptions&o,std::span<const native::GapParticle>endpoints,std::size_t cap=100000){
    std::vector<const native::GapParticle*>particles;particles.reserve(endpoints.size());for(const auto&p:endpoints)particles.push_back(&p);
    return charge_of_particles(g,o,particles,cap);
}
inline bool same_net_charge(const FreeCharge&a,const FreeCharge&b){
    if(a.net_C.size()!=b.net_C.size()||a.gross_C.size()!=a.net_C.size()||b.gross_C.size()!=b.net_C.size())return false;
    for(std::size_t n=0;n<a.net_C.size();++n)if(std::abs((a.net_C[n]-b.net_C[n]).wide())>
        256*std::numeric_limits<double>::epsilon()*((long double)a.gross_C[n]+b.gross_C[n]))return false;
    return true;
}

// Deposit-rounding support of the endpoint continuity guard (see README and
// field_model.hpp endpoint_tolerance_rule). In exact arithmetic the CIC point
// deposits and the cut-at-cell-planes, two-point-Gauss path current conserve
// charge exactly per node; the residual is rounding only. For carrier p with
// macro-charge Q_p, cell counts c_d (C = sum c_d) and path R_p = sum_d |dx_d|/h_d
// (<= 0.6 by the motion displacement guard), with eps = 2^-52:
//   start and end point deposits (cell coordinate (x-o)/h: <= eps (c_d+1);
//   1-a: eps/2; double-word products O(eps^2))            <= 2 eps Q_p (C+7)
//   path current: rounded macrocharge, product/quotient roundings, quadrature
//   coordinates (<= eps(2 c_d+3) per factor), cut parameters (misassigned
//   weights O(eps(3c_d+4))), endpoint representation          <= eps Q_p R_p (16C+67)
//   accumulation of M_n carriers' contributions into plain binary64 edge sums
//   (recursive summation, <= 8 terms per carrier and edge, 6 edges per node)
//                                                            <= 24 eps M_n sum Q_p R_p
// The bound is applied at every node of the carrier's cell range: the cells of
// the start, the end and every path segment lie in [floor(min(a,b)-delta),
// floor(max(a,b)+delta)] per axis, delta = eps(3c_d+8) cells covering the
// coordinate rounding. Nodes outside every carrier's range get no absolute term.
class DepositSupportBuilder {
public:
    static DepositSupport build(const Mesh&g,const native::GapPicOptions&o,std::span<const native::ArenaRateEntry>census,
                                std::span<const native::GapParticle>endpoints){
        constexpr long double eps=std::numeric_limits<double>::epsilon();
        need(census.size()==endpoints.size(),"deposit support census/endpoint count");
        std::vector<const native::ArenaRateEntry*>rows;rows.reserve(census.size());for(const auto&r:census)rows.push_back(&r);
        std::sort(rows.begin(),rows.end(),[](auto a,auto b){return a->native.id<b->native.id;});
        const std::array<int,3>counts{g.nx,g.ny,g.nz};const std::array<double,3>spacing{g.dx,g.dy,g.dz},origin{-g.side/2,-g.side/2,g.zmin};
        const long double C=counts[0]+counts[1]+counts[2];
        std::vector<long double>s1(g.nn),s2(g.nn);std::vector<std::uint64_t>m(g.nn);
        for(std::size_t i=0;i<rows.size();++i){const auto&p=rows[i]->native;const auto&e=endpoints[i];
            need(p.id==e.id&&(i==0||rows[i-1]->native.id<p.id),"deposit support carrier identity/order");
            const long double Q=magnitude(Words(p.weight)*p.charge_C);need(std::isfinite(double(Q)),"deposit support charge operand");
            if(Q==0)continue;
            const native::PicVector start{p.position_m[0]-o.width_m/2,p.position_m[1]-o.gap_m/2,p.position_m[2]};
            const native::PicVector end{e.position_m[0]-o.width_m/2,e.position_m[1]-o.gap_m/2,e.position_m[2]};
            std::array<int,3>lo{},hi{};long double R=0;
            for(int d=0;d<3;++d){const double a=(start[d]-origin[d])/spacing[d],b=(end[d]-origin[d])/spacing[d];
                need(std::isfinite(a)&&std::isfinite(b),"deposit support nonfinite coordinate");
                const long double delta=eps*(3*counts[d]+8);R+=std::abs((long double)b-a)+2*delta;
                lo[d]=std::clamp(int(std::floor((long double)std::min(a,b)-delta)),0,counts[d]-1);
                hi[d]=std::clamp(int(std::floor((long double)std::max(a,b)+delta)),0,counts[d]-1);}
            need(R<=1,"deposit support path exceeds the motion displacement guard");
            const long double t1=Q*(2*C+14+(16*C+67)*R),t2=Q*R;
            for(int k=lo[2];k<=hi[2]+1;++k)for(int j=lo[1];j<=hi[1]+1;++j)for(int ii=lo[0];ii<=hi[0]+1;++ii){
                const int n=g.node(ii,j,k);s1[n]+=t1;s2[n]+=t2;++m[n];}
        }
        DepositSupport out;out.absolute_C_.assign(g.nn,0.);out.carriers_=rows.size();
        for(int n=0;n<g.nn;++n){if(!m[n])continue;
            // 1e-9 relative margin covers the rounding of this bound's own sums.
            const long double bound=eps*(s1[n]+24*(long double)m[n]*s2[n])*(1+1e-9L);
            double b=double(bound);if((long double)b<bound)b=std::nextafter(b,INFINITY);out.absolute_C_[n]=b;}
        return out;
    }
};
} // namespace hall_corrective_v2
