#pragma once
// Canonical carrier charge (unchanged arithmetic from 8dddec0 census_charge.hpp),
// the v2 deposit-rounding support of the endpoint continuity guard, and (v3) the
// bank transfers of kinetic sinks and surface births with the free charge of
// newborn carriers.
#include "field_model.hpp"
#include "motion_receipt.hpp"
#include "sink_ledger.hpp"

namespace hall_corrective_v3 {
namespace native=raikou::mfem_backend;
inline std::string Model::motion_mesh_sha256(const native::GapPicOptions&o)const{
    return ::hall_corrective_v3::motion_mesh_sha256(mesh_,o,mesh_json_);
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
        return finish(g,s1,s2,m,rows.size());
    }
    // v3: every carrier of the step with its actual start and end point (native
    // frame): census survivors (endpoint), kinetic sinks (exact native hit point)
    // and newborns (birth point -> endpoint or hit point). Same per-carrier bound.
    static DepositSupport build(const Mesh&g,const native::GapPicOptions&o,std::span<const native::ArenaRateEntry>census,
                                const hall_driver::MotionReceipt&r){
        if(r.sinks().empty()&&r.births().empty())return build(g,o,census,r.endpoints());
        constexpr long double eps=std::numeric_limits<double>::epsilon();
        need(census.size()==r.endpoints().size()+r.sinks().size(),"deposit support census/endpoint/sink count");
        std::vector<const native::ArenaRateEntry*>rows;rows.reserve(census.size());for(const auto&x:census)rows.push_back(&x);
        std::sort(rows.begin(),rows.end(),[](auto a,auto b){return a->native.id<b->native.id;});
        const std::array<int,3>counts{g.nx,g.ny,g.nz};const std::array<double,3>spacing{g.dx,g.dy,g.dz},origin{-g.side/2,-g.side/2,g.zmin};
        const long double C=counts[0]+counts[1]+counts[2];
        std::vector<long double>s1(g.nn),s2(g.nn);std::vector<std::uint64_t>m(g.nn);
        auto one=[&](const native::PicVector&a_native,const native::PicVector&b_native,long double Q){
            if(Q==0)return;
            const native::PicVector start{a_native[0]-o.width_m/2,a_native[1]-o.gap_m/2,a_native[2]};
            const native::PicVector end{b_native[0]-o.width_m/2,b_native[1]-o.gap_m/2,b_native[2]};
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
        };
        std::size_t ie=0,is=0;
        for(std::size_t i=0;i<rows.size();++i){const auto&p=rows[i]->native;
            need(i==0||rows[i-1]->native.id<p.id,"deposit support carrier identity/order");
            const long double Q=magnitude(Words(p.weight)*p.charge_C);need(std::isfinite(double(Q)),"deposit support charge operand");
            if(is<r.sinks().size()&&r.sinks()[is].native_id==p.id){one(p.position_m,r.sinks()[is].hit_point_m,Q);++is;}
            else{need(ie<r.endpoints().size()&&r.endpoints()[ie].id==p.id,"deposit support carrier identity/order");one(p.position_m,r.endpoints()[ie].position_m,Q);++ie;}}
        return finish(g,s1,s2,m,rows.size());
    }
private:
    static DepositSupport finish(const Mesh&g,const std::vector<long double>&s1,const std::vector<long double>&s2,const std::vector<std::uint64_t>&m,std::size_t carriers){
        constexpr long double eps=std::numeric_limits<double>::epsilon();
        DepositSupport out;out.absolute_C_.assign(g.nn,0.);out.carriers_=carriers;
        for(int n=0;n<g.nn;++n){if(!m[n])continue;
            // 1e-9 relative margin covers the rounding of this bound's own sums.
            const long double bound=eps*(s1[n]+24*(long double)m[n]*s2[n])*(1+1e-9L);
            double b=double(bound);if((long double)b<bound)b=std::nextafter(b,INFINITY);out.absolute_C_[n]=b;}
        return out;
    }
};

// v3 bank transfer of one step: collected sinks add their macro-charge to the
// conductor bank and escaped sinks to the exterior bank (Words CIC deposit at the
// exact native hit point, the operations of charge_of_endpoints, so the bank
// matches the truncated path deposit); a charged newborn on a foil face takes its
// macro-charge from that foil's conductor bank at its birth point. Order: sinks,
// then newborns, each ascending native ID. CIC weight that rounding puts on a node
// outside the bank support (cell coordinate of a plane a few ulps off its node
// plane) stays in the bank as residue, admitted only up to 1e-12 of the macro-charge
// per deposit and bounded per node in ChargeLedger::residue_C.
inline BankTransfer bank_transfer(const Mesh&g,const native::GapPicOptions&o,const hall_driver::MotionReceipt&r,const std::vector<bool>&exterior_support){
    BankTransfer t;if(r.sinks().empty()&&r.births().empty())return t;
    t.empty=false;FreeCharge cond{WordVec(g.nn),Vec(g.nn)},ext{WordVec(g.nn),Vec(g.nn)};t.residue_C.assign(g.nn,0.);
    auto field=[&](const native::PicVector&x){auto f=x;f[0]-=o.width_m/2;f[1]-=o.gap_m/2;return f;};
    auto conductor=[&](int n){return g.node_body[n]>=0||g.short_node[n];};
    auto deposit=[&](bool to_conductor,const native::PicVector&x,Words q){
        FreeCharge one{WordVec(g.nn),Vec(g.nn)};deposit_charge(g,one,field(x),q);auto&bank=to_conductor?cond:ext;
        // the 8 CIC nodes of x
        const long double Q=magnitude(q);
        for(int n=0;n<g.nn;++n){if(one.net_C[n].hi==0&&one.net_C[n].lo==0)continue;
            bank.net_C[n]=bank.net_C[n]+one.net_C[n];bank.gross_C[n]+=one.gross_C[n];
            const bool supported=to_conductor?conductor(n):bool(exterior_support[n]);
            if(!supported){const long double m=magnitude(one.net_C[n]);
                need(m<=1e-12L*Q,to_conductor?"kinetic sink/newborn conductor bank deposit outside the conductor nodes":
                                             "kinetic sink exterior bank deposit outside the exterior support");
                const long double x2=(long double)t.residue_C[n]+m;double v=double(x2);if((long double)v<x2)v=std::nextafter(v,INFINITY);t.residue_C[n]=v;}}
    };
    for(const auto&k:r.sinks()){if(k.charge_C==0)continue;const auto q=Words(k.weight)*k.charge_C;
        need(canonical(q)&&(q.hi!=0||q.lo!=0),"kinetic sink unresolved macro-charge");deposit(k.collected,k.hit_point_m,q);}
    for(const auto&b:r.births()){if(b.charge_C==0)continue;
        need(b.origin>=0,"charged interior newborn needs the source path current and charge impulses (not admitted by this Hall version)");
        const auto q=Words(b.weight)*b.charge_C;need(canonical(q)&&(q.hi!=0||q.lo!=0),"newborn unresolved macro-charge");deposit(true,b.start_m,-q);}
    t.conductor_C=std::move(cond.net_C);t.exterior_C=std::move(ext.net_C);t.gross_C.resize(g.nn);
    for(int n=0;n<g.nn;++n){const long double x=(long double)cond.gross_C[n]+ext.gross_C[n];double v=double(x);if((long double)v<x)v=std::nextafter(v,INFINITY);t.gross_C[n]=v;}
    return t;
}
// Census of the next step: the survivors (residents and newborns), ascending native ID.
inline std::vector<native::ArenaRateEntry> next_census(std::span<const native::ArenaRateEntry>census,const hall_driver::MotionReceipt&r){
    std::vector<const native::ArenaRateEntry*>rows;for(const auto&x:census)rows.push_back(&x);
    std::sort(rows.begin(),rows.end(),[](auto a,auto b){return a->native.id<b->native.id;});
    std::vector<native::ArenaRateEntry>out;out.reserve(r.endpoints().size());std::size_t ie=0;
    for(const auto*row:rows){if(ie<r.endpoints().size()&&r.endpoints()[ie].id==row->native.id){auto x=*row;x.native=r.endpoints()[ie++];out.push_back(std::move(x));}}
    need(ie==r.endpoints().size(),"next census endpoint order");
    return out;
}
} // namespace hall_corrective_v3
