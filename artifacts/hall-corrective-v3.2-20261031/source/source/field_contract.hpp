#pragma once
// Hall corrective field contract, v3. Derived from Hall corrective v2 (itself
// derived from hall-field-correction-20261030/field_contract.hpp, commit 8dddec0);
// the field arithmetic, the compensated words, the absolute Gauss/divB admission
// and the C19 drive are unchanged. v2: identities are canonical binary encodings.
// v3: the exterior charge bank may also hold charge on the kinetic escape
// support (the target-end node planes, see Model::exterior_support), where
// escaped carriers are banked; BankTransfer carries the per-step conductor and
// exterior bank changes of kinetic sinks and surface births.
// v3.2: the charge ledger carries the admitted Gauss ledger (per-node signed
// accumulation of every accepted step's admitted Gauss-increment mismatch minus
// its admitted endpoint continuity defect, compensated words); the absolute Gauss
// admission measures the defect relative to it (admitted_gauss_ledger_rule). The
// ledger is part of the field/charge state digest and of the checkpoint.
#include "hall-finite-foil-em-20261029/axial_fields.hpp"
#include "encoding.hpp"
#include <iomanip>

namespace hall_driver {std::string total_field_sha256(const std::vector<double>&E,const std::vector<double>&B);}
namespace hall_corrective_v3 {
using axial::Vec;using axial::Mesh;using Json=raikou::Json;
inline constexpr char api[]="hall-compensated-field-owner-v3";
inline constexpr char c19_sha256[]="2c54d2dbaec6b6c8fcf2e45b5abd837aa72999f5fd14217f530277ffd5949076";
inline void need(bool ok,const char*message){hall::require(ok,message);}
inline std::string sha256(const std::string&s){return bytes_sha256(s);}
inline bool digest(const std::string&s){return s.size()==64&&s.find_first_not_of("0123456789abcdef")==std::string::npos;}
// Error-free two-sum/product building blocks; no change to native Count2.
// These words represent electromagnetic field/charge arithmetic, not inventory.
struct Words {
    double hi=0,lo=0;
    Words()=default;explicit Words(double value):hi(value){}
    Words(double a,double b){auto pair=axial::two_sum(a,b);hi=pair.first;lo=pair.second;}
    double rounded()const{return hi+lo;}
    long double wide()const{return (long double)hi+lo;}
};
inline Words operator+(Words a,Words b){auto first=axial::two_sum(a.hi,b.hi);
    return Words(first.first,first.second+a.lo+b.lo);}
inline Words operator-(Words a){return Words(-a.hi,-a.lo);}
inline Words operator-(Words a,Words b){return a+(-b);}
inline Words operator*(Words a,double b){const double p=a.hi*b;
    return Words(p,std::fma(a.hi,b,-p)+a.lo*b);}
inline Words operator/(Words a,double b){const double q=a.hi/b;const auto r=a-Words(q)*b;
    return Words(q,(r.hi+r.lo)/b);}
inline long double magnitude(Words a){return std::abs((long double)a.hi)+std::abs((long double)a.lo);}
using WordVec=std::vector<Words>;
static_assert(sizeof(Words)==16&&std::is_standard_layout_v<Words>,"word pairs encode as consecutive binary64 values");
inline bool canonical(Words a){const Words n(a.hi,a.lo);return std::isfinite(a.hi)&&std::isfinite(a.lo)&&n.hi==a.hi&&n.lo==a.lo;}
// Word vectors encode as count + (hi, lo) binary64 pairs in index order.
inline void encode(Digest&d,const WordVec&v){d.u64(v.size());for(const auto&x:v){d.f64(x.hi);d.f64(x.lo);}}
inline WordVec words(const Vec&hi,const Vec&lo){need(hi.size()==lo.size(),"field word dimensions");WordVec v;v.reserve(hi.size());
    for(std::size_t i=0;i<hi.size();++i){Words x;x.hi=hi[i];x.lo=lo[i];need(canonical(x),"noncanonical field words");v.push_back(x);}return v;}
inline WordVec incidence(const Mesh&g,const WordVec&v){need(v.size()==std::size_t(g.ne),"current dimensions");WordVec out(g.nn);
    for(int e=0;e<g.ne;++e){const auto q=v[e]/g.edges[e].length;out[g.edges[e].a]=out[g.edges[e].a]-q;out[g.edges[e].b]=out[g.edges[e].b]+q;}return out;}
inline WordVec incidence(const Mesh&g,const Vec&v){WordVec w;w.reserve(v.size());for(double x:v)w.emplace_back(x);return incidence(g,w);}
inline WordVec field_charge(const Mesh&g,const axial::State&s){auto e=words(s.E,s.E_low);for(int i=0;i<g.ne;++i)e[i]=e[i]*(-g.M[i]);return incidence(g,e);}
struct FreeCharge {WordVec net_C;Vec gross_C;};
struct ChargeLedger {
    FreeCharge particles;WordVec conductor_C,exterior_C;
    std::string conductor_declaration,exterior_declaration;
    // v3: bank deposit residue bound (C) per node: cumulative magnitude of the
    // rounding-level CIC weight of sink/newborn bank deposits that fell on nodes
    // outside the bank support (empty = zero everywhere, the v2 ledger).
    Vec residue_C;
    // v3.2: admitted Gauss ledger L (C) per node, signed compensated words:
    // L_n += delta_n - c_n for every accepted step, delta_n = Q1-Q0-(div plasma+div metal+div outer)dt
    // (admitted by the Gauss increment guard), c_n = after+banks-before-div(plasma)dt (admitted by the
    // endpoint continuity guard; zero in field-only steps). Empty = zero everywhere (bootstrap).
    WordVec gauss_admitted_C;
};
// Guard values of an admission. v3.2: max_gauss_C is the checked quantity max_n |D_n - L_n|
// (D_n the absolute Gauss defect, L_n the admitted Gauss ledger; equal to the v3.1 value where L = 0),
// max_gauss_bound_C the unchanged bound expression (maximum over nodes); the remaining fields are
// reported values: max_n |D_n|, max_n |L_n|, max_n |D_n-L_n|/bound_n (the per-node guard ratio) and
// max_n |D_n|/(bound_n+|L_n|) (implied by the guard).
struct Constraints {double max_gauss_C=0,max_gauss_bound_C=0,max_divB_T_m=0,max_divB_bound_T_m=0;
    double max_gauss_defect_C=0,max_gauss_ledger_C=0,max_gauss_node_ratio=0,max_gauss_defect_over_bound_plus_ledger=0;};
// v3.2 rule text (part of the configuration identity).
inline constexpr char admitted_gauss_ledger_rule[]=
    "admitted-gauss-ledger-v1: per node n the state carries L_n (compensated words, zero at bootstrap); every accepted "
    "step adds delta_n - c_n with delta_n = (Q1-Q0)_n - (div plasma + div metal + div outer)_n dt (the per-node value "
    "bounded by the Gauss increment guard) and c_n = (after + bank conductor + bank exterior - before)_n - div(plasma)_n dt "
    "(the per-node value bounded by the endpoint continuity guard); the absolute Gauss admission requires "
    "|D_n - L_n| <= bound_n with D_n = Q_n - free_n - conductor_n - exterior_n and the unchanged "
    "bound_n = epsilon*max(1e-30,2e-13*initial_vacuum_rhs_L2) + 256 eps (activity_n + gross_n + |conductor_n| + |exterior_n|); "
    "this implies |D_n| <= bound_n + |L_n|";
inline bool zero_words(const WordVec&v){for(const auto&x:v)if(x.hi!=0||x.lo!=0)return false;return true;}
// Per-step bank changes of kinetic sinks (collected -> conductor, escaped ->
// exterior) and surface births (conductor -> free). gross_C bounds the
// magnitude of every deposit (endpoint continuity scale). empty: no sink and
// no birth in the step; nothing is added to the banks (v2 arithmetic).
struct BankTransfer {WordVec conductor_C,exterior_C;Vec gross_C,residue_C;bool empty=true;};
inline void free_charge_guard(const Mesh&g,const FreeCharge&q){need(q.net_C.size()==std::size_t(g.nn)&&q.gross_C.size()==q.net_C.size(),"free charge shape");
    for(int n=0;n<g.nn;++n)need(canonical(q.net_C[n])&&std::isfinite(q.gross_C[n])&&q.gross_C[n]>=0&&
        std::abs(q.net_C[n].wide())<=q.gross_C[n]*(1+32*std::numeric_limits<double>::epsilon()),"free charge/gross operand invalid");}
// exterior_support: nodes where exterior charge may be held (EM exterior nodes
// and the kinetic escape planes). nullptr = v2 rule (EM exterior nodes only).
inline Constraints admit_constraints(const Mesh&g,const axial::State&s,const ChargeLedger&q,double initial_solve_bound_C,
                                     const std::vector<bool>*exterior_support=nullptr){
    need(std::isfinite(initial_solve_bound_C)&&initial_solve_bound_C>=0,"immutable initial solver bound");
    need(s.E.size()==std::size_t(g.ne)&&s.B.size()==std::size_t(g.nb)&&s.E_low.size()==s.E.size()&&s.B_low.size()==s.B.size(),"field state shape");
    need(std::isfinite(s.time_s)&&s.time_s>=0&&s.revision<UINT64_MAX,"field clock/revision");
    const auto E=words(s.E,s.E_low),B=words(s.B,s.B_low);free_charge_guard(g,q.particles);
    need(q.conductor_C.size()==std::size_t(g.nn)&&q.exterior_C.size()==std::size_t(g.nn)&&
        digest(q.conductor_declaration)&&digest(q.exterior_declaration),"explicit conductor/exterior charge declarations required");
    need(q.gauss_admitted_C.empty()||q.gauss_admitted_C.size()==std::size_t(g.nn),"admitted Gauss ledger shape");
    std::vector<long double>activity(g.nn);for(int e=0;e<g.ne;++e){const auto&edge=g.edges[e];
        need(!g.fixed_E[e]||(E[e].hi==0&&E[e].lo==0),"nonzero PEC tangential field");
        const long double a=magnitude(E[e])*g.M[e]/edge.length;activity[edge.a]+=a;activity[edge.b]+=a;}
    const auto Q=field_charge(g,s);Constraints r;
    for(int n=0;n<g.nn;++n){need(canonical(q.conductor_C[n])&&canonical(q.exterior_C[n]),"noncanonical owned charge");
        const long double residue=q.residue_C.empty()?0.L:(long double)q.residue_C[n];
        need((g.node_body[n]>=0||g.short_node[n])||(q.conductor_C[n].hi==0&&q.conductor_C[n].lo==0)||magnitude(q.conductor_C[n])<=residue,
            "conductor charge outside declared material/short");
        need((exterior_support?bool((*exterior_support)[n]):bool(g.outer_node[n]))||(q.exterior_C[n].hi==0&&q.exterior_C[n].lo==0)||
             magnitude(q.exterior_C[n])<=residue,
            exterior_support?"exterior charge outside the EM boundary and the kinetic escape planes":"exterior charge outside EM boundary; kinetic sinks unsupported");
        const Words ledger=q.gauss_admitted_C.empty()?Words():q.gauss_admitted_C[n];
        need(canonical(ledger),"noncanonical admitted Gauss ledger");
        const auto defect=Q[n]-q.particles.net_C[n]-q.conductor_C[n]-q.exterior_C[n];
        const auto centered=defect-ledger;   // v3.2: defect relative to the admitted Gauss ledger
        const long double bound=initial_solve_bound_C+256*std::numeric_limits<double>::epsilon()*
            (activity[n]+q.particles.gross_C[n]+magnitude(q.conductor_C[n])+magnitude(q.exterior_C[n]));
        if(!(std::abs(centered.wide())<=bound)){std::ostringstream why;why<<std::setprecision(17)
            <<"absolute Gauss constraint disagrees with actual declared charge: node="<<n
            <<" defect_C="<<double(centered.wide())<<" bound_C="<<double(bound)<<" field_activity_C="<<double(activity[n]);
            if(ledger.hi!=0||ledger.lo!=0)why<<" admitted_gauss_ledger_C="<<double(ledger.wide())<<" uncorrected_defect_C="<<double(defect.wide());
            throw std::runtime_error(why.str());}
        const long double size=std::abs(centered.wide()),raw=std::abs(defect.wide()),led=std::abs(ledger.wide());
        r.max_gauss_C=std::max(r.max_gauss_C,double(size));r.max_gauss_bound_C=std::max(r.max_gauss_bound_C,double(bound));
        r.max_gauss_defect_C=std::max(r.max_gauss_defect_C,double(raw));r.max_gauss_ledger_C=std::max(r.max_gauss_ledger_C,double(led));
        if(bound>0){r.max_gauss_node_ratio=std::max(r.max_gauss_node_ratio,double(size/bound));
            r.max_gauss_defect_over_bound_plus_ledger=std::max(r.max_gauss_defect_over_bound_plus_ledger,double(raw/(bound+led)));}}
    for(int k=0;k<g.nz;++k)for(int j=0;j<g.ny;++j)for(int i=0;i<g.nx;++i){Words d;long double scale=0;
        for(auto [a,b,h]:std::array<std::tuple<int,int,double>,3>{{{g.bx(i+1,j,k),g.bx(i,j,k),g.dx},{g.by(i,j+1,k),g.by(i,j,k),g.dy},{g.bz(i,j,k+1),g.bz(i,j,k),g.dz}}}){
            d=d+(B[a]-B[b])/h;scale+=(magnitude(B[a])+magnitude(B[b]))/h;}
        const long double bound=128*std::numeric_limits<double>::epsilon()*scale;
        need(std::abs(d.wide())<=bound,"absolute divB constraint failed");
        r.max_divB_T_m=std::max(r.max_divB_T_m,double(std::abs(d.wide())));r.max_divB_bound_T_m=std::max(r.max_divB_bound_T_m,double(bound));}
    return r;
}
class CurrentDrive {
    std::vector<std::pair<double,double>>knots_;
public:
    explicit CurrentDrive(const std::filesystem::path&path){std::ifstream in(path,std::ios::binary);
        need(bool(in),"missing sealed C19 current");std::string bytes((std::istreambuf_iterator<char>(in)),{});
        need(sha256(bytes)==c19_sha256,"C19 bytes differ from selected current-only drive");
        std::istringstream all(bytes);std::string line;std::getline(all,line);
        while(std::getline(all,line)){std::istringstream row(line);std::string cell;std::vector<std::string>v;
            while(std::getline(row,cell,','))v.push_back(cell);need(v.size()==12,"C19 row shape");
            double t=std::stod(v[4]),I=std::stod(v[8]);need(std::isfinite(t)&&std::isfinite(I)&&(knots_.empty()||t>knots_.back().first),"C19 finite monotone support");knots_.emplace_back(t,I);}
        need(knots_.size()>1,"C19 support empty");}
    double start()const{return knots_.front().first;}double end()const{return knots_.back().first;}
    // v3.1 (diagnostic): spacing of the C19 knots around t (controller receipt).
    double knot_interval(double t)const{need(std::isfinite(t)&&t>=start()&&t<=end(),"C19 query outside original support");
        auto high=std::upper_bound(knots_.begin(),knots_.end(),t,[](double x,const auto&p){return x<p.first;});
        if(high==knots_.end())--high;if(high==knots_.begin())++high;return high->first-(high-1)->first;}
    double current(double t)const{need(std::isfinite(t)&&t>=start()&&t<=end(),"C19 query outside original support");if(t==end())return knots_.back().second;
        auto high=std::upper_bound(knots_.begin(),knots_.end(),t,[](double x,const auto&p){return x<p.first;});auto low=high-1;
        return low->second+(t-low->first)/(high->first-low->first)*(high->second-low->second);}
};
struct MeshSpec {int transverse=45,axial=4,air=1;double width_m=.009;};
// Complete binary record of the immutable mesh: every scalar parameter, every
// operator array (M, H, Y, VE, VB, VN), the edge/curl/cell incidence and the
// material, short, PEC, gauge and exterior masks.
inline void encode_mesh(Digest&d,const Mesh&g){
    d.text("hall-v2-axial-mesh-v1");
    for(int v:{g.nx,g.ny,g.nz,g.nz_foil,g.kshort,g.kfeed,g.nn,g.ne,g.nb,g.nex,g.ney,g.nbx,g.nby})d.i64(v);
    for(double v:{g.side,g.depth,g.dx,g.dy,g.dz,g.zmin})d.f64(v);
    for(const Vec*v:{&g.M,&g.H,&g.Y,&g.VE,&g.VB,&g.VN})d.f64s(*v);
    d.u64(g.edges.size());for(const auto&e:g.edges){d.i64(e.a);d.i64(e.b);d.i64(e.axis);d.f64(e.length);}
    d.u64(g.C.size());for(const auto&c:g.C){d.i64(c.count);for(int k=0;k<4;++k){d.i64(c.index[k]);d.f64(c.coefficient[k]);}}
    d.u64(g.cells.size());for(const auto&c:g.cells){for(int e:c.edges)d.i64(e);d.i64(c.i);d.i64(c.j);d.i64(c.k);d.i64(c.body);d.f64(c.volume);}
    d.i32s_as_i64(g.node_body);d.bools(g.short_node);d.bools(g.fixed_E);d.bools(g.fixed_A);d.bools(g.outer_node);
}
// Field state: both stored words of E and B, clock and revision.
inline void encode_state(Digest&d,const axial::State&s){
    d.text("hall-v2-field-words-v1");d.f64s(s.E);d.f64s(s.B);d.f64s(s.E_low);d.f64s(s.B_low);d.f64(s.time_s);d.u64(s.revision);
}
// v3.2: the admitted Gauss ledger follows (always nn words; an empty ledger encodes as zeros, so the
// identity does not depend on the in-memory representation of a zero ledger).
inline void encode_ledger(Digest&d,const ChargeLedger&q){
    d.text("hall-v3-charge-ledger-v1");encode(d,q.particles.net_C);d.f64s(q.particles.gross_C);
    encode(d,q.conductor_C);encode(d,q.exterior_C);d.text(q.conductor_declaration);d.text(q.exterior_declaration);d.f64s(q.residue_C);
    d.text("hall-v3.2-admitted-gauss-ledger-v1");
    if(q.gauss_admitted_C.empty()){d.u64(q.particles.net_C.size());for(std::size_t i=0;i<q.particles.net_C.size();++i){d.f64(0.);d.f64(0.);}}
    else encode(d,q.gauss_admitted_C);
}
// JSON records are used only for the human-readable checkpoint export.
inline Json record(const WordVec&v){Json out=Json::array();for(auto x:v){Json pair=Json::array();pair.push_back(x.hi);pair.push_back(x.lo);out.push_back(std::move(pair));}return out;}
inline Json record(const ChargeLedger&q){Json j=Json::object();j["free_C"]=record(q.particles.net_C);j["gross_free_C"]=q.particles.gross_C;
    j["conductor_C"]=record(q.conductor_C);j["exterior_C"]=record(q.exterior_C);
    j["conductor_declaration"]=q.conductor_declaration;j["exterior_declaration"]=q.exterior_declaration;return j;}
inline Json state_record(const axial::State&s){Json j=Json::object();j["E"]=s.E;j["B"]=s.B;j["E_low"]=s.E_low;j["B_low"]=s.B_low;
    j["time_s"]=s.time_s;j["revision"]=s.revision;return j;}
// Clears nested JSON bottom-up so that its destructor never allocates.
inline void drain(Json&j)noexcept{
    if(j.is_object()){auto&o=j.get_ref<Json::object_t&>();for(auto&item:o)drain(item.second);o.clear();}
    else if(j.is_array()){auto&a=j.get_ref<Json::array_t&>();for(auto&item:a)drain(item);a.clear();}
}
class Midpoint {
    Vec E_,B_;std::string seal_;
    mutable std::string total_;   // v3.1: total-field digest of this immutable object, computed once on demand
    static std::string role_seal(const Vec&E,const Vec&B){Digest d("authoritative-compensated-midpoint-v2");d.f64s(E);d.f64s(B);return d.hex();}
public:
    Midpoint(Vec E,Vec B):E_(std::move(E)),B_(std::move(B)),seal_(role_seal(E_,B_)){}
    const Vec&E()const noexcept{return E_;}const Vec&B()const noexcept{return B_;}const std::string&seal()const noexcept{return seal_;}
    // hall_driver::total_field_sha256(E(), B()) (unchanged record, same value), cached per object.
    const std::string&total_field_sha256()const{if(total_.empty())total_=hall_driver::total_field_sha256(E_,B_);return total_;}
};
// v3.1: a total-field digest computed by the caller from these exact vector objects in the same call
// sequence; motion admission then skips its recomputation (same value, same refusal otherwise).
struct VerifiedField {const Vec*E=nullptr;const Vec*B=nullptr;std::string digest;};
} // namespace hall_corrective_v3
