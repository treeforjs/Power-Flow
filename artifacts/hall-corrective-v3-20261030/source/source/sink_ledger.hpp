#pragma once
// Per-boundary kinetic sink ledger and surface birth ledger (Hall corrective v3).
// Counts are exact integers; weights and the elementary-charge totals are exact
// (ExactParticleCount, the sealed exact-count header of native candidate003,
// SHA-256 46e01fe2...); charge and kinetic energy are compensated double-word
// sums of the per-carrier values (charge: exact two-product w*q per carrier).
#include "field_contract.hpp"
#include "motion_receipt.hpp"
#include "raikou/exact_particle_count.h"
#include <array>

namespace hall_corrective_v3 {
using ExactCount=raikou::mfem_backend::ExactParticleCount;
inline constexpr const char* sink_boundary_names[5]={"collected/lower-foil","collected/upper-foil","escaped/outer-xy",
                                                    "collected/short-bridge","escaped/axial"};
inline constexpr const char* conductor_names[3]={"lower-foil","upper-foil","short-bridge"};
inline ExactCount exact_sum(const ExactCount&a,const ExactCount&b){auto r=a;for(double w:b.positive_carriers())r=r.credited(w);return r;}
inline void encode_exact(Digest&d,const ExactCount&c){const auto l=c.positive_carriers();d.u64(l.size());for(double w:l)d.f64(w);}
inline Json exact_json(const ExactCount&c){Json a=Json::array();for(double w:c.positive_carriers())a.push_back(w);return a;}
struct ChargeTally {
    std::uint64_t count=0,non_elementary=0;ExactCount weight,positive_e,negative_e,kinetic_exact;Words charge_C,kinetic_J;
    void add(double w,double q,double weighted_kinetic){
        ++count;weight=weight.credited(w);const auto qw=Words(w)*q;charge_C=charge_C+qw;kinetic_J=kinetic_J+Words(weighted_kinetic);
        need(std::isfinite(weighted_kinetic)&&weighted_kinetic>=0,"ledger kinetic energy operand");kinetic_exact=kinetic_exact.credited(weighted_kinetic);
        const double z=q/raikou::mfem_backend::pic_e,r=std::nearbyint(z);
        if(q!=0&&z==r&&std::abs(r)<=8&&r*raikou::mfem_backend::pic_e==q){for(int k=0;k<int(std::abs(r));++k){if(r>0)positive_e=positive_e.credited(w);else negative_e=negative_e.credited(w);}}
        else if(q!=0)++non_elementary;
    }
    void add(const ChargeTally&o){count+=o.count;non_elementary+=o.non_elementary;weight=exact_sum(weight,o.weight);kinetic_exact=exact_sum(kinetic_exact,o.kinetic_exact);
        positive_e=exact_sum(positive_e,o.positive_e);negative_e=exact_sum(negative_e,o.negative_e);charge_C=charge_C+o.charge_C;kinetic_J=kinetic_J+o.kinetic_J;}
    void encode(Digest&d)const{d.u64(count);d.u64(non_elementary);encode_exact(d,weight);encode_exact(d,positive_e);encode_exact(d,negative_e);encode_exact(d,kinetic_exact);
        d.f64(charge_C.hi);d.f64(charge_C.lo);d.f64(kinetic_J.hi);d.f64(kinetic_J.lo);}
    Json json()const{Json j=Json::object();j["count"]=count;j["weight_exact_carriers"]=exact_json(weight);j["weight"]=weight.zero()?0.:weight.rounded_for_rate();
        j["positive_elementary_charges_exact_carriers"]=exact_json(positive_e);j["negative_elementary_charges_exact_carriers"]=exact_json(negative_e);
        j["non_elementary_carriers"]=non_elementary;j["kinetic_J_exact_carriers"]=exact_json(kinetic_exact);Json c=Json::array();c.push_back(charge_C.hi);c.push_back(charge_C.lo);j["charge_C_words"]=std::move(c);
        Json k=Json::array();k.push_back(kinetic_J.hi);k.push_back(kinetic_J.lo);j["kinetic_J_words"]=std::move(k);return j;}
};
struct SinkLedger {
    std::array<ChargeTally,5> boundary;
    static SinkLedger of(std::span<const SinkRecord> sinks){SinkLedger l;for(const auto&s:sinks){
        need(s.boundary>=0&&s.boundary<5,"sink ledger boundary");l.boundary[s.boundary].add(s.weight,s.charge_C,s.weighted_kinetic_J);}return l;}
    void add(const SinkLedger&o){for(int b=0;b<5;++b)boundary[b].add(o.boundary[b]);}
    std::uint64_t count()const{std::uint64_t n=0;for(const auto&b:boundary)n+=b.count;return n;}
    Words collected_charge_C()const{return boundary[0].charge_C+boundary[1].charge_C+boundary[3].charge_C;}
    Words escaped_charge_C()const{return boundary[2].charge_C+boundary[4].charge_C;}
    Words collected_kinetic_J()const{return boundary[0].kinetic_J+boundary[1].kinetic_J+boundary[3].kinetic_J;}
    Words escaped_kinetic_J()const{return boundary[2].kinetic_J+boundary[4].kinetic_J;}
    void encode(Digest&d)const{d.text("hall-v3-sink-ledger-v1");for(const auto&b:boundary)b.encode(d);}
    std::string sha256()const{Digest d("hall-v3-sink-ledger-digest-v1");encode(d);return d.hex();}
    Json json()const{Json j=Json::object();for(int b=0;b<5;++b)j[sink_boundary_names[b]]=boundary[b].json();return j;}
};
// Newborn ledger per origin (lower-foil face, upper-foil face, interior point):
// count, exact weight, charge, and the birth kinetic energy fl(w KE) supplied by the source.
inline constexpr const char* newborn_origin_names[3]={"lower-foil-face","upper-foil-face","interior"};
struct BirthLedger {
    std::array<ChargeTally,3> origin;
    static BirthLedger of(std::span<const NewbornRecord> births){BirthLedger l;for(const auto&b:births){const int k=b.origin<0?2:b.origin;
        need(k>=0&&k<3,"newborn ledger origin");l.origin[k].add(b.weight,b.charge_C,b.birth_weighted_kinetic_J);}return l;}
    void add(const BirthLedger&o){for(int c=0;c<3;++c)origin[c].add(o.origin[c]);}
    std::uint64_t count()const{std::uint64_t n=0;for(const auto&c:origin)n+=c.count;return n;}
    Words kinetic_J()const{return origin[0].kinetic_J+origin[1].kinetic_J+origin[2].kinetic_J;}
    void encode(Digest&d)const{d.text("hall-v3-newborn-ledger-v1");for(const auto&c:origin)c.encode(d);}
    Json json()const{Json j=Json::object();for(int c=0;c<3;++c)j[newborn_origin_names[c]]=origin[c].json();return j;}
};
} // namespace hall_corrective_v3
