#pragma once
// Heat/source owner v7 numerical admission: a copy of the owner v6.1 admission
// (pic-electrode-heat-source-v6.1-20261031, sealed) in namespace electrode_heat_v7, so that v6.1
// and v7 types never meet in overload resolution when both are compiled into one program.
// Unchanged v6.1 text below (v6 admission of pic-electrode-heat-source-v6-20261030).
// Refusal texts, tolerances and the binary-rescaling test are unchanged. The only
// change is the speed of detail::Encoding: integers, doubles and double arrays are
// appended as little-endian 8-byte words in bulk instead of byte by byte; the byte
// sequence, and therefore every digest, is identical to v6.
#include <openssl/sha.h>
#include <array>
#include <cstring>
#include <type_traits>
#include <bit>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <new>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
namespace electrode_heat_v7 {

// T: largest admitted relative error of a scalar checked against an exact
// reference value (exposure expected count and its wall-heat diagnostic).
inline constexpr double relative_tolerance=0x1p-40;
// Bound on the error of the extended-precision admission comparison itself.
// A scalar is admitted only when its computed relative error is at most
// relative_tolerance-evaluation_bound, so its exact relative error is at most T.
inline constexpr double evaluation_bound=0x1p-58;
inline constexpr double admission_bound=relative_tolerance-evaluation_bound;
static_assert(std::numeric_limits<double>::is_iec559&&std::numeric_limits<double>::digits==53);
static_assert(std::numeric_limits<long double>::digits>=64&&std::numeric_limits<long double>::min_exponent<=-16000&&std::numeric_limits<long double>::max_exponent>=16000,
              "x87 extended long double required for exact-reference admission");

enum class RefusalKind {
 weighted_representation, // a weighted packet value is not the exact binary rescaling of its normal-range evaluation
 sealed_guard_at_weight,  // a sealed numerical guard failed only at the actual weight
 exposure_scalar,         // expected exposure count/diagnostic not represented within T of its finite-log value
 release_scalar,          // neutral release count not representable
 stock_drain,             // rounded release exceeds the exact accepted stock
 interval_bound,          // v7: step longer than the source-owner interval bound (API-v7.txt 3.2)
 surface_domain           // v7: surface state outside a law's domain (water above 3000 K, emission above T_max)
};
inline const char* kind_name(RefusalKind k){
 switch(k){
  case RefusalKind::weighted_representation:return "weighted-representation";
  case RefusalKind::sealed_guard_at_weight:return "sealed-guard-at-weight";
  case RefusalKind::exposure_scalar:return "exposure-scalar";
  case RefusalKind::release_scalar:return "release-scalar";
  case RefusalKind::stock_drain:return "stock-drain";
  case RefusalKind::interval_bound:return "interval-bound";
  case RefusalKind::surface_domain:return "surface-domain";
 }
 return "unknown";
}

// Typed scoped refusal. quantity names the unresolved quantity, log_abs_scope is
// log(|exact magnitude|) of that quantity, identity names the patch, stock and/or
// canonical particle, relative_error is the observed relative error (NaN when not
// applicable, +inf for a value lost entirely) and tolerance the admitted bound
// (0 for the exact rescaling test). For weighted_representation refusals,
// observed is the value obtained at the actual weight (or the rounded rescaling)
// and reference*2^binary_exponent the exact value it had to equal.
struct Refusal:std::invalid_argument {
 const std::string quantity;const double log_abs_scope;const std::string identity;const RefusalKind kind;
 const double relative_error,tolerance;const std::string detail;
 const double observed,reference;const int binary_exponent;
 Refusal(std::string q,double scope,std::string id,RefusalKind k,double rel=std::numeric_limits<double>::quiet_NaN(),double tol=0,std::string more={},
         double observed_value=std::numeric_limits<double>::quiet_NaN(),double reference_value=std::numeric_limits<double>::quiet_NaN(),int exponent=0)
  :std::invalid_argument(message(q,scope,id,k,rel,tol,more)),quantity(std::move(q)),log_abs_scope(scope),identity(std::move(id)),kind(k),relative_error(rel),tolerance(tol),detail(std::move(more)),
   observed(observed_value),reference(reference_value),binary_exponent(exponent){}
 static std::string message(const std::string&q,double x,const std::string&id,RefusalKind k,double rel,double tol,const std::string&more){
  std::ostringstream s;s<<std::setprecision(17)<<q<<" unresolved positive/nonzero contribution; log(abs(scope))="<<x<<"; kind="<<kind_name(k)<<"; identity="<<id<<"; relative_error="<<rel<<"; tolerance="<<tol;
  if(!more.empty())s<<"; "<<more;
  s<<"; whole proposal refused";return s.str();}
};

namespace detail {
struct Encoding {
 std::string bytes;
 static_assert(std::endian::native==std::endian::little,"v6.1 bulk encoding requires a little-endian host (v6 byte order)");
 void u(std::uint64_t x){char b[8];std::memcpy(b,&x,8);bytes.append(b,8);}
 void d(double x){u(std::bit_cast<std::uint64_t>(x));}
 void s(const std::string&x){u(x.size());bytes+=x;}
 // Contiguous double ranges are appended in one copy (same bytes as the v6 per-word loop).
 template<class V>void doubles(const V&v){u(v.size());
  if constexpr(requires{v.data();}&&std::is_same_v<std::remove_cv_t<std::remove_reference_t<decltype(*v.data())>>,double>){if(v.size())bytes.append(reinterpret_cast<const char*>(v.data()),v.size()*sizeof(double));}
  else for(double x:v)d(x);}
 std::string sha()const{unsigned char digest[SHA256_DIGEST_LENGTH];SHA256(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size(),digest);std::ostringstream s;for(auto c:digest)s<<std::hex<<std::setfill('0')<<std::setw(2)<<unsigned(c);return s.str();}
};
inline std::string hex(double x){std::ostringstream s;s<<std::hexfloat<<x;return s.str();}

// v4 zero/nonfinite admission of a weighted product, now with identity.
inline void require_product(double factor,double weight,double actual,const std::string&q,const std::string&id){
 if(!std::isfinite(actual))throw Refusal(q,std::log(std::abs(factor))+std::log(weight),id,RefusalKind::weighted_representation,std::numeric_limits<double>::quiet_NaN());
 if(factor!=0&&actual==0)throw Refusal(q,std::log(std::abs(factor))+std::log(weight),id,RefusalKind::weighted_representation,INFINITY);
}
inline void require_log_product(double log_abs,double actual,const std::string&q,const std::string&id){
 if(!std::isfinite(actual)||actual==0)throw Refusal(q,log_abs,id,RefusalKind::weighted_representation,actual==0?INFINITY:std::numeric_limits<double>::quiet_NaN());
}

// w = tilde*2^e with tilde in [1,2): the binade-normalized weight.
struct Binade {double w=1,tilde=1;int e=0;};
inline Binade binade(double w){
 if(!(std::isfinite(w)&&w>0))throw std::invalid_argument("v6 binade requires a positive finite weight");
 const int e=std::ilogb(w);return {w,std::scalbn(w,-e),e};
}
// out=x*2^e; true iff that value is a finite binary64 number exactly (no rounding).
inline bool rescale(double x,int e,double&out){out=std::scalbn(x,e);return std::isfinite(out)&&std::scalbn(out,-e)==x;}
// |actual-x*2^e|/|x*2^e| in extended precision (report only).
inline double scaled_relative_error(double actual,double x,int e){
 if(!std::isfinite(actual)||!std::isfinite(x))return std::numeric_limits<double>::quiet_NaN();
 if(x==0)return actual==0?0:INFINITY;
 const long double exact=std::ldexp(static_cast<long double>(x),e);
 return double(std::fabs((static_cast<long double>(actual)-exact)/exact));
}

// Quantity names for sealed numerical guards that can fail only at the actual
// weight once the binade-normalized evaluation (identical inputs, weight in
// [1,2)) has passed. The sealed text is pinned by the prerequisite hashes.
inline std::string guard_quantity(const std::string&m){
 struct G{const char*prefix;const char*quantity;};
 static constexpr G table[]={
  {"electrode heat: source incremental aggregate energy failed","source incremental aggregate energy closure J"},
  {"Micro-layer aggregate energy closure failed","micro-layer aggregate energy closure J"},
  {"electrode heat: source path/common-stage work unresolved","source path work closure J"},
  {"electrode heat: actual handoff remainder work/continuity failed","ordinary handoff work/continuity closure"},
  {"Invalid emitted phase-space or energy ledger","emitted energy ledger J"},
  {"electrode heat: packet","neutral packet kinematics"}};
 for(const auto&g:table)if(m.rfind(g.prefix,0)==0)return g.quantity;
 return "weighted sealed numerical guard";
}

// Exact binary-rescaling admission of the weighted values of one evaluation.
// unit: the same quantity at weight 1 (v4 diagnostic, used for the log scope);
// tilde: the same quantity at weight w~=w*2^-e in [1,2); actual: at weight w.
// predicted(): before the actual-weight evaluation, w~-value*2^e must be a finite
// binary64 number exactly. consistent(): the actual value must equal it exactly.
// In the normal range every operation of the weighted kernels is homogeneous in
// the weight, so both tests pass for every weight whose evaluation neither
// underflows nor overflows; any rounding caused by the weight's magnitude fails.
class Weighted {
 Binade b_;std::string identity_;
 [[noreturn]] void refuse(const std::string&q,double unit,double tilde,double observed,const std::string&more)const{
  throw Refusal(q,scope(unit,tilde),identity_,RefusalKind::weighted_representation,scaled_relative_error(observed,tilde,b_.e),0,more,observed,tilde,b_.e);}
public:
 Weighted(double w,std::string identity):b_(binade(w)),identity_(std::move(identity)){}
 const Binade&scale()const{return b_;}double w()const{return b_.w;}double tilde()const{return b_.tilde;}bool trivial()const{return b_.e==0;}
 const std::string&identity()const{return identity_;}
 double scope(double unit,double tilde)const{
  if(std::isfinite(unit)&&unit!=0)return std::log(std::abs(unit))+std::log(b_.w);
  if(std::isfinite(tilde)&&tilde!=0)return std::log(std::abs(tilde))-std::log(b_.tilde)+std::log(b_.w);
  return std::log(b_.w);}
 void predicted(const std::string&q,double unit,double tilde,const std::string&more={})const{double p;if(!rescale(tilde,b_.e,p))refuse(q,unit,tilde,p,more);}
 void consistent(const std::string&q,double unit,double tilde,double actual,const std::string&more={})const{double p;if(!(rescale(tilde,b_.e,p)&&actual==p))refuse(q,unit,tilde,actual,more);}
 template<class V>void predicted_array(const std::string&q,const V&tilde)const{
  for(std::size_t i=0;i<tilde.size();++i){double p;if(!rescale(tilde[i],b_.e,p))refuse(q,std::numeric_limits<double>::quiet_NaN(),tilde[i],p,"entry "+std::to_string(i)+" of "+std::to_string(tilde.size()));}}
 template<class V>void consistent_array(const std::string&q,const V&tilde,const V&actual)const{
  if(tilde.size()!=actual.size())throw Refusal(q,std::log(b_.w),identity_,RefusalKind::weighted_representation,std::numeric_limits<double>::quiet_NaN(),0,"array dimensions differ between weights");
  for(std::size_t i=0;i<tilde.size();++i){double p;if(!(rescale(tilde[i],b_.e,p)&&actual[i]==p))refuse(q,std::numeric_limits<double>::quiet_NaN(),tilde[i],actual[i],"entry "+std::to_string(i)+" of "+std::to_string(tilde.size()));}}
 // A weight-independent value (geometry, time, identity) must be identical at both weights.
 void same(const std::string&q,bool equal)const{if(!equal)throw Refusal(q,std::log(b_.w),identity_,RefusalKind::weighted_representation,std::numeric_limits<double>::quiet_NaN(),0,"weight-independent value differs between weights");}
 // Run the actual-weight sealed evaluation after its binade-normalized evaluation passed:
 // a sealed numerical guard failing now fails only because of the weight's magnitude.
 // With a trivial binade (w in [1,2)) no separate evaluation exists and nothing is converted.
 template<class F>auto at_weight(F&&f,double unit_scale)const->decltype(f()){
  if(trivial())return f();
  try{return f();}
  catch(const Refusal&){throw;}
  catch(const std::bad_alloc&){throw;}
  catch(const std::exception&x){
   const double s=std::isfinite(unit_scale)&&unit_scale>0?unit_scale:1;
   throw Refusal(guard_quantity(x.what()),std::log(s)+std::log(b_.w),identity_,RefusalKind::sealed_guard_at_weight,std::numeric_limits<double>::quiet_NaN(),0,std::string("sealed guard: ")+x.what());}
 }
};
}
}
