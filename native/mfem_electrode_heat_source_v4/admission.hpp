#pragma once
#include <openssl/sha.h>
#include <bit>
#include <stdexcept>
#include <sstream>
#include <iomanip>
#include <limits>
#include <cmath>
namespace electrode_heat_v4 {
struct Refusal:std::invalid_argument {
 const std::string quantity;const double log_abs_scope;
 Refusal(std::string q,double scope):std::invalid_argument(message(q,scope)),quantity(std::move(q)),log_abs_scope(scope){}
 static std::string message(const std::string&q,double x){std::ostringstream s;s<<std::setprecision(17)<<q<<" unresolved positive/nonzero contribution; log(abs(scope))="<<x<<"; whole proposal refused";return s.str();}
};
namespace detail {
struct Encoding {
 std::string bytes;
 void u(std::uint64_t x){for(int i=0;i<8;++i)bytes.push_back(char((x>>(8*i))&255));}
 void d(double x){u(std::bit_cast<std::uint64_t>(x));}
 void s(const std::string&x){u(x.size());bytes+=x;}
 template<class V>void doubles(const V&v){u(v.size());for(double x:v)d(x);}
 std::string sha()const{unsigned char digest[SHA256_DIGEST_LENGTH];SHA256(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size(),digest);std::ostringstream s;for(auto c:digest)s<<std::hex<<std::setfill('0')<<std::setw(2)<<unsigned(c);return s.str();}
};
inline void require_product(double factor,double weight,double actual,const std::string&q){
 if(!std::isfinite(actual))throw Refusal(q,std::log(std::abs(factor))+std::log(weight));
 if(factor!=0&&actual==0)throw Refusal(q,std::log(std::abs(factor))+std::log(weight));
}
inline void require_log_product(double log_abs,double actual,const std::string&q){if(!std::isfinite(actual)||actual==0)throw Refusal(q,log_abs);}
}
}
