#pragma once
#include <openssl/bn.h>
#include <bit>
#include <cmath>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include "raikou/json.hpp"
namespace metal_transfer {
inline void need(bool x,const char* text){if(!x)throw std::invalid_argument(text);}
// Exact represented-binary currency, including products. The integer/exponent
// is canonical restart state; floating projections are diagnostics only.
class D {
 using BN=std::unique_ptr<BIGNUM,decltype(&BN_free)>;
 BN n_{BN_new(),BN_free};int e_=0;
 void normalize(){need(bool(n_),"BIGNUM allocation");if(BN_is_zero(n_.get())){e_=0;BN_set_negative(n_.get(),0);return;}while(!BN_is_odd(n_.get())){need(BN_rshift1(n_.get(),n_.get()),"normalize");++e_;}}
 public:
 D(){need(bool(n_),"BIGNUM allocation");BN_zero(n_.get());}
 D(const D& x):n_(BN_dup(x.n_.get()),BN_free),e_(x.e_){need(bool(n_),"BIGNUM copy");}
 D(D&&) noexcept=default;D& operator=(D&&) noexcept=default;
 D& operator=(const D& x){if(this!=&x){D t(x);*this=std::move(t);}return *this;}
 explicit D(double x):D(){need(std::isfinite(x),"nonfinite exact input");const auto bits=std::bit_cast<std::uint64_t>(x);const auto eb=(bits>>52)&2047;const auto m=(bits&((1ULL<<52)-1))|(eb?(1ULL<<52):0);need(BN_set_word(n_.get(),m),"set word");e_=eb?int(eb)-1023-52:-1074;BN_set_negative(n_.get(),(bits>>63)&&m);normalize();}
 friend D operator+(const D&a,const D&b){D z;z.e_=std::min(a.e_,b.e_);D x(a),y(b);need(BN_lshift(x.n_.get(),x.n_.get(),a.e_-z.e_)&&BN_lshift(y.n_.get(),y.n_.get(),b.e_-z.e_)&&BN_add(z.n_.get(),x.n_.get(),y.n_.get()),"exact add");z.normalize();return z;}
 friend D operator-(const D&a,const D&b){D y(b);if(!BN_is_zero(y.n_.get()))BN_set_negative(y.n_.get(),!BN_is_negative(y.n_.get()));return a+y;}
 friend D operator*(const D&a,const D&b){D z;std::unique_ptr<BN_CTX,decltype(&BN_CTX_free)>ctx(BN_CTX_new(),BN_CTX_free);need(bool(ctx)&&BN_mul(z.n_.get(),a.n_.get(),b.n_.get(),ctx.get()),"exact product");z.e_=a.e_+b.e_;z.normalize();return z;}
 int compare(const D& b)const{D diff=*this-b;return BN_is_zero(diff.n_.get())?0:BN_is_negative(diff.n_.get())?-1:1;}
 long double ld()const{if(BN_is_zero(n_.get()))return 0;D x(*this);const bool neg=BN_is_negative(x.n_.get());BN_set_negative(x.n_.get(),0);int shift=std::max(0,BN_num_bits(x.n_.get())-62);need(BN_rshift(x.n_.get(),x.n_.get(),shift),"project");long double v=std::ldexp(static_cast<long double>(BN_get_word(x.n_.get())),e_+shift);need(std::isfinite(v),"exact projection overflow");return neg?-v:v;}
 double value()const{double v=static_cast<double>(ld());need(std::isfinite(v),"double projection overflow");return v;}
 nlohmann::json json()const{char* p=BN_bn2hex(n_.get());need(p,"hex");std::string s(p);OPENSSL_free(p);return {{"hex",s},{"exponent",e_}};}
 static D read(const nlohmann::json& j){need(j.is_object()&&j.size()==2&&j.at("hex").is_string()&&j.at("exponent").is_number_integer(),"exact restart schema");std::string s=j.at("hex");need(!s.empty()&&s.size()<8192,"exact restart length");for(std::size_t i=0;i<s.size();++i)need((i==0&&s[i]=='-')||(s[i]>='0'&&s[i]<='9')||(s[i]>='A'&&s[i]<='F'),"exact restart digit");D z;BIGNUM* p=nullptr;need(BN_hex2bn(&p,s.c_str())==int(s.size()),"exact hex parse");z.n_.reset(p);z.e_=j.at("exponent");need(z.e_>=-8192&&z.e_<=8192,"exact exponent range");z.normalize();need(z.json()==j,"noncanonical exact restart");return z;}
};
} // namespace metal_transfer
