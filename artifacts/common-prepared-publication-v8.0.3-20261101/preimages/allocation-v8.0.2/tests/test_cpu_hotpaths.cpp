#include "water_parent_bridge.hpp"
#include "exact_product_error.hpp"
#include <bit>
#include <cfenv>
#include <cstring>
#include <iostream>
#include <random>
#include <typeinfo>

namespace w=raikou::pic_water;
using Json=nlohmann::json;
static std::uint64_t products=0,fast=0,fallback=0,queries=0,refusals=0,zero=0;
static void demand(bool b,const std::string& s){if(!b)throw std::runtime_error(s);}
static bool same(long double a,long double b){
    if(std::isnan(a)&&std::isnan(b))return true; // NaN classification, not payload, is the caller contract.
    if constexpr(LDBL_MANT_DIG==64&&sizeof(long double)==16)return std::memcmp(&a,&b,10)==0;
    return a==b&&std::signbit(a)==std::signbit(b);
}
static void product(long double a,long double b,bool wrong=false){
    volatile long double aa=a,bb=b;
    long double x=aa*bb;
    if(wrong&&std::isfinite(x))x=std::nextafter(x,std::numeric_limits<long double>::infinity());
    const bool shortcut=pf_exact::product_error_fast_domain(a,b,x);
    const long double expected=std::fma(a,b,-x),actual=pf_exact::product_error(a,b,x);
    demand(same(expected,actual),"exact product error differs");
    ++products;shortcut?++fast:++fallback;
}
static void product_checks(){
    const long double inf=std::numeric_limits<long double>::infinity();
    std::vector<long double> edges={0,-0.0L,1,-1,3,0x1p-4000L,0x1p4000L,
       std::nextafter(0x1p-4000L,0.0L),std::nextafter(0x1p-4000L,inf),
       std::nextafter(0x1p4000L,0.0L),std::nextafter(0x1p4000L,inf),
       std::numeric_limits<long double>::min(),std::numeric_limits<long double>::denorm_min(),
       std::numeric_limits<long double>::max(),inf,-inf,std::numeric_limits<long double>::quiet_NaN(),
       std::nextafter(1.0L,0.0L),std::nextafter(1.0L,inf),299792458,1e-300L,1e300L};
    const int original=std::fegetround();
    for(int rounding:{FE_TONEAREST,FE_DOWNWARD,FE_UPWARD,FE_TOWARDZERO}){
        demand(std::fesetround(rounding)==0,"rounding mode setup failed");
        if(rounding!=FE_TONEAREST)demand(!pf_exact::product_error_environment(),"nonnearest fast path enabled");
        for(auto a:edges)for(auto b:edges){product(a,b);product(-a,b);product(a,b,true);}
        demand(std::fegetround()==rounding,"helper changed rounding mode");
    }
    std::fesetround(FE_TONEAREST);
#if defined(__x86_64__) || defined(__i386__)
    unsigned short old;__asm__ volatile("fnstcw %0":"=m"(old));
    for(unsigned short precision:{0x0000u,0x0200u}){
        unsigned short control=(old&~0x0300u)|precision;
        __asm__ volatile("fldcw %0"::"m"(control));
        demand(!pf_exact::product_error_environment(),"reduced precision fast path enabled");
        product(1.234567890123456789L,3.456789012345678901L);
    }
    __asm__ volatile("fldcw %0"::"m"(old));
#endif
    std::mt19937_64 rng(20261031);
    for(std::uint64_t i=0;i<2000000;++i){
        auto operand=[&](){
            const long double sig=1+std::ldexp(static_cast<long double>(rng()>>1),-63);
            const int exponent=static_cast<int>(rng()%8201)-4100;
            return std::ldexp((rng()&1)?sig:-sig,exponent);
        };
        auto a=operand(),b=operand();
        if(i%4==0){a=static_cast<double>(a);b=static_cast<double>(b);}
        if(i%4==1)b=std::nextafter(a,inf)-a;
        product(a,b,i%31==0);
    }
    std::fesetround(original);
    demand(fast>100000&&fallback>10000,"missing fast/fallback coverage");
}
struct Outcome {std::uint64_t bits=0;std::string type,text;bool operator==(const Outcome&)const=default;};
template<class F>static Outcome query(F&& f,double E){
    try{return {std::bit_cast<std::uint64_t>(f(E)),"value",""};}
    catch(const std::exception& e){return {0,typeid(e).name(),e.what()};}
}
static Json& channel(Json& data){for(auto& c:data.at("channels"))if(c.at("id")==w::parent_channel)return c;throw std::runtime_error("no parent table");}
static void sigma_case(const std::string& name,const Json& raw,bool expect_cache){
    auto d=std::make_shared<const w::Data>(raw);
    w::ParentBridge bridge(d,w::parent_model,std::string(64,'a'),std::string(64,'b'),18.0106*1.66053906660e-27);
    auto cache=w::ParentBridge::parent_sigma_cache(*d);
    demand(bool(cache)==expect_cache,name+": unexpected cache admission");
    const auto actual=bridge.channel().sigma_m2;
    const auto original=[&](double E){const double s=d->observable_value_SI(w::parent_channel,E);if(s>0)d->require_MCC_event_kernel(w::parent_channel,w::parent_model,E);return s;};
    const auto& p=d->channel_metadata(w::parent_channel);
    const double inf=std::numeric_limits<double>::infinity();
    std::vector<double> energies={-inf,-1,-0.0,0,inf,std::numeric_limits<double>::quiet_NaN()};
    auto neighbors=[&](double E){energies.push_back(std::nextafter(E,-inf));energies.push_back(E);energies.push_back(std::nextafter(E,inf));};
    for(double E:p.at("support_eV").get<std::vector<double>>())neighbors(E);
    if(!p.at("zero_below_eV").is_null())neighbors(p.at("zero_below_eV").get<double>());
    if(p.at("model")=="linear_table")for(double E:p.at("energy_eV").get<std::vector<double>>())neighbors(E);
    std::mt19937_64 random(121);
    for(int i=0;i<10000;++i)energies.push_back(1100.*static_cast<double>(random()>>11)/9007199254740992.);
    for(double E:energies){
        auto a=query(actual,E),b=query(original,E);
        demand(a==b,name+": value/refusal/type/text differs at "+std::to_string(E)+": "+a.text+" / "+b.text);
        ++queries;if(a.type!="value")++refusals;else if((a.bits&0x7fffffffffffffffull)==0)++zero;
    }
}
static void sigma_checks(const std::string& path){
    auto original=w::Data::read(path).metadata();
    sigma_case("pinned",original,true);
    auto j=original;channel(j)["zero_below_eV"]=nullptr;sigma_case("no threshold",j,true);
    j=original;channel(j)["zero_below_eV"]=13.5;sigma_case("threshold at first point",j,true);
    j=original;for(auto& v:channel(j)["values"])v=0.;sigma_case("zero table",j,true);
    j=original;j["conditional_event_models"][0]["data_support_eV"]={20.,900.};sigma_case("narrow kernel domain",j,true);
    j=original;j["conditional_event_models"][0].erase("products");sigma_case("missing products",j,true);
    j=original;j["conditional_event_models"][0]["products"]=nullptr;sigma_case("null products",j,true);
    j=original;j["conditional_event_models"][0].erase("data_support_eV");sigma_case("missing kernel domain fallback",j,false);
    j=original;j["conditional_event_models"][0]["data_support_eV"]={"bad",1000.};sigma_case("malformed domain fallback",j,false);
    j=original;
    for(const auto& c:original.at("channels"))if(c.at("model")=="nist_BEB"){
        auto p=c;p["id"]=w::parent_channel;channel(j)=p;break;
    }
    demand(channel(j).at("model")=="nist_BEB","BEB fallback fixture unavailable");
    sigma_case("BEB fallback",j,false);
    demand(refusals>1000&&zero>10000,"missing cross-section refusal/zero coverage");
}
int main(int argc,char**argv){
    try{demand(argc==2,"dataset argument required");product_checks();sigma_checks(argv[1]);
        std::cout<<Json{{"passed",true},{"product_queries",products},{"product_fast",fast},{"product_fallback",fallback},
          {"sigma_queries",queries},{"sigma_refusals",refusals},{"sigma_zero",zero},{"sigma_fixture_count",10},
          {"comparison","exact finite/zero bits; exception dynamic type+text; NaN classification"}}.dump(2)<<'\n';return 0;
    }catch(const std::exception& e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}
}
