// Source owner v6.1: exact fast evaluation of the sealed thermal-field emission kernel.
// See fast_emission.hpp. The helpers in the anonymous namespace below are verbatim copies
// of the sealed kernel's internal helpers (heat-v2 binding source-kernel
// src/surface_emission.cpp): same expressions, same operation order, same libm calls.
// Compile with the sealed kernel's flags (-O2 -ffp-contract=off -fno-fast-math).
#include "fast_emission.hpp"
#include <algorithm>
#include <array>
#include <bit>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <deque>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>
#include <immintrin.h>

namespace electrode_heat_v6_1::emission {
namespace se = surface_emission;
namespace {
using se::eps0;
using se::kb;
using se::me;
using se::planck;
using se::qe;
// ---- verbatim copies of the sealed kernel's internal constants and helpers
constexpr double pi=3.1415926535897932384626433832795;
constexpr double hbar=planck/(2*pi),kb_eV=kb/qe,image_c=qe/(16*pi*eps0);
void require(bool okay,const char* message){if(!okay)throw std::invalid_argument(message);}
void config(const se::Config& c) {
    require(std::isfinite(c.work_function_eV) && c.work_function_eV>=2 && c.work_function_eV<=6,"Explicit work function outside numerical domain2..6eV");
    require(c.energy_panels>=16 && c.energy_panels<=1024 && c.action_nodes>=16 && c.action_nodes<=256 && c.action_nodes%2==0,
            "Bounded quadrature order required");
    require(std::isfinite(c.tail_widths) && c.tail_widths>=32 && c.tail_widths<=80,"Tail width outside32..80");
}
// Runtime evaluation as in the sealed object (legendre is a separate, non-inlined function there).
[[gnu::noinline,gnu::noipa]] std::pair<std::vector<double>,std::vector<double>> legendre(int n) {
    std::vector<double> x(n),w(n);
    for(int i=0;i<(n+1)/2;++i) {
        double z=std::cos(pi*(i+.75)/(n+.5)),derivative=0;
        for(int j=0;j<100;++j) {
            double p=1,q=0;
            for(int k=1;k<=n;++k){const double next=((2*k-1)*z*p-(k-1)*q)/k;q=p;p=next;}
            derivative=n*(z*p-q)/(z*z-1);const double next=z-p/derivative;
            if(std::abs(next-z)<2e-15){z=next;break;}z=next;
        }
        // Re-evaluate derivative at converged root.
        double p=1,q=0;for(int k=1;k<=n;++k){const double next=((2*k-1)*z*p-(k-1)*q)/k;q=p;p=next;}
        derivative=n*(z*p-q)/(z*z-1);
        x[i]=-z;x[n-1-i]=z;w[i]=w[n-1-i]=2/((1-z*z)*derivative*derivative);
    }
    return {x,w};
}
double softplus(double x){if(x>40)return x+std::log1p(std::exp(-x));return std::log1p(std::exp(x));}
double log_transmission(double G){return -softplus(G);}
double log_supply(double E,double kT) {
    if(kT==0)return E<0?std::log(-E):-std::numeric_limits<double>::infinity();
    const double eta=-E/kT;
    return std::log(kT)+(eta< -36?eta:std::log(softplus(eta)));
}
double fermi1_negative(double eta) {
    const double z=1/(1+std::exp(-eta)),L=softplus(eta);double series=0,power=z;
    for(int k=1;k<=100;++k){const double term=power/(double(k)*k);series+=term;if(term<1e-17*std::max(series,1e-300))break;power*=z;}
    return series+.5*L*L;
}
double mean_transverse(double E,double kT) {
    if(kT==0)return -.5*E;
    const double eta=-E/kT;
    if(eta< -36)return kT;
    const double L=softplus(eta),f1=eta<=0?fermi1_negative(eta):.5*eta*eta+pi*pi/6-fermi1_negative(-eta);
    return kT*f1/L;
}

// ---- memoized rules and action tables (computed by the code above, once per order)
struct Rule {std::vector<double> x,w;};
struct ActionTable {
    std::vector<double> cosine,sine2,weights,ws;   // ws[i] = weights[i]*sine2[i] (the sealed term numerator)
    double ws_lo=0,ws_hi=0;                        // bounds of the real sum of ws
    bool bounds_ok=false;
};
std::mutex cache_mutex;
const Rule& rule(int n) {
    static std::map<int,std::unique_ptr<const Rule>> cache;
    std::lock_guard<std::mutex> lock(cache_mutex);
    auto& slot=cache[n];
    if(!slot){auto [x,w]=legendre(n);slot=std::make_unique<const Rule>(Rule{std::move(x),std::move(w)});}
    return *slot;
}
const ActionTable& action_table(int n) {
    static std::map<int,std::unique_ptr<const ActionTable>> cache;
    {std::lock_guard<std::mutex> lock(cache_mutex);if(auto it=cache.find(n);it!=cache.end())return *it->second;}
    const auto [x,w]=legendre(n);
    auto t=std::make_unique<ActionTable>();
    for(std::size_t i=0;i<x.size();++i){const double angle=pi*(x[i]+1)/2;t->cosine.push_back(std::cos(angle));t->sine2.push_back(std::pow(std::sin(angle),2));t->weights.push_back(w[i]*pi/2);}
    long double sum=0;bool ok=true;
    for(std::size_t i=0;i<x.size();++i){const double v=t->weights[i]*t->sine2[i];t->ws.push_back(v);ok=ok&&std::isfinite(v)&&v>=0&&(v==0||v>=1e-200)&&t->cosine[i]<=1&&t->cosine[i]>=-1;sum+=v;}
    t->ws_lo=double(sum)*(1-0x1p-40);t->ws_hi=double(sum)*(1+0x1p-40);
    t->bounds_ok=ok&&std::isfinite(t->ws_hi)&&t->ws_lo>0;
    std::lock_guard<std::mutex> lock(cache_mutex);
    auto& slot=cache[n];
    if(!slot)slot=std::move(t);
    return *slot;
}

// ---- barrier (sealed Barrier members; tables from the cache)
struct Barrier {
    double F,phi,lowering,top,top_slope;
    const ActionTable* table=nullptr;
    Barrier(double field,const se::Config& c):F(field),phi(c.work_function_eV),lowering(2*std::sqrt(image_c*field)),top(phi-lowering) {
        const double z=std::sqrt(image_c/F),omega=std::sqrt(2*qe*image_c/(me*z*z*z));
        top_slope=2*pi*qe/(hbar*omega);
        table=&action_table(c.action_nodes);
    }
};
struct Chord {double mid,half;};
inline Chord chord(const Barrier& b,double E) {
    const double A=b.phi-E,disc=std::sqrt((A-b.lowering)*(A+b.lowering));
    const double z1=2*image_c/(A+disc),z2=(A+disc)/(2*b.F),mid=(z1+z2)/2,half=(z2-z1)/2;
    return {mid,half};
}
inline double action_from_integral(const Barrier& b,double half,long double integral) {
    const double G=2*std::sqrt(2*me*qe)/hbar*std::sqrt(b.F)*half*half*double(integral);
    require(std::isfinite(G) && G>=0,"Nonfinite Schottky-Nordheim action");return G;
}
// Sealed gamow, scalar.
double gamow(const Barrier& b,double E) {
    if(E>=b.top)return -b.top_slope*(E-b.top);
    const auto [mid,half]=chord(b,E);
    const auto& t=*b.table;
    long double integral=0;
    for(std::size_t i=0;i<t.ws.size();++i)integral+=t.ws[i]/std::sqrt(mid+half*t.cosine[i]);
    return action_from_integral(b,half,integral);
}
// Binary64 action terms ws[i]/sqrt(mid+half*cos[i]) with vector IEEE operations.
[[gnu::target("avx2")]] void terms_avx2(const double* ws,const double* cs,std::size_t n,double mid,double half,double* out) {
    const __m256d m=_mm256_set1_pd(mid),h=_mm256_set1_pd(half);
    std::size_t i=0;
    for(;i+4<=n;i+=4)_mm256_storeu_pd(out+i,_mm256_div_pd(_mm256_loadu_pd(ws+i),_mm256_sqrt_pd(_mm256_add_pd(m,_mm256_mul_pd(h,_mm256_loadu_pd(cs+i))))));
    for(;i<n;++i)out[i]=ws[i]/std::sqrt(mid+half*cs[i]);
}
void terms_sse2(const double* ws,const double* cs,std::size_t n,double mid,double half,double* out) {
    const __m128d m=_mm_set1_pd(mid),h=_mm_set1_pd(half);
    std::size_t i=0;
    for(;i+2<=n;i+=2)_mm_storeu_pd(out+i,_mm_div_pd(_mm_loadu_pd(ws+i),_mm_sqrt_pd(_mm_add_pd(m,_mm_mul_pd(h,_mm_loadu_pd(cs+i))))));
    for(;i<n;++i)out[i]=ws[i]/std::sqrt(mid+half*cs[i]);
}
bool have_avx2(){static const bool v=(__builtin_cpu_init(),__builtin_cpu_supports("avx2")!=0);return v;}
// Sealed gamow of up to four below-top energies; each long-double sum in the sealed term order.
void gamow4(const Barrier& b,const double* E,int count,double* G) {
    const auto& t=*b.table;const std::size_t n=t.ws.size();
    alignas(32) double terms[4][256];
    double half[4]={0,0,0,0};
    for(int k=0;k<count;++k){
        const auto c=chord(b,E[k]);half[k]=c.half;
        if(have_avx2())terms_avx2(t.ws.data(),t.cosine.data(),n,c.mid,c.half,terms[k]);
        else terms_sse2(t.ws.data(),t.cosine.data(),n,c.mid,c.half,terms[k]);
    }
    long double s0=0,s1=0,s2=0,s3=0;
    if(count==4)for(std::size_t i=0;i<n;++i){s0+=terms[0][i];s1+=terms[1][i];s2+=terms[2][i];s3+=terms[3][i];}
    else for(std::size_t i=0;i<n;++i){s0+=terms[0][i];if(count>1)s1+=terms[1][i];if(count>2)s2+=terms[2][i];}
    const long double s[4]={s0,s1,s2,s3};
    for(int k=0;k<count;++k)G[k]=action_from_integral(b,half[k],s[k]);
}
// Rigorous upper bound of a below-top node's log weight, or +inf when it cannot be shown
// that the node's sealed action G is finite, nonnegative and above 40.
double upper_log_weight(const Barrier& b,double E,double pw_up,double supply_up) {
    const auto& t=*b.table;
    if(!t.bounds_ok)return INFINITY;
    const auto [mid,half]=chord(b,E);
    if(!(std::isfinite(mid)&&std::isfinite(half)))return INFINITY;
    const double xa=mid+half,xb=mid-half,xmax=std::max(xa,xb),xmin=std::min(xa,xb);
    if(!(std::isfinite(xmax)&&xmin>=DBL_MIN))return INFINITY;
    const double rmax=std::sqrt(xmax),rmin=std::sqrt(xmin);
    // K: the first three factors of the sealed product, rounded as there.
    const double K=2*std::sqrt(2*me*qe)/hbar*std::sqrt(b.F)*half*half;
    const double glo=K*(t.ws_lo/rmax)*(1-0x1p-40),ghi=K*(t.ws_hi/rmin)*(1+0x1p-40);
    if(!(std::isfinite(K)&&K>=0&&std::isfinite(ghi)&&ghi<1e300&&glo>40))return INFINITY;
    return pw_up+supply_up-glo+1e-6;
}

// ---- configuration keys, statistics, retarding memo, self-check
using Key=std::tuple<std::uint64_t,int,int,std::uint64_t>;
Key key_of(const se::Config& c){return {std::bit_cast<std::uint64_t>(c.work_function_eV),c.energy_panels,c.action_nodes,std::bit_cast<std::uint64_t>(c.tail_widths)};}
std::mutex state_mutex;
Statistics stats;
bool fast_enabled=true;
struct ConfigState {bool checked=false,fast=false;bool runtime_richardson=false;};
std::map<Key,ConfigState> configs;
struct RetardEntry {std::uint64_t T;Key key;std::shared_ptr<const se::Distribution> d;};
std::deque<RetardEntry> retarding;
constexpr std::size_t retarding_capacity=512;

double richardson_log(bool runtime) {
    if(runtime){volatile double r=se::ideal_richardson_A_m2_K2();return std::log(r/(kb_eV*kb_eV));}
    return std::log(4*pi*me*qe*kb*kb/(planck*planck*planck)/(kb_eV*kb_eV));
}

// Extracting field, T > 0: the sealed quadrature with exactly-zero nodes dropped.
se::Distribution extracting(double T,double F,const se::Config& c,bool runtime_richardson,Statistics& local) {
    const double lowering=se::schottky_lowering_eV(F);
    se::Distribution d{};d.config=c;d.T_K=T;d.extraction_field_V_m=F;d.barrier_lowering_eV=lowering;
    const double kT=kb_eV*T;
    const double decay=F/(1.5*6.830890e9*std::sqrt(c.work_function_eV));
    const double low=-std::max({2.,c.tail_widths*decay,c.tail_widths*kT});
    const double high=c.work_function_eV+c.tail_widths*kT;
    const Barrier barrier(F,c);
    d.fermi_decay_width_eV=-2e-5/(gamow(barrier,1e-5)-gamow(barrier,-1e-5));
    std::vector<double> breaks{low,high};
    {
        const double top=c.work_function_eV-lowering,top_width=1/barrier.top_slope;
        for(double scale:{-48.,-8.,-1.,0.,1.,8.,48.})for(double value:{scale*decay,top+scale*top_width})
            if(value>low && value<high)breaks.push_back(value);
    }
    std::sort(breaks.begin(),breaks.end());breaks.erase(std::unique(breaks.begin(),breaks.end()),breaks.end());
    const auto& r8=rule(8);const auto& x=r8.x;const auto& w=r8.w;
    const std::size_t Q=x.size(),P=std::size_t(c.energy_panels),N=(breaks.size()-1)*P*Q;
    const double w_max=*std::max_element(w.begin(),w.end());
    // state: 0 exact, 1 pending (bounded), 2 dropped
    thread_local std::vector<double> energy,logw,upper;thread_local std::vector<unsigned char> state;
    energy.assign(N,0);logw.assign(N,-INFINITY);upper.assign(N,INFINITY);state.assign(N,0);
    // pw(i) = log((b-a)/2*w[q]) and supply(i), exactly as sealed, for the node index i.
    auto exact_weight=[&](std::size_t i,double logD){
        const std::size_t q=i%Q,cell=(i/Q)%P,segment=i/(Q*P)+1;
        const double a=breaks[segment-1]+(breaks[segment]-breaks[segment-1])*int(cell)/c.energy_panels;
        const double b=breaks[segment-1]+(breaks[segment]-breaks[segment-1])*(int(cell)+1)/c.energy_panels;
        const double supply=log_supply(energy[i],kT);
        require(std::isfinite(supply),"v6.1 fast emission: nonfinite supply at positive temperature");
        return std::log((b-a)/2*w[q])+supply+logD;
    };
    double M=-INFINITY;std::size_t pending=0;
    for(std::size_t segment=1;segment<breaks.size();++segment)for(int cell=0;cell<c.energy_panels;++cell) {
        const double a=breaks[segment-1]+(breaks[segment]-breaks[segment-1])*cell/c.energy_panels;
        const double b=breaks[segment-1]+(breaks[segment]-breaks[segment-1])*(cell+1)/c.energy_panels;
        bool have_cell_bound=false;double pw_up=0,supply_up=0;
        for(std::size_t q=0;q<Q;++q) {
            const std::size_t i=((segment-1)*P+std::size_t(cell))*Q+q;
            const double E=(a+b)/2+(b-a)/2*x[q];
            energy[i]=E;
            if(E>=barrier.top){logw[i]=exact_weight(i,log_transmission(-barrier.top_slope*(E-barrier.top)));M=std::max(M,logw[i]);continue;}
            if(!have_cell_bound){
                have_cell_bound=true;
                pw_up=std::log((b-a)/2*w_max)+1e-9;
                const double e_lo=a-(std::abs(a)+std::abs(b))*0x1p-50-1e-300;
                supply_up=log_supply(e_lo,kT)+1e-9;
            }
            upper[i]=upper_log_weight(barrier,E,pw_up,supply_up);
            state[i]=1;++pending;
        }
    }
    local.nodes+=N;
    if(pending) {
        // Seed M with the exact weight of the pending node with the largest bound.
        std::size_t best=N;
        for(std::size_t i=0;i<N;++i)if(state[i]==1&&(best==N||upper[i]>upper[best]))best=i;
        double G=0;gamow4(barrier,&energy[best],1,&G);
        logw[best]=exact_weight(best,log_transmission(G));state[best]=0;M=std::max(M,logw[best]);
        const double tau=M-800;
        thread_local std::vector<std::size_t> todo;todo.clear();
        for(std::size_t i=0;i<N;++i)if(state[i]==1){if(upper[i]<tau)state[i]=2;else todo.push_back(i);}
        for(std::size_t k=0;k<todo.size();k+=4) {
            const int count=int(std::min<std::size_t>(4,todo.size()-k));
            double Es[4]={0,0,0,0},Gs[4]={0,0,0,0};
            for(int j=0;j<count;++j)Es[j]=energy[todo[k+j]];
            gamow4(barrier,Es,count,Gs);
            for(int j=0;j<count;++j){const std::size_t i=todo[k+j];logw[i]=exact_weight(i,log_transmission(Gs[j]));state[i]=0;}
        }
    }
    double largest=-INFINITY;bool any=false;
    for(std::size_t i=0;i<N;++i)if(state[i]==0){largest=any?std::max(largest,logw[i]):logw[i];any=true;}
    // Every dropped node is below largest-800 (tau <= largest-800).
    for(std::size_t i=0;i<N;++i)if(state[i]==2){require(upper[i]<largest-800,"v6.1 fast emission: dropped node bound violated");++local.nodes_bounded_zero;}
    require(N>0,"Empty numerical emission quadrature");
    long double sum=0;
    for(std::size_t i=0;i<N;++i)if(state[i]==0){sum+=std::exp(logw[i]-largest);++local.nodes_exact;}
    const double logIntegral=largest+std::log(double(sum));
    d.log_current_A_m2=richardson_log(runtime_richardson)+logIntegral;d.has_log_current=true;
    d.current_A_m2=std::exp(d.log_current_A_m2);d.current_underflow=d.current_A_m2==0;
    require(std::isfinite(d.current_A_m2),"Emission current overflow");
    long double cdf=0,moment=0;
    for(std::size_t i=0;i<N;++i) {
        if(state[i]!=0)continue;
        const double p=std::exp(logw[i]-logIntegral);if(p==0)continue;
        const double transverse=mean_transverse(energy[i],kT);cdf+=p;moment+=static_cast<long double>(p)*(energy[i]+transverse);
        d.nodes.push_back({energy[i],p,double(cdf),transverse});
    }
    require(!d.nodes.empty() && std::abs(double(cdf)-1)<1e-10,"Spectrum normalization failed");d.nodes.back().cdf=1;
    d.mean_metal_energy_eV=double(moment);d.wall_heat_W_m2=-d.current_A_m2*d.mean_metal_energy_eV;
    require(std::isfinite(d.mean_metal_energy_eV) && std::isfinite(d.wall_heat_W_m2),"Emission energy moment overflow");return d;
}
bool same(double a,double b){return std::bit_cast<std::uint64_t>(a)==std::bit_cast<std::uint64_t>(b);}
// Self-check of one configuration against the sealed evaluate (reference fields inside its domain).
ConfigState check_configuration(const se::Config& c) {
    ConfigState s;s.checked=true;
    Statistics scratch;
    for(bool runtime:{false,true}) {
        bool ok=true;int compared=0;
        for(double T:{300.,1000.})for(double F:{3e2,2e6,3e7}) {
            se::Distribution reference;
            try{reference=se::evaluate(T,F,c);}catch(const std::exception&){continue;}
            try{ok=ok&&same_distribution(extracting(T,F,c,runtime,scratch),reference);}catch(const std::exception&){ok=false;}
            ++compared;
        }
        if(ok&&compared>0){s.fast=true;s.runtime_richardson=runtime;return s;}
        if(same(richardson_log(false),richardson_log(true)))break;
    }
    return s;
}
}  // namespace

bool same_distribution(const se::Distribution& a,const se::Distribution& b) noexcept {
    const auto& x=a.config;const auto& y=b.config;
    if(!(same(x.work_function_eV,y.work_function_eV)&&x.energy_panels==y.energy_panels&&x.action_nodes==y.action_nodes&&same(x.tail_widths,y.tail_widths)))return false;
    if(!(same(a.T_K,b.T_K)&&same(a.extraction_field_V_m,b.extraction_field_V_m)&&same(a.current_A_m2,b.current_A_m2)&&same(a.log_current_A_m2,b.log_current_A_m2)&&
         a.has_log_current==b.has_log_current&&a.current_underflow==b.current_underflow&&same(a.mean_metal_energy_eV,b.mean_metal_energy_eV)&&
         same(a.wall_heat_W_m2,b.wall_heat_W_m2)&&same(a.barrier_lowering_eV,b.barrier_lowering_eV)&&same(a.fermi_decay_width_eV,b.fermi_decay_width_eV)))return false;
    if(a.nodes.size()!=b.nodes.size())return false;
    for(std::size_t i=0;i<a.nodes.size();++i){const auto& p=a.nodes[i];const auto& q=b.nodes[i];
        if(!(same(p.normal_energy_eV,q.normal_energy_eV)&&same(p.probability,q.probability)&&same(p.cdf,q.cdf)&&same(p.mean_transverse_eV,q.mean_transverse_eV)))return false;}
    return true;
}
se::Distribution evaluate_sealed(double T,double F,const se::Config& c){return se::evaluate(T,F,c);}
se::Distribution evaluate(double T,double F,const se::Config& c) {
    // Sealed domain checks first, in the sealed order (same refusals).
    config(c);require(std::isfinite(T) && (T==0 || (T>=1 && T<=3000)),"Numerical free-electron law domain is T=0 or1..3000K; no extrapolation/clamp");
    require(std::isfinite(F) && std::abs(F)<=5e9 && (F<=0 || F>=1),"Signed normal-field numerical domain exceeded");
    const double lowering=F>0?se::schottky_lowering_eV(F):0;
    require(lowering<c.work_function_eV*.98,"Barrier collapse outside admitted numerical domain");
    const Key key=key_of(c);
    bool enabled,runtime=false;
    {
        std::lock_guard<std::mutex> lock(state_mutex);
        ++stats.calls;enabled=fast_enabled&&T>0;
    }
    if(!enabled){std::lock_guard<std::mutex> lock(state_mutex);++stats.sealed_calls;}
    if(!enabled)return se::evaluate(T,F,c);
    if(F<=0) {
        const std::uint64_t Tb=std::bit_cast<std::uint64_t>(T);
        std::shared_ptr<const se::Distribution> hit;
        {std::lock_guard<std::mutex> lock(state_mutex);
         for(const auto& e:retarding)if(e.T==Tb&&e.key==key){hit=e.d;break;}
         if(hit)++stats.retarding_reused;else ++stats.retarding_sealed;}
        if(hit){se::Distribution d=*hit;d.extraction_field_V_m=F;return d;}
        auto d=se::evaluate(T,F,c);
        auto stored=std::make_shared<se::Distribution>(d);stored->extraction_field_V_m=0;
        std::lock_guard<std::mutex> lock(state_mutex);
        retarding.push_back({Tb,key,std::move(stored)});
        if(retarding.size()>retarding_capacity)retarding.pop_front();
        return d;
    }
    ConfigState s;
    {std::lock_guard<std::mutex> lock(state_mutex);if(auto it=configs.find(key);it!=configs.end())s=it->second;}
    if(!s.checked){
        s=check_configuration(c);
        std::lock_guard<std::mutex> lock(state_mutex);
        configs[key]=s;++stats.self_checks;if(!s.fast)++stats.sealed_fallback_configs;
    }
    if(!s.fast){std::lock_guard<std::mutex> lock(state_mutex);++stats.sealed_calls;return se::evaluate(T,F,c);}
    runtime=s.runtime_richardson;
    Statistics local;
    auto d=extracting(T,F,c,runtime,local);
    std::lock_guard<std::mutex> lock(state_mutex);
    ++stats.extracting_fast;stats.nodes+=local.nodes;stats.nodes_exact+=local.nodes_exact;stats.nodes_bounded_zero+=local.nodes_bounded_zero;
    return d;
}
Statistics statistics() noexcept {std::lock_guard<std::mutex> lock(state_mutex);return stats;}
void reset_statistics() noexcept {std::lock_guard<std::mutex> lock(state_mutex);stats=Statistics{};}
void set_fast_path(bool enabled) noexcept {std::lock_guard<std::mutex> lock(state_mutex);fast_enabled=enabled;}
}  // namespace electrode_heat_v6_1::emission
