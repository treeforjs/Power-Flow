#pragma once
#include "native_heating_fields.hpp"
#include <numeric>
#include <limits>
namespace validation_native {
constexpr double eps=8.8541878188e-12;
struct Voltage {
    std::vector<double> start,end,value;
    int current=0;
    explicit Voltage(const std::filesystem::path& path) {
        std::ifstream in(path);std::string line;std::getline(in,line);
        while(std::getline(in,line)){auto r=split(line);start.push_back(std::stod(r.at(0)));end.push_back(std::stod(r.at(1)));value.push_back(std::stod(r.at(2)));}
        if(start.empty())throw std::runtime_error("Missing pulsed voltage intervals");
        for(std::size_t i=0;i<start.size();++i)if(!std::isfinite(start[i]) || !std::isfinite(end[i]) || !std::isfinite(value[i]) || end[i]<=start[i] || value[i]>=0 || value[i] < -1e6 || (i && std::abs(start[i]-end[i-1])>1e-21))throw std::runtime_error("Invalid pulsed voltage intervals");
    }
    double at(double t) {
        if(t<start.front()-1e-20 || t>end.back()+1e-20)throw std::runtime_error("Voltage extrapolation prohibited");
        const double lookup=t+32*std::numeric_limits<double>::epsilon()*std::max(std::abs(t),1e-12);
        current=std::clamp(int(std::upper_bound(start.begin(),start.end(),lookup)-start.begin())-1,0,int(start.size())-1);
        return value[current];
    }
};
struct Thermal {
    std::vector<double> times;
    std::array<std::vector<double>,4> temperatures;
    std::array<std::uint64_t,4> sizes{};
    explicit Thermal(const std::filesystem::path& path) {
        std::ifstream in(path,std::ios::binary);std::uint64_t magic=0,n=0;
        in.read(reinterpret_cast<char*>(&magic),8);in.read(reinterpret_cast<char*>(&n),8);
        in.read(reinterpret_cast<char*>(sizes.data()),32);
        if(!in || magic!=0x544845524d303031ULL || n<2 || n>10000)throw std::runtime_error("Invalid thermal history header");
        times.resize(n);in.read(reinterpret_cast<char*>(times.data()),8*n);
        for(int f=0;f<4;++f){if(sizes[f]<1 || sizes[f]>4096)throw std::runtime_error("Invalid thermal patch count");temperatures[f].resize(n*sizes[f]);in.read(reinterpret_cast<char*>(temperatures[f].data()),8*n*sizes[f]);}
        if(!in || in.peek()!=std::char_traits<char>::eof())throw std::runtime_error("Truncated thermal history");
        for(std::size_t i=0;i<times.size();++i)if(!std::isfinite(times[i]) || (i && times[i]<=times[i-1]))throw std::runtime_error("Invalid thermal clock");
        for(const auto& v:temperatures)for(double T:v)if(!std::isfinite(T) || T<=0 || T>1e6)throw std::runtime_error("Invalid absolute temperature");
    }
    std::array<std::vector<double>,4> at(double t) const {
        if(t<times.front()-1e-20 || t>times.back()+1e-20)throw std::runtime_error("Temperature extrapolation prohibited");
        auto k=std::clamp(int(std::upper_bound(times.begin(),times.end(),t)-times.begin())-1,0,int(times.size())-2);
        double q=std::clamp((t-times[k])/(times[k+1]-times[k]),0.,1.);std::array<std::vector<double>,4> result;
        for(int f=0;f<4;++f){result[f].resize(sizes[f]);for(std::size_t i=0;i<sizes[f];++i)result[f][i]=(1-q)*temperatures[f][k*sizes[f]+i]+q*temperatures[f][(k+1)*sizes[f]+i];}
        return result;
    }
};
struct Fields {
    std::filesystem::path directory;
    std::vector<double> times;
    std::vector<double> bounds;
    std::unique_ptr<NativeField> lower,upper;
    int current=-1;double fraction=0;
    explicit Fields(const std::filesystem::path& path):directory(path) {
        std::ifstream in(path/"em_field_snapshots.csv");std::string line;std::getline(in,line);
        while(std::getline(in,line))times.push_back(std::stod(split(line).at(0)));
        if(times.size()<2)throw std::runtime_error("Missing field history");bounds.resize(times.size(),-1);
    }
    static double maximum(NativeField& f) {
        double value=0;mfem::Vector g(2);
        // Q1 gradient norm attains its maximum at an element corner.
        for(int e=0;e<f.mesh.GetNE();++e)for(double x:{0.,1.})for(double y:{0.,1.}) {
            mfem::IntegrationPoint ip;ip.Set2(x,y);auto* tr=f.mesh.GetElementTransformation(e);tr->SetIntPoint(&ip);f.az->GetGradient(*tr,g);value=std::max(value,g.Norml2());
        }
        return value;
    }
    void set(double t) {
        if(t<times.front()-1e-20 || t>times.back()+1e-20)throw std::runtime_error("Field extrapolation prohibited");
        const double lookup=t+32*std::numeric_limits<double>::epsilon()*std::max(std::abs(t),1e-12);
        int k=std::clamp(int(std::upper_bound(times.begin(),times.end(),lookup)-times.begin())-1,0,int(times.size())-2);
        if(k!=current) {
            if(k==current+1 && upper)lower=std::move(upper);else lower=std::make_unique<NativeField>(directory,times[k]);
            upper=std::make_unique<NativeField>(directory,times[k+1]);current=k;
            if(bounds[k]<0)bounds[k]=maximum(*lower);if(bounds[k+1]<0)bounds[k+1]=maximum(*upper);
        }
        fraction=std::clamp((t-times[k])/(times[k+1]-times[k]),0.,1.);
    }
    double bound() const {return std::max(bounds[current],bounds[current+1]);}
    PicVector value(double x,double y,double t,bool magnetic) {
        set(t);auto a=magnetic?lower->magnetic(x,y):lower->electric(x,y);auto b=magnetic?upper->magnetic(x,y):upper->electric(x,y);
        for(int i=0;i<3;++i)a[i]=(1-fraction)*a[i]+fraction*b[i];return a;
    }
};
struct Reservoir {std::vector<double> integral,weighted_temperature;double span=0;};
template<class Pic,class Field,class Bias>
auto allowed_step(Pic& pic,Field& fields,Bias& voltage,const Thermal& thermal,
 const std::array<Reservoir,4>& reservoir,const GapPicOptions& o,double t,double end,int next_frame) {
            double dt=pic.particles().empty()?1e-12:std::min(1e-14,.05*pic_me/(pic_e*std::max(fields.bound(),1e-20)));
            if(!pic.particles().empty()) {
                std::vector<double> count(o.nx*o.ny,0.);double vmax=0;
                for(const auto& p:pic.particles()){int i=std::clamp(int(p.position_m[0]/25e-6),0,o.nx-1),j=std::clamp(int(p.position_m[1]/25e-6),0,o.ny-1);count[j*o.nx+i]+=p.weight;double u=std::hypot(p.proper_velocity_m_s[0],p.proper_velocity_m_s[1],p.proper_velocity_m_s[2]);vmax=std::max(vmax,u/std::sqrt(1+u*u/(pic_c*pic_c)));}
                double density=*std::max_element(count.begin(),count.end())/(25e-6*25e-6*o.depth_m);
                dt=std::min(dt,.1/std::sqrt(std::max(density,1e-100)*pic_e*pic_e/(eps*pic_me)));
                if(vmax>0)dt=std::min(dt,.2*25e-6/vmax);
            }
            dt=std::min({dt,end-t,fields.times[fields.current+1]-t,voltage.end[voltage.current]-t,std::max(1e-22,next_frame*1e-9-t)});
            if(dt<=1e-23 || !std::isfinite(dt))throw std::runtime_error("Nonprogressing time step");
            auto T=thermal.at(t+.5*dt);double jmax=0;
            for(const auto& v:T)for(double vT:v)jmax=std::max(jmax,pic_richardson_supply_A_m2(vT,4.3,1201732.289489034));
            // Source-induced density also bounds the new explicit plasma step.
            double source_bound=std::cbrt(.01*eps*pic_me*25e-6/(std::max(jmax,1e-100)*pic_e));
            if(dt>source_bound){dt=source_bound;T=thermal.at(t+.5*dt);}
            if(pic.particles().empty()) {
                bool will_emit=false;
                for(int f=0;f<4;++f) {
                    double sum=std::accumulate(reservoir[f].integral.begin(),reservoir[f].integral.end(),0.);
                    for(double vT:T[f])sum+=pic_richardson_supply_A_m2(vT,4.3,1201732.289489034)*dt;
                    if(sum*(f<2?.005:.0002)/T[f].size()*o.depth_m/pic_e>=1.)will_emit=true;
                }
                if(will_emit) {
                    dt=std::min(dt,std::min(1e-14,.05*pic_me/(pic_e*std::max(fields.bound(),1e-20))));
                    T=thermal.at(t+.5*dt);
                }
            }

return std::pair<double,std::array<std::vector<double>,4>>{dt,std::move(T)};
}
template<class Source>
bool accumulate_source(Source& pic,std::array<Reservoir,4>& reservoir,
 std::array<std::vector<double>,4>& T,double dt,const GapPicOptions& o,std::uint64_t step) {
            bool injected=false;
            for(int f=0;f<4;++f) {
                auto& pool=reservoir[f];pool.span+=dt;
                for(std::size_t i=0;i<T[f].size();++i){double amount=pic_richardson_supply_A_m2(T[f][i],4.3,1201732.289489034)*dt;pool.integral[i]+=amount;pool.weighted_temperature[i]+=amount*T[f][i];}
                double ne=std::accumulate(pool.integral.begin(),pool.integral.end(),0.)*(f<2?.005:.0002)/T[f].size()*o.depth_m/pic_e;
                if(ne<1.)continue;
                auto temp=T[f],J=pool.integral;
                for(std::size_t i=0;i<J.size();++i){if(J[i]>0)temp[i]=pool.weighted_temperature[i]/J[i];J[i]/=pool.span;}
                const GapSurface side=f==0?GapSurface::Back:f==1?GapSurface::Gap:f==2?GapSurface::Left:GapSurface::Right;
                pic.inject_surface_profile(GapElectrode::Upper,side,temp,J,pool.span,32,1000000ULL*step+f*100);
                if(f>=2){std::reverse(temp.begin(),temp.end());std::reverse(J.begin(),J.end());}
                pic.inject_surface_profile(GapElectrode::Lower,side,temp,J,pool.span,32,1000000ULL*step+f*100+50);
                std::fill(pool.integral.begin(),pool.integral.end(),0.);std::fill(pool.weighted_temperature.begin(),pool.weighted_temperature.end(),0.);pool.span=0;injected=true;
            }
return injected;
}
} // namespace validation_native
