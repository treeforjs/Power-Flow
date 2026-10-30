#pragma once
// Benchmark-only independent reference, extracted unchanged from the qualified
// CPU-v2 driver. Not a new application provider or force implementation.
#include "native_heating_fields.hpp"
namespace native_reference {
struct Voltage {
    std::vector<double> start,end,value;
    int current=0;
    double last_time=std::numeric_limits<double>::quiet_NaN();
    explicit Voltage(const std::filesystem::path& path) {
        std::ifstream in(path);std::string line;std::getline(in,line);
        while(std::getline(in,line)){auto r=split(line);start.push_back(std::stod(r.at(0)));end.push_back(std::stod(r.at(1)));value.push_back(std::stod(r.at(2)));}
        if(start.empty())throw std::runtime_error("Missing pulsed voltage intervals");
        for(std::size_t i=0;i<start.size();++i)if(!std::isfinite(start[i]) || !std::isfinite(end[i]) || !std::isfinite(value[i]) || end[i]<=start[i] || value[i]>=0 || value[i] < -1e6 || (i && std::abs(start[i]-end[i-1])>1e-21))throw std::runtime_error("Invalid pulsed voltage intervals");
    }
    double at(double t) {
        if(t==last_time) return value[current];
        if(t<start.front()-1e-20 || t>end.back()+1e-20)throw std::runtime_error("Voltage extrapolation prohibited");
        const double lookup=t+32*std::numeric_limits<double>::epsilon()*std::max(std::abs(t),1e-12);
        current=std::clamp(int(std::upper_bound(start.begin(),start.end(),lookup)-start.begin())-1,0,int(start.size())-1);
        last_time=t;
        return value[current];
    }
};
struct Fields {
    std::filesystem::path directory;
    std::vector<double> times;
    std::vector<double> bounds;
    std::unique_ptr<NativeField> lower,upper;
    int current=-1;double fraction=0;
    double last_time=std::numeric_limits<double>::quiet_NaN();
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
        if(t==last_time) return;
        if(t<times.front()-1e-20 || t>times.back()+1e-20)throw std::runtime_error("Field extrapolation prohibited");
        const double lookup=t+32*std::numeric_limits<double>::epsilon()*std::max(std::abs(t),1e-12);
        int k=std::clamp(int(std::upper_bound(times.begin(),times.end(),lookup)-times.begin())-1,0,int(times.size())-2);
        if(k!=current) {
            if(k==current+1 && upper)lower=std::move(upper);else lower=std::make_unique<NativeField>(directory,times[k]);
            upper=std::make_unique<NativeField>(directory,times[k+1]);current=k;
            if(bounds[k]<0)bounds[k]=maximum(*lower);if(bounds[k+1]<0)bounds[k+1]=maximum(*upper);
        }
        fraction=std::clamp((t-times[k])/(times[k+1]-times[k]),0.,1.);
        last_time=t;
    }
    double bound() const {return std::max(bounds[current],bounds[current+1]);}
    PicVector value(double x,double y,double t,bool magnetic) {
        set(t);auto a=magnetic?lower->magnetic(x,y):lower->electric(x,y);auto b=magnetic?upper->magnetic(x,y):upper->electric(x,y);
        for(int i=0;i<3;++i)a[i]=(1-fraction)*a[i]+fraction*b[i];return a;
    }
};
}
