#pragma once
#include "raikou/mfem_gap_pic.h"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <memory>
#include <cmath>
#include <limits>
using namespace raikou::mfem_backend;
std::vector<std::string> split(std::string line) {
    std::vector<std::string> out;std::istringstream in(line);std::string v;
    while(std::getline(in,v,',')) out.push_back(v);
    return out;
}
struct NativeField {
    mfem::Mesh mesh;
    std::vector<double> xs,ys;
    std::vector<int> elements;
    std::unique_ptr<mfem::GridFunction> az,rate;
    std::unique_ptr<Q1ScalarGather> az_gather,rate_gather;
    double last_x=std::numeric_limits<double>::quiet_NaN(),last_y=last_x;
    int last_element=-1;
    mfem::IntegrationPoint last_point;
    explicit NativeField(const std::filesystem::path& dir,double time):mesh((dir/"em_mesh.mesh").c_str(),1,1) {
        std::ifstream axes(dir/"em_axes.csv");std::string line;std::getline(axes,line);
        while(std::getline(axes,line)) {auto r=split(line);(r.at(0)=="x"?xs:ys).push_back(std::stod(r.at(2)));}
        if(xs.size()<2 || ys.size()<2) throw std::runtime_error("Missing native field axes");
        std::ifstream index(dir/"em_field_snapshots.csv");std::getline(index,line);
        double best=1e100;std::string a,b;
        while(std::getline(index,line)) {auto r=split(line);double error=std::abs(std::stod(r.at(0))-time);
            if(error<best){best=error;a=r.at(1);b=r.at(2);}}
        if(best>1e-17) throw std::runtime_error("Require an exact native snapshot at test start");
        std::ifstream ia(dir/a),ir(dir/b);az=std::make_unique<mfem::GridFunction>(&mesh,ia);rate=std::make_unique<mfem::GridFunction>(&mesh,ir);
        az_gather=std::make_unique<Q1ScalarGather>(*az);
        rate_gather=std::make_unique<Q1ScalarGather>(*rate);
        const int nx=static_cast<int>(xs.size())-1;elements.resize(nx*(ys.size()-1));
        for(int e=0;e<mesh.GetNE();++e) {
            mfem::Vector p(2);mesh.GetElementTransformation(e)->Transform(mfem::Geometries.GetCenter(mfem::Geometry::SQUARE),p);
            elements.at(bin(ys,p[1])*nx+bin(xs,p[0]))=e;
        }
    }
    static int bin(const std::vector<double>& a,double x) {
        if(x<a.front() || x>a.back()) throw std::runtime_error("Particle left native field coverage");
        return std::clamp(int(std::upper_bound(a.begin(),a.end(),x)-a.begin())-1,0,int(a.size())-2);
    }
    void locate(double px,double py) {
        if(last_element>=0 && px==last_x && py==last_y) return;
        // PIC box coordinates to original centered heating coordinates.
        const double x=px-.0045,y=py-.0045;
        const int i=bin(xs,x),j=bin(ys,y),e=elements.at(j*(xs.size()-1)+i);
        last_point.Set2((x-xs[i])/(xs[i+1]-xs[i]),(y-ys[j])/(ys[j+1]-ys[j]));
        last_x=px;last_y=py;last_element=e;
    }
    PicVector magnetic(double px,double py) {
        locate(px,py);
        const auto g=az_gather->gradient(last_element,last_point);
        return {g[1],-g[0],0};
    }
    PicVector electric(double px,double py) {
        locate(px,py);
        return {0,0,-rate_gather->value(last_element,last_point)};
    }
};
