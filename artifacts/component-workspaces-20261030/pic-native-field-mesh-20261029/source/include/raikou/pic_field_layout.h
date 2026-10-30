#pragma once
#include <mfem.hpp>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>
#include <utility>

namespace raikou::mfem_backend {
struct PicKineticPoint { double x,y,z; };
struct PicElectricPoint { double x,y; };
struct PicKineticCellPoint { int i,j;double xi,eta; };
struct PicNativeGridPoint { double x,y; };
inline PicElectricPoint electric_point(PicKineticPoint p) { return {p.x,p.y}; }
// Only raw native-grid sampling uses this conversion. Existing Fields::value
// already performs its own offset internally and must receive PIC coordinates.
inline PicNativeGridPoint native_grid_point(PicKineticPoint p) { return {p.x-.0045,p.y-.0045}; }

struct PicFieldMeshOptions {
    int margin_cells=0;
    double growth_ratio=1.05;
    // Both zero selects fixed-growth nested axes. Positive values request a
    // fitted physical full extent; every original core node remains unchanged.
    double full_width_m=0,full_height_m=0;
};

class PicFieldLayout {
    static void require(bool condition,const char* message) {
        if(!condition)throw std::invalid_argument(message);
    }
    static std::vector<double> extend(const std::vector<double>& core,double width,
        int count,double ratio,double requested_full) {
        if(!count)return core;
        const double first=width/(core.size()-1);
        const auto length=[&](double r){double sum=0,h=first;for(int i=0;i<count;++i){sum+=h;h*=r;}return sum;};
        double extent=0;
        if(requested_full>0){
            require(requested_full>width,"Field extent must contain the kinetic core");
            extent=(requested_full-width)/2;
            require(extent>=first*count&&length(1.15)>=extent,"Fitted exterior requires more margin cells");
            double low=1,high=1.15;
            for(int k=0;k<100;++k){const double mid=.5*(low+high);if(length(mid)<extent)low=mid;else high=mid;}
            ratio=.5*(low+high);
        }
        std::vector<double> distance(count+1,0.);double h=first;
        for(int k=0;k<count;++k){distance[k+1]=distance[k]+h;h*=ratio;}
        if(requested_full>0)distance.back()=extent;
        std::vector<double> axis;axis.reserve(core.size()+2*count);
        for(int k=count;k>0;--k)axis.push_back(-distance[k]);
        axis.insert(axis.end(),core.begin(),core.end());
        for(int k=1;k<=count;++k)axis.push_back(width+distance[k]);
        for(std::size_t k=1;k<axis.size();++k)require(std::isfinite(axis[k])&&axis[k]>axis[k-1],"Invalid field axis");
        return axis;
    }
    static std::pair<int,double> axis_locate(const std::vector<double>& axis,
        int core_cells,int margin,double width,double point) {
        require(std::isfinite(point)&&point>=axis.front()&&point<=axis.back(),"Point outside electric field mesh");
        if(point>=0&&point<=width){
            const double s=point/width*core_cells;const int i=std::min(int(s),core_cells-1);
            return {i+margin,s-i}; // exact original core arithmetic
        }
        const int i=std::clamp(int(std::upper_bound(axis.begin(),axis.end(),point)-axis.begin())-1,0,int(axis.size()-2));
        return {i,(point-axis[i])/(axis[i+1]-axis[i])};
    }
public:
    const int core_nx,core_ny,margin;
    const double core_width,core_height;
    const PicFieldMeshOptions configuration;
    std::vector<double> xs,ys;
    PicFieldLayout(int nx,int ny,double width,double height,PicFieldMeshOptions options={})
        :core_nx(nx),core_ny(ny),margin(options.margin_cells),core_width(width),core_height(height),configuration(options) {
        require(nx>0&&ny>0&&std::isfinite(width)&&width>0&&std::isfinite(height)&&height>0,"Invalid electric core");
        require(margin>=0&&margin<=300&&std::isfinite(options.growth_ratio)&&options.growth_ratio>=1&&options.growth_ratio<=1.15,"Invalid bounded electric exterior");
        require(std::isfinite(options.full_width_m)&&std::isfinite(options.full_height_m)&&options.full_width_m>=0&&options.full_height_m>=0,"Invalid field extent");
        require((options.full_width_m==0)==(options.full_height_m==0),"Specify both fitted field extents");
        require(margin>0||(options.full_width_m==0&&options.full_height_m==0),"Uniform field requires zero extent override");
        require(std::uint64_t(nx+2*margin+1)*(ny+2*margin+1)<=1010000,"Electric field node ceiling exceeded");
        auto core=mfem::Mesh::MakeCartesian2D(nx,ny,mfem::Element::QUADRILATERAL,true,width,height,false);
        std::vector<double>x(nx+1),y(ny+1);mfem::IntegrationPoint p;p.Set2(0,0);mfem::Vector v(2);
        for(int i=0;i<nx;++i){core.GetElementTransformation(i)->Transform(p,v);x[i]=v[0];}
        for(int j=0;j<ny;++j){core.GetElementTransformation(j*nx)->Transform(p,v);y[j]=v[1];}
        p.Set2(1,1);core.GetElementTransformation(nx*ny-1)->Transform(p,v);x.back()=v[0];y.back()=v[1];
        xs=extend(x,width,margin,options.growth_ratio,options.full_width_m);ys=extend(y,height,margin,options.growth_ratio,options.full_height_m);
    }
    int nx() const {return int(xs.size())-1;}
    int ny() const {return int(ys.size())-1;}
    mfem::Mesh make_mesh() const {
        if(!margin)return mfem::Mesh::MakeCartesian2D(core_nx,core_ny,mfem::Element::QUADRILATERAL,true,core_width,core_height,false);
        auto mesh=mfem::Mesh::MakeCartesian2D(nx(),ny(),mfem::Element::QUADRILATERAL,true,nx(),ny(),false);
        for(int n=0;n<mesh.GetNV();++n){auto* p=mesh.GetVertex(n);const int i=std::lround(p[0]),j=std::lround(p[1]);p[0]=xs.at(i);p[1]=ys.at(j);}
        return mesh;
    }
    int locate(PicElectricPoint p,mfem::IntegrationPoint& ip) const {
        const auto x=axis_locate(xs,core_nx,margin,core_width,p.x),y=axis_locate(ys,core_ny,margin,core_height,p.y);
        ip.Set2(x.second,y.second);return y.first*nx()+x.first;
    }
};
} // namespace raikou::mfem_backend
