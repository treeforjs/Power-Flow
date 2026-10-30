#pragma once
#include "batch_snapshot.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace electric_exterior {
// These frames cannot be implicitly substituted for one another.
struct ParticlePoint { double x,y,z; };
struct ElectricPoint { double x,y; };
struct NativePoint { double x,y; };
inline ElectricPoint electric_point(ParticlePoint p){return {p.x,p.y};}
inline NativePoint native_point(ParticlePoint p){return {p.x-.0045,p.y-.0045};}
struct LocalCoordinate { int cell;double xi; };

class Axis {
    std::vector<double> core_;
    double width_;
    int margin_;
public:
    std::vector<double> nodes;
    double growth=1,first_width=0,last_width=0;
    Axis(std::vector<double> core,double width,double extent,int count,double fixed_growth=0):core_(std::move(core)),width_(width),margin_(count) {
        snapshot_io::require(core_.size()>=2 && core_.front()==0 && core_.back()==width && count>=0,"Invalid core axis");
        for(std::size_t i=1;i<core_.size();++i)snapshot_io::require(core_[i]>core_[i-1],"Nonincreasing core axis");
        if(count==0){snapshot_io::require(extent==0,"Zero margin requires captured uniform extent");nodes=core_;return;}
        snapshot_io::require(count<=300 && std::isfinite(extent) && extent>0,"Invalid bounded exterior margin");
        first_width=width/(core_.size()-1);
        snapshot_io::require(extent>=count*first_width,"Exterior cannot coarsen below near-core cell size");
        const auto length=[&](double ratio){double sum=0,h=first_width;for(int i=0;i<count;++i){sum+=h;h*=ratio;}return sum;};
        if(fixed_growth){snapshot_io::require(fixed_growth>=1&&fixed_growth<=1.15,"Invalid fixed growth");growth=fixed_growth;}
        else {double low=1,high=1.15;
            snapshot_io::require(length(high)>=extent,"Requested margin needs more cells to keep growth<=1.15");
            for(int i=0;i<100;++i){double mid=.5*(low+high);if(length(mid)<extent)low=mid;else high=mid;}
            growth=.5*(low+high);
        }
        std::vector<double> distance(count+1,0.);double h=first_width;
        for(int i=0;i<count;++i){distance[i+1]=distance[i]+h;h*=growth;}
        // Only far endpoint is corrected for bisection roundoff. The original
        // near-core nodes and the first outer cell are never rescaled.
        if(!fixed_growth)distance.back()=extent;
        last_width=distance.back()-distance[count-1];
        nodes.reserve(core_.size()+2*count);
        for(int i=count;i>0;--i)nodes.push_back(-distance[i]);
        nodes.insert(nodes.end(),core_.begin(),core_.end());
        for(int i=1;i<=count;++i)nodes.push_back(width+distance[i]);
        for(std::size_t i=1;i<nodes.size();++i)snapshot_io::require(std::isfinite(nodes[i])&&nodes[i]>nodes[i-1],"Invalid graded axis");
    }
    LocalCoordinate locate(double x)const {
        snapshot_io::require(std::isfinite(x)&&x>=nodes.front()&&x<=nodes.back(),"Electric point outside exported axis");
        // Preserve the exact sealed normalized arithmetic on the entire original
        // particle region. Its coordinates and source weights do not change
        // when remote vacuum cells are appended.
        if(x>=0 && x<=width_) {
            const double s=x/width_*(core_.size()-1);
            const int i=std::min(int(s),int(core_.size()-2));return {i+margin_,s-i};
        }
        const int i=std::clamp(int(std::upper_bound(nodes.begin(),nodes.end(),x)-nodes.begin())-1,0,int(nodes.size()-2));
        return {i,(x-nodes[i])/(nodes[i+1]-nodes[i])};
    }
    int cells()const{return int(nodes.size())-1;}
    int margin_cells()const{return margin_;}
};

struct LocatedPoint {int element;mfem::IntegrationPoint reference;};
class Layout {
    static double nested_width(const pic_batch::HostSnapshot&s,int count){
        snapshot_io::require(count>=1&&count<=120,"Nested exterior margin must be 1..120 cells");
        double extent=0,h=s.parameters.width/s.parameters.self_nx;
        for(int i=0;i<count;++i){extent+=h;h*=1.05;}return s.parameters.width+2*extent;
    }
    static std::vector<double> core_x(const pic_batch::HostSnapshot&s) {
        std::vector<double>x(s.parameters.self_nx+1);for(int i=0;i<s.parameters.self_nx;++i)x[i]=s.self_cells[i].x[0];x.back()=s.self_cells[s.parameters.self_nx-1].x[1];return x;
    }
    static std::vector<double> core_y(const pic_batch::HostSnapshot&s) {
        std::vector<double>y(s.parameters.self_ny+1);for(int i=0;i<s.parameters.self_ny;++i)y[i]=s.self_cells[i*s.parameters.self_nx].y[0];y.back()=s.self_cells[(s.parameters.self_ny-1)*s.parameters.self_nx].y[3];return y;
    }
public:
    const double core_width,core_height,full_width;
    const int margin_count;
    Axis x,y;
    Layout(const pic_batch::HostSnapshot&s,int width_mm,int margin):
      core_width(s.parameters.width),core_height(s.parameters.height),
      full_width(width_mm==0?nested_width(s,margin):(width_mm==9?s.parameters.width:width_mm*.001)),margin_count(margin),
      x(core_x(s),core_width,(full_width-core_width)/2,margin,width_mm==0?1.05:0),y(core_y(s),core_height,(full_width-core_height)/2,margin,width_mm==0?1.05:0) {
        snapshot_io::require(width_mm==0||width_mm==9||width_mm==13||width_mm==17||width_mm==25||width_mm==40||width_mm==80||width_mm==160,"Undeclared static domain");
        snapshot_io::require((width_mm==9)==(margin==0),"Captured limit must have zero exterior cells");
        snapshot_io::require(std::uint64_t(x.cells()+1)*(y.cells()+1)<=1010000,"Node ceiling exceeded");
    }
    mfem::Mesh mesh()const {
        if(margin_count==0)return mfem::Mesh::MakeCartesian2D(x.cells(),y.cells(),mfem::Element::QUADRILATERAL,true,core_width,core_height,false);
        auto m=mfem::Mesh::MakeCartesian2D(x.cells(),y.cells(),mfem::Element::QUADRILATERAL,true,x.cells(),y.cells(),false);
        for(int v=0;v<m.GetNV();++v){auto*p=m.GetVertex(v);int i=std::lround(p[0]),j=std::lround(p[1]);snapshot_io::require(i>=0&&i<=x.cells()&&j>=0&&j<=y.cells(),"Index mesh vertex outside axis");p[0]=x.nodes[i];p[1]=y.nodes[j];}
        return m;
    }
    LocatedPoint locate(ElectricPoint p)const {
        auto a=x.locate(p.x),b=y.locate(p.y);mfem::IntegrationPoint ip;ip.Set2(a.xi,b.xi);return {b.cell*x.cells()+a.cell,ip};
    }
    LocatedPoint locate(ParticlePoint p)const{return locate(electric_point(p));}
    void validate_mesh(mfem::Mesh&mesh)const {
        snapshot_io::require(mesh.GetNE()==x.cells()*y.cells(),"Mesh cell count mismatch");
        for(int j=0;j<y.cells();++j)for(int i=0;i<x.cells();++i){
            mfem::Vector center(2);mesh.GetElementTransformation(j*x.cells()+i)->Transform(mfem::Geometries.GetCenter(mfem::Geometry::SQUARE),center);
            snapshot_io::require(center[0]>x.nodes[i]&&center[0]<x.nodes[i+1]&&center[1]>y.nodes[j]&&center[1]<y.nodes[j+1],"Mesh/axis element order mismatch");
        }
    }
};

struct ScalarSample {double value,dx,dy;};
inline void deposit_point(const Layout&layout,mfem::FiniteElementSpace&space,
                          mfem::Vector&rhs,ParticlePoint point,double charge_per_m) {
    snapshot_io::require(std::isfinite(charge_per_m),"Nonfinite point charge");
    const auto located=layout.locate(point);const auto*fe=space.GetFE(located.element);
    mfem::Array<int>dofs;mfem::Vector shape(fe->GetDof());fe->CalcShape(located.reference,shape);
    space.GetElementDofs(located.element,dofs);rhs.AddElementVector(dofs,shape*=charge_per_m);
}
inline ScalarSample gather(const Layout&layout,const mfem::GridFunction&field,ElectricPoint point) {
    const auto located=layout.locate(point);auto*tr=field.FESpace()->GetMesh()->GetElementTransformation(located.element);tr->SetIntPoint(&located.reference);
    mfem::Vector gradient(2);field.GetGradient(*tr,gradient);const double value=field.GetValue(located.element,located.reference);
    snapshot_io::require(std::isfinite(value)&&std::isfinite(gradient[0])&&std::isfinite(gradient[1]),"Nonfinite graded scalar sample");
    return {value,gradient[0],gradient[1]};
}
} // namespace electric_exterior
