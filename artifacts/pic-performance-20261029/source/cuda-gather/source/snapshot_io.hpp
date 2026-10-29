#pragma once
#include "gather_push_core.hpp"
#include "raikou/mfem_gap_pic.h"
#include <filesystem>
#include <fstream>
#include <set>
#include <stdexcept>
#include <vector>
#include <cmath>
namespace snapshot_io {
namespace gp=gather_push;
using raikou::mfem_backend::PicVector;
inline void require(bool ok,const char* message) {if(!ok)throw std::runtime_error(message);}
template<class T> void write(const std::filesystem::path& path,const std::vector<T>& data) {
    std::ofstream f(path,std::ios::binary);f.write(reinterpret_cast<const char*>(data.data()),data.size()*sizeof(T));
    require(bool(f),"Cannot write fixture");
}
template<class T> std::vector<T> read(const std::filesystem::path& path) {
    std::ifstream f(path,std::ios::binary|std::ios::ate);require(bool(f),"Cannot read fixture");
    const auto size=f.tellg();require(size>=0 && size%sizeof(T)==0,"Fixture has invalid size");
    std::vector<T> data(static_cast<std::size_t>(size)/sizeof(T));f.seekg(0);
    f.read(reinterpret_cast<char*>(data.data()),size);require(bool(f),"Truncated fixture");return data;
}
inline gp::Vec vec(PicVector v) {return {{v[0],v[1],v[2]}};}
inline PicVector vec(gp::Vec v) {return {v.v[0],v.v[1],v.v[2]};}
inline std::vector<double> values(const mfem::GridFunction& field) {
    return {field.HostRead(),field.HostRead()+field.Size()};
}
inline std::vector<gp::Cell> export_grid(const mfem::GridFunction& field,const std::vector<int>& elements) {
    const auto* space=field.FESpace();auto* mesh=space->GetMesh();
    require(space->GetVDim()==1 && mesh->Dimension()==2 && mesh->SpaceDimension()==2 && !mesh->GetNodes(),
            "Only straight scalar 2D Q1 meshes are supported");
    require(elements.size()==static_cast<std::size_t>(mesh->GetNE()),"Grid cell count mismatch");
    std::set<int> seen;std::vector<gp::Cell> cells(elements.size());
    const double cx[4]={0,1,1,0},cy[4]={0,0,1,1};
    for(std::size_t i=0;i<elements.size();++i) {
        const int e=elements[i];require(e>=0 && e<mesh->GetNE() && seen.insert(e).second,"Duplicate or absent native cell");
        const auto* fe=space->GetFE(e);mfem::Array<int> dofs;mfem::DofTransformation transform;
        const auto* q1=dynamic_cast<const mfem::H1_QuadrilateralElement*>(fe);
        require(q1 && q1->GetBasisType()==mfem::BasisType::GaussLobatto && fe->GetOrder()==1 && fe->GetDof()==4 &&
                fe->GetMapType()==mfem::FiniteElement::VALUE,"Unsupported finite element");
        space->GetElementDofs(e,dofs,transform);
        require(dofs.Size()==4 && transform.IsIdentity(),"Unsupported DOF transform");
        for(int k=0;k<4;++k) {
            const auto& ip=fe->GetNodes().IntPoint(k);
            require(ip.x==cx[k] && ip.y==cy[k],"Q1 local ordering changed");
            require(dofs[k]>=0 && dofs[k]<field.Size(),"Invalid scalar DOF");
            cells[i].dof[k]=dofs[k];mfem::Vector point(2);
            mesh->GetElementTransformation(e)->Transform(ip,point);
            cells[i].x[k]=point[0];cells[i].y[k]=point[1];
        }
    }
    for(double p:values(field))require(std::isfinite(p),"Nonfinite snapshot");
    return cells;
}

} // namespace snapshot_io
