#include "raikou/q1_scalar_gather.h"
#include <cmath>
#include <cstring>
#include <iostream>
#include <stdexcept>

using raikou::mfem_backend::Q1ScalarGather;
namespace {
void check(bool pass,const char* message) {
    if(!pass) throw std::runtime_error(message);
}
void identical(double a,double b) {
    check(std::memcmp(&a,&b,sizeof(a))==0,"Q1 gather changed the MFEM result bits");
}
void dof_renumbering_invalidation() {
    auto mesh=mfem::Mesh::MakeCartesian2D(3,2,mfem::Element::QUADRILATERAL,true,2.,3.,false);
    mfem::H1_FECollection collection(1,2);
    mfem::FiniteElementSpace space(&mesh,&collection);
    mfem::GridFunction field(&space);
    for(int i=0;i<field.Size();++i) field[i]=i+1.;
    Q1ScalarGather gather(field);
    check(gather.available(),"Q1 gather initially cached");
    const long mesh_sequence=mesh.GetSequence(),space_sequence=space.GetSequence();
    // MFEM does not advance either sequence when it renumbers element DOFs.
    gather.invalidate();
    space.ReorderElementToDofTable();
    check(mesh.GetSequence()==mesh_sequence && space.GetSequence()==space_sequence,
          "DOF renumbering regression must exercise unchanged sequence counters");
    check(!gather.available(),"explicit invalidation rejects renumbered connectivity");
    gather.invalidate(); // Idempotent, including an already-disabled cache.
    Q1ScalarGather rebuilt(field);
    check(rebuilt.available(),"a new gather caches the renumbered connectivity");
    for(int e=0;e<mesh.GetNE();++e) {
        for(double x:{0.,.23,1.}) for(double y:{0.,.61,1.}) {
            mfem::IntegrationPoint point;point.Set2(x,y);
            const double value=field.GetValue(e,point);
            identical(gather.value(e,point),value);
            identical(rebuilt.value(e,point),value);
            auto* transformation=mesh.GetElementTransformation(e);
            transformation->SetIntPoint(&point);
            mfem::Vector reference(2);field.GetGradient(*transformation,reference);
            const auto fallback=gather.gradient(e,point),cached=rebuilt.gradient(e,point);
            for(int i=0;i<2;++i) {
                identical(fallback[i],reference[i]);identical(cached[i],reference[i]);
            }
        }
    }
}
void compare(mfem::Mesh& mesh,int order,bool supported) {
    mfem::H1_FECollection collection(order,2);
    mfem::FiniteElementSpace space(&mesh,&collection);
    mfem::GridFunction field(&space);
    for(int i=0;i<field.Size();++i) field[i]=std::sin(17.*i)+i*.13;
    Q1ScalarGather gather(field);
    check(gather.available()==supported,"Q1 gather structural guard");
    for(int e=0;e<mesh.GetNE();++e) {
        for(double x:{0.,.001,.13,.5,.999,1.}) for(double y:{0.,.001,.37,.5,.999,1.}) {
            mfem::IntegrationPoint point;point.Set2(x,y);
            identical(gather.value(e,point),field.GetValue(e,point));
            auto* transformation=mesh.GetElementTransformation(e);
            transformation->SetIntPoint(&point);
            mfem::Vector reference(2);field.GetGradient(*transformation,reference);
            const auto actual=gather.gradient(e,point);
            identical(actual[0],reference[0]);identical(actual[1],reference[1]);
        }
    }
    // Connectivity becomes stale after a mesh/space update, so use MFEM again.
    mesh.UniformRefinement();space.Update();field.Update();field=3.;
    check(!gather.available(),"Q1 gather rejects stale connectivity");
    mfem::IntegrationPoint point;point.Set2(.23,.61);
    identical(gather.value(0,point),field.GetValue(0,point));
    auto* transformation=mesh.GetElementTransformation(0);
    transformation->SetIntPoint(&point);
    mfem::Vector reference(2);field.GetGradient(*transformation,reference);
    const auto actual=gather.gradient(0,point);
    identical(actual[0],reference[0]);identical(actual[1],reference[1]);
}
}
int main() {
    try {
        dof_renumbering_invalidation();
        auto rectangles=mfem::Mesh::MakeCartesian2D(7,9,mfem::Element::QUADRILATERAL,true,.009,.013,false);
        // Nonuniform rectangle spacing, with coordinates far from the origin.
        for(int v=0;v<rectangles.GetNV();++v) {
            auto* p=rectangles.GetVertex(v);p[0]=2.+p[0]*p[0]*100.;p[1]=-3.+p[1]*p[1]*80.;
        }
        compare(rectangles,1,true);
        auto distorted=mfem::Mesh::MakeCartesian2D(3,4,mfem::Element::QUADRILATERAL,true,.009,.009,false);
        for(int v=0;v<distorted.GetNV();++v) {
            auto* p=distorted.GetVertex(v);p[0]+=p[1]*.3;p[1]+=p[0]*p[0];
        }
        compare(distorted,1,true);
        auto quadratics=mfem::Mesh::MakeCartesian2D(2,3,mfem::Element::QUADRILATERAL,true,.009,.009,false);
        compare(quadratics,2,false);
        auto triangles=mfem::Mesh::MakeCartesian2D(2,3,mfem::Element::TRIANGLE,true,.009,.009,false);
        compare(triangles,1,false);
        std::cout<<"Q1 scalar gathers, explicit DOF invalidation and fallbacks match MFEM bit for bit\n";
    } catch(const std::exception& error) {std::cerr<<error.what()<<'\n';return 1;}
}
