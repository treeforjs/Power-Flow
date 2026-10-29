#include "raikou/q1_scalar_gather.h"
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
using raikou::mfem_backend::Q1ScalarGather;
int main(int argc,char** argv) {
    if(argc!=2)return 2;
    auto mesh=mfem::Mesh::MakeCartesian2D(360,360,mfem::Element::QUADRILATERAL,true,.009,.009,false);
    mfem::H1_FECollection collection(1,2);mfem::FiniteElementSpace space(&mesh,&collection);
    mfem::GridFunction field(&space);
    for(int i=0;i<field.Size();++i)field[i]=std::sin(i*.2);
    Q1ScalarGather gather(field);
    if(!gather.available())return 3;
    std::ofstream out(argv[1]);out<<std::setprecision(17)<<"[\n";
    bool first=true;
    for(int repetition=0;repetition<5;++repetition)for(int which=0;which<2;++which) {
        const bool cached=(repetition+which)%2;
        const auto start=std::chrono::steady_clock::now();double checksum=0;
        for(int n=0;n<2000000;++n) {
            const int element=(n*31)%mesh.GetNE();
            mfem::IntegrationPoint point;point.Set2(double(n%997)/997,double(n%991)/991);
            if(cached) {const auto g=gather.gradient(element,point);checksum+=g[0]+g[1];}
            else {auto* tr=mesh.GetElementTransformation(element);tr->SetIntPoint(&point);mfem::Vector g(2);field.GetGradient(*tr,g);checksum+=g[0]+g[1];}
        }
        const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        if(!first)out<<",\n";first=false;
        out<<"{\"case\":\"360x360_Q1_gradient\",\"label\":\""<<(cached?"candidate":"baseline")<<"\",\"samples\":2000000,\"seconds\":"<<seconds<<",\"checksum\":"<<checksum<<'}';
    }
    out<<"\n]\n";
}
