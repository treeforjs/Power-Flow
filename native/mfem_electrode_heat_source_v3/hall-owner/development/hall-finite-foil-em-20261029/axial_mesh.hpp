#pragma once
#include "foil_em.hpp"

namespace axial {
using foil::Vec;using foil::V3;using foil::dot;using foil::norm;using foil::maximum;
using foil::subtract;using foil::finite;using foil::eps;using foil::mu;using foil::Z0;
struct Edge {int a,b,axis;double length;};
struct Cell {std::array<int,12> edges;int i,j,k,body;double volume;};
// SI 3D Yee complex. The physical foil and ideal short are interior to the
// terminated air box. No out-of-plane derivative is discarded.
struct Mesh {
    int nx,ny,nz,nz_foil,kshort,kfeed,nn,ne,nb,nex,ney,nbx,nby;
    double side,depth=.025,dx,dy,dz,zmin;
    Vec M,H,Y,VE,VB,VN;std::vector<Edge>edges;std::vector<foil::CurlRow>C;
    std::vector<Cell>cells;std::vector<int>node_body;
    std::vector<bool>short_node,fixed_E,fixed_A,outer_node;
    Mesh(int n=45,int axial_cells=4,double width=.009,int air_cells=1)
      :nx(n),ny(n),nz(axial_cells+2*air_cells),nz_foil(axial_cells),kshort(air_cells),
       kfeed(air_cells+axial_cells),nn((n+1)*(n+1)*(nz+1)),
       ne(n*(n+1)*(nz+1)+(n+1)*n*(nz+1)+(n+1)*(n+1)*nz),
       nb((n+1)*n*nz+n*(n+1)*nz+n*n*(nz+1)),
       nex(n*(n+1)*(nz+1)),ney((n+1)*n*(nz+1)),
       nbx((n+1)*n*nz),nby(n*(n+1)*nz),side(width),dx(width/n),dy(width/n),
       dz(depth/axial_cells),zmin(-depth/2-air_cells*dz),
       M(ne),H(nb),Y(ne),VE(ne),VB(nb),VN(nn),edges(ne),C(nb),
       node_body(nn,-1),short_node(nn),fixed_E(ne),fixed_A(ne),outer_node(nn) {
        hall::require(n>=9&&n<=500&&axial_cells>=2&&air_cells>=1&&width>.0054,"invalid axial mesh");
        for(double x:{-.0027,-.0025,.0025,.0027})
            hall::require(std::abs((x+side/2)/dx-std::round((x+side/2)/dx))<1e-8,"unaligned physical foil");
        for(int k=0;k<=nz;++k)for(int j=0;j<=ny;++j)for(int i=0;i<=nx;++i){
            int n0=node(i,j,k);double x=i*dx-side/2,y=j*dy-side/2;
            VN[n0]=dual(i,nx,dx)*dual(j,ny,dy)*dual(k,nz,dz);
            outer_node[n0]=i==0||i==nx||j==0||j==ny||k==0||k==nz;
            bool inx=std::abs(x)<=.0025+1e-12;
            if(inx&&k>=kshort&&k<=kfeed){
                if(y>=-.0027-1e-12&&y<=-.0025+1e-12)node_body[n0]=0;
                if(y>=.0025-1e-12&&y<=.0027+1e-12)node_body[n0]=1;
            }
            short_node[n0]=k==kshort&&inx&&std::abs(y)<=.0027+1e-12;
        }
        auto edge=[&](int e,int a,int b,int axis,double length,double volume,bool outer,double boundary_area){
            edges[e]={a,b,axis,length};VE[e]=volume;M[e]=eps*volume;
            Y[e]=boundary_area/Z0;fixed_A[e]=outer;
            fixed_E[e]=axis!=2&&short_node[a]&&short_node[b];
        };
        for(int k=0;k<=nz;++k)for(int j=0;j<=ny;++j)for(int i=0;i<nx;++i){
            double ay=dual(j,ny,dy),az=dual(k,nz,dz);
            edge(ex(i,j,k),node(i,j,k),node(i+1,j,k),0,dx,dx*ay*az,
                 j==0||j==ny||k==0||k==nz,(j==0||j==ny?dx*az:0)+(k==0||k==nz?dx*ay:0));
        }
        for(int k=0;k<=nz;++k)for(int j=0;j<ny;++j)for(int i=0;i<=nx;++i){
            double ax=dual(i,nx,dx),az=dual(k,nz,dz);
            edge(ey(i,j,k),node(i,j,k),node(i,j+1,k),1,dy,ax*dy*az,
                 i==0||i==nx||k==0||k==nz,(i==0||i==nx?dy*az:0)+(k==0||k==nz?ax*dy:0));
        }
        for(int k=0;k<nz;++k)for(int j=0;j<=ny;++j)for(int i=0;i<=nx;++i){
            double ax=dual(i,nx,dx),ay=dual(j,ny,dy);
            edge(ez(i,j,k),node(i,j,k),node(i,j,k+1),2,dz,ax*ay*dz,
                 i==0||i==nx||j==0||j==ny,(i==0||i==nx?ay*dz:0)+(j==0||j==ny?ax*dz:0));
        }
        for(int k=0;k<nz;++k)for(int j=0;j<ny;++j)for(int i=0;i<=nx;++i){
            int b=bx(i,j,k);VB[b]=dual(i,nx,dx)*dy*dz;
            C[b]={{ez(i,j+1,k),ez(i,j,k),ey(i,j,k+1),ey(i,j,k)},{1/dy,-1/dy,-1/dz,1/dz},4};
        }
        for(int k=0;k<nz;++k)for(int j=0;j<=ny;++j)for(int i=0;i<nx;++i){
            int b=by(i,j,k);VB[b]=dx*dual(j,ny,dy)*dz;
            C[b]={{ex(i,j,k+1),ex(i,j,k),ez(i+1,j,k),ez(i,j,k)},{1/dz,-1/dz,-1/dx,1/dx},4};
        }
        for(int k=0;k<=nz;++k)for(int j=0;j<ny;++j)for(int i=0;i<nx;++i){
            int b=bz(i,j,k);VB[b]=dx*dy*dual(k,nz,dz);
            C[b]={{ey(i+1,j,k),ey(i,j,k),ex(i,j+1,k),ex(i,j,k)},{1/dx,-1/dx,-1/dy,1/dy},4};
        }
        for(int b=0;b<nb;++b)H[b]=VB[b]/mu;
        for(int k=0;k<nz;++k)for(int j=0;j<ny;++j)for(int i=0;i<nx;++i){
            double x=(i+.5)*dx-side/2,y=(j+.5)*dy-side/2;int body=-1;
            if(k>=kshort&&k<kfeed&&std::abs(x)<.0025&&std::abs(y)>.0025&&std::abs(y)<.0027)body=y>0?1:0;
            cells.push_back({{ex(i,j,k),ex(i,j+1,k),ex(i,j,k+1),ex(i,j+1,k+1),
                              ey(i,j,k),ey(i+1,j,k),ey(i,j,k+1),ey(i+1,j,k+1),
                              ez(i,j,k),ez(i+1,j,k),ez(i,j+1,k),ez(i+1,j+1,k)},i,j,k,body,dx*dy*dz});
        }
    }
    static double dual(int i,int n,double h){return (i==0||i==n)?.5*h:h;}
    int node(int i,int j,int k)const{return (k*(ny+1)+j)*(nx+1)+i;}
    int ex(int i,int j,int k)const{return (k*(ny+1)+j)*nx+i;}
    int ey(int i,int j,int k)const{return nex+(k*ny+j)*(nx+1)+i;}
    int ez(int i,int j,int k)const{return nex+ney+(k*(ny+1)+j)*(nx+1)+i;}
    int bx(int i,int j,int k)const{return (k*ny+j)*(nx+1)+i;}
    int by(int i,int j,int k)const{return nbx+(k*(ny+1)+j)*nx+i;}
    int bz(int i,int j,int k)const{return nbx+nby+(k*ny+j)*nx+i;}
    V3 node_position(int n)const{int i=n%(nx+1),j=(n/(nx+1))%(ny+1),k=n/((nx+1)*(ny+1));return {i*dx-side/2,j*dy-side/2,zmin+k*dz};}
    Vec gradient(const Vec&p)const{Vec v(ne);for(int e=0;e<ne;++e)v[e]=(p[edges[e].b]-p[edges[e].a])/edges[e].length;return v;}
    Vec gradient_transpose(const Vec&v)const{Vec p(nn);for(int e=0;e<ne;++e){double q=v[e]/edges[e].length;p[edges[e].a]-=q;p[edges[e].b]+=q;}return p;}
    Vec curl(const Vec&E)const{Vec B(nb);for(int b=0;b<nb;++b)for(int k=0;k<4;++k)B[b]+=C[b].coefficient[k]*E[C[b].index[k]];return B;}
    Vec curl_transpose(const Vec&B)const{Vec E(ne);for(int b=0;b<nb;++b)for(int k=0;k<4;++k)E[C[b].index[k]]+=C[b].coefficient[k]*B[b];return E;}
    Vec curl_adjoint(Vec B)const{for(int b=0;b<nb;++b)B[b]*=H[b];return curl_transpose(B);}
    Vec magnetic_divergence(const Vec&B)const{Vec p(nx*ny*nz);for(int k=0;k<nz;++k)for(int j=0;j<ny;++j)for(int i=0;i<nx;++i)p[(k*ny+j)*nx+i]=(B[bx(i+1,j,k)]-B[bx(i,j,k)])/dx+(B[by(i,j+1,k)]-B[by(i,j,k)])/dy+(B[bz(i,j,k+1)]-B[bz(i,j,k)])/dz;return p;}
    Vec field_charge(Vec E)const{for(int e=0;e<ne;++e)E[e]*=-M[e];return gradient_transpose(E);}
    Vec conduction(const Vec&sigma_cells)const{
        hall::require(sigma_cells.size()==cells.size(),"cell material binding mismatch");Vec s(ne);
        for(std::size_t c=0;c<cells.size();++c){hall::require(std::isfinite(sigma_cells[c])&&sigma_cells[c]>=0&&
            (cells[c].body>=0||sigma_cells[c]==0),"nonfinite or vacuum conductivity");
            for(int e:cells[c].edges)s[e]+=.25*cells[c].volume*sigma_cells[c];}
        return s;
    }
    Vec cell_joule_W(const Vec&sigma_cells,const Vec&E)const{Vec p(cells.size());for(std::size_t c=0;c<cells.size();++c)for(int e:cells[c].edges)p[c]+=.25*cells[c].volume*sigma_cells[c]*E[e]*E[e];return p;}
    Vec feed_shape()const{
        Vec f(ne);double weight=0;
        int i0=std::lround((side/2-.0025)/dx),i1=std::lround((side/2+.0025)/dx);
        int j0=std::lround((side/2+.0025)/dy),j1=std::lround((side/2+.0027)/dy);
        for(int j=j0;j<=j1;++j)for(int i=i0;i<=i1;++i){
            double w=((i==i0||i==i1)?.5:1)*((j==j0||j==j1)?.5:1)*dx*dy/(.005*.0002);
            weight+=w;for(int y=ny-j;y<j;++y)f[ey(i,y,kfeed)]-=w*dy;
        }
        hall::require(std::abs(weight-1)<=2e-14,"feed contact quadrature normalization");return f;
    }
    // Positive means inflow into nodes above the mid-gap and between a chosen
    // axial cut and feed. This exact incidence measurement includes transverse
    // breakdown, external faces, displacement and local accumulation.
    double upper_region_inflow(const Vec&load,int axial_cut)const{
        Vec membership(nn);for(int k=axial_cut+1;k<=kfeed;++k)for(int j=0;j<=ny;++j)
            if(j*dy-side/2>0)for(int i=0;i<=nx;++i)membership[node(i,j,k)]=1;
        return dot(gradient(membership),load);
    }
    Json json()const{return {{"schema","axial-finite-foil-yee-v1"},{"nx",nx},{"ny",ny},{"nz",nz},{"foil_axial_cells",nz_foil},{"air_cells_per_end",kshort},{"outer_z_min_m",zmin},{"outer_z_max_m",zmin+nz*dz},{"width_m",side},{"foil_depth_m",depth},{"dx_m",dx},{"dy_m",dy},{"dz_m",dz},{"short_z_m",-depth/2},{"feed_z_m",depth/2},{"field_dimensions",3},{"z_invariant",false},{"short","ideal PEC sheet; no short heating"},{"feed","distributed impressed branch upper to lower at positive-z end"},{"exterior_converged",false}};}
};
} // namespace axial
