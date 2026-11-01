#pragma once
#include "../hall-coupled-transport-20261029/coupled_pic.hpp"
#include <functional>
#include <numeric>
#include <set>

namespace foil {
using Vec=std::vector<double>;using V3=std::array<double,3>;
constexpr double eps=hall::eps,mu=1/(eps*pic_c*pic_c),Z0=mu*pic_c;
inline double dot(const Vec&a,const Vec&b){long double s=0;for(std::size_t i=0;i<a.size();++i)s+=(long double)a[i]*b[i];return double(s);}
inline double norm(const Vec&a){return std::sqrt(dot(a,a));}
inline double maximum(const Vec&a){double s=0;for(double x:a)s=std::max(s,std::abs(x));return s;}
inline Vec subtract(Vec a,const Vec&b){for(std::size_t i=0;i<a.size();++i)a[i]-=b[i];return a;}
inline bool finite(const Vec&v){return std::all_of(v.begin(),v.end(),[](double x){return std::isfinite(x);});}
struct CurlRow {std::array<int,4>index{};std::array<double,4>coefficient{};int count=0;};
struct Mesh {
    int nx,ny,nn,nex,ney,ne,nb;double width,height,dx,dy,depth;
    bool foils;std::vector<int>owner;std::vector<bool>fixed_E,outer_node;
    Vec volume_E,volume_B,M,H,Y;std::vector<CurlRow>C;
    Mesh(int n,double side=.009,double represented_depth=.025,bool conductors=true):nx(n),ny(n),nn((n+1)*(n+1)),nex(n*(n+1)),ney(n*(n+1)),ne(nex+ney+nn),nb(ney+nex+n*n),width(side),height(side),dx(side/n),dy(side/n),depth(represented_depth),foils(conductors),owner(nn,-1),fixed_E(ne),outer_node(nn),volume_E(ne),volume_B(nb),M(ne),H(nb),Y(ne),C(nb){
        hall::require(n>=4&&n<=1000&&side>.0054&&depth>0,"invalid finite foil mesh");
        if(foils)for(double pos:{-.0027,-.0025,.0025,.0027})hall::require(std::abs((pos+side/2)/dx-std::round((pos+side/2)/dx))<1e-8,"foil boundaries must align with the resolved mesh");
        for(int j=0;j<=ny;++j)for(int i=0;i<=nx;++i){const int p=node(i,j);double x=i*dx-width/2,y=j*dy-height/2;outer_node[p]=i==0||j==0||i==nx||j==ny;
            if(foils&&x>=-.0025-1e-12&&x<=.0025+1e-12){if(y>=-.0027-1e-12&&y<=-.0025+1e-12)owner[p]=0;if(y>=.0025-1e-12&&y<=.0027+1e-12)owner[p]=1;}
            volume_E[ez(i,j)]=dual_x(i)*dual_y(j)*depth;fixed_E[ez(i,j)]=owner[p]>=0;
            if(i==0||i==nx)Y[ez(i,j)]+=dual_y(j)*depth/Z0;if(j==0||j==ny)Y[ez(i,j)]+=dual_x(i)*depth/Z0;
        }
        for(int j=0;j<=ny;++j)for(int i=0;i<nx;++i){int e=ex(i,j);volume_E[e]=dx*dual_y(j)*depth;fixed_E[e]=owner[node(i,j)]>=0&&owner[node(i,j)]==owner[node(i+1,j)];if(j==0||j==ny)Y[e]=dx*depth/Z0;
            int b=by(i,j);volume_B[b]=volume_E[e];C[b]={{ez(i+1,j),ez(i,j),0,0},{-1/dx,1/dx,0,0},2};}
        for(int j=0;j<ny;++j)for(int i=0;i<=nx;++i){int e=ey(i,j);volume_E[e]=dual_x(i)*dy*depth;fixed_E[e]=owner[node(i,j)]>=0&&owner[node(i,j)]==owner[node(i,j+1)];if(i==0||i==nx)Y[e]=dy*depth/Z0;
            int b=bx(i,j);volume_B[b]=volume_E[e];C[b]={{ez(i,j+1),ez(i,j),0,0},{1/dy,-1/dy,0,0},2};}
        for(int j=0;j<ny;++j)for(int i=0;i<nx;++i){int b=bz(i,j);volume_B[b]=dx*dy*depth;C[b]={{ey(i+1,j),ey(i,j),ex(i,j+1),ex(i,j)},{1/dx,-1/dx,-1/dy,1/dy},4};}
        for(int e=0;e<ne;++e)M[e]=eps*volume_E[e];for(int b=0;b<nb;++b)H[b]=volume_B[b]/mu;
    }
    int node(int i,int j)const{return j*(nx+1)+i;}
    int ex(int i,int j)const{return j*nx+i;}
    int ey(int i,int j)const{return nex+j*(nx+1)+i;}
    int ez(int i,int j)const{return nex+ney+node(i,j);}
    int bx(int i,int j)const{return j*(nx+1)+i;}
    int by(int i,int j)const{return ney+j*nx+i;}
    int bz(int i,int j)const{return ney+nex+j*nx+i;}
    double dual_x(int i)const{return (i==0||i==nx)?.5*dx:dx;}
    double dual_y(int j)const{return (j==0||j==ny)?.5*dy:dy;}
    Vec gradient(const Vec&phi)const{Vec E(ne);for(int j=0;j<=ny;++j)for(int i=0;i<nx;++i)E[ex(i,j)]=(phi[node(i+1,j)]-phi[node(i,j)])/dx;for(int j=0;j<ny;++j)for(int i=0;i<=nx;++i)E[ey(i,j)]=(phi[node(i,j+1)]-phi[node(i,j)])/dy;return E;}
    Vec gradient_transpose(const Vec&v)const{Vec p(nn);for(int j=0;j<=ny;++j)for(int i=0;i<nx;++i){double q=v[ex(i,j)]/dx;p[node(i,j)]-=q;p[node(i+1,j)]+=q;}for(int j=0;j<ny;++j)for(int i=0;i<=nx;++i){double q=v[ey(i,j)]/dy;p[node(i,j)]-=q;p[node(i,j+1)]+=q;}return p;}
    Vec curl(const Vec&E)const{Vec B(nb);for(int b=0;b<nb;++b)for(int k=0;k<C[b].count;++k)B[b]+=C[b].coefficient[k]*E[C[b].index[k]];return B;}
    Vec curl_transpose(const Vec&B)const{Vec E(ne);for(int b=0;b<nb;++b)for(int k=0;k<C[b].count;++k)E[C[b].index[k]]+=C[b].coefficient[k]*B[b];return E;}
    Vec curl_adjoint(const Vec&B)const{Vec v=B;for(int b=0;b<nb;++b)v[b]*=H[b];return curl_transpose(v);}
    Vec magnetic_divergence(const Vec&B)const{Vec r(nx*ny);for(int j=0;j<ny;++j)for(int i=0;i<nx;++i)r[j*nx+i]=(B[bx(i+1,j)]-B[bx(i,j)])/dx+(B[by(i,j+1)]-B[by(i,j)])/dy;return r;}
    Vec field_charge(const Vec&E)const{Vec v=E;for(int e=0;e<ne;++e)v[e]*=-M[e];return gradient_transpose(v);}
    Json json()const{return {{"schema","finite-foil-compatible-mesh-v1"},{"nx",nx},{"ny",ny},{"width_m",width},{"height_m",height},{"depth_m",depth},{"dx_m",dx},{"dy_m",dy},{"foil_width_m",.005},{"foil_thickness_m",.0002},{"gap_m",.005},{"geometry_fixed",true},{"field_dimensions",2},{"velocity_dimensions",3},{"z_invariant",true},{"axial_instabilities_resolved",false}};}
};

// Preconditioned CG with independently recomputed residual and at most two
// residual-replacement corrections. The final bound is unchanged by refinement.
inline Vec solve(const std::function<Vec(const Vec&)>&apply,const Vec&b,const Vec&diagonal,int maxit=3000){
    Vec x(b.size()),r=b;const double target=std::max(1e-30,2e-13*norm(b));
    for(int correction=0;correction<3;++correction){
        Vec z(r.size());for(std::size_t i=0;i<r.size();++i)z[i]=r[i]/diagonal[i];Vec p=z;double rz=dot(r,z);
        for(int iteration=0;iteration<maxit&&norm(r)>target*.2;++iteration){auto Ap=apply(p);double denominator=dot(p,Ap);hall::require(denominator>0&&std::isfinite(denominator),"non-SPD field operator");double alpha=rz/denominator;for(std::size_t i=0;i<r.size();++i){x[i]+=alpha*p[i];r[i]-=alpha*Ap[i];z[i]=r[i]/diagonal[i];}double nrz=dot(r,z),beta=nrz/rz;for(std::size_t i=0;i<r.size();++i)p[i]=z[i]+beta*p[i];rz=nrz;}
        r=subtract(b,apply(x));if(norm(r)<=target){hall::require(finite(x),"nonfinite field solve");return x;}
    }
    throw std::runtime_error("field solve independent residual failed");
}
struct State {Vec E,B;double time_s=7.4978805549199997e-10,voltage_V=0;std::uint64_t revision=0;};
struct FieldReceipt {
    double field_change_J=0,particle_work_J=0,outward_boundary_energy_J=0,background_supply_work_J=0,circuit_supply_work_J=0,circuit_charge_C=0,balance_J=0,voltage_error_V=0;
    Json json()const{return {{"field_change_J",field_change_J},{"work_to_particles_J",particle_work_J},{"outward_boundary_energy_J",outward_boundary_energy_J},{"background_supply_work_J",background_supply_work_J},{"circuit_supply_work_J",circuit_supply_work_J},{"circuit_charge_C",circuit_charge_C},{"balance_J",balance_J},{"voltage_error_V",voltage_error_V},{"energy_identity","dU + work_to_particles + outward_boundary - background_supply - circuit_supply = 0"}};}
};
class Fields {
    const Mesh&g;double dt;Vec diagonal,sqrtM,circuit_response;double circuit_denominator=0;
    Vec apply(const Vec&y)const{Vec E(y.size());for(int e=0;e<g.ne;++e)if(!g.fixed_E[e])E[e]=y[e]/sqrtM[e];auto cc=g.curl_adjoint(g.curl(E));Vec r(y.size());for(int e=0;e<g.ne;++e)r[e]=g.fixed_E[e]?y[e]:y[e]+(.5*dt*g.Y[e]*E[e]+.25*dt*dt*cc[e])/sqrtM[e];return r;}
    Vec invert(const Vec&rhs)const{Vec b=rhs;for(int e=0;e<g.ne;++e)b[e]=g.fixed_E[e]?0:b[e]/sqrtM[e];auto y=solve([&](const Vec&v){return apply(v);},b,diagonal,100);for(int e=0;e<g.ne;++e)y[e]/=sqrtM[e];return y;}
public:
    Vec weighting_potential,port_shape,background_E,background_B,background_force;
    double capacitance_F=0,dc_current_A=0;
    Fields(const Mesh&mesh,double step,double voltage,double dc_current):g(mesh),dt(step),diagonal(g.ne,1),sqrtM(g.ne),weighting_potential(g.nn),port_shape(g.ne),background_E(g.ne),background_B(g.nb),background_force(g.ne),dc_current_A(dc_current){
        hall::require(dt>0&&dt<=1e-14,"10fs full-field guard");
        for(int e=0;e<g.ne;++e){sqrtM[e]=std::sqrt(g.M[e]);if(!g.fixed_E[e])diagonal[e]+=.5*dt*g.Y[e]/g.M[e];}
        for(int b=0;b<g.nb;++b)for(int k=0;k<g.C[b].count;++k){int e=g.C[b].index[k];if(!g.fixed_E[e])diagonal[e]+=.25*dt*dt*g.H[b]*g.C[b].coefficient[k]*g.C[b].coefficient[k]/g.M[e];}
        if(g.foils){
            auto lap=[&](const Vec&p){auto v=g.gradient(p);for(int e=0;e<g.ne;++e)v[e]*=g.M[e];return g.gradient_transpose(v);};
            Vec known(g.nn),diag(g.nn);for(int p=0;p<g.nn;++p)if(g.owner[p]==1)known[p]=1;
            for(int j=0;j<=g.ny;++j)for(int i=0;i<g.nx;++i){double k=g.M[g.ex(i,j)]/(g.dx*g.dx);diag[g.node(i,j)]+=k;diag[g.node(i+1,j)]+=k;}
            for(int j=0;j<g.ny;++j)for(int i=0;i<=g.nx;++i){double k=g.M[g.ey(i,j)]/(g.dy*g.dy);diag[g.node(i,j)]+=k;diag[g.node(i,j+1)]+=k;}
            auto rhs=lap(known);for(int p=0;p<g.nn;++p){rhs[p]=g.owner[p]>=0?0:-rhs[p];if(g.owner[p]>=0)diag[p]=1;}
            auto x=solve([&](const Vec&v){auto w=v;for(int p=0;p<g.nn;++p)if(g.owner[p]>=0)w[p]=0;auto a=lap(w);for(int p=0;p<g.nn;++p)if(g.owner[p]>=0)a[p]=v[p];return a;},rhs,diag);
            for(int p=0;p<g.nn;++p)weighting_potential[p]=known[p]+x[p];port_shape=g.gradient(weighting_potential);for(double&v:port_shape)v=-v;
            Vec rhs_port(g.ne);for(int e=0;e<g.ne;++e){background_E[e]=voltage*port_shape[e];rhs_port[e]=g.M[e]*port_shape[e];capacitance_F+=g.M[e]*port_shape[e]*port_shape[e];}
            circuit_response=invert(rhs_port);for(int e=0;e<g.ne;++e)circuit_denominator+=port_shape[e]*g.M[e]*circuit_response[e];hall::require(circuit_denominator>0,"invalid circuit capacitance");
        }
        if(dc_current!=0){
            hall::require(g.foils,"DC support requires foil geometry");Vec J(g.nn),diag(g.nn);
            for(int j=0;j<g.ny;++j)for(int i=0;i<g.nx;++i){double x=(i+.5)*g.dx-g.width/2,y=(j+.5)*g.dy-g.height/2;if(std::abs(x)<.0025&&std::abs(y)>.0025&&std::abs(y)<.0027){double current=(y>0?1:-1)*dc_current/(.005*.0002)*g.dx*g.dy*g.depth/4;for(int dj:{0,1})for(int di:{0,1})J[g.node(i+di,j+dj)]+=current;}}
            auto lap=[&](const Vec&p){Vec A(g.ne);for(int n=0;n<g.nn;++n)A[g.nex+g.ney+n]=p[n];auto r=g.curl_adjoint(g.curl(A));return Vec(r.begin()+g.nex+g.ney,r.end());};
            for(int b=0;b<g.ney+g.nex;++b)for(int k=0;k<g.C[b].count;++k){int n=g.C[b].index[k]-g.nex-g.ney;diag[n]+=g.H[b]*g.C[b].coefficient[k]*g.C[b].coefficient[k];}
            for(int n=0;n<g.nn;++n)if(g.outer_node[n]){J[n]=0;diag[n]=1;}
            auto az=solve([&](const Vec&p){auto x=p;for(int n=0;n<g.nn;++n)if(g.outer_node[n])x[n]=0;auto y=lap(x);for(int n=0;n<g.nn;++n)if(g.outer_node[n])y[n]=p[n];return y;},J,diag);
            Vec A(g.ne);for(int n=0;n<g.nn;++n)A[g.nex+g.ney+n]=az[n];background_B=g.curl(A);
        }
        auto source=g.curl_adjoint(background_B);for(int e=0;e<g.ne;++e)background_force[e]=g.Y[e]*background_E[e]-source[e];
    }
    State initial(double voltage)const{return {background_E,background_B,7.4978805549199997e-10,voltage,0};}
    std::pair<State,FieldReceipt>advance(const State&old,const Vec&J,double voltage)const{
        hall::require(int(J.size())==g.ne&&finite(J)&&finite(old.E)&&finite(old.B),"invalid field/current state");
        auto curl=g.curl(old.E);Vec Bcenter=old.B;for(int b=0;b<g.nb;++b)Bcenter[b]-=.5*dt*curl[b];auto adj=g.curl_adjoint(Bcenter);Vec rhs(g.ne);
        for(int e=0;e<g.ne;++e)rhs[e]=dt*(adj[e]-J[e]-g.Y[e]*old.E[e]+background_force[e]);auto delta=invert(rhs);
        double circuit=0;if(g.foils){double predicted=0;for(int e=0;e<g.ne;++e)predicted+=port_shape[e]*g.M[e]*delta[e];circuit=(capacitance_F*(voltage-old.voltage_V)-predicted)/circuit_denominator;for(int e=0;e<g.ne;++e)delta[e]+=circuit*circuit_response[e];}
        State next=old;Vec mid(g.ne);for(int e=0;e<g.ne;++e){next.E[e]=old.E[e]+delta[e];mid[e]=.5*(next.E[e]+old.E[e]);hall::require(!g.fixed_E[e]||next.E[e]==0,"PEC tangent field changed");}
        auto Cmid=g.curl(mid);for(int b=0;b<g.nb;++b)next.B[b]=old.B[b]-dt*Cmid[b];
        FieldReceipt receipt;for(int e=0;e<g.ne;++e){receipt.field_change_J+=g.M[e]*mid[e]*(next.E[e]-old.E[e]);receipt.particle_work_J+=dt*mid[e]*J[e];receipt.outward_boundary_energy_J+=dt*g.Y[e]*mid[e]*mid[e];receipt.background_supply_work_J+=dt*mid[e]*background_force[e];receipt.circuit_supply_work_J+=circuit*mid[e]*g.M[e]*port_shape[e];}
        for(int b=0;b<g.nb;++b)receipt.field_change_J+=g.H[b]*.5*(next.B[b]+old.B[b])*(next.B[b]-old.B[b]);
        receipt.circuit_charge_C=circuit*capacitance_F;double projected=0;for(int e=0;e<g.ne;++e)projected+=port_shape[e]*g.M[e]*next.E[e];if(g.foils)receipt.voltage_error_V=projected/capacitance_F-voltage;
        receipt.balance_J=receipt.field_change_J+receipt.particle_work_J+receipt.outward_boundary_energy_J-receipt.background_supply_work_J-receipt.circuit_supply_work_J;
        next.time_s=old.time_s+dt;next.voltage_V=voltage;next.revision=old.revision+1;return {std::move(next),receipt};
    }
    double energy(const State&s)const{double U=0;for(int e=0;e<g.ne;++e)U+=.5*g.M[e]*s.E[e]*s.E[e];for(int b=0;b<g.nb;++b)U+=.5*g.H[b]*s.B[b]*s.B[b];return U;}
};
}
