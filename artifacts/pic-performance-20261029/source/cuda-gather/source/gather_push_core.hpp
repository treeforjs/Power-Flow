#pragma once
// Isolated research prototype. No source, deposition, drift, or solver work.
// One arithmetic implementation is compiled on host and device. In an eventual
// integration, pic_boris must call boris() below after retaining its exceptions;
// do not maintain a second force implementation in the application backend.
#if defined(__CUDACC__)
#define PIC_HD __host__ __device__
#else
#define PIC_HD
#endif
#if !defined(__CUDACC_RTC__)
#include <cmath>
#endif

namespace gather_push {
struct Vec { double v[3]{}; };
struct Cell { int dof[4]; double x[4], y[4]; };
// Cell local order is LL, LR, UR, UL, explicitly checked against MFEM nodes.
struct Grid { const Cell* cells; int nx, ny; };
struct Inputs { double x, y, z; Vec u; double qm; };
struct Outputs { Vec E, B, u; unsigned long long status; };
struct Tables {
    Grid self, native;
    const double *phi, *psi, *az0, *az1, *rate0, *rate1;
    const double *xs, *ys;
    double width, height, depth, short_z, fraction, voltage, dt;
};
enum Status : unsigned long long {
    success=0, invalid_input=1, outside_domain=2, singular_geometry=3,
    nonfinite_field=4, boris_overflow=5
};
PIC_HD inline bool finite(double x) {
#if defined(__CUDACC_RTC__)
    return isfinite(x);
#else
    return std::isfinite(x);
#endif
}
PIC_HD inline double root(double x) {
#if defined(__CUDACC_RTC__)
    return sqrt(x);
#else
    return std::sqrt(x);
#endif
}
PIC_HD inline bool finite(Vec x) {
    return finite(x.v[0]) && finite(x.v[1]) && finite(x.v[2]);
}
PIC_HD inline double norm2(Vec x) {
    return x.v[0]*x.v[0]+x.v[1]*x.v[1]+x.v[2]*x.v[2];
}
PIC_HD inline double gamma(Vec u) {
    constexpr double c=299792458.;
    return root(1+norm2(u)/(c*c));
}
PIC_HD inline Vec cross(Vec a, Vec b) {
    return {{a.v[1]*b.v[2]-a.v[2]*b.v[1],
             a.v[2]*b.v[0]-a.v[0]*b.v[2],
             a.v[0]*b.v[1]-a.v[1]*b.v[0]}};
}
// Arithmetic and ordering transcribed from the sealed pic_boris. Proper u,
// not ordinary v; macro weight does not enter q/m. Host tests compare result
// bits with the original function. Device uses --fmad=false, FP64 throughout.
PIC_HD inline unsigned long long boris(Vec u,double qm,Vec E,Vec B,double dt,Vec& result) {
    if(!finite(u) || !finite(E) || !finite(B) || !finite(qm) || !finite(dt) || dt<=0)
        return invalid_input;
    Vec minus{},t{},prime{},plus{},s{};
    for(int i=0;i<3;++i) minus.v[i]=u.v[i]+.5*dt*qm*E.v[i];
    for(int i=0;i<3;++i) t.v[i]=.5*dt*qm*B.v[i]/gamma(minus);
    const auto mt=cross(minus,t);
    for(int i=0;i<3;++i) {prime.v[i]=minus.v[i]+mt.v[i];s.v[i]=2*t.v[i]/(1+norm2(t));}
    const auto ps=cross(prime,s);
    for(int i=0;i<3;++i) {plus.v[i]=minus.v[i]+ps.v[i];result.v[i]=plus.v[i]+.5*dt*qm*E.v[i];}
    return finite(result) && finite(norm2(result)) ? success : boris_overflow;
}
struct Sample { double value, dx, dy; unsigned long long status; };
// Degree-one specialization of MFEM4.10 Poly_1D::Basis::Eval (Barycentric).
// Retain its separate value/derivative arithmetic, including divisions. The
// direct polynomial (1-x,x) has observably different rounding near cancellation.
// Derived from MFEM, copyright2010-2026 LLNS, BSD3-Clause; see MFEM-LICENSE.
PIC_HD inline void basis1d(double y,double* u,double* d=nullptr) {
    const int k=y>=.5 ? 1 : 0;
    const double lk=k ? y : y-1,l=lk*(k ? y-1 : y);
    if(!d) {
        if(k) {u[0]=l*(-1)/y;u[1]=lk;} else {u[0]=lk*(-1);u[1]=l/(y-1);}
        return;
    }
    const double si=1/(k ? y : y-1);
    if(k) {u[0]=l*si*(-1);u[1]=lk;} else {u[0]=lk*(-1);u[1]=l*si;}
    const double lp=l*si+lk;
    if(k) {d[0]=(lp*(-1)-u[0])/y;d[1]=si*u[1];}
    else {d[0]=si*u[0];d[1]=(lp-u[1])/(y-1);}
}
PIC_HD inline Sample sample(Grid grid,const double* field,int element,double x,double y) {
    if(!grid.cells || !field || grid.nx<1 || grid.nx>4096 || grid.ny<1 || grid.ny>4096 ||
       element<0 || element>=grid.nx*grid.ny || !finite(x) || !finite(y) || x<0 || x>1 || y<0 || y>1)
        return {0,0,0,outside_domain};
    const auto& cell=grid.cells[element];
    double vx[2],vy[2],gx[2],gy[2],ddx[2],ddy[2];
    basis1d(x,vx);basis1d(y,vy);basis1d(x,gx,ddx);basis1d(y,gy,ddy);
    const double s[4]={vx[0]*vy[0],vx[1]*vy[0],vx[1]*vy[1],vx[0]*vy[1]};
    const double dx[4]={ddx[0]*gy[0],ddx[1]*gy[0],ddx[1]*gy[1],ddx[0]*gy[1]},
                 dy[4]={gx[0]*ddy[0],gx[1]*ddy[0],gx[1]*ddy[1],gx[0]*ddy[1]};
    // Mesh geometry uses MFEM's BiLinear2DFiniteElement, not its H1 basis.
    const double jx[4]={-1+y,1-y,y,-y},jy[4]={-1+x,-x,x,1-x};
    double value=0,rx=0,ry=0,j00=0,j01=0,j10=0,j11=0;
    for(int k=0;k<4;++k) {
        const double p=field[cell.dof[k]];
        value+=s[k]*p;rx+=dx[k]*p;ry+=dy[k]*p;
        j00+=jx[k]*cell.x[k];j01+=jy[k]*cell.x[k];
        j10+=jx[k]*cell.y[k];j11+=jy[k]*cell.y[k];
    }
    const double determinant=j00*j11-j01*j10;
    if(!finite(determinant) || determinant<=0) return {0,0,0,singular_geometry};
    // General bilinear physical map, no constant-J or uniform-native shortcut.
    const double inverse=1/determinant;
    const double a=j11*inverse,b=(-j01)*inverse,c=(-j10)*inverse,d=j00*inverse;
    double px=0,py=0;px+=rx*a;px+=ry*c;py+=rx*b;py+=ry*d;
    Sample result{value,px,py,success};
    if(!finite(value) || !finite(result.dx) || !finite(result.dy)) result.status=nonfinite_field;
    return result;
}
PIC_HD inline int bin(const double* axes,int cells,double x) {
    if(!finite(x) || x<axes[0] || x>axes[cells]) return -1;
    // std::upper_bound semantics, including an exact internal grid edge.
    int first=0,count=cells+1;
    while(count>0) {
        const int step=count/2,mid=first+step;
        if(!(x<axes[mid])) {first=mid+1;count-=step+1;} else count=step;
    }
    const int index=first-1;
    return index<0 ? 0 : (index>=cells ? cells-1 : index);
}
PIC_HD inline Outputs evaluate(const Tables& t,Inputs p) {
    Outputs out{};
    if(!finite(p.x) || !finite(p.y) || !finite(p.z) || !finite(p.u) || !finite(p.qm) ||
       !finite(t.dt) || t.dt<=0 || !finite(t.width) || t.width<=0 || !finite(t.height) || t.height<=0 ||
       !finite(t.depth) || t.depth<=0 || !finite(t.short_z) || !finite(t.voltage) || !finite(t.fraction) ||
       t.fraction<0 || t.fraction>1 || t.self.nx<1 || t.self.nx>4096 || t.self.ny<1 || t.self.ny>4096 ||
       t.native.nx<1 || t.native.nx>4096 || t.native.ny<1 || t.native.ny>4096 || !t.xs || !t.ys)
       {out.status=invalid_input;return out;}
    if(p.x<0 || p.x>t.width || p.y<0 || p.y>t.height) {out.status=outside_domain;return out;}
    // Same operation order as GapPic::locate; upper endpoint stays in last cell.
    const double sx=p.x/t.width*t.self.nx,sy=p.y/t.height*t.self.ny;
    int ix=static_cast<int>(sx),iy=static_cast<int>(sy);
    if(ix>=t.self.nx)ix=t.self.nx-1;
    if(iy>=t.self.ny)iy=t.self.ny-1;
    const int e=iy*t.self.nx+ix;
    const auto self=sample(t.self,t.phi,e,sx-ix,sy-iy);
    const auto unit=sample(t.self,t.psi,e,sx-ix,sy-iy);
    if(self.status || unit.status) {out.status=self.status ? self.status : unit.status;return out;}
    // NativeField's translation is part of the existing model.
    const double px=p.x-.0045,py=p.y-.0045;
    const int ni=bin(t.xs,t.native.nx,px),nj=bin(t.ys,t.native.ny,py);
    if(ni<0 || nj<0) {out.status=outside_domain;return out;}
    const int ne=nj*t.native.nx+ni;
    const double nx=(px-t.xs[ni])/(t.xs[ni+1]-t.xs[ni]),
                 ny=(py-t.ys[nj])/(t.ys[nj+1]-t.ys[nj]);
    const auto a0=sample(t.native,t.az0,ne,nx,ny),a1=sample(t.native,t.az1,ne,nx,ny),
               r0=sample(t.native,t.rate0,ne,nx,ny),r1=sample(t.native,t.rate1,ne,nx,ny);
    if(a0.status || a1.status || r0.status || r1.status) {out.status=nonfinite_field;return out;}
    out.E={{-self.dx,-self.dy,0}};
    // Preserve the callback sequence: self E, native E, axial E correction.
    out.E.v[0]+=0;out.E.v[1]+=0;
    out.E.v[2]+=(1-t.fraction)*(-r0.value)+t.fraction*(-r1.value);
    const double g=(p.z-t.short_z)/t.depth;
    Vec correction{{-unit.dx,-unit.dy,0}};
    for(int k=0;k<2;++k)correction.v[k]*=t.voltage*(2*g-1);
    correction.v[2]=-2*t.voltage/t.depth*(unit.value-.5);
    for(int k=0;k<3;++k)out.E.v[k]+=correction.v[k];
    // The magnetic callback REPLACES uniform B in the original advance loop.
    out.B={{(1-t.fraction)*a0.dy+t.fraction*a1.dy,
            (1-t.fraction)*(-a0.dx)+t.fraction*(-a1.dx),0}};
    out.status=boris(p.u,p.qm,out.E,out.B,t.dt,out.u);
    return out;
}
} // namespace gather_push
