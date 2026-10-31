// Hall vacuum port load from the sealed Hall axial field operators (hall-finite-foil-em-20261029, the 8dddec0-era
// axial::Mesh/dc_initial/Fields used unchanged by Hall corrective v3/v3.1; read-only include).
//   mode dc   : DC magnetostatic inductance L_dc = 2 W_B / I^2 and DC resistance R_dc = V_T / I from dc_initial.
//   mode ramp : Fields::prepare (dt = 10 fs, the Hall step) with an impressed current I0 + k (t - t0), plasma = 0,
//               uniform 300 K conductivity; prints t, I, V_T and the step ledger every `every` steps.
// Usage: hall_port_load dc SIGMA_S_PER_M I0_A [N_TRANSVERSE WIDTH_M]   (default 45 cells over 9 mm: the Hall mesh)
//        hall_port_load ramp SIGMA I0_A DIDT_A_PER_S NSTEPS EVERY
#include "hall-finite-foil-em-20261029/axial_fields.hpp"
#include <cstdio>
#include <cstdlib>
#include <string>
int main(int argc,char**argv){
    if(argc<4){std::fprintf(stderr,"usage\n");return 2;}
    const std::string mode=argv[1];const double sigma=std::atof(argv[2]),I0=std::atof(argv[3]);
    const int nt=(mode=="dc"&&argc>5)?std::atoi(argv[4]):45;const double width=(mode=="dc"&&argc>5)?std::atof(argv[5]):.009;
    axial::Mesh g(nt,4,width,1);axial::Vec s(g.cells.size());for(std::size_t c=0;c<s.size();++c)if(g.cells[c].body>=0)s[c]=sigma;
    axial::Fields f(g,s,I0);
    long double WB=0,WE=0;for(int b=0;b<g.nb;++b)WB+=.5L*g.H[b]*f.initial.B[b]*f.initial.B[b];
    for(int e=0;e<g.ne;++e)WE+=.5L*g.M[e]*f.initial.E[e]*f.initial.E[e];
    std::printf("mesh nx=%d ny=%d nz=%d ne=%d nb=%d dx=%.6e dz=%.6e\n",g.nx,g.ny,g.nz,g.ne,g.nb,g.dx,g.dz);
    std::printf("dc I0=%.17g V_T=%.17g R_dc=%.17g W_B=%.17Lg L_dc=%.17Lg W_E=%.17Lg C_dc=%.17Lg\n",I0,f.initial.terminal_V,
        f.initial.terminal_V/I0,WB,2*WB/(I0*I0),WE,2*WE/(f.initial.terminal_V*f.initial.terminal_V));
    if(mode=="dc")return 0;
    const double k=std::atof(argv[4]);const long n=std::atol(argv[5]),every=std::atol(argv[6]);const double dt=1e-14,t0=7.4978805549199997e-10;
    auto st=f.start(t0);axial::Vec plasma(g.ne);
    long double Wc=0,J=0,Lo=0,Fd=0;
    for(long i=0;i<n;++i){const double tm=st.time_s+.5*dt;const double I=I0+k*(tm-t0);
        auto [next,r]=f.prepare(st,plasma,s,dt,I);Wc+=r.circuit_work_J;J+=r.joule_J;Lo+=r.outer_loss_J;Fd+=r.field_delta_J;st=std::move(next);
        if((i+1)%every==0)std::printf("step %ld t=%.17g I=%.17g V_T=%.17g circuit_J=%.10Lg joule_J=%.10Lg outer_loss_J=%.10Lg field_delta_J=%.10Lg balance=%.3e\n",
            i+1,tm,I,r.terminal_V,Wc,J,Lo,Fd,r.balance_J);
        std::fflush(stdout);}
    return 0;}
