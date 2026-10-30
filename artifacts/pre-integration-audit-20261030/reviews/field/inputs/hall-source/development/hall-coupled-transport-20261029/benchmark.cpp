#include "coupled_pic.hpp"
#include <algorithm>
#include <chrono>
#include <complex>
#include <iomanip>
#include <iostream>
#include <numbers>
using namespace hall;
namespace fs=std::filesystem;
constexpr double pi=std::numbers::pi;
Json rows=Json::array();fs::path output;
void emit(Json j){rows.push_back(j);std::cout<<j.dump()<<'\n';}

void uniform_drift(){
    const double n=1e14,E=1e5,Bmag=10.,M=proton_mass/pic_me;
    const double T=20*2*pi/(pic_e/pic_me*Bmag);
    for(int sign:{-1,1})for(int refine:{1,2,4}){
        int steps=static_cast<int>(std::ceil(T/1e-14))*refine;double dt=T/steps;
        std::array<PicVector,2>u{},x{};std::array<double,2>work{},defect{},activity{};
        const double B=sign*Bmag;
        for(int it=0;it<steps;++it)for(int s=0;s<2;++s){
            const double q=s?pic_e:-pic_e,m=s?proton_mass:pic_me,qm=q/m;
            const auto next=pic_boris(u[s],qm,{E,0,0},{0,0,B},dt);
            auto vm=velocity(u[s]),vp=velocity(next);for(int k=0;k<3;++k)x[s][k]+=.5*dt*(vm[k]+vp[k]);
            PicVector minus=u[s],plus=next;minus[0]+=.5*dt*qm*E;plus[0]-=.5*dt*qm*E;
            double w=q*E*.5*dt*((minus[0]+u[s][0])/(gamma(minus)+gamma(u[s]))+(next[0]+plus[0])/(gamma(next)+gamma(plus)));
            GapParticle p;p.mass_kg=m;p.proper_velocity_m_s=u[s];double k0=pic_kinetic_energy_J(p);p.proper_velocity_m_s=next;
            work[s]+=w;activity[s]+=std::abs(w);defect[s]+=pic_kinetic_energy_J(p)-k0-w;u[s]=next;
        }
        std::array<double,2>jy{},expected{},mean{};
        for(int s=0;s<2;++s){const double q=s?pic_e:-pic_e,m=s?proton_mass:pic_me,a=q*B/m*T;
            mean[s]=x[s][1]/T;jy[s]=q*n*mean[s];expected[s]=q*n*(-E/B)*(1-std::sin(a)/a);}
        double relative=(jy[0]+jy[1]-expected[0]-expected[1])/(n*pic_e*E/Bmag);
        emit({{"family","uniform_drift"},{"B_T",B},{"dt_s",dt},{"steps",steps},{"duration_s",T},{"density_cm3",n/1e6},{"mass_ratio",M},{"mean_vy_m_s",mean},{"Jy_A_m2",jy},{"expected_Jy_A_m2",expected},{"relative_total_current_error",relative},{"work_J_per_particle",work},{"kick_work_defect_J_per_particle",defect},{"work_activity_J_per_particle",activity},{"ion_gyro_angle_over_window",pic_e*Bmag/proton_mass*T}});
    }
    // Common balanced E cross B drift: equal densities give zero net current.
    double defect=0;for(int s=0;s<2;++s){const double v=-E/Bmag,g=1/std::sqrt(1-v*v/(pic_c*pic_c));PicVector u{0,g*v,0};auto v1=velocity(pic_boris(u,(s?pic_e/proton_mass:-pic_e/pic_me),{E,0,0},{0,0,Bmag},1e-14));defect+=std::abs(v1[1]-v);}
    emit({{"family","balanced_common_drift"},{"velocity_error_m_s",defect},{"net_current_expected_A_m2",0}});
}

struct WaveExact {std::array<double,2>xi,vy,vz;};
WaveExact exact_wave(double t,double M,double b,double delta,double wp){
    // Independent analytic matrix-function solution for the linear cold mode.
    const double a=1+b*b,d=1/M+b*b/(M*M),off=-1.,low=-1/M;
    const double det=a*d-off*low;
    std::array<double,2>eq{d*b*b*delta/det,-low*b*b*delta/det};
    const double disc=std::sqrt((a-d)*(a-d)+4*off*low),l1=.5*(a+d+disc),l2=.5*(a+d-disc);
    std::array<double,2>z{delta-eq[0],-eq[1]},Az{a*z[0]+off*z[1],low*z[0]+d*z[1]};
    WaveExact w{eq,{0,0},{0,0}};
    for(int s=0;s<2;++s){double p1=(Az[s]-l2*z[s])/(l1-l2),p2=z[s]-p1;
        w.xi[s]+=p1*std::cos(std::sqrt(l1)*wp*t)+p2*std::cos(std::sqrt(l2)*wp*t);
        w.vy[s]=-wp*(p1*std::sqrt(l1)*std::sin(std::sqrt(l1)*wp*t)+p2*std::sqrt(l2)*std::sin(std::sqrt(l2)*wp*t));
        w.vz[s]=-(s?b*wp/M:-b*wp)*(w.xi[s]-(s?0:delta));}
    return w;
}

// Independent Q1 path audit, constructed from the two physical half-drifts.
struct PathAudit {
    GapPicOptions o;std::array<std::vector<double>,2>dq,h,v;
    explicit PathAudit(GapPicOptions opts):o(opts){for(int s=0;s<2;++s){dq[s].resize((o.nx+1)*(o.ny+1));h[s].resize(o.nx*(o.ny+1));v[s].resize((o.nx+1)*o.ny);}}
    void node(int s,PicVector p,double q){double x=p[0]/o.width_m*o.nx,y=p[1]/o.gap_m*o.ny;int i=std::min(o.nx-1,int(x)),j=std::min(o.ny-1,int(y));x-=i;y-=j;dq[s][j*(o.nx+1)+i]+=q*(1-x)*(1-y);dq[s][j*(o.nx+1)+i+1]+=q*x*(1-y);dq[s][(j+1)*(o.nx+1)+i]+=q*(1-x)*y;dq[s][(j+1)*(o.nx+1)+i+1]+=q*x*y;}
    void segment(int s,PicVector a,PicVector b,double q){
        const double hx=o.width_m/o.nx,hy=o.gap_m/o.ny;std::vector<double>cuts{0,1};
        for(int d=0;d<2;++d){double spacing=d?hy:hx;if(b[d]!=a[d])for(int k=int(std::floor(std::min(a[d],b[d])/spacing))+1;k*spacing<std::max(a[d],b[d]);++k){double t=(k*spacing-a[d])/(b[d]-a[d]);if(t>0&&t<1)cuts.push_back(t);}}
        std::sort(cuts.begin(),cuts.end());
        for(std::size_t k=1;k<cuts.size();++k){double t=.5*(cuts[k]+cuts[k-1]);double x=(a[0]+t*(b[0]-a[0]))/hx,y=(a[1]+t*(b[1]-a[1]))/hy;int i=std::min(o.nx-1,int(x)),j=std::min(o.ny-1,int(y));x-=i;y-=j;
            const double dx=q*(b[0]-a[0])*(cuts[k]-cuts[k-1])/hx,dy=q*(b[1]-a[1])*(cuts[k]-cuts[k-1])/hy;
            h[s][j*o.nx+i]+=dx*(1-y);h[s][(j+1)*o.nx+i]+=dx*y;v[s][j*(o.nx+1)+i]+=dy*(1-x);v[s][j*(o.nx+1)+i+1]+=dy*x;}
    }
    std::pair<double,double>finish(const GapPic&p){double continuity=0,parity=0;
        for(int s=0;s<2;++s){for(int j=0;j<=o.ny;++j)for(int i=0;i<o.nx;++i){double q=h[s][j*o.nx+i];dq[s][j*(o.nx+1)+i]+=q;dq[s][j*(o.nx+1)+i+1]-=q;}
            for(int j=0;j<o.ny;++j)for(int i=0;i<=o.nx;++i){double q=v[s][j*(o.nx+1)+i];dq[s][j*(o.nx+1)+i]+=q;dq[s][(j+1)*(o.nx+1)+i]-=q;}for(double q:dq[s])continuity=std::max(continuity,std::abs(q));}
        for(std::size_t i=0;i<h[0].size();++i)parity=std::max(parity,std::abs(h[0][i]+h[1][i]-p.horizontal_transfers_C()[i]));
        for(std::size_t i=0;i<v[0].size();++i)parity=std::max(parity,std::abs(v[0][i]+v[1][i]-p.vertical_transfers_C()[i]));return {continuity,parity};
    }
};

void native_wave(int ny,int ppc,int steps,double M,bool restart){
    constexpr double wp=1e12,b=1,delta=1e-9,T=4e-12;
    const double n=wp*wp*eps*pic_me/(pic_e*pic_e),k=2*pi/.001;
    GapPicOptions o;o.nx=2;o.ny=ny;o.magnetic_field_T={b*wp*pic_me/pic_e,0,0};
    CoupledPic c(o,M);std::vector<GapParticle>ps;std::vector<int>sp;std::vector<PicVector>refs;
    const double volume=o.width_m*o.gap_m*o.depth_m,weight=n*volume/(o.nx*ny*ppc),dt=T/steps,Q=2*n*pic_e*volume;
    for(int j=0;j<ny*ppc;++j)for(int i=0;i<o.nx;++i)for(int s=0;s<2;++s){
        PicVector ref{(i+.5)*o.width_m/o.nx,(j+.5)*o.gap_m/(ny*ppc),0};GapParticle p;p.position_m=ref;p.position_m[1]+=(s?0:delta*std::sin(k*ref[1]));p.weight=weight;p.charge_C=s?pic_e:-pic_e;p.mass_kg=s?M*pic_me:pic_me;ps.push_back(p);sp.push_back(s);refs.push_back(ref);}
    c.authored("manufactured-initial-neutral-plasma",ps,sp,refs);c.native().solve_fields();
    const double U0=c.native().diagnostics().field_energy_J;
    mfem::FunctionCoefficient phi([&](const mfem::Vector&x){return n*pic_e*delta/(eps*k)*(std::cos(k*x[1])-1);});
    const double phierr=c.native().potential().ComputeL2Error(phi)/(n*pic_e*delta/(eps*k)*std::sqrt(o.width_m*o.gap_m));
    double energy=0,gauss=0,continuity=0,parity=0,charge=0,residual=0,restart_position=0,restart_velocity=0;
    std::unique_ptr<CoupledPic>resumed;Json history=Json::array();
    for(int it=0;it<steps;++it){
        auto before=c.native().particles();c.advance(dt);const auto&after=c.native().particles();require(before.size()==after.size(),"unexpected boundary loss in wave");PathAudit audit(o);
        for(std::size_t j=0;j<before.size();++j){const auto&a=before[j];const auto&z=after[j];const int s=c.tags().at(a.id).species;auto mid=a.position_m,vel=velocity(a.proper_velocity_m_s);for(int d=0;d<3;++d)mid[d]+=.5*dt*vel[d];audit.node(s,a.position_m,-a.charge_C*a.weight);audit.node(s,z.position_m,z.charge_C*z.weight);audit.segment(s,a.position_m,mid,a.charge_C*a.weight);audit.segment(s,mid,z.position_m,a.charge_C*a.weight);}
        auto [cont,par]=audit.finish(c.native());continuity=std::max(continuity,cont/Q);parity=std::max(parity,par/Q);
        const auto d=c.native().diagnostics();energy=std::max(energy,std::abs(d.field_energy_J+d.particle_kinetic_energy_J-U0)/U0);gauss=std::max(gauss,std::abs(d.gauss_defect_C)/Q);charge=std::max(charge,std::abs(d.particle_charge_defect_C)/Q);residual=std::max(residual,d.poisson_residual_C_per_m*o.depth_m/Q);
        if(resumed){resumed->advance(dt);const auto&r=resumed->native().particles();for(std::size_t j=0;j<r.size();++j)for(int d0=0;d0<3;++d0){restart_position=std::max(restart_position,std::abs(r[j].position_m[d0]-after[j].position_m[d0]));restart_velocity=std::max(restart_velocity,std::abs(r[j].proper_velocity_m_s[d0]-after[j].proper_velocity_m_s[d0]));}}
        if(restart&&it==steps/2-1){c.save(output/"wave-checkpoint");resumed=std::make_unique<CoupledPic>(o,M);resumed->load(output/"wave-checkpoint");}
        if(it%(std::max(1,steps/100))==0||it==steps-1){std::array<double,2>xi{},vy{},vz{},ke{};
            for(const auto&p:after){const auto&t=c.tags().at(p.id);const double f=2*p.weight*std::sin(k*t.reference[1])/(n*volume);xi[t.species]+=f*(p.position_m[1]-t.reference[1]);auto v0=velocity(p.proper_velocity_m_s);vy[t.species]+=f*v0[1];vz[t.species]+=f*v0[2];ke[t.species]+=p.weight*pic_kinetic_energy_J(p);}
            auto ex=exact_wave(d.time_s,M,b,delta,wp);history.push_back({{"t_s",d.time_s},{"xi_m",xi},{"vy_m_s",vy},{"vz_m_s",vz},{"exact_xi_m",ex.xi},{"exact_vy_m_s",ex.vy},{"exact_vz_m_s",ex.vz},{"K_J",ke},{"U_J",d.field_energy_J}});}
    }
    auto last=history.back();auto ex=exact_wave(T,M,b,delta,wp);std::array<double,2>err{},jy{},jz{};
    for(int s=0;s<2;++s){double x=last["xi_m"][s];err[s]=std::abs(x-ex.xi[s])/(s?delta/M:delta);jy[s]=(s?1:-1)*pic_e*n*last["vy_m_s"][s].get<double>();jz[s]=(s?1:-1)*pic_e*n*last["vz_m_s"][s].get<double>();}
    std::string name="wave-n"+std::to_string(ny)+"-ppc"+std::to_string(ppc)+"-steps"+std::to_string(steps)+"-M"+std::to_string(int(M));std::ofstream(output/(name+".json"))<<history.dump(2)<<'\n';
    double workdef=c.work_abs_defect[0]+c.work_abs_defect[1];
    std::array<double,2>max_displacement{},max_current{};
    for(const auto&frame:history)for(int s=0;s<2;++s){max_displacement[s]=std::max(max_displacement[s],std::abs(frame["xi_m"][s].get<double>()-frame["exact_xi_m"][s].get<double>())/(s?delta/M:delta));for(const auto&key:{std::pair{"vy_m_s","exact_vy_m_s"},std::pair{"vz_m_s","exact_vz_m_s"}})max_current[s]=std::max(max_current[s],std::abs(frame[key.first][s].get<double>()-frame[key.second][s].get<double>())/(wp*(s?delta/M:delta)));}
    double work_restart=0;if(resumed)for(int s=0;s<2;++s)work_restart=std::max(work_restart,std::abs(c.work[s]-resumed->work[s]));
    emit({{"family","native_wave"},{"name",name},{"ny",ny},{"ppc_y",ppc},{"steps",steps},{"dt_s",dt},{"duration_s",T},{"mass_ratio",M},{"density_cm3",n/1e6},{"B_T",o.magnetic_field_T[0]},{"delta_m",delta},{"initial_phi_relative_L2_error",phierr},{"final_displacement_error_scaled",err},{"max_displacement_error_scaled",max_displacement},{"max_velocity_error_scaled",max_current},{"final_xi_m",last["xi_m"]},{"expected_xi_m",ex.xi},{"final_Jy_A_m2",jy},{"final_Jz_A_m2",jz},{"max_energy_relative_defect",energy},{"max_gauss_relative_defect",gauss},{"max_species_continuity_relative_defect",continuity},{"max_native_path_parity_relative_defect",parity},{"max_charge_relative_defect",charge},{"max_poisson_relative_residual",residual},{"work_J",c.work},{"kick_work_relative_defect",workdef/U0},{"restart_tested",restart},{"restart_position_error_m",restart_position},{"restart_velocity_error_m_s",restart_velocity},{"restart_work_error_J",work_restart},{"U0_J",U0}});
}

void transactions(){
    GapPicOptions o;o.nx=3;o.ny=4;CoupledPic p(o,25);p.native().solve_fields();GapParticle e;e.position_m={.0004,.0005,0};e.weight=1000;auto ion=e;ion.charge_C=pic_e;ion.mass_kg=25*pic_me;
    bool reject=false;auto invalid=ion;invalid.position_m[0]=-1;
    try{p.authored("invalid-second-product",{e,invalid},{0,1},{e.position_m,ion.position_m});}catch(const std::exception&){reject=true;}
    require(reject&&p.native().particles().empty(),"batch rollback");auto E=p.native().electric_field(.0004,.0005);require(dot(E,E)==0,"rejected batch retained solved fields");
    p.authored("manufactured-pair",{e,ion},{0,1},{e.position_m,ion.position_m});
    bool stale=false;try{static_cast<void>(p.native().electric_field(.0004,.0005));}catch(const std::logic_error&){stale=true;}
    p.native().solve_fields();const auto d=p.native().diagnostics();require(d.field_energy_J==0&&d.particle_charge_C==0,"collocated pair net field");
    bool duplicate=false;try{p.authored("manufactured-pair",{e,ion},{0,1},{e.position_m,ion.position_m});}catch(const std::exception&){duplicate=true;}
    GapPic partial(o);partial.add_particle(e);partial.solve_fields();const double partial_energy=partial.diagnostics().field_energy_J;
    p.save(output/"pair-checkpoint");CoupledPic r(o,25);r.load(output/"pair-checkpoint");
    bool species_reject=false;try{CoupledPic wrong(o,26);wrong.load(output/"pair-checkpoint");}catch(const std::exception&){species_reject=true;}
    bool guard=false;try{p.advance(1.1e-14);}catch(const std::exception&){guard=true;}
    bool corrupt_rejected=false;
    fs::copy(output/"pair-checkpoint",output/"corrupt-checkpoint",fs::copy_options::recursive);
    {std::ofstream f(output/"corrupt-checkpoint/native.bin",std::ios::app|std::ios::binary);f.put('x');}
    try{CoupledPic corrupt(o,25);corrupt.load(output/"corrupt-checkpoint");}catch(const std::exception&){corrupt_rejected=true;}
    GapPic walls(o);auto we=e,wi=ion;we.position_m[1]=1e-9;wi.position_m[1]=o.gap_m-1e-9;we.proper_velocity_m_s={0,-1e6,0};wi.proper_velocity_m_s={0,1e6,0};
    const auto wall_ids=walls.add_authored_batch({we,wi});walls.advance(1e-14);const auto wd=walls.diagnostics();
    const double initialK=we.weight*pic_kinetic_energy_J(we)+wi.weight*pic_kinetic_energy_J(wi);
    emit({{"family","signed_walls"},{"electron_lower_charge_C",wd.lower_absorbed_charge_C},{"ion_upper_charge_C",wd.upper_absorbed_charge_C},{"charge_defect_C",wd.particle_charge_defect_C},{"relative_kinetic_throughput_defect",std::abs(wd.lower_absorbed_energy_J+wd.upper_absorbed_energy_J-initialK)/initialK},{"resident_count",walls.particles().size()}});
    emit({{"family","transactions"},{"corrupt_checkpoint_rejected",corrupt_rejected},{"single_product_spurious_field_J",partial_energy},{"invalid_second_product_rollback",reject},{"field_cache_invalidated",stale},{"collocated_pair_field_J",d.field_energy_J},{"duplicate_event_rejected",duplicate},{"wrong_species_restart_rejected",species_reject},{"step_cap_rejected",guard},{"first_ID",r.native().particles()[0].id},{"second_ID",r.native().particles()[1].id},{"charge_C",d.particle_charge_C},{"source_owner","authored reservoir; no reaction or water rate"}});
}

using C=std::complex<double>;using State=std::array<C,4>;
double dispersion(double w,double k,double b,double M){return w*w-k*k-w*(1/(w-b)+1/(M*(w+b/M)));}
double root(double k,double b,double M){double lo=0,hi=b*(1-1e-12);for(int it=0;it<100;++it){double mid=.5*(lo+hi);if(dispersion(mid,k,b,M)>0)hi=mid;else lo=mid;}return .5*(lo+hi);}
std::array<std::array<C,4>,4> maxwell(double kh,double b,double M,int h,bool current,bool faraday){
    std::array<std::array<C,4>,4>a{};a[0][1]=h*kh;if(current){a[0][2]=1;a[0][3]=-1;}if(faraday)a[1][0]=-h*kh;
    a[2][0]=-1;a[2][2]=C(0,-h*b);a[3][0]=1/M;a[3][3]=C(0,h*b/M);return a;
}
mfem::DenseMatrix midpoint(const std::array<std::array<C,4>,4>&a,double dt){
    mfem::DenseMatrix L(8),R(8),S(8);L=0.;R=0.;
    for(int i=0;i<4;++i)for(int j=0;j<4;++j){double r=a[i][j].real(),im=a[i][j].imag();L(i,j)=-dt*r/2;L(i,j+4)=dt*im/2;L(i+4,j)=-dt*im/2;L(i+4,j+4)=-dt*r/2;R(i,j)=dt*r/2;R(i,j+4)=-dt*im/2;R(i+4,j)=dt*im/2;R(i+4,j+4)=dt*r/2;}
    for(int i=0;i<8;++i){L(i,i)+=1;R(i,i)+=1;}mfem::DenseMatrixInverse inverse(L);mfem::Vector column(8),answer(8);
    for(int j=0;j<8;++j){for(int i=0;i<8;++i)column[i]=R(i,j);inverse.Mult(column,answer);for(int i=0;i<8;++i)S(i,j)=answer[i];}return S;
}
double energy(const State&s,double M){return .5*(std::norm(s[0])+std::norm(s[1])+std::norm(s[2])+M*std::norm(s[3]));}
State step(const mfem::DenseMatrix&S,const State&x){mfem::Vector a(8),b(8);for(int i=0;i<4;++i){a[i]=x[i].real();a[i+4]=x[i].imag();}S.Mult(a,b);State y;for(int i=0;i<4;++i)y[i]={b[i],b[i+4]};return y;}
void em_wave(int N,double dt,int polarity=1,bool current=true,bool faraday=true){
    const double k=.2,b=.5,M=proton_mass/pic_me,kh=k*std::sin(pi/N)/(pi/N),w=root(kh,b,M),wc=root(k,b,M);
    const int h=polarity;const double B=polarity*b,amplitude=1e-6;
    State x{amplitude,C(0,-h*kh/w)*amplitude,C(0,-1/(w-h*B))*amplitude,C(0,1/(M*(w+h*B/M)))*amplitude};
    const auto initial=x;const double U=energy(x,M);auto S=midpoint(maxwell(kh,B,M,h,current,faraday),dt);
    const int steps=static_cast<int>(std::ceil(2*2*pi/(w*dt)));double phase=0,error=0,work=0,workdef=0;Json history=Json::array();
    for(int i=0;i<steps;++i){auto y=step(S,x);phase-=std::arg(y[0]/x[0]);error=std::max(error,std::abs(energy(y,M)-U)/U);
        const C E=.5*(x[0]+y[0]),J=.5*(x[3]+y[3]-x[2]-y[2]);const double W=dt*std::real(std::conj(E)*J);work+=W;
        const double dk=.5*(std::norm(y[2])-std::norm(x[2])+M*(std::norm(y[3])-std::norm(x[3])));workdef=std::max(workdef,std::abs(dk-W)/U);
        x=y;if(i%std::max(1,steps/100)==0||i==steps-1)history.push_back({{"t_omega_pe",(i+1)*dt},{"E_real",x[0].real()},{"E_imag",x[0].imag()},{"B_real",x[1].real()},{"B_imag",x[1].imag()},{"total_energy_normalized",energy(x,M)/U}});
    }
    const double measured=phase/(steps*dt),numerical=2/dt*std::atan(w*dt/2);
    const double hall=.5*(b*k*k+std::sqrt(b*b*k*k*k*k+4*b*b*k*k/M));
    const auto reference=initial[0]*std::exp(C(0,-w*steps*dt));double mismatch=std::abs(x[0]-reference)/amplitude;
    std::string name="em-N"+std::to_string(N)+"-dt"+std::to_string(dt)+"-pol"+std::to_string(polarity)+"-J"+std::to_string(current)+"-F"+std::to_string(faraday);std::ofstream(output/(name+".json"))<<history.dump(2)<<'\n';
    emit({{"family","em_wave"},{"name",name},{"N",N},{"dt_omega_pe",dt},{"dt_s",dt/1e12},{"steps",steps},{"B0_T",B*1e12*pic_me/pic_e},{"density_cm3",1e24*eps*pic_me/(pic_e*pic_e)/1e6},{"current_feedback",current},{"faraday_feedback",faraday},{"helicity",h},{"omega_measured_over_pe",measured},{"omega_continuum_over_pe",wc},{"omega_spatial_over_pe",w},{"omega_midpoint_over_pe",numerical},{"omega_Hall_MHD_over_pe",hall},{"dispersion_relative_error",std::abs(measured-wc)/wc},{"midpoint_relative_error",std::abs(measured-numerical)/numerical},{"temporal_relative_error",std::abs(measured-w)/w},{"max_energy_relative_defect",error},{"max_midpoint_work_relative_defect",workdef},{"endpoint_E_relative_mismatch",mismatch},{"electron_ion_relative_current_amplitude_normalized",std::abs(initial[3]-initial[2])},{"linear_uniform_density_continuity_defect",0},{"delta_B_over_B0",std::abs(initial[1])/b}});
}
int main(int argc,char**argv){try{
    require(argc==2,"supply fresh output directory");output=argv[1];fs::create_directories(output);
    const auto start=std::chrono::steady_clock::now();transactions();uniform_drift();
    for(int n:{16,32,64})native_wave(n,4,400,25,false);
    for(int steps:{800,1600})native_wave(64,4,steps,25,steps==800);
    for(int ppc:{2,8})native_wave(32,ppc,400,25,false);
    native_wave(32,4,400,proton_mass/pic_me,false);
    for(int N:{16,32,64,128})em_wave(N,.01);
    for(double dt:{.005,.0025})em_wave(64,dt);
    em_wave(64,.01,-1);em_wave(64,.01,1,false,true);em_wave(64,.01,1,true,false);
    Json limits=Json::array();for(double k:{.025,.05,.1,.2}){double b=.5,M=proton_mass/pic_me,w=root(k,b,M),h=.5*(b*k*k+std::sqrt(b*b*k*k*k*k+4*b*b*k*k/M));limits.push_back({{"k_de",k},{"omega_over_pe",w},{"Hall_omega_over_pe",h},{"relative_difference",std::abs(w-h)/w}});}emit({{"family","hall_limit"},{"rows",limits}});
    emit({{"family","execution"},{"elapsed_s",std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()},{"CPU_only",true}});
    std::ofstream(output/"results.json")<<rows.dump(2)<<'\n';return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
