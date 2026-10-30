#pragma once
#include "axial_mesh.hpp"

namespace axial {
inline Vec scalar_dirichlet(const Mesh&g,const Vec&weights,const Vec&rhs,
                            const std::vector<bool>&fixed) {
    Vec diagonal(g.nn),b=rhs;
    for(int e=0;e<g.ne;++e){double q=weights[e]/std::pow(g.edges[e].length,2);
        diagonal[g.edges[e].a]+=q;diagonal[g.edges[e].b]+=q;}
    for(int n=0;n<g.nn;++n)if(fixed[n]){diagonal[n]=1;b[n]=0;}
    auto apply=[&](Vec v){for(int n=0;n<g.nn;++n)if(fixed[n])v[n]=0;
        auto d=g.gradient(v);for(int e=0;e<g.ne;++e)d[e]*=weights[e];auto a=g.gradient_transpose(d);
        // Fixed unknowns are identically zero in this solve.
        for(int n=0;n<g.nn;++n)if(fixed[n])a[n]=0;return a;};
    return foil::solve(apply,b,diagonal,4000);
}
struct Initial {
    Vec E,B,potential,feed,background_force,solid_current,short_current;
    double current_A,terminal_V,dc_continuity_A=0,magnetic_residual_A_m=0;
};
inline Initial dc_initial(const Mesh&g,const Vec&sigma_cells,double current) {
    hall::require(std::isfinite(current),"invalid supported DC start");
    Initial out;out.current_A=current;out.feed=g.feed_shape();auto S=g.conduction(sigma_cells);
    Vec impressed=out.feed;for(double&v:impressed)v*=current;
    auto rhs=g.gradient_transpose(impressed);std::vector<bool>fixed(g.nn);
    for(int n=0;n<g.nn;++n)fixed[n]=g.node_body[n]<0||g.short_node[n];
    // The declared uniform-temperature DC start has an exact affine metal
    // potential. Evaluate it directly to avoid roundoff-sized transverse modes
    // from a highly anisotropic iterative conductor solve. Admit it with the
    // SAME independently recomputed 2e-13 residual bound as the iterative path.
    double uniform_sigma=0;bool uniform=true;
    for(std::size_t c=0;c<g.cells.size();++c)if(g.cells[c].body>=0){
        hall::require(sigma_cells[c]>0,"initial DC metal must conduct");
        if(uniform_sigma==0)uniform_sigma=sigma_cells[c];else if(sigma_cells[c]!=uniform_sigma)uniform=false;}
    Vec phi(g.nn);
    if(uniform){
        double gradient=current/(.005*.0002*uniform_sigma);
        for(int k=g.kshort;k<=g.kfeed;++k)for(int j=0;j<=g.ny;++j)for(int i=0;i<=g.nx;++i){
            int n=g.node(i,j,k);if(g.node_body[n]>=0)phi[n]=(g.node_body[n]? -1:1)*gradient*(k-g.kshort)*g.dz;}
        auto applied=g.gradient(phi);for(int e=0;e<g.ne;++e)applied[e]*=S[e];auto residual=g.gradient_transpose(applied);Vec b=rhs;
        for(int n=0;n<g.nn;++n){residual[n]=fixed[n]?0:residual[n]-rhs[n];if(fixed[n])b[n]=0;}
        hall::require(norm(residual)<=std::max(1e-30,2e-13*norm(b)),"analytic DC conductor independent residual failed");
    }else try{phi=scalar_dirichlet(g,S,rhs,fixed);}catch(const std::exception&e){throw std::runtime_error(std::string("DC conductor scalar: ")+e.what());}
    // The same scalar potential extends through vacuum. All three derivatives
    // are retained, so transverse/axial scalar curl terms cancel identically.
    for(int n=0;n<g.nn;++n)fixed[n]=g.node_body[n]>=0||g.short_node[n];
    auto known=g.gradient(phi);for(int e=0;e<g.ne;++e)known[e]*=g.VE[e];
    rhs=g.gradient_transpose(known);for(double&v:rhs)v=-v;
    Vec vacuum;try{vacuum=scalar_dirichlet(g,g.VE,rhs,fixed);}catch(const std::exception&e){throw std::runtime_error(std::string("DC vacuum scalar: ")+e.what());}
    for(int n=0;n<g.nn;++n)phi[n]+=vacuum[n];out.potential=phi;
    out.E=g.gradient(phi);for(double&v:out.E)v=-v;
    out.terminal_V=-dot(out.E,out.feed);
    out.solid_current.resize(g.ne);Vec closed(g.ne);
    for(int e=0;e<g.ne;++e)closed[e]=(out.solid_current[e]=S[e]*out.E[e])+impressed[e];
    auto residual=g.gradient_transpose(closed);out.short_current.resize(g.ne);
    if(uniform){
        // Explicit uniform-DC return scenario: the short sheet carries the
        // reflected feed-contact quadrature in the opposite direction. An ideal
        // PEC sheet does not specify a resistive current profile by itself.
        for(int j=0;j<g.ny;++j)for(int i=0;i<=g.nx;++i)
            out.short_current[g.ey(i,j,g.kshort)]=-impressed[g.ey(i,j,g.kfeed)];
        auto div=g.gradient_transpose(out.short_current);Vec b(g.nn),error(g.nn);
        for(int n=0;n<g.nn;++n)if(g.short_node[n]){b[n]=residual[n];error[n]=residual[n]+div[n];}
        hall::require(norm(error)<=std::max(1e-30,2e-13*norm(b)),"analytic DC short independent residual failed");
    }else{
        Vec short_weight(g.ne);for(int e=0;e<g.ne;++e)if(g.fixed_E[e])short_weight[e]=g.VE[e];
        bool anchored=false;for(int n=0;n<g.nn;++n){fixed[n]=!g.short_node[n];
            if(!anchored&&g.short_node[n]){fixed[n]=true;anchored=true;}rhs[n]=-residual[n];}
        Vec short_phi;try{short_phi=scalar_dirichlet(g,short_weight,rhs,fixed);}catch(const std::exception&e){throw std::runtime_error(std::string("DC short scalar: ")+e.what());}
        out.short_current=g.gradient(short_phi);for(int e=0;e<g.ne;++e)out.short_current[e]*=short_weight[e];
    }
    for(int e=0;e<g.ne;++e)closed[e]+=out.short_current[e];
    out.dc_continuity_A=maximum(g.gradient_transpose(closed));
    hall::require(out.dc_continuity_A<=2e-10*std::max(1.,std::abs(current)),"initial feed/solid/short KCL");

    // Magnetostatic vector potential: tangential A=0 on the outer air box,
    // Interior Coulomb penalty G^T V_E A with scalar gauge test functions zero
    // on the outer boundary, compatible with essential tangential A. Penalizing
    // outer nodal divergence as well would add a spurious harmonic gauge force.
    // This gauge is initialization metadata only;
    // the evolved variables are total E and B and never a sampled dA/dt.
    Vec diagonal(g.ne);
    for(int b=0;b<g.nb;++b)for(int q=0;q<4;++q){int e=g.C[b].index[q];
        diagonal[e]+=g.H[b]*std::pow(g.C[b].coefficient[q],2);}
    for(int e=0;e<g.ne;++e){const auto&v=g.edges[e];
        diagonal[e]+=g.VE[e]*g.VE[e]/(mu*v.length*v.length)*
            ((g.outer_node[v.a]?0:1/g.VN[v.a])+(g.outer_node[v.b]?0:1/g.VN[v.b]));
        if(g.fixed_A[e]){diagonal[e]=1;closed[e]=0;}}
    auto apply=[&](Vec A){for(int e=0;e<g.ne;++e)if(g.fixed_A[e])A[e]=0;
        auto K=g.curl_adjoint(g.curl(A));for(int e=0;e<g.ne;++e)A[e]*=g.VE[e];
        auto div=g.gradient_transpose(A);for(int n=0;n<g.nn;++n)div[n]=g.outer_node[n]?0:div[n]/(mu*g.VN[n]);auto grad=g.gradient(div);
        for(int e=0;e<g.ne;++e)K[e]=g.fixed_A[e]?0:K[e]+g.VE[e]*grad[e];return K;};
    Vec A;try{A=foil::solve(apply,closed,diagonal,5000);}catch(const std::exception&e){throw std::runtime_error(std::string("DC Coulomb vector: ")+e.what());}out.B=g.curl(A);
    auto ampere=g.curl_adjoint(out.B);out.background_force.resize(g.ne);
    for(int e=0;e<g.ne;++e){
        double defect=out.solid_current[e]+impressed[e]+out.short_current[e]-ampere[e];
        if(!g.fixed_A[e])out.magnetic_residual_A_m=std::max(out.magnetic_residual_A_m,std::abs(defect));
        // Only the terminated outer boundary maintains the initial DC field.
        if(g.fixed_A[e])out.background_force[e]=g.Y[e]*out.E[e]+defect;
    }
    hall::require(out.magnetic_residual_A_m<=2e-9*std::max(1.,maximum(closed)),"initial interior Ampere residual");
    return out;
}
struct State {
    Vec E,B;double time_s=0;std::uint64_t revision=0;
    // Canonical high+low field values retain changes smaller than an ULP of the
    // initial tesla-scale DC field. These are stored state, not energy repairs.
    Vec E_low,B_low;
};
inline std::pair<double,double> two_sum(double a,double b){double s=a+b,v=s-a;return {s,(a-(s-v))+(b-v)};}
inline void add_field(double old,double low,double increment,double&next,double&next_low){
    auto first=two_sum(old,increment);auto last=two_sum(first.first,low+first.second);next=last.first;next_low=last.second;
}
struct Receipt {
    double field_delta_J=0,plasma_work_J=0,joule_J=0,circuit_work_J=0,outer_loss_J=0,outer_supply_J=0;
    double balance_J=0,terminal_V=0,current_A=0,ampere_residual_A_m=0,kcl_residual_A=0,gauss_increment_C=0;
    double region_conduction_A=0,region_plasma_A=0,region_displacement_A=0,region_boundary_A=0;
    Vec cell_joule_W,solid_current,short_current,mid_E;
    Json json()const{return {{"field_delta_J",field_delta_J},{"plasma_work_J",plasma_work_J},{"joule_J",joule_J},
        {"circuit_work_J",circuit_work_J},{"outer_loss_J",outer_loss_J},{"outer_supply_J",outer_supply_J},
        {"balance_J",balance_J},{"terminal_V",terminal_V},{"current_A",current_A},
        {"ampere_residual_A_m",ampere_residual_A_m},{"kcl_residual_A",kcl_residual_A},{"gauss_increment_C",gauss_increment_C},
        {"region_conduction_A",region_conduction_A},{"region_plasma_A",region_plasma_A},
        {"region_displacement_A",region_displacement_A},{"region_boundary_A",region_boundary_A}};}
};
class Fields {
    const Mesh&g;
public:
    Initial initial;
    Fields(const Mesh&mesh,const Vec&sigma_cells,double initial_current):g(mesh),initial(dc_initial(mesh,sigma_cells,initial_current)){}
    State start(double t)const{return {initial.E,initial.B,t,0,Vec(g.ne),Vec(g.nb)};}
    std::pair<State,Receipt> prepare(const State&old,const Vec&plasma,const Vec&sigma_cells,double dt,double current)const{
        hall::require(dt>0&&dt<=1e-14&&std::isfinite(dt)&&std::isfinite(current)&&
                      int(plasma.size())==g.ne&&finite(plasma),"invalid axial common step");
        hall::require(int(old.E.size())==g.ne&&int(old.B.size())==g.nb&&int(old.E_low.size())==g.ne&&int(old.B_low.size())==g.nb&&
            finite(old.E)&&finite(old.B)&&finite(old.E_low)&&finite(old.B_low),"invalid axial field state");
        auto S=g.conduction(sigma_cells);Vec diagonal(g.ne,1),sqrtM(g.ne);
        for(int e=0;e<g.ne;++e){sqrtM[e]=std::sqrt(g.M[e]);if(!g.fixed_E[e])diagonal[e]+=.5*dt*(S[e]+g.Y[e])/g.M[e];}
        for(int b=0;b<g.nb;++b)for(int q=0;q<4;++q){int e=g.C[b].index[q];
            if(!g.fixed_E[e])diagonal[e]+=.25*dt*dt*g.H[b]*std::pow(g.C[b].coefficient[q],2)/g.M[e];}
        auto apply=[&](const Vec&v){Vec E(g.ne);for(int e=0;e<g.ne;++e)if(!g.fixed_E[e])E[e]=v[e]/sqrtM[e];
            auto K=g.curl_adjoint(g.curl(E));Vec result(g.ne);
            for(int e=0;e<g.ne;++e)result[e]=g.fixed_E[e]?v[e]:v[e]+(.5*dt*(S[e]+g.Y[e])*E[e]+.25*dt*dt*K[e])/sqrtM[e];return result;};
        auto curl=g.curl(old.E);Vec mid_B=old.B;
        for(int b=0;b<g.nb;++b)mid_B[b]-=.5*dt*curl[b];auto ampere=g.curl_adjoint(mid_B);Vec rhs(g.ne);
        for(int e=0;e<g.ne;++e)if(!g.fixed_E[e])rhs[e]=dt*(ampere[e]-plasma[e]-current*initial.feed[e]
            -(S[e]+g.Y[e])*old.E[e]+initial.background_force[e])/sqrtM[e];
        auto increment=foil::solve(apply,rhs,diagonal,200);State next=old;Receipt r;
        r.current_A=current;r.mid_E.resize(g.ne);r.solid_current.resize(g.ne);r.short_current.resize(g.ne);
        for(int e=0;e<g.ne;++e){add_field(old.E[e],old.E_low[e],increment[e]/sqrtM[e],next.E[e],next.E_low[e]);
            r.mid_E[e]=old.E[e]+(old.E_low[e]+.5*increment[e]/sqrtM[e]);
            hall::require(!g.fixed_E[e]||next.E[e]==0,"ideal short field changed");}
        auto cmid=g.curl(r.mid_E);for(int b=0;b<g.nb;++b)add_field(old.B[b],old.B_low[b],-dt*cmid[b],next.B[b],next.B_low[b]);
        for(int b=0;b<g.nb;++b)mid_B[b]=.5*(old.B[b]+next.B[b]);ampere=g.curl_adjoint(mid_B);
        Vec displacement(g.ne),boundary(g.ne),total(g.ne);
        for(int e=0;e<g.ne;++e){r.solid_current[e]=S[e]*r.mid_E[e];
            const long double delta=(long double)next.E[e]-old.E[e]+next.E_low[e]-old.E_low[e];
            const long double midpoint=(long double)old.E[e]+old.E_low[e]+.5L*delta;
            displacement[e]=g.M[e]*double(delta)/dt;
            boundary[e]=g.Y[e]*r.mid_E[e]-initial.background_force[e];
            total[e]=r.solid_current[e]+plasma[e]+current*initial.feed[e]+boundary[e]+displacement[e];
            double residue=ampere[e]-total[e];
            if(g.fixed_E[e]){r.short_current[e]=residue;total[e]+=residue;}
            else r.ampere_residual_A_m=std::max(r.ampere_residual_A_m,std::abs(residue));
            r.field_delta_J+=double(g.M[e]*midpoint*delta);
            r.plasma_work_J+=dt*r.mid_E[e]*plasma[e];r.joule_J+=dt*r.mid_E[e]*r.solid_current[e];
            r.circuit_work_J-=dt*current*initial.feed[e]*r.mid_E[e];
            r.outer_loss_J+=dt*g.Y[e]*r.mid_E[e]*r.mid_E[e];r.outer_supply_J+=dt*initial.background_force[e]*r.mid_E[e];
        }
        for(int b=0;b<g.nb;++b){const long double delta=(long double)next.B[b]-old.B[b]+next.B_low[b]-old.B_low[b];
            const long double midpoint=(long double)old.B[b]+old.B_low[b]+.5L*delta;r.field_delta_J+=double(g.H[b]*midpoint*delta);}
        r.terminal_V=-dot(initial.feed,r.mid_E);
        r.cell_joule_W=g.cell_joule_W(sigma_cells,r.mid_E);
        r.balance_J=r.field_delta_J+r.plasma_work_J+r.joule_J+r.outer_loss_J-r.outer_supply_J-r.circuit_work_J;
        r.kcl_residual_A=maximum(g.gradient_transpose(total));
        auto charge0=g.field_charge(old.E),charge1=g.field_charge(next.E),low0=g.field_charge(old.E_low),low1=g.field_charge(next.E_low);Vec physical_current=total;
        for(int e=0;e<g.ne;++e)physical_current[e]-=displacement[e];
        auto continuity=g.gradient_transpose(physical_current);double charge_activity=0;
        for(int n=0;n<g.nn;++n){double change=charge1[n]-charge0[n]+low1[n]-low0[n];charge_activity=std::max(charge_activity,std::abs(change));
            r.gauss_increment_C=std::max(r.gauss_increment_C,std::abs(change-dt*continuity[n]));}
        hall::require(r.gauss_increment_C<=1e-22+2e-10*charge_activity,"axial full-domain Gauss increment");
        int cut=g.kshort+g.nz_foil/2-1;
        r.region_conduction_A=g.upper_region_inflow(r.solid_current,cut);
        r.region_plasma_A=g.upper_region_inflow(plasma,cut);
        r.region_displacement_A=g.upper_region_inflow(displacement,cut);
        r.region_boundary_A=g.upper_region_inflow(boundary,cut);
        double activity=std::abs(r.field_delta_J)+std::abs(r.plasma_work_J)+std::abs(r.joule_J)+
            std::abs(r.outer_loss_J)+std::abs(r.outer_supply_J)+std::abs(r.circuit_work_J);
        if(!(std::abs(r.balance_J)<=1e-20+1e-10*activity)){
            std::ostringstream message;message<<std::setprecision(17)<<"axial field/current/heat work balance: dt="<<dt
                <<" revision="<<old.revision<<" defect_J="<<r.balance_J<<" allowed_J="<<(1e-20+1e-10*activity)
                <<" field_delta_J="<<r.field_delta_J<<" Joule_J="<<r.joule_J<<" circuit_J="<<r.circuit_work_J;
            throw std::runtime_error(message.str());}
        next.time_s+=dt;++next.revision;return {std::move(next),std::move(r)};
    }
    double energy(const State&s)const{long double value=0;for(int e=0;e<g.ne;++e){long double v=(long double)s.E[e]+s.E_low[e];value+=.5L*g.M[e]*v*v;}
        for(int b=0;b<g.nb;++b){long double v=(long double)s.B[b]+s.B_low[b];value+=.5L*g.H[b]*v*v;}return double(value);}
};
} // namespace axial
