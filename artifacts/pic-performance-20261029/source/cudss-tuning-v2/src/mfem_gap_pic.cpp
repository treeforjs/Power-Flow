#include "raikou/mfem_gap_pic.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <map>
#include <memory>
#include <fstream>
#include <filesystem>
#include <cstring>
#include <stdexcept>
#include <cstdlib>
#include "bench_profile.hpp"
#include "cudss_poisson.hpp"

#ifdef MFEM_USE_MPI
namespace raikou_pic_parallel {
using namespace mfem;
#include "mfem_pic_distributed.hpp"
}
#endif
namespace raikou::mfem_backend {
namespace {
constexpr double epsilon0 = 8.8541878188e-12;
constexpr double kb = 1.380649e-23;
double norm2(const PicVector& a) { return a[0]*a[0]+a[1]*a[1]+a[2]*a[2]; }
double gamma(const PicVector& u) { return std::sqrt(1+norm2(u)/(pic_c*pic_c)); }
bool finite(const PicVector& a) {
    return std::all_of(a.begin(), a.end(), [](double v) { return std::isfinite(v); });
}
PicVector cross(const PicVector& a, const PicVector& b) {
    return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
}
GapPicOptions checked(GapPicOptions o) {
    if (o.nx < 1 || o.ny < 1 || o.nx > 4096 || o.ny > 4096 ||
        !std::isfinite(o.width_m) || o.width_m <= 0 ||
        !std::isfinite(o.gap_m) || o.gap_m <= 0 ||
        !std::isfinite(o.depth_m) || o.depth_m <= 0 ||
        !std::isfinite(o.lower_potential_V) || !std::isfinite(o.upper_potential_V) ||
        !finite(o.magnetic_field_T))
        throw std::invalid_argument("Gap PIC requires finite geometry/fields and positive mesh sizes");
    if(o.finite_foils) {
        for(const auto& f:{o.lower_foil,o.upper_foil}) {
            if(!(0<f.left_m && f.left_m<f.right_m && f.right_m<o.width_m &&
                 0<f.bottom_m && f.bottom_m<f.top_m && f.top_m<o.gap_m))
                throw std::invalid_argument("Finite foils must lie strictly inside the vacuum box");
            for(double x:{f.left_m,f.right_m})
                if(std::abs(x/o.width_m*o.nx-std::round(x/o.width_m*o.nx))>1e-9)
                    throw std::invalid_argument("Foil edges must align to the mesh");
            for(double y:{f.bottom_m,f.top_m})
                if(std::abs(y/o.gap_m*o.ny-std::round(y/o.gap_m*o.ny))>1e-9)
                    throw std::invalid_argument("Foil surfaces must align to the mesh");
        }
        if(o.lower_foil.top_m>=o.upper_foil.bottom_m)
            throw std::invalid_argument("Finite electrodes need a positive intervening gap");
    }
    if(o.finite_axial_extent && (!o.finite_foils || !std::isfinite(o.short_z_m) || !std::isfinite(o.feed_z_m) ||
       o.feed_z_m<=o.short_z_m || std::abs(o.feed_z_m-o.short_z_m-o.depth_m)>1e-12*o.depth_m))
        throw std::invalid_argument("Finite axial bounds must span the extrusion depth");
    return o;
}
void check_particle(const GapParticle& p) {
    if (!finite(p.position_m) || !finite(p.proper_velocity_m_s) ||
        !std::isfinite(norm2(p.proper_velocity_m_s)) ||
        !std::isfinite(p.charge_C) || !std::isfinite(p.mass_kg) || p.mass_kg <= 0 ||
        !std::isfinite(p.weight) || p.weight <= 0 ||
        !std::isfinite(p.charge_C*p.weight) ||
        !std::isfinite(p.charge_C/p.mass_kg) ||
        !std::isfinite(pic_kinetic_energy_J(p)*p.weight))
        throw std::invalid_argument("Gap PIC particle is nonfinite or has invalid mass/weight");
}
}

void pic_parallel_workers() {
#ifdef MFEM_USE_MPI
    raikou_pic_parallel::distributed_linear::Workers();
#else
    throw std::runtime_error("MPI fields unavailable");
#endif
}
void pic_parallel_finish() {
#ifdef MFEM_USE_MPI
    raikou_pic_parallel::distributed_linear::Send(raikou_pic_parallel::distributed_linear::Finish,0);
#endif
}
void pic_parallel_record(const std::string& directory) {
#ifdef MFEM_USE_MPI
    raikou_pic_parallel::distributed_linear::output=directory;
    raikou_pic_parallel::distributed_linear::Record();
#else
    (void)directory;
#endif
}

struct PicPoissonState {
    mfem::ConstantCoefficient epsilon;
    mfem::BilinearForm form;
    std::unique_ptr<mfem::SparseMatrix> full_matrix;
    mfem::GridFunction xfield,yfield;
#ifdef MFEM_USE_MPI
    std::unique_ptr<raikou_pic_parallel::LinearSystem> distributed;
    const mfem::SparseMatrix* address=nullptr;
#endif
    explicit PicPoissonState(mfem::FiniteElementSpace& space):epsilon(epsilon0),form(&space),xfield(&space),yfield(&space) {
        form.AddDomainIntegrator(new mfem::DiffusionIntegrator(epsilon));
        form.Assemble();form.Finalize();full_matrix=std::make_unique<mfem::SparseMatrix>(form.SpMat());
        mfem::FunctionCoefficient xcoord([](const mfem::Vector& x){return x[0];});
        mfem::FunctionCoefficient ycoord([](const mfem::Vector& x){return x[1];});
        xfield.ProjectCoefficient(xcoord);yfield.ProjectCoefficient(ycoord);
    }
};
GapPic::~GapPic()=default;
void GapPic::release_field_solver(){poisson_.reset();}

double pic_kinetic_energy_J(const GapParticle& p) {
    return p.mass_kg*norm2(p.proper_velocity_m_s)/(gamma(p.proper_velocity_m_s)+1);
}

double pic_richardson_supply_A_m2(double T,double phi,double A) {
    if (!std::isfinite(T) || T<=0 || !std::isfinite(phi) || phi<=0 ||
        !std::isfinite(A) || A<=0)
        throw std::invalid_argument("Richardson supply requires explicit positive finite T, work function and prefactor");
    const double J=std::exp(std::log(A)+2*std::log(T)-phi*pic_e/(kb*T));
    if (!std::isfinite(J)) throw std::runtime_error("Richardson supply overflow");
    return J;
}

PicVector pic_boris(const PicVector& u, double qm, const PicVector& E,
                    const PicVector& B, double dt) {
    if (!finite(u) || !finite(E) || !finite(B) || !std::isfinite(qm) ||
        !std::isfinite(dt) || dt <= 0)
        throw std::invalid_argument("Boris update requires finite data and positive dt");
    PicVector minus{}, t{}, prime{}, plus{}, s{}, result{};
    for (int i=0; i<3; ++i) minus[i] = u[i] + .5*dt*qm*E[i];
    for (int i=0; i<3; ++i) t[i] = .5*dt*qm*B[i]/gamma(minus);
    const auto mt = cross(minus,t);
    for (int i=0; i<3; ++i) { prime[i]=minus[i]+mt[i]; s[i]=2*t[i]/(1+norm2(t)); }
    const auto ps = cross(prime,s);
    for (int i=0; i<3; ++i) { plus[i]=minus[i]+ps[i]; result[i]=plus[i]+.5*dt*qm*E[i]; }
    if (!finite(result) || !std::isfinite(norm2(result)))
        throw std::runtime_error("Boris update overflow");
    return result;
}

GapPic::GapPic(GapPicOptions options)
    : options_(checked(options)),
      mesh_(mfem::Mesh::MakeCartesian2D(options_.nx, options_.ny,
          mfem::Element::QUADRILATERAL, true, options_.width_m, options_.gap_m, false)),
      collection_(1,2), space_(&mesh_, &collection_), potential_(&space_) {
    potential_=0;
    horizontal_transfers_.resize(options_.nx*(options_.ny+1));
    vertical_transfers_.resize((options_.nx+1)*options_.ny);
}

int GapPic::locate(double x, double y, mfem::IntegrationPoint& point) const {
    if (!std::isfinite(x) || !std::isfinite(y) || x<0 || x>options_.width_m || y<0 || y>options_.gap_m)
        throw std::invalid_argument("Gap PIC point is outside rectangular domain");
    const double sx=x/options_.width_m*options_.nx, sy=y/options_.gap_m*options_.ny;
    const int ix=std::min(static_cast<int>(sx),options_.nx-1);
    const int iy=std::min(static_cast<int>(sy),options_.ny-1);
    point.Set2(sx-ix,sy-iy);
    return iy*options_.nx+ix;
}

void GapPic::add_particle(const GapParticle& p) {
    check_particle(p);
    mfem::IntegrationPoint point;
    static_cast<void>(locate(p.position_m[0],p.position_m[1],point));
    if(inside_conductor(p.position_m[0],p.position_m[1]))
        throw std::invalid_argument("Cannot place a particle inside a conductor");
    if(options_.finite_axial_extent && (p.position_m[2]<options_.short_z_m || p.position_m[2]>options_.feed_z_m))
        throw std::invalid_argument("Authored particle outside axial bounds");
    auto authored=p;
    authored.id=next_particle_id_++;
    particles_.push_back(authored);
    const double charge=p.charge_C*p.weight;
    // Authored populations have separate ownership from surface emission.
    ledger_.authored_charge_C += charge;
    ledger_.authored_kinetic_energy_J += pic_kinetic_energy_J(p)*p.weight;
    fields_current_=false;
}

void GapPic::inject_thermal_electrons(GapElectrode electrode, double T, double J,
                                     double dt, int count, std::uint64_t seed, GapSurface surface, double start_fraction, double end_fraction) {
    if ((electrode!=GapElectrode::Lower && electrode!=GapElectrode::Upper) ||
        !std::isfinite(T) || T<=0 || T>1e6 || !std::isfinite(J) || J<0 ||
        !std::isfinite(dt) || dt<=0 || count<1 ||
        !std::isfinite(start_fraction) || !std::isfinite(end_fraction) ||
        start_fraction<0 || start_fraction>=end_fraction || end_fraction>1)
        throw std::invalid_argument("Invalid thermal electron boundary source");
    double xlo=0,xhi=options_.width_m,ylo=electrode==GapElectrode::Lower ? 0 : options_.gap_m,yhi=ylo;
    int normal_axis=1; double normal_sign=electrode==GapElectrode::Lower ? 1 : -1;
    if(options_.finite_foils) {
        const auto f=electrode==GapElectrode::Lower ? options_.lower_foil : options_.upper_foil;
        xlo=f.left_m; xhi=f.right_m;
        switch(surface) {
        case GapSurface::Gap: ylo=electrode==GapElectrode::Lower ? f.top_m : f.bottom_m; yhi=ylo; break;
        case GapSurface::Back: ylo=electrode==GapElectrode::Lower ? f.bottom_m : f.top_m; yhi=ylo; normal_sign*=-1; break;
        case GapSurface::Left: xhi=xlo; ylo=f.bottom_m; yhi=f.top_m; normal_axis=0; normal_sign=-1; break;
        case GapSurface::Right: xlo=xhi; ylo=f.bottom_m; yhi=f.top_m; normal_axis=0; normal_sign=1; break;
        default: throw std::invalid_argument("Unknown foil surface");
        }
    } else if(surface!=GapSurface::Gap) throw std::invalid_argument("Flat-plate benchmark only has gap faces");
    const double dx=xhi-xlo,dy=yhi-ylo;
    xhi=xlo+end_fraction*dx; yhi=ylo+end_fraction*dy;
    xlo+=start_fraction*dx; ylo+=start_fraction*dy;
    if (J==0) return;
    const double length=normal_axis==1 ? xhi-xlo : yhi-ylo;
    const double weight=J*length*options_.depth_m*dt/(pic_e*count);
    if (!std::isfinite(weight) || weight<=0)
        throw std::invalid_argument("Thermal source weight overflow/underflow");
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> uniform(0,1);
    std::normal_distribution<double> normal(0,std::sqrt(kb*T/pic_me));
    std::vector<GapParticle> incoming;
    incoming.reserve(static_cast<std::size_t>(count));
    double energy=0;
    for (int n=0; n<count; ++n) {
        GapParticle p;
        p.id=next_particle_id_+static_cast<std::uint64_t>(n);
        const double along=uniform(rng);
        p.position_m={xlo+along*(xhi-xlo),ylo+along*(yhi-ylo),0};
        if(options_.finite_axial_extent)p.position_m[2]=options_.short_z_m+uniform(rng)*options_.depth_m;
        const double draw=std::max(uniform(rng),std::numeric_limits<double>::min());
        PicVector velocity={normal(rng),std::sqrt(-2*kb*T/pic_me*std::log(draw)),normal(rng)};
        if(normal_axis==0) std::swap(velocity[0],velocity[1]);
        velocity[normal_axis]*=normal_sign;
        const double v2=norm2(velocity);
        if (v2>=pic_c*pic_c) throw std::runtime_error("Nonrelativistic thermal sample exceeded c");
        for (int i=0;i<3;++i) p.proper_velocity_m_s[i]=velocity[i]/std::sqrt(1-v2/(pic_c*pic_c));
        p.weight=weight;
        check_particle(p);
        energy+=pic_kinetic_energy_J(p)*p.weight;
        incoming.push_back(p);
    }
    particles_.insert(particles_.end(),incoming.begin(),incoming.end());
    next_particle_id_+=static_cast<std::uint64_t>(count);
    const double charge=-pic_e*weight*count;
    if (electrode==GapElectrode::Lower) ledger_.lower_injected_charge_C+=charge;
    else ledger_.upper_injected_charge_C+=charge;
    ledger_.injected_energy_J+=energy;
    fields_current_=false;
}

void GapPic::inject_surface_profile(GapElectrode electrode, GapSurface surface,
    const std::vector<double>& T, const std::vector<double>& J,
    double dt, int count, std::uint64_t seed) {
    if(T.empty() || T.size()!=J.size() || count<1 || !std::isfinite(dt) || dt<=0)
        throw std::invalid_argument("Invalid surface profile");
    double total=0;
    for(std::size_t i=0;i<T.size();++i) {
        if(!std::isfinite(T[i]) || T[i]<=0 || T[i]>1e6 || !std::isfinite(J[i]) || J[i]<0)
            throw std::invalid_argument("Invalid surface profile patch");
        total+=J[i];
    }
    if(!std::isfinite(total)) throw std::invalid_argument("Surface source overflow");
    if(total==0) return;
    std::mt19937_64 rng(seed);
    std::discrete_distribution<std::size_t> patch(J.begin(),J.end());
    for(int n=0;n<count;++n) {
        const auto i=patch(rng);
        inject_thermal_electrons(electrode,T[i],total/count,dt,1,rng(),surface,
            double(i)/T.size(),double(i+1)/T.size());
    }
}

void GapPic::set_external_fields(std::function<PicVector(double,double,double)> electric,
                                std::function<PicVector(double,double,double)> magnetic) {
    external_electric_=std::move(electric); external_magnetic_=std::move(magnetic);
}

void GapPic::set_external_electric_3d(std::function<PicVector(double,double,double,double)> correction) {
    external_electric_3d_=std::move(correction);
}

void GapPic::set_potentials(double lower, double upper) {
    if (!std::isfinite(lower) || !std::isfinite(upper))
        throw std::invalid_argument("Electrode potentials must be finite");
    options_.lower_potential_V=lower;
    options_.upper_potential_V=upper;
    fields_current_=false;
}

void GapPic::set_empty_potentials(double lower,double upper) {
    if(!particles_.empty() || !fields_current_ || !std::isfinite(lower) || !std::isfinite(upper) ||
       options_.upper_potential_V==options_.lower_potential_V)
        throw std::invalid_argument("Empty potential rescaling requires solved vacuum with nonzero prior bias");
    const double ratio=(upper-lower)/(options_.upper_potential_V-options_.lower_potential_V);
    for(int i=0;i<potential_.Size();++i)potential_[i]=lower+ratio*(potential_[i]-options_.lower_potential_V);
    ledger_.lower_field_charge_C*=ratio;ledger_.upper_field_charge_C*=ratio;
    ledger_.gauss_defect_C*=ratio;ledger_.field_energy_J*=ratio*ratio;
    ledger_.poisson_residual_C_per_m*=std::abs(ratio);
    options_.lower_potential_V=lower;options_.upper_potential_V=upper;
}

void GapPic::inject_richardson_electrons(GapElectrode electrode,double T,double phi,
                                        double A,double dt,int count,std::uint64_t seed,GapSurface surface) {
    inject_thermal_electrons(electrode,T,pic_richardson_supply_A_m2(T,phi,A),dt,count,seed,surface);
}

bool GapPic::inside_conductor(double x,double y) const {
    if(!options_.finite_foils) return false;
    for(const auto& f:{options_.lower_foil,options_.upper_foil})
        if(x>f.left_m && x<f.right_m && y>f.bottom_m && y<f.top_m) return true;
    return false;
}

void GapPic::solve_fields() {
 pic_bench::Timer timer("solve_fields_total");
    if(!poisson_)poisson_=std::make_unique<PicPoissonState>(space_);
    auto& form=poisson_->form;
    mfem::LinearForm rhs(&space_); rhs=0;
    mfem::Array<int> dofs;
    mfem::Vector shape;
    double charge=0;
    for (const auto& p : particles_) {
        mfem::IntegrationPoint point;
        const int element=locate(p.position_m[0],p.position_m[1],point);
        const auto* fe=space_.GetFE(element);
        shape.SetSize(fe->GetDof()); fe->CalcShape(point,shape);
        space_.GetElementDofs(element,dofs);
        const double q=p.charge_C*p.weight/options_.depth_m;
        rhs.AddElementVector(dofs,shape *= q);
        charge+=q*options_.depth_m;
    }
    const double deposited=rhs.Sum()*options_.depth_m;
    mfem::Array<int> marker(mesh_.bdr_attributes.Max()); marker=0;
    mfem::Array<int> vertices;
    for (int b=0;b<mesh_.GetNBE();++b) {
        mesh_.GetBdrElementVertices(b,vertices);
        const double y0=mesh_.GetVertex(vertices[0])[1], y1=mesh_.GetVertex(vertices[1])[1];
        if (y0==y1 && (y0==0 || y0==options_.gap_m)) marker[mesh_.GetBdrAttribute(b)-1]=1;
    }
    mfem::Array<int> essential;
    space_.GetEssentialTrueDofs(marker,essential);
    mfem::FunctionCoefficient voltage([this](const mfem::Vector& x) {
        return options_.lower_potential_V+(options_.upper_potential_V-options_.lower_potential_V)*x[1]/options_.gap_m;
    });
    mfem::GridFunction trial(&space_); trial=potential_; trial.ProjectBdrCoefficient(voltage,marker);
    auto& xfield=poisson_->xfield;auto& yfield=poisson_->yfield;
    std::vector<int> owner(static_cast<std::size_t>(trial.Size()),-1);
    if(options_.finite_foils) {
        essential.SetSize(0);
        const double tol=1e-12*std::max(options_.width_m,options_.gap_m);
        for(int i=0;i<trial.Size();++i) {
            int who=0;
            for(const auto& f:{options_.lower_foil,options_.upper_foil}) {
                if(xfield[i]>=f.left_m-tol && xfield[i]<=f.right_m+tol &&
                   yfield[i]>=f.bottom_m-tol && yfield[i]<=f.top_m+tol) {
                    owner[i]=who; essential.Append(i);
                    trial[i]=who==0 ? options_.lower_potential_V : options_.upper_potential_V;
                }
                ++who;
            }
        }
    }
    mfem::Vector full_rhs(rhs);
    auto& full_matrix=*poisson_->full_matrix;
    mfem::OperatorPtr matrix;
    mfem::Vector x,b;
    form.FormLinearSystem(essential,trial,rhs,matrix,x,b);
    auto* sparse=dynamic_cast<mfem::SparseMatrix*>(matrix.Ptr());
    if (!sparse) throw std::runtime_error("Gap PIC expected sparse Poisson matrix");
#ifdef MFEM_USE_MPI
    if(options_.distributed_fields) {
        if(!poisson_->distributed){poisson_->distributed=std::make_unique<raikou_pic_parallel::LinearSystem>(*sparse);poisson_->address=sparse;}
        if(poisson_->address!=sparse)throw std::runtime_error("Cached Poisson operator changed identity");
        poisson_->distributed->Solve(b,x);
    } else
#endif
    {
    mfem::GSSmoother preconditioner(*sparse);
    mfem::CGSolver solver;
    solver.SetOperator(*matrix); solver.SetPreconditioner(preconditioner);
    solver.SetRelTol(1e-12); solver.SetAbsTol(1e-30); solver.SetMaxIter(4000); solver.SetPrintLevel(-1);
    solver.iterative_mode=true;
    solver.Mult(b,x);
    if (!solver.GetConverged()) throw std::runtime_error("Gap PIC Poisson solve failed");
    }
    mfem::Vector residual(b.Size()); matrix->Mult(x,residual); residual-=b;
    form.RecoverFEMSolution(x,rhs,trial);
    mfem::Vector reaction(trial.Size()); full_matrix.Mult(trial,reaction); reaction-=full_rhs;
    // Find electrode nodal charges through Q1 nodal coordinate fields.
    double lower=0,upper=0;
    for (int i=0;i<essential.Size();++i) {
        const int d=essential[i];
        if (options_.finite_foils ? owner[d]==0 : yfield[d]<.5*options_.gap_m) lower+=reaction[d]*options_.depth_m;
        else upper+=reaction[d]*options_.depth_m;
    }
    mfem::Vector action(trial.Size()); full_matrix.Mult(trial,action);
    ledger_.particle_charge_C=charge;
    ledger_.deposited_charge_C=deposited;
    ledger_.lower_field_charge_C=lower; ledger_.upper_field_charge_C=upper;
    ledger_.gauss_defect_C=charge+lower+upper;
    ledger_.field_energy_J=.5*(trial*action)*options_.depth_m;
    ledger_.poisson_residual_C_per_m=residual.Norml2();
    potential_=trial;
    fields_current_=true;
}

PicVector GapPic::electric_field(double x, double y) {
    if (!fields_current_) throw std::logic_error("Solve gap fields before gathering");
    mfem::IntegrationPoint point;
    const int element=locate(x,y,point);
    auto* transformation=mesh_.GetElementTransformation(element);
    transformation->SetIntPoint(&point);
    mfem::Vector gradient(2); potential_.GetGradient(*transformation,gradient);
    return {-gradient[0],-gradient[1],0};
}

void GapPic::drift(double dt) {
    const int nx=options_.nx,ny=options_.ny;
    const double hx=options_.width_m/nx,hy=options_.gap_m/ny;
    std::vector<double> defect(static_cast<std::size_t>((nx+1)*(ny+1)),0);
    const auto nodal=[&](double x,double y,double q) {
        mfem::IntegrationPoint ip;
        const int e=locate(x,y,ip),ix=e%nx,iy=e/nx;
        defect[iy*(nx+1)+ix]+=q*(1-ip.x)*(1-ip.y);
        defect[iy*(nx+1)+ix+1]+=q*ip.x*(1-ip.y);
        defect[(iy+1)*(nx+1)+ix]+=q*(1-ip.x)*ip.y;
        defect[(iy+1)*(nx+1)+ix+1]+=q*ip.x*ip.y;
    };
    const auto folded=[&](double x) {
        if(options_.finite_foils) return std::clamp(x,0.,options_.width_m);
        double v=std::fmod(x,2*options_.width_m);
        if (v<0) v+=2*options_.width_m;
        return v>options_.width_m ? 2*options_.width_m-v : v;
    };
    const auto edge_segment=[&](double x0,double y0,double x1,double y1,double q) {
        mfem::IntegrationPoint ip;
        const int e=locate(.5*(x0+x1),.5*(y0+y1),ip),ix=e%nx,iy=e/nx;
        const double xi=.5*(x0+x1)/hx-ix,eta=.5*(y0+y1)/hy-iy;
        const double dx=(x1-x0)/hx,dy=(y1-y0)/hy;
        const double bottom=q*dx*(1-eta),top=q*dx*eta,left=q*dy*(1-xi),right=q*dy*xi;
        horizontal_transfers_[iy*nx+ix]+=bottom;
        horizontal_transfers_[(iy+1)*nx+ix]+=top;
        vertical_transfers_[iy*(nx+1)+ix]+=left;
        vertical_transfers_[iy*(nx+1)+ix+1]+=right;
        defect[iy*(nx+1)+ix]+=bottom+left;
        defect[iy*(nx+1)+ix+1]+=-bottom+right;
        defect[(iy+1)*(nx+1)+ix]+=top-left;
        defect[(iy+1)*(nx+1)+ix+1]+=-top-right;
    };
    std::vector<GapParticle> remaining;
    remaining.reserve(particles_.size());
    for (auto p : particles_) {
        const double g=gamma(p.proper_velocity_m_s);
        const double vx=p.proper_velocity_m_s[0]/g,vy=p.proper_velocity_m_s[1]/g;
        const double x0=p.position_m[0],y0=p.position_m[1],q=p.charge_C*p.weight;
        double flight=dt;
        bool absorbed=false;
        int hit=-1;
        if(options_.finite_foils) {
            const auto exit=[&](double t) { if(t>=0 && t<=flight) {flight=t;absorbed=true;hit=2;} };
            if(vx<0) exit(-x0/vx);
            if(vx>0) exit((options_.width_m-x0)/vx);
            if(vy<0) exit(-y0/vy);
            if(vy>0) exit((options_.gap_m-y0)/vy);
            int who=0;
            for(const auto& f:{options_.lower_foil,options_.upper_foil}) {
                double enter=-std::numeric_limits<double>::infinity(),leave=std::numeric_limits<double>::infinity();
                bool possible=true;
                const auto slab=[&](double pos,double v,double low,double high) {
                    if(v==0) {if(pos<low || pos>high) possible=false; return;}
                    double a=(low-pos)/v,b=(high-pos)/v;
                    if(a>b) std::swap(a,b);
                    enter=std::max(enter,a); leave=std::min(leave,b);
                };
                slab(x0,vx,f.left_m,f.right_m); slab(y0,vy,f.bottom_m,f.top_m);
                if(possible && leave>0 && enter<=leave && std::max(enter,0.)<=flight) {
                    flight=std::max(enter,0.); absorbed=true; hit=who;
                }
                ++who;
            }
        } else {
            if (vy<0 && -y0/vy<=flight) { flight=-y0/vy; absorbed=true;hit=0; }
            if (vy>0 && (options_.gap_m-y0)/vy<=flight) { flight=(options_.gap_m-y0)/vy; absorbed=true;hit=1; }
        }
        if(options_.finite_axial_extent) {
            const double vz=p.proper_velocity_m_s[2]/g;
            const auto axial=[&](double t,bool short_end) {
                if(t<0 || t>flight)return;
                const double x=x0+vx*t,y=y0+vy*t;
                const bool bridge=short_end && x>=options_.lower_foil.left_m && x<=options_.lower_foil.right_m &&
                    y>=options_.lower_foil.top_m && y<=options_.upper_foil.bottom_m;
                flight=t;absorbed=true;hit=bridge?3:4;
            };
            if(vz<0)axial((options_.short_z_m-p.position_m[2])/vz,true);
            if(vz>0)axial((options_.feed_z_m-p.position_m[2])/vz,false);
        }
        std::vector<double> times{0,flight};
        const auto split=[&](double pos,double velocity,double spacing) {
            if (velocity==0) return;
            const double finish=pos+velocity*flight;
            const double low=std::min(pos,finish)/spacing,high=std::max(pos,finish)/spacing;
            if (!std::isfinite(low) || !std::isfinite(high) || high-low>100000)
                throw std::runtime_error("Gap PIC drift crosses too many cells; reduce time step");
            for (double k=std::floor(low)+1;k<high;++k) {
                const double t=(k*spacing-pos)/velocity;
                if (t>0 && t<flight) times.push_back(t);
            }
        };
        split(x0,vx,hx); split(y0,vy,hy);
        std::sort(times.begin(),times.end());
        nodal(x0,y0,-q);
        for (std::size_t k=1;k<times.size();++k) {
            const double t0=times[k-1],t1=times[k];
            const double a=std::clamp(y0+vy*t0,0.,options_.gap_m);
            const double b=std::clamp(y0+vy*t1,0.,options_.gap_m);
            edge_segment(folded(x0+vx*t0),a,folded(x0+vx*t1),b,q);
        }
        double raw=x0+vx*flight;
        p.position_m[0]=folded(raw);
        p.position_m[1]=std::clamp(y0+vy*flight,0.,options_.gap_m);
        p.position_m[2]+=flight*p.proper_velocity_m_s[2]/g;
        if (!finite(p.position_m)) throw std::runtime_error("Gap PIC drift overflow");
        // End charge is either still in the domain or belongs to the explicit
        // boundary sink. Both enter continuity with the same endpoint shape.
        nodal(p.position_m[0],p.position_m[1],q);
        if (absorbed) {
            const double energy=pic_kinetic_energy_J(p)*p.weight;
            if (hit==0) { ledger_.lower_absorbed_charge_C+=q; ledger_.lower_absorbed_energy_J+=energy; }
            else if(hit==1) { ledger_.upper_absorbed_charge_C+=q; ledger_.upper_absorbed_energy_J+=energy; }
            else if(hit==3){ledger_.short_absorbed_charge_C+=q;ledger_.short_absorbed_energy_J+=energy;}
            else { ledger_.escaped_charge_C+=q; ledger_.escaped_kinetic_energy_J+=energy;if(hit==4)ledger_.axial_escaped_charge_C+=q; }
            continue;
        }
        raw=std::fmod(raw,2*options_.width_m);
        if (raw<0) raw+=2*options_.width_m;
        if (!options_.finite_foils && raw>options_.width_m) p.proper_velocity_m_s[0]*=-1;
        remaining.push_back(p);
    }
    for (const double value:defect) ledger_.maximum_nodal_continuity_defect_C=
        std::max(ledger_.maximum_nodal_continuity_defect_C,std::abs(value));
    particles_=std::move(remaining);
    fields_current_=false;
}

void GapPic::advance(double dt) {
 pic_bench::Timer timer("advance_total");
    if (!std::isfinite(dt) || dt<=0 || !std::isfinite(ledger_.time_s+dt))
        throw std::invalid_argument("Gap PIC advance requires positive finite dt");
    std::fill(horizontal_transfers_.begin(),horizontal_transfers_.end(),0);
    std::fill(vertical_transfers_.begin(),vertical_transfers_.end(),0);
    // Midpoint position splitting; full-time velocities. Absorption is resolved
    // at each half drift. Wall interaction timing must be checked by refinement.
    drift(.5*dt); solve_fields();
    for (auto& p : particles_) {
        auto E=electric_field(p.position_m[0],p.position_m[1]);
        auto B=options_.magnetic_field_T;
        if(external_electric_) {
            auto extra=external_electric_(p.position_m[0],p.position_m[1],ledger_.time_s+.5*dt);
            for(int i=0;i<3;++i) E[i]+=extra[i];
        }
        if(external_electric_3d_) {
            auto extra=external_electric_3d_(p.position_m[0],p.position_m[1],p.position_m[2],ledger_.time_s+.5*dt);
            for(int i=0;i<3;++i)E[i]+=extra[i];
        }
        if(external_magnetic_) B=external_magnetic_(p.position_m[0],p.position_m[1],ledger_.time_s+.5*dt);
        p.proper_velocity_m_s=pic_boris(p.proper_velocity_m_s,p.charge_C/p.mass_kg,E,B,dt);
    }
    drift(.5*dt); solve_fields();
    ledger_.time_s+=dt;
}

void GapPic::advance_empty(double dt) {
    if(!particles_.empty() || !fields_current_ || !std::isfinite(dt) || dt<=0 || !std::isfinite(ledger_.time_s+dt))
        throw std::invalid_argument("Empty advance requires solved vacuum and a finite positive interval");
    ledger_.time_s+=dt;
    std::fill(horizontal_transfers_.begin(),horizontal_transfers_.end(),0);
    std::fill(vertical_transfers_.begin(),vertical_transfers_.end(),0);
}

void GapPic::save_checkpoint(const std::string& path) const {
    if(!fields_current_)throw std::logic_error("Checkpoint requires solved fields");
    std::ofstream f(path+".tmp",std::ios::binary);
    const auto put=[&](const auto& value){f.write(reinterpret_cast<const char*>(&value),sizeof(value));};
    const std::uint64_t magic=0x5241494b50494332ULL,n=particles_.size();
    put(magic);put(options_);put(ledger_);put(next_particle_id_);put(n);
    for(const auto& particle:particles_)put(particle);
    f.close();if(!f)throw std::runtime_error("Checkpoint write failed");
    std::filesystem::rename(path+".tmp",path);
}

void GapPic::load_checkpoint(const std::string& path) {
    if(!particles_.empty() || ledger_.time_s!=0)throw std::logic_error("Restart requires fresh PIC state");
    std::ifstream f(path,std::ios::binary);
    const auto get=[&](auto& value){f.read(reinterpret_cast<char*>(&value),sizeof(value));};
    std::uint64_t magic=0,n=0,id=0;GapPicOptions saved;GapPicDiagnostics ledger;
    get(magic);get(saved);get(ledger);get(id);get(n);
    if(!f || magic!=0x5241494b50494332ULL || n>2000000 || saved.nx!=options_.nx || saved.ny!=options_.ny ||
       saved.width_m!=options_.width_m || saved.gap_m!=options_.gap_m || saved.depth_m!=options_.depth_m ||
       saved.lower_potential_V!=options_.lower_potential_V || saved.upper_potential_V!=options_.upper_potential_V ||
       saved.finite_foils!=options_.finite_foils || saved.finite_axial_extent!=options_.finite_axial_extent ||
       saved.short_z_m!=options_.short_z_m || saved.feed_z_m!=options_.feed_z_m || saved.distributed_fields!=options_.distributed_fields ||
       saved.magnetic_field_T!=options_.magnetic_field_T ||
       std::memcmp(&saved.lower_foil,&options_.lower_foil,sizeof(GapFoil)) ||
       std::memcmp(&saved.upper_foil,&options_.upper_foil,sizeof(GapFoil)) ||
       !std::isfinite(ledger.time_s) || ledger.time_s<0 || id<1)
        throw std::runtime_error("Checkpoint incompatible or invalid");
    std::vector<GapParticle> recovered(n);
    for(auto& particle:recovered){get(particle);check_particle(particle);mfem::IntegrationPoint ip;static_cast<void>(locate(particle.position_m[0],particle.position_m[1],ip));if(inside_conductor(particle.position_m[0],particle.position_m[1]) || (options_.finite_axial_extent && (particle.position_m[2]<options_.short_z_m || particle.position_m[2]>options_.feed_z_m)))throw std::runtime_error("Checkpoint particle outside vacuum bounds");if(particle.id>=id)throw std::runtime_error("Checkpoint particle identity invalid");}
    if(!f || f.peek()!=std::char_traits<char>::eof())throw std::runtime_error("Checkpoint truncated or trailing bytes");
    particles_=std::move(recovered);ledger_=ledger;next_particle_id_=id;solve_fields();
}

GapPicDiagnostics GapPic::diagnostics() const {
    if (!fields_current_) throw std::logic_error("Solve gap fields before diagnostics");
    auto result=ledger_;
    for (const auto& p : particles_) result.particle_kinetic_energy_J+=pic_kinetic_energy_J(p)*p.weight;
    result.particle_charge_defect_C=result.particle_charge_C+
        result.lower_absorbed_charge_C+result.upper_absorbed_charge_C+result.short_absorbed_charge_C-
        result.lower_injected_charge_C-result.upper_injected_charge_C-result.authored_charge_C+result.escaped_charge_C;
    return result;
}
} // namespace raikou::mfem_backend
