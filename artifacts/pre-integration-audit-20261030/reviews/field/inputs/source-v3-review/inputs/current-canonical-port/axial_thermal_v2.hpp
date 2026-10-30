#pragma once
#include "../hall-finite-foil-em-20261029/axial_fields.hpp"
#include "bindings/heat-v2/ss304_heat_material_v2.hpp"

namespace hall_driver {
using axial::Mesh;using axial::Vec;using axial::Receipt;
class CurrentTrace {
    std::vector<std::pair<double,double>> knots;
public:
    explicit CurrentTrace(const std::string&path){
        std::ifstream file(path);hall::require(bool(file),"missing sealed C19 trace");std::string line;
        std::getline(file,line);hall::require(line.find("time_s")!=std::string::npos&&line.find("current_A")!=std::string::npos,"C19 header mismatch");
        while(std::getline(file,line)){std::istringstream row(line);std::vector<std::string>cells;std::string value;
            while(std::getline(row,value,','))cells.push_back(value);hall::require(cells.size()==12,"C19 row schema");
            double time=std::stod(cells[4]),current=std::stod(cells[8]);
            hall::require(std::isfinite(time)&&std::isfinite(current)&&(knots.empty()||time>knots.back().first),"invalid C19 knot");
            knots.emplace_back(time,current);
        }
        hall::require(knots.size()>1,"empty C19 trace");
    }
    double start_s()const{return knots.front().first;}
    double end_s()const{return knots.back().first;}
    double current(double t)const{
        hall::require(std::isfinite(t)&&t>=start_s()&&t<=end_s(),"C19 query outside support; no extrapolation");
        if(t==end_s())return knots.back().second;
        auto upper=std::upper_bound(knots.begin(),knots.end(),t,[](double t,const auto&v){return t<v.first;});
        auto lower=upper-1;double f=(t-lower->first)/(upper->first-lower->first);
        return lower->second+f*(upper->second-lower->second);
    }
};
struct ThermalTrial {
    std::vector<electrode_heat_v2::Trial> bodies;
    double admitted_joule_J=0,thermal_delta_J=0,numerical_residual_J=0,roundoff_bound_J=0;
};
class Thermal {
    const Mesh&g;int nx,ny;double initial_time_s;
    std::vector<std::unique_ptr<electrode_heat_v2::Coupler>> owners;
    struct JouleCoefficient:mfem::Coefficient {
        const Mesh&g;const Vec&power;int k,body;
        JouleCoefficient(const Mesh&mesh,const Vec&values,int slice,int which):g(mesh),power(values),k(slice),body(which){}
        mfem::real_t Eval(mfem::ElementTransformation&tr,const mfem::IntegrationPoint&ip)override{
            mfem::Vector x(2);tr.Transform(ip,x);int i=std::floor((x[0]+g.side/2)/g.dx),j=std::floor((x[1]+g.side/2)/g.dy);
            hall::require(i>=0&&i<g.nx&&j>=0&&j<g.ny,"thermal quadrature outside EM mesh");
            int c=(k*g.ny+j)*g.nx+i;hall::require(g.cells[c].body==body,"Joule transfer body mismatch");
            return power.at(c)/g.cells[c].volume;
        }
    };
    std::vector<std::unique_ptr<electrode_heat_v2::Coupler>> make_owners()const{
        std::vector<std::unique_ptr<electrode_heat_v2::Coupler>> result;
        auto law=electrode_heat_v2::conditional_ss304_v2(ss304v2::Qualification::evaluated_solid_1500,
            ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation,8000.,electrode_heat_v2::material_headers_sha256);
        for(int k=g.kshort;k<g.kfeed;++k)for(int body=0;body<2;++body){
            electrode_heat_v2::Geometry geometry;geometry.nx=nx;geometry.ny=ny;geometry.x0=-.0025;
            geometry.y0=body?.0025:-.0027;geometry.width=.005;geometry.thickness=.0002;geometry.height=g.dz;
            geometry.revision="axial-foil-yee-heat-v2/"+std::to_string(g.nx)+"/"+std::to_string(g.nz_foil);
            geometry.electrode="foil-"+std::to_string(body)+"-slice-"+std::to_string(k);
            result.push_back(std::make_unique<electrode_heat_v2::Coupler>(geometry,law,std::vector<electrode_heat_v2::Patch>{},300.,electrode_heat_v2::source_kernel_sha256,initial_time_s));
        }
        return result;
    }
public:
    // Explicit evaluated-solid scenario: 298.2..1500 K, rho_ref=8000 kg/m3,
    // fixed reference volume and Cp-as-Cv. No mechanics or axial diffusion.
    Thermal(const Mesh&mesh,double start):g(mesh),nx(std::lround(.005/g.dx)),ny(std::lround(.0002/g.dy)),initial_time_s(start),owners(make_owners()){}
    static Json material_binding(){return {{"namespace","electrode_heat_v2"},
        {"material_headers_sha256",electrode_heat_v2::material_headers_sha256},
        {"source_kernel_sha256",electrode_heat_v2::source_kernel_sha256},
        {"qualification","evaluated_solid_1500"},{"minimum_temperature_K",298.2},
        {"maximum_temperature_K",1500.},{"rho_ref_kg_m3",8000.},
        {"volume_convention","fixed_reference_volume_cp_as_cv_approximation"},
        {"axial_diffusion",false},{"source_patch_count",0}};}
    std::vector<electrode_heat_v2::State> states()const{
        std::vector<electrode_heat_v2::State> result;
        for(const auto&owner:owners)result.push_back(owner->state());return result;
    }
    Vec conductivity()const{
        Vec result(g.cells.size());
        for(std::size_t c=0;c<g.cells.size();++c){const auto&cell=g.cells[c];if(cell.body<0)continue;
            result[c]=owners.at(2*(cell.k-g.kshort)+cell.body)->electrical_conductivity_at(
                (cell.i+.5)*g.dx-g.side/2,(cell.j+.5)*g.dy-g.side/2);
        }
        return result;
    }
    double time_s()const{return owners.front()->state().time;}
    std::uint64_t revision()const{return owners.front()->state().revision;}
    double max_temperature_K()const{double maximum=300;for(const auto&owner:owners)for(double T:owner->state().T)maximum=std::max(maximum,T);return maximum;}
    double energy_J()const{double sum=0;for(const auto&owner:owners)sum+=owner->stored_energy_J();return sum;}
    ThermalTrial prepare(const Receipt&field,double dt,const std::string&field_state_hash){
        hall::require(std::isfinite(dt)&&dt>0&&dt<=1e-14&&time_s()+dt>time_s(),"v2 heat common interval guard");
        hall::require(field.cell_joule_W.size()==g.cells.size()&&!field_state_hash.empty(),"v2 heat field binding dimensions/identity");
        hall::require(std::isfinite(field.joule_J)&&field.joule_J>=0,"v2 heat invalid Joule operand");
        for(std::size_t c=0;c<g.cells.size();++c)
            hall::require(std::isfinite(field.cell_joule_W[c])&&field.cell_joule_W[c]>=0&&
                (g.cells[c].body>=0||field.cell_joule_W[c]==0),"v2 heat Joule power outside metal or invalid");
        ThermalTrial trial;
        for(int k=g.kshort;k<g.kfeed;++k)for(int body=0;body<2;++body){
            auto&owner=*owners.at(2*(k-g.kshort)+body);
            hall::require(owner.state().time==time_s()&&owner.state().revision==revision(),"thermal owners lost common clock");
            JouleCoefficient power(g,field.cell_joule_W,k,body);
            auto part=owner.prepare(owner.make_epoch(dt,field_state_hash),power,{},{});const auto&b=part.budget();
            trial.admitted_joule_J+=b.joule_J;trial.thermal_delta_J+=b.thermal_delta_J;
            trial.numerical_residual_J+=b.residual_J;trial.roundoff_bound_J+=b.roundoff_bound_J;
            trial.bodies.push_back(std::move(part));
        }
        hall::require(std::abs(trial.admitted_joule_J-field.joule_J)<=1e-22+2e-12*std::abs(field.joule_J),
                      "field-to-heat transfer is not the same Joule operand");
        return trial;
    }
    bool preflight(const ThermalTrial&trial)const noexcept{
        if(trial.bodies.size()!=owners.size())return false;
        for(std::size_t i=0;i<owners.size();++i)if(!owners[i]->preflight(trial.bodies[i]))return false;
        return true;
    }
    bool commit(ThermalTrial&&trial)noexcept{
        if(!preflight(trial))return false;
        for(std::size_t i=0;i<owners.size();++i)if(!owners[i]->commit(std::move(trial.bodies[i])))std::terminate();
        return true;
    }
    void save(const std::filesystem::path&dir)const{
        hall::require(!std::filesystem::exists(dir),"thermal checkpoint destination exists");std::filesystem::create_directory(dir);
        for(std::size_t i=0;i<owners.size();++i)owners[i]->save((dir/(std::to_string(i)+".bin")).string());
    }
    void load(const std::filesystem::path&dir){
        hall::require(revision()==0,"v2 heat restart requires fresh aggregate");
        auto candidate=make_owners();
        for(std::size_t i=0;i<candidate.size();++i)candidate[i]->load((dir/(std::to_string(i)+".bin")).string());
        const auto&first=candidate.front()->state();
        for(const auto&owner:candidate)hall::require(owner->state().revision==first.revision&&
            owner->state().time==first.time,"v2 heat checkpoint slice clocks differ");
        // A corrupt later slice must not leave earlier accepted slices restored.
        owners.swap(candidate);
    }
};
} // namespace hall_driver
