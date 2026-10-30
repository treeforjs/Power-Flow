#pragma once
#include "axial_fields.hpp"
#include "thermal-owner/ss304_heat_material.hpp"
#include "thermal-owner/material/source_seal.hpp"

namespace axial {
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
    std::vector<electrode_heat::Trial> bodies;
    double admitted_joule_J=0,thermal_delta_J=0,numerical_residual_J=0,roundoff_bound_J=0;
};
class Thermal {
    const Mesh&g;int nx,ny;
    std::vector<std::unique_ptr<electrode_heat::Coupler>> owners;
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
public:
    Thermal(const Mesh&mesh,double start):g(mesh),nx(std::lround(.005/g.dx)),ny(std::lround(.0002/g.dy)){
        auto law=electrode_heat::conditional_ss304(ss304::Qualification::evaluated_solid_1500,
            ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation,8000.,ss304_source_seal_sha256);
        for(int k=g.kshort;k<g.kfeed;++k)for(int body=0;body<2;++body){
            electrode_heat::Geometry geometry;geometry.nx=nx;geometry.ny=ny;geometry.x0=-.0025;
            geometry.y0=body?.0025:-.0027;geometry.width=.005;geometry.thickness=.0002;geometry.height=g.dz;
            geometry.revision="axial-foil-yee-v1/"+std::to_string(g.nx)+"/"+std::to_string(g.nz_foil);
            geometry.electrode="foil-"+std::to_string(body)+"-slice-"+std::to_string(k);
            owners.push_back(std::make_unique<electrode_heat::Coupler>(geometry,law,std::vector<electrode_heat::Patch>{},300.,"joule-only-no-inventory-v1",start));
        }
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
    void load(const std::filesystem::path&dir){for(std::size_t i=0;i<owners.size();++i)owners[i]->load((dir/(std::to_string(i)+".bin")).string());}
};
} // namespace axial
