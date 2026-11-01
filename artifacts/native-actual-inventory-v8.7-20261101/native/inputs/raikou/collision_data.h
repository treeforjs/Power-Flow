#pragma once
// Strict observational collision data. This layer never evolves particles/F0,
// debits a parent, creates products, deposits current or commits a coupled step.
#include "raikou/electron_boltzmann.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace raikou::kinetics {

class CollisionCoverageError : public std::runtime_error {
public:
    explicit CollisionCoverageError(const std::string& text) : std::runtime_error(text) {}
};

struct CollisionEnergyWeight { double energy_eV{0}, physical_electrons{0}; };
struct CollisionCoverage {
    double total_weight{0}, covered_weight{0}, uncovered_weight{0};
    // The denominator is ALL represented electrons, including uncovered ones.
    // This is a known-support contribution, never a complete rate if uncovered>0.
    double known_support_coefficient_m3_s{0};
    std::string observable;
    bool complete{false};
};

class StrictCollisionDataset {
    Json data_;
    static void require(bool condition,const std::string& text) {
        if (!condition) throw std::invalid_argument("Collision data: "+text);
    }
    static double number(const Json& value) {
        require(value.is_number(),"expected a number");
        const double x=value.get<double>();require(std::isfinite(x),"nonfinite data");return x;
    }
    [[nodiscard]] const Json& record(const std::string& id) const {
        for (const auto& p:data_.at("channels")) if (p.at("id")==id) return p;
        throw std::invalid_argument("Unknown collision observable: "+id);
    }
    static double speed(double eV) {
        constexpr double rest_eV=510998.95069,c=299792458.;
        const double q=eV/rest_eV;
        return c*std::sqrt(q*(q+2))/(q+1);
    }
    static double beb(const Json& p,double energy) {
        double sum=0;
        const double a=p.at("bohr_radius_m").get<double>(),R=p.at("rydberg_eV").get<double>();
        for (const auto& orbital:p.at("orbitals")) {
            const double B=orbital.at("B_eV").get<double>();
            if (energy<=B) continue;
            const double t=energy/B,u=orbital.at("U_eV").get<double>()/B;
            const double Q=orbital.at("Q").get<double>(),N=orbital.at("N").get<double>();
            const double S=4*std::numbers::pi*a*a*N*(R/B)*(R/B),l=std::log(t);
            sum+=S/(t+u+1)*(Q*l/2*(1-1/(t*t))+(2-Q)*(1-1/t-l/(t+1)));
        }
        return sum;
    }
public:
    explicit StrictCollisionDataset(Json data) : data_(std::move(data)) {
        require(data_.at("schema")=="raikou-strict-gas-collision-data-v1","unsupported schema");
        require(data_.at("phase")=="gas"&&data_.at("target")=="H2O","gas-phase H2O required");
        require(!data_.at("full_fragmentation_ready").get<bool>(),"this dataset is not full fragmentation");
        require(data_.at("channels").is_array()&&!data_.at("channels").empty(),"empty channel inventory");
        std::vector<std::string> ids;
        for (const auto& p:data_.at("channels")) {
            const auto id=p.at("id").get<std::string>();
            require(!id.empty()&&std::find(ids.begin(),ids.end(),id)==ids.end(),"duplicate channel ID");ids.push_back(id);
            require(p.at("phase")=="gas"&&p.at("target")=="H2O","wrong target phase");
            require(p.at("outside_support")=="error","non-strict extrapolation policy");
            require(data_.at("sources").contains(p.at("source").get<std::string>()),"missing source provenance");
            const auto& domain=p.at("support_eV");require(domain.size()==2,"bad support domain");
            const double lo=number(domain[0]),hi=number(domain[1]);require(lo>=0&&hi>lo,"invalid energy support");
            require(p.at("unit_SI")=="m2"||p.at("unit_SI")=="eV m2","unsupported dimensions");
            if (!p.at("zero_below_eV").is_null()) require(number(p.at("zero_below_eV"))>=0,"negative threshold");
            const auto model=p.at("model").get<std::string>();
            if (model=="linear_table") {
                const auto& e=p.at("energy_eV");const auto& y=p.at("values");
                require(e.is_array()&&y.is_array()&&e.size()==y.size()&&e.size()>=2,"invalid table dimensions");
                require(number(p.at("ordinate_scale_to_SI"))>0,"invalid unit conversion");
                require(number(e.front())==lo&&number(e.back())==hi,"support/table mismatch");
                double previous=-1;
                for (std::size_t i=0;i<e.size();++i) {
                    const double x=number(e[i]);require(x>previous&&x>=0,"energy grid not strictly increasing");
                    require(number(y[i])>=0,"negative cross section");previous=x;
                }
                if (!p.at("zero_below_eV").is_null()) require(number(p.at("zero_below_eV"))<=lo,"threshold conflicts with data onset");
            } else if (model=="nist_BEB") {
                require(hi<=10000&&p.at("unit_SI")=="m2","nonrelativistic BEB domain exceeded");
                require(number(p.at("bohr_radius_m"))>0&&number(p.at("rydberg_eV"))>0,"bad BEB constants");
                require(p.at("orbitals").size()==4,"frozen valence model requires four orbitals");
                for (const auto& o:p.at("orbitals")) {
                    require(number(o.at("B_eV"))>0&&number(o.at("U_eV"))>0&&number(o.at("N"))>0,"bad orbital");
                    const auto q=number(o.at("Q"));require(q>=0&&q<=2,"invalid dipole constant");
                }
            } else require(false,"unknown model");
        }
    }
    static StrictCollisionDataset read(const std::string& path) {
        std::ifstream in(path);require(bool(in),"cannot open dataset");Json j;in>>j;
        in>>std::ws;require(in.eof(),"trailing content after dataset");
        return StrictCollisionDataset(std::move(j));
    }
    [[nodiscard]] const Json& metadata() const { return data_; }
    [[nodiscard]] const Json& channel_metadata(const std::string& id) const { return record(id); }
    [[nodiscard]] double observable_value_SI(const std::string& id,double energy_eV) const {
        require(std::isfinite(energy_eV)&&energy_eV>=0,"incident energy must be finite nonnegative eV");
        const auto& p=record(id);const auto& domain=p.at("support_eV");
        const double lo=domain[0].get<double>(),hi=domain[1].get<double>();
        if (energy_eV>hi) throw CollisionCoverageError(id+": above supported incident energy; no endpoint hold");
        if (!p.at("zero_below_eV").is_null()&&energy_eV<=p.at("zero_below_eV").get<double>()) return 0;
        if (energy_eV<lo) throw CollisionCoverageError(id+": threshold-to-first-point or low-energy coverage gap");
        if (p.at("model")=="nist_BEB") return beb(p,energy_eV);
        const auto e=p.at("energy_eV").get<std::vector<double>>();const auto y=p.at("values").get<std::vector<double>>();
        const double scale=p.at("ordinate_scale_to_SI").get<double>();
        if (energy_eV==e.back()) return y.back()*scale;
        const auto k=static_cast<std::size_t>(std::upper_bound(e.begin(),e.end(),energy_eV)-e.begin()-1);
        const double q=(energy_eV-e[k])/(e[k+1]-e[k]);return ((1-q)*y[k]+q*y[k+1])*scale;
    }
    [[nodiscard]] CollisionCoverage coverage(const std::string& id,const std::vector<CollisionEnergyWeight>& population) const {
        const auto& p=record(id);require(p.at("unit_SI")=="m2","rate query requires area, not an energy-loss moment");
        CollisionCoverage out;out.observable=p.at("observable").get<std::string>();
        double numerator=0;
        for (const auto& sample:population) {
            require(std::isfinite(sample.physical_electrons)&&sample.physical_electrons>=0,"invalid physical weight");
            require(std::isfinite(sample.energy_eV)&&sample.energy_eV>=0,"invalid incident energy");
            if (sample.physical_electrons==0) continue;
            out.total_weight+=sample.physical_electrons;
            try {
                numerator+=sample.physical_electrons*observable_value_SI(id,sample.energy_eV)*speed(sample.energy_eV);
                out.covered_weight+=sample.physical_electrons;
            } catch (const CollisionCoverageError&) { out.uncovered_weight+=sample.physical_electrons; }
        }
        require(out.total_weight>0,"empty weighted electron population");
        require(std::isfinite(numerator)&&std::isfinite(out.total_weight),"rate accumulation overflow");
        out.known_support_coefficient_m3_s=numerator/out.total_weight;out.complete=out.uncovered_weight==0;return out;
    }
    [[nodiscard]] double strict_observable_rate_m3_s(const std::string& id,const std::vector<CollisionEnergyWeight>& population) const {
        const auto evidence=coverage(id,population);
        if (!evidence.complete) throw CollisionCoverageError(id+": incomplete weighted energy support; no tail renormalization");
        return evidence.known_support_coefficient_m3_s;
    }
    // This gate protects event multiplicity but does NOT supply an event sampler.
    // Yield-only and aggregate data cannot become an exclusive parent debit.
    void require_exclusive_products(const std::string& id) const {
        const auto& p=record(id);
        if (p.at("exclusive_products").is_null()) throw CollisionCoverageError(id+": observable yield has no exclusive product contract");
    }
    [[nodiscard]] Json kernel_status(const std::string& id,const std::string& selected_model="") const {
        const auto& p=record(id);Json available=Json::array();
        if (data_.contains("conditional_event_models")) for (const auto& model:data_.at("conditional_event_models")) {
            if (model.at("channel")==id) {
                available.push_back(model);
                if (!selected_model.empty()&&model.at("id")==selected_model)
                    return {{"status","conditional_available"},{"reason","explicit source-rate/kinematic comparison model; not measured differential validation"},{"model",model}};
            }
        }
        return {{"status",selected_model.empty()?"not_selected":"unavailable"},
            {"reason",p.at("exclusive_products").is_null()?"inclusive/aggregate observable has no exclusive product mapping":"no implicitly selected secondary/internal/angular/recoil closure"},
            {"conditional_models",available}};
    }
    void require_MCC_event_kernel(const std::string& id,const std::string& selected_model="",
                                  double incident_eV=std::numeric_limits<double>::quiet_NaN()) const {
        const auto status=kernel_status(id,selected_model);
        if (status.at("status")!="conditional_available") throw CollisionCoverageError(id+": "+status.at("reason").get<std::string>());
        if (!status.at("model").contains("products")||!status.at("model").at("products").is_object())
            throw CollisionCoverageError(id+": selected kernel has no explicit conditional product mapping");
        const auto& interval=status.at("model").at("data_support_eV");
        if (!std::isfinite(incident_eV)||incident_eV<interval[0].get<double>()||incident_eV>interval[1].get<double>())
            throw CollisionCoverageError(id+": selected conditional kernel outside its explicit incident-energy domain");
        if (!(observable_value_SI(id,incident_eV)>0)) throw CollisionCoverageError(id+": no positive supported event cross section");
    }
    [[nodiscard]] std::vector<CollisionProcess> reduced_transport_processes(double minimum_eV,double maximum_eV) const {
        require(minimum_eV>=.001&&maximum_eV<=100&&maximum_eV>minimum_eV,"reduced transport support is0.001..100eV");
        std::vector<CollisionProcess> result;
        for (const std::string id:{"H2O:elastic_momentum","H2O:excite_v010"}) {
            const auto& p=record(id);const bool elastic=p.at("collision_process_admission")=="elastic_momentum_only";
            require(elastic||p.at("collision_process_admission")=="one_way_bend","unadmitted reduced process");
            CollisionProcess q;q.label=id;q.target="H2O";q.angular=CollisionProcess::Angular::Isotropic;
            q.kind=elastic?CollisionKind::Elastic:CollisionKind::Excitation;
            if (elastic) q.mass_ratio=9.1093837139e-31/(18.0106*1.66053906660e-27);
            else { q.product="H2O_v010";q.threshold_ev=.198;q.group=CollisionProcess::Group::Vibrational; }
            const double first=elastic?minimum_eV:std::max(minimum_eV,.198);
            if (first>=maximum_eV) continue;
            q.energy_ev.push_back(first);
            for (const double e:p.at("energy_eV").get<std::vector<double>>()) if (e>first&&e<maximum_eV) q.energy_ev.push_back(e);
            if (maximum_eV>first) q.energy_ev.push_back(maximum_eV);
            require(q.energy_ev.size()>=2,"reduced native process needs at least two points");
            for (const auto e:q.energy_ev) q.cross_section_m2.push_back(observable_value_SI(id,e));
            result.push_back(std::move(q));
        }
        return result;
    }
};
} // namespace raikou::kinetics
