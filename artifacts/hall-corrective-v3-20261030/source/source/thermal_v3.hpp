#pragma once
// Hall thermal owner v3: the single electrode thermal owner of the coupled step
// (temperature, sigma(T), Joule heating) with external surface loads.
//
// Derived from hall-canonical-driver-20261029/axial_thermal_v2.hpp (commit
// 8dddec0, sealed, unchanged). Unchanged: one heat-v2 Coupler per foil and foil
// axial cell, the same geometry, material law (evaluated solid 298.2..1500 K),
// Joule coefficient and per-slice prepare/preflight/commit/save/load. With no
// declared source patches and no loads the arithmetic is the sealed one, so the
// heat state is bit-identical to Hall v2.
// Changes:
//  1. ThermalDeclaration: heat-only source patches (desorption disabled, zero
//     inventory) per foil slice, so that a source owner can read this owner's
//     accepted surface temperature (patch mean or point value) through a const
//     view of the slice Coupler instead of running a private heat engine. Heat-only
//     patches do not change the thermal solve; the declaration enters the
//     material binding and the configuration identity.
//  2. prepare() takes surface loads (heat-v2 Transfers routed to one foil slice):
//     Hall's own collected-carrier wall heat and external loads (emission
//     cooling, desorption binding/kinetic/internal energy, ion impact, certified
//     source exchanges). Exactly-once accounting: producer and energy-owner
//     watermarks of the sealed Coupler (strictly increasing sequences, no
//     duplicate within a step), reserved Hall producer names, no external Joule
//     channel, and the per-step closure thermal_delta = Joule + wall + external
//     within the sealed roundoff bound.
//  3. Typed scoped refusals (ThermalRefusal) instead of an untyped engine error or
//     a process abort: every load is checked for finiteness and against the
//     enthalpy capacity of its slice before the engine runs (a load larger than
//     the whole material range, or a net heating/cooling that no temperature
//     field inside 298.2..1500 K can absorb, can never be admitted); engine
//     errors after that are rethrown typed with the slice identity.
#include "field_contract.hpp"
#include "hall-canonical-driver-20261029/axial_thermal_v2.hpp"
#include <map>
#include <span>

namespace hall_corrective_v3 {
namespace eh=electrode_heat_v2;

// Typed, scoped thermal refusal. Every refusal leaves all accepted state unchanged.
struct ThermalRefusal:std::runtime_error {
    std::string kind,quantity,identity;double observed=0,bound=0,log_abs_scope=0;
    ThermalRefusal(std::string k,std::string q,std::string id,double obs,double b,const std::string&detail)
        :std::runtime_error("thermal refusal ["+k+"] "+q+" at "+id+": "+detail),kind(std::move(k)),quantity(std::move(q)),identity(std::move(id)),
         observed(obs),bound(b),log_abs_scope(obs!=0&&std::isfinite(obs)?std::log(std::abs(obs)):std::numeric_limits<double>::infinity()){}
};

// One surface energy transfer into one foil slice (body 0 lower, 1 upper; slice =
// axial cell index k, kshort <= k < kfeed). The Transfer's (x, y) is the point on
// the foil cross-section boundary in the field frame (the slice Coupler's frame).
struct SurfaceLoad {int body=-1,slice=-1;std::string channel;eh::Transfer transfer;};
inline constexpr char hall_load_prefix[]="hall-v3/";
inline constexpr char wall_heat_producer[]="hall-v3/kinetic-sink-wall-heat";
inline constexpr char wall_heat_owner[]="hall-v3/kinetic-sink-energy";
inline constexpr char wall_heat_channel[]="kinetic-sink-collected-wall-heat";
// Energy deposited into the solid by one transfer (same partition as the sealed Coupler).
inline double transfer_energy_J(const eh::Transfer&x){
    if(x.certified_source_exchange)return x.source_wall_J;
    const double boundary=x.complete_nottingham?x.nottingham_wall_J:x.incoming_K_J-x.outgoing_K_J-x.reference_out_J;
    return boundary+x.solid_collision_J+x.authored_heat_J;
}
inline bool mentions_joule(std::string s){for(char&c:s)c=char(std::tolower(static_cast<unsigned char>(c)));return s.find("joule")!=std::string::npos;}
inline void encode_transfer(Digest&d,const eh::Transfer&x){
    d.text("hall-v3-surface-transfer-v1");for(const std::string*s:{&x.producer,&x.state_id,&x.parent_id,&x.channel_id,&x.energy_owner})d.text(*s);
    d.u64(x.sequence);d.u64(x.owner_sequence);d.f64(x.event_time_s);d.boolean(x.surface);
    for(double v:{x.x,x.y,x.incoming_K_J,x.outgoing_K_J,x.reference_out_J,x.solid_collision_J,x.authored_heat_J})d.f64(v);
    d.boolean(x.complete_nottingham);d.f64(x.nottingham_wall_J);d.f64(x.work_function_reference_part_J);d.boolean(x.certified_source_exchange);
    for(double v:{x.source_wall_J,x.source_particle_K_change_J,x.source_electron_reference_change_J,x.source_image_energy_change_J,x.source_field_work_J})d.f64(v);
}

// Heat-only source patches declared on foil slices (key: body, axial cell k).
struct SlicePatches {std::string source_hash;std::vector<eh::Patch> patches;};
struct ThermalDeclaration {
    std::map<std::pair<int,int>,SlicePatches> slices;
    bool empty()const noexcept{return slices.empty();}
    std::string sha256()const{
        Digest d("hall-v3-thermal-declaration-v1");d.u64(slices.size());
        for(const auto&[key,s]:slices){d.i64(key.first);d.i64(key.second);d.text(s.source_hash);d.u64(s.patches.size());
            for(const auto&p:s.patches){d.text(p.id);d.i64(int(p.face));d.f64(p.begin);d.f64(p.end);d.boolean(p.desorption_enabled);
                d.f64(p.initial_particles);d.f64(p.initial_low_particles);d.text(p.source_hash);}}
        return d.hex();
    }
};

struct ThermalTrialV3 {
    std::vector<eh::Trial> bodies;
    double admitted_joule_J=0,thermal_delta_J=0,numerical_residual_J=0,roundoff_bound_J=0;
    // Exactly-once load accounting of this step.
    double wall_heat_J=0;            // Hall's collected-carrier wall heat (incoming kinetic energy)
    double external_J=0,external_heating_J=0,external_cooling_J=0;
    std::uint64_t loads=0,external_loads=0;
    std::map<std::string,Words> channel_J;   // per load channel (compensated sum)
    std::string loads_sha256;                 // canonical record of every applied load in application order
};

class ThermalV3 {
    const Mesh&g;int nx,ny;double initial_time_s;ThermalDeclaration declaration_;eh::Material law_;
    std::vector<std::unique_ptr<eh::Coupler>> owners;
    struct JouleCoefficient:mfem::Coefficient {   // sealed axial_thermal_v2.hpp JouleCoefficient, unchanged
        const Mesh&g;const Vec&power;int k,body;
        JouleCoefficient(const Mesh&mesh,const Vec&values,int slice,int which):g(mesh),power(values),k(slice),body(which){}
        mfem::real_t Eval(mfem::ElementTransformation&tr,const mfem::IntegrationPoint&ip)override{
            mfem::Vector x(2);tr.Transform(ip,x);int i=std::floor((x[0]+g.side/2)/g.dx),j=std::floor((x[1]+g.side/2)/g.dy);
            hall::require(i>=0&&i<g.nx&&j>=0&&j<g.ny,"thermal quadrature outside EM mesh");
            int c=(k*g.ny+j)*g.nx+i;hall::require(g.cells[c].body==body,"Joule transfer body mismatch");
            return power.at(c)/g.cells[c].volume;
        }
    };
    static eh::Material material_law(){return eh::conditional_ss304_v2(ss304v2::Qualification::evaluated_solid_1500,
        ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation,8000.,eh::material_headers_sha256);}
    std::vector<std::unique_ptr<eh::Coupler>> make_owners()const{
        std::vector<std::unique_ptr<eh::Coupler>> result;
        for(int k=g.kshort;k<g.kfeed;++k)for(int body=0;body<2;++body){
            eh::Geometry geometry;geometry.nx=nx;geometry.ny=ny;geometry.x0=-.0025;
            geometry.y0=body?.0025:-.0027;geometry.width=.005;geometry.thickness=.0002;geometry.height=g.dz;
            geometry.revision="axial-foil-yee-heat-v2/"+std::to_string(g.nx)+"/"+std::to_string(g.nz_foil);
            geometry.electrode="foil-"+std::to_string(body)+"-slice-"+std::to_string(k);
            auto it=declaration_.slices.find({body,k});
            if(it==declaration_.slices.end())
                result.push_back(std::make_unique<eh::Coupler>(geometry,material_law(),std::vector<eh::Patch>{},300.,eh::source_kernel_sha256,initial_time_s));
            else result.push_back(std::make_unique<eh::Coupler>(geometry,material_law(),it->second.patches,300.,it->second.source_hash,initial_time_s));
        }
        return result;
    }
    static ThermalDeclaration checked(const Mesh&g,ThermalDeclaration d){
        for(auto&[key,s]:d.slices){need(key.first>=0&&key.first<2&&key.second>=g.kshort&&key.second<g.kfeed,"thermal declaration slice outside the foils");
            need(!s.source_hash.empty()&&!s.patches.empty()&&s.patches.size()<=256,"thermal declaration requires a source hash and 1..256 patches");
            for(const auto&p:s.patches)need(!p.desorption_enabled&&p.initial_particles==0&&p.initial_low_particles==0&&!p.coverage_temkin&&
                p.source_hash==s.source_hash,"thermal declaration patches are heat-only (no desorption, zero inventory)");}
        return d;
    }
    std::string slice_text(std::size_t i)const{return "foil-"+std::to_string(i%2)+"-slice-"+std::to_string(g.kshort+int(i/2));}
public:
    // Explicit evaluated-solid scenario: 298.2..1500 K, rho_ref=8000 kg/m3,
    // fixed reference volume and Cp-as-Cv. No mechanics or axial diffusion.
    ThermalV3(const Mesh&mesh,double start,ThermalDeclaration declaration={})
        :g(mesh),nx(std::lround(.005/g.dx)),ny(std::lround(.0002/g.dy)),initial_time_s(start),
         declaration_(checked(mesh,std::move(declaration))),law_(material_law()),owners(make_owners()){}
    // Material binding: the sealed record when no patches are declared (same digest
    // as hall_driver::Thermal::material_binding), else with the declared patch count
    // and the declaration digest.
    static Json material_binding(const ThermalDeclaration&d){
        Json j=hall_driver::Thermal::material_binding();
        if(!d.empty()){std::size_t n=0;for(const auto&[k,s]:d.slices)n+=s.patches.size();
            j["source_patch_count"]=n;j["heat_only_source_patches_sha256"]=d.sha256();}
        return j;
    }
    const ThermalDeclaration& declaration()const noexcept{return declaration_;}
    std::vector<eh::State> states()const{std::vector<eh::State> r;for(const auto&o:owners)r.push_back(o->state());return r;}
    std::size_t slices()const noexcept{return owners.size();}
    std::size_t index(int body,int k)const{need(body>=0&&body<2&&k>=g.kshort&&k<g.kfeed,"thermal slice outside the foils");return std::size_t(2*(k-g.kshort)+body);}
    // Const view of one slice owner (accepted state, temperatures, patches) for a
    // source owner that consumes this owner's temperature.
    const eh::Coupler& slice(int body,int k)const{return *owners.at(index(body,k));}
    Vec conductivity()const{
        Vec result(g.cells.size());
        for(std::size_t c=0;c<g.cells.size();++c){const auto&cell=g.cells[c];if(cell.body<0)continue;
            result[c]=owners.at(2*(cell.k-g.kshort)+cell.body)->electrical_conductivity_at((cell.i+.5)*g.dx-g.side/2,(cell.j+.5)*g.dy-g.side/2);}
        return result;
    }
    double time_s()const{return owners.front()->state().time;}
    std::uint64_t revision()const{return owners.front()->state().revision;}
    double max_temperature_K()const{double m=300;for(const auto&o:owners)for(double T:o->state().T)m=std::max(m,T);return m;}
    double energy_J()const{double s=0;for(const auto&o:owners)s+=o->stored_energy_J();return s;}
    double slice_volume_m3()const{return .005*.0002*g.dz;}
    // Enthalpy capacity bounds of one slice at its accepted state (J).
    struct Capacity {double heating_room_J=0,cooling_room_J=0,full_range_J=0;};
    Capacity capacity(std::size_t i)const{
        const auto&T=owners.at(i)->state().T;const double lo=law_.knots.front().T,hi=law_.knots.back().T,V=slice_volume_m3();
        double tmin=hi,tmax=lo;for(double t:T){tmin=std::min(tmin,t);tmax=std::max(tmax,t);}
        return {V*(law_.energy_reference(hi)-law_.energy_reference(tmin)),V*(law_.energy_reference(tmax)-law_.energy_reference(lo)),
                V*(law_.energy_reference(hi)-law_.energy_reference(lo))};
    }
    // Typed admission of one step's loads (and the Joule power) against the
    // capacity of every slice; runs before the engine.
    void admit_loads(const axial::Receipt&field,double dt,std::span<const SurfaceLoad>loads)const{
        std::vector<long double>in(owners.size()),out(owners.size());
        for(std::size_t c=0;c<g.cells.size();++c){const auto&cell=g.cells[c];if(cell.body<0)continue;
            in[2*(cell.k-g.kshort)+cell.body]+=(long double)field.cell_joule_W[c]*dt;}
        for(const auto&l:loads){const std::size_t i=index(l.body,l.slice);const auto&x=l.transfer;
            for(double v:{x.incoming_K_J,x.outgoing_K_J,x.reference_out_J,x.solid_collision_J,x.authored_heat_J,x.nottingham_wall_J,
                          x.work_function_reference_part_J,x.source_wall_J,x.source_particle_K_change_J,x.source_electron_reference_change_J,
                          x.source_image_energy_change_J,x.source_field_work_J,x.x,x.y,x.event_time_s})
                if(!std::isfinite(v))throw ThermalRefusal("nonfinite-load","surface load operand",slice_text(i)+";producer="+x.producer+";channel="+l.channel,
                    v,0,"a surface load operand is not finite; no engine call");
            const double e=transfer_energy_J(x);const auto cap=capacity(i);
            if(!(std::abs(e)<=cap.full_range_J))throw ThermalRefusal("load-exceeds-material-range","surface load energy J",
                slice_text(i)+";producer="+x.producer+";sequence="+std::to_string(x.sequence)+";channel="+l.channel,e,cap.full_range_J,
                "one load exceeds the enthalpy of the whole slice between 298.2 K and 1500 K; no engine call");
            if(e>0)in[i]+=e;else out[i]-=e;}
        for(std::size_t i=0;i<owners.size();++i){const auto cap=capacity(i);const long double net=in[i]-out[i];
            if(!std::isfinite(double(net)))throw ThermalRefusal("nonfinite-load","net slice energy J",slice_text(i),double(net),0,"net slice energy not finite");
            if(net>(long double)cap.heating_room_J)throw ThermalRefusal("net-heating-exceeds-capacity","net slice energy J",slice_text(i),double(net),
                cap.heating_room_J,"Joule plus surface heating exceeds the enthalpy room up to 1500 K of the coldest accepted node; no engine call");
            if(-net>(long double)cap.cooling_room_J)throw ThermalRefusal("net-cooling-exceeds-capacity","net slice energy J",slice_text(i),double(net),
                -cap.cooling_room_J,"surface cooling exceeds the enthalpy room down to 298.2 K of the hottest accepted node; no engine call");}
    }
    // Field receipt (Joule power) and loads -> one thermal trial. Loads of a slice
    // are applied by that slice's sealed Coupler (watermarks, deposit point checks).
    ThermalTrialV3 prepare(const axial::Receipt&field,double dt,const std::string&field_state_hash,std::span<const SurfaceLoad>loads={}){
        hall::require(std::isfinite(dt)&&dt>0&&dt<=1e-14&&time_s()+dt>time_s(),"v2 heat common interval guard");
        hall::require(field.cell_joule_W.size()==g.cells.size()&&!field_state_hash.empty(),"v2 heat field binding dimensions/identity");
        hall::require(std::isfinite(field.joule_J)&&field.joule_J>=0,"v2 heat invalid Joule operand");
        for(std::size_t c=0;c<g.cells.size();++c)
            hall::require(std::isfinite(field.cell_joule_W[c])&&field.cell_joule_W[c]>=0&&
                (g.cells[c].body>=0||field.cell_joule_W[c]==0),"v2 heat Joule power outside metal or invalid");
        ThermalTrialV3 trial;std::vector<std::vector<eh::Transfer>> routed(owners.size());
        if(!loads.empty()){
            for(const auto&l:loads){const auto&x=l.transfer;
                const bool own=x.producer==wall_heat_producer&&x.energy_owner==wall_heat_owner&&l.channel==wall_heat_channel;
                if(!own&&(x.producer.rfind(hall_load_prefix,0)==0||x.energy_owner.rfind(hall_load_prefix,0)==0))
                    throw ThermalRefusal("reserved-producer","surface load identity",x.producer+"/"+x.energy_owner,0,0,"external load uses a Hall-reserved producer or energy owner");
                if(mentions_joule(l.channel)||mentions_joule(x.producer)||mentions_joule(x.energy_owner)||mentions_joule(x.channel_id))
                    throw ThermalRefusal("external-joule","surface load channel",x.producer+"/"+l.channel,transfer_energy_J(x),0,
                        "Joule heat is owned by the Hall field/heat proposal; an external Joule load would count it twice");
                if(!x.surface)throw ThermalRefusal("volume-load","surface load","producer="+x.producer,transfer_energy_J(x),0,"external loads are surface transfers");
                if(l.channel.empty())throw ThermalRefusal("unnamed-channel","surface load channel","producer="+x.producer,transfer_energy_J(x),0,"load channel required");
                (void)index(l.body,l.slice);}
            admit_loads(field,dt,loads);
            Digest d("hall-v3-surface-load-list-v1");d.u64(loads.size());
            for(const auto&l:loads){const auto&x=l.transfer;routed[index(l.body,l.slice)].push_back(x);
                const double e=transfer_energy_J(x);++trial.loads;auto&w=trial.channel_J[l.channel];w=w+Words(e);
                if(x.producer==wall_heat_producer)trial.wall_heat_J+=e;
                else{++trial.external_loads;trial.external_J+=e;if(e>0)trial.external_heating_J+=e;else trial.external_cooling_J-=e;}
                d.i64(l.body);d.i64(l.slice);d.text(l.channel);encode_transfer(d,x);}
            trial.loads_sha256=d.hex();
        }
        for(int k=g.kshort;k<g.kfeed;++k)for(int body=0;body<2;++body){
            const std::size_t i=std::size_t(2*(k-g.kshort)+body);auto&owner=*owners.at(i);
            hall::require(owner.state().time==time_s()&&owner.state().revision==revision(),"thermal owners lost common clock");
            JouleCoefficient power(g,field.cell_joule_W,k,body);
            try{auto part=owner.prepare(owner.make_epoch(dt,field_state_hash),power,routed[i],std::vector<std::vector<eh::Draw>>(owner.state().inventory.size()));
                const auto&b=part.budget();
                trial.admitted_joule_J+=b.joule_J;trial.thermal_delta_J+=b.thermal_delta_J;
                trial.numerical_residual_J+=b.residual_J;trial.roundoff_bound_J+=b.roundoff_bound_J;
                trial.bodies.push_back(std::move(part));}
            catch(const std::bad_alloc&){throw;}
            catch(const ThermalRefusal&){throw;}
            catch(const std::exception&e){
                if(loads.empty())throw;   // v2 behavior (no load): the sealed message unchanged
                throw ThermalRefusal("thermal-engine","slice enthalpy step",slice_text(i),0,0,e.what());}
        }
        hall::require(std::abs(trial.admitted_joule_J-field.joule_J)<=1e-22+2e-12*std::abs(field.joule_J),
                      "field-to-heat transfer is not the same Joule operand");
        if(trial.loads){const double expect=trial.admitted_joule_J+trial.wall_heat_J+trial.external_J,
            scale=std::abs(trial.admitted_joule_J)+std::abs(trial.wall_heat_J)+trial.external_heating_J+trial.external_cooling_J+std::abs(trial.thermal_delta_J);
            if(!(std::abs(trial.thermal_delta_J-expect)<=trial.roundoff_bound_J+2e-12*scale))
                throw ThermalRefusal("load-accounting","thermal delta minus Joule and loads J","all-slices",trial.thermal_delta_J-expect,
                    trial.roundoff_bound_J,"exactly-once closure failed: thermal change differs from Joule + wall heat + external loads");}
        return trial;
    }
    bool preflight(const ThermalTrialV3&trial)const noexcept{
        if(trial.bodies.size()!=owners.size())return false;
        for(std::size_t i=0;i<owners.size();++i)if(!owners[i]->preflight(trial.bodies[i]))return false;
        return true;
    }
    bool commit(ThermalTrialV3&&trial)noexcept{
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
        for(const auto&owner:candidate)hall::require(owner->state().revision==first.revision&&owner->state().time==first.time,"v2 heat checkpoint slice clocks differ");
        owners.swap(candidate);
    }
};
} // namespace hall_corrective_v3
