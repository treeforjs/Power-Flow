#pragma once
#include "raikou/mfem_gap_pic.h"
#include "raikou/json.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
using namespace raikou::mfem_backend;
using Json=raikou::Json;
namespace hall {
constexpr double eps=8.8541878188e-12, proton_mass=1.67262192595e-27;
inline double dot(PicVector a,PicVector b){double s=0;for(int k=0;k<3;++k)s+=a[k]*b[k];return s;}
inline double gamma(PicVector u){return std::sqrt(1+dot(u,u)/(pic_c*pic_c));}
inline PicVector velocity(PicVector u){const double g=gamma(u);for(double&v:u)v/=g;return u;}
inline void require(bool b,const char* text){if(!b)throw std::runtime_error(text);}
inline std::uint64_t file_checksum(const std::filesystem::path& p){
    std::ifstream f(p,std::ios::binary);require(bool(f),"checkpoint missing");
    std::uint64_t h=14695981039346656037ULL;char c;
    while(f.get(c)){h^=static_cast<unsigned char>(c);h*=1099511628211ULL;}
    require(f.eof(),"checkpoint read error");return h;
}
struct Tag {int species; std::string event; PicVector reference;};
// Prototype adapter: one particle population owner, immutable two-species table,
// explicit authored events, and a separately versioned restart envelope.
class CoupledPic {
    GapPic pic_;
    double mass_ratio_;
    std::map<std::uint64_t,Tag> tags_;
    std::set<std::string> events_;
    std::uint64_t next_=1;
    void observe(const GapParticle&p,const PicVector&E,const PicVector&,const PicVector&u1,double dt){
        const int s=tags_.at(p.id).species;
        PicVector minus=p.proper_velocity_m_s,plus=u1,D{};
        for(int k=0;k<3;++k){const double a=.5*dt*p.charge_C/p.mass_kg*E[k];minus[k]+=a;plus[k]-=a;}
        for(int k=0;k<3;++k)D[k]=.5*dt*((minus[k]+p.proper_velocity_m_s[k])/(gamma(minus)+gamma(p.proper_velocity_m_s))+(u1[k]+plus[k])/(gamma(u1)+gamma(plus)));
        const double w=p.weight*p.charge_C*dot(E,D);
        auto q=p;q.proper_velocity_m_s=u1;
        const double dk=(pic_kinetic_energy_J(q)-pic_kinetic_energy_J(p))*p.weight;
        work[s]+=w;work_defect[s]+=dk-w;work_abs_defect[s]+=std::abs(dk-w);work_activity[s]+=std::abs(w);
    }
    void attach(){pic_.set_push_observer([this](const auto&p,const auto&E,const auto&B,const auto&u,double dt){observe(p,E,B,u,dt);});}
public:
    std::array<double,2> work{},work_defect{},work_abs_defect{},work_activity{};
    explicit CoupledPic(GapPicOptions o,double ratio):pic_(o),mass_ratio_(ratio){
        require(!o.distributed_fields,"prototype requires serial CPU fields");
        require(std::isfinite(ratio)&&ratio>0,"invalid ion mass");attach();
    }
    GapPic& native(){return pic_;}
    const auto& tags()const{return tags_;}
    double ratio()const{return mass_ratio_;}
    Json registry()const{return Json::array({{{"id",0},{"name","electron"},{"charge_C",-pic_e},{"mass_kg",pic_me}},{{"id",1},{"name","manufactured_positive_ion"},{"charge_C",pic_e},{"mass_kg",pic_me*mass_ratio_}}});}
    void authored(const std::string&event,const std::vector<GapParticle>&p,const std::vector<int>&species,const std::vector<PicVector>&reference){
        require(!event.empty()&&!events_.count(event),"duplicate or missing authored event");
        require(p.size()==species.size()&&p.size()==reference.size()&&!p.empty(),"authored metadata size");
        auto tags=tags_;auto events=events_;events.insert(event);auto next=next_;
        for(std::size_t i=0;i<p.size();++i){
            require(species[i]==0||species[i]==1,"unknown species");
            require(p[i].charge_C==(species[i]?pic_e:-pic_e)&&p[i].mass_kg==(species[i]?pic_me*mass_ratio_:pic_me),"species charge/mass mismatch");
            for(double v:reference[i])require(std::isfinite(v),"invalid reference position");
            tags.emplace(next++,Tag{species[i],event,reference[i]});
        }
        const auto ids=pic_.add_authored_batch(p);
        // Sole insertion owner: IDs are allocated once, before any field update.
        require(ids.front()==next_&&ids.back()==next-1,"adapter/native ID mismatch");
        tags_.swap(tags);events_.swap(events);next_=next;
    }
    void guard(double dt){
        const auto&o=pic_.options();const double h=std::min(o.width_m/o.nx,o.gap_m/o.ny);
        require(dt>0&&dt<=1e-14,"10 fs populated step cap");
        std::vector<double>w2(o.nx*o.ny,0.);const double volume=o.width_m*o.gap_m*o.depth_m/(o.nx*o.ny);
        for(const auto&p:pic_.particles()){
            require(std::abs(p.charge_C/p.mass_kg)*std::sqrt(dot(o.magnetic_field_T,o.magnetic_field_T))*dt<=.05,"gyro angle guard");
            require(std::sqrt(dot(velocity(p.proper_velocity_m_s),velocity(p.proper_velocity_m_s)))*dt<=.2*h,"displacement guard");
            const int i=std::min(o.nx-1,int(p.position_m[0]/o.width_m*o.nx)),j=std::min(o.ny-1,int(p.position_m[1]/o.gap_m*o.ny));
            w2[j*o.nx+i]+=p.weight*p.charge_C*p.charge_C/(eps*p.mass_kg*volume);
        }
        for(double x:w2)require(dt*std::sqrt(x)<=.1,"sum-species plasma frequency guard");
    }
    void advance(double dt){guard(dt);pic_.advance(dt);}
    void save(const std::filesystem::path&dir){
        require(!std::filesystem::exists(dir),"checkpoint destination exists");
        auto tmp=dir;tmp+=".staging";require(std::filesystem::create_directory(tmp),"checkpoint staging exists");
        pic_.save_checkpoint((tmp/"native.bin").string());Json tags=Json::array();
        for(const auto&[id,t]:tags_)tags.push_back({{"id",id},{"species",t.species},{"event",t.event},{"reference",t.reference}});
        Json j={{"format","hall-authored-pic-envelope-v1"},{"model","collisionless-2D3V-Poisson-prescribed-B"},{"registry",registry()},{"tags",tags},{"events",events_},{"next",next_},{"work_J",work},{"work_defect_J",work_defect},{"work_abs_defect_J",work_abs_defect},{"work_activity_J",work_activity},{"native_checksum_fnv1a64",file_checksum(tmp/"native.bin")},{"executable_checksum_fnv1a64",file_checksum("/proc/self/exe")}};
        std::ofstream f(tmp/"species.json");f<<j.dump(2)<<'\n';f.close();require(bool(f),"metadata write failed");
        std::filesystem::rename(tmp,dir);
    }
    void load(const std::filesystem::path&dir){
        require(tags_.empty()&&pic_.particles().empty(),"restart requires fresh adapter");
        Json j;std::ifstream(dir/"species.json")>>j;
        require(j.at("format")=="hall-authored-pic-envelope-v1"&&j.at("model")=="collisionless-2D3V-Poisson-prescribed-B"&&j.at("registry")==registry(),"restart species/model mismatch");
        require(j.at("executable_checksum_fnv1a64").get<std::uint64_t>()==file_checksum("/proc/self/exe"),"restart requires same executable");
        require(j.at("native_checksum_fnv1a64").get<std::uint64_t>()==file_checksum(dir/"native.bin"),"restart native checksum mismatch");
        auto events=j.at("events").get<std::set<std::string>>();std::map<std::uint64_t,Tag>tags;
        for(const auto&t:j.at("tags")){
            Tag tag{t.at("species"),t.at("event"),t.at("reference").get<PicVector>()};
            require((tag.species==0||tag.species==1)&&events.count(tag.event),"invalid restart tag");
            for(double v:tag.reference)require(std::isfinite(v),"invalid restart reference");
            require(tags.emplace(t.at("id"),tag).second,"duplicate restart ID");
        }
        pic_.load_checkpoint((dir/"native.bin").string());
        const auto next=j.at("next").get<std::uint64_t>();
        require(next>0&&tags.size()==next-1,"restart ID sequence mismatch");
        std::uint64_t expected=1;for(const auto&[id,tag]:tags)require(id==expected++,"nonconsecutive restart lineage");
        for(const auto&p:pic_.particles()){
            const int s=tags.at(p.id).species;
            require(p.charge_C==(s?pic_e:-pic_e)&&p.mass_kg==(s?pic_me*mass_ratio_:pic_me)&&p.id<next,"restart particle/species mismatch");
        }
        tags_=std::move(tags);events_=std::move(events);next_=next;
        work=j.at("work_J").get<std::array<double,2>>();work_defect=j.at("work_defect_J").get<std::array<double,2>>();work_abs_defect=j.at("work_abs_defect_J").get<std::array<double,2>>();work_activity=j.at("work_activity_J").get<std::array<double,2>>();
        for(const auto&values:{work,work_defect,work_abs_defect,work_activity})for(double v:values)require(std::isfinite(v),"restart nonfinite energy ledger");
        attach();
    }
};
}
