#include "water_network_scope.hpp"
#include <openssl/evp.h>
#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>

namespace water_scope {
namespace {
constexpr std::uint32_t none = 0xffffffffu;
void need(bool ok, const std::string& why) { if (!ok) throw std::runtime_error("water scope: " + why); }
void digest_valid(const std::string& d) {
    need(d.size()==64 && d.find_first_not_of("0123456789abcdef")==std::string::npos,"invalid identity digest");
}
void append_u64(std::string& s, std::uint64_t n) { for (int i=0;i<8;++i) s.push_back(char((n>>(8*i))&255)); }
void append_text(std::string& s, const std::string& v) { append_u64(s,v.size()); s+=v; }
std::uint64_t read_le(std::span<const unsigned char> b, std::size_t off, std::size_t n) {
    need(n<=8 && off<=b.size() && n<=b.size()-off,"truncated table section");
    std::uint64_t v=0; for (std::size_t i=0;i<n;++i) v|=std::uint64_t(b[off+i])<<(8*i); return v;
}
std::vector<InventoryRow> canonical(std::span<const InventoryRow> rows) {
    std::vector<InventoryRow> result(rows.begin(),rows.end());
    for (const auto& row:result) need(std::isfinite(row.physical_count) && row.physical_count>=0,"invalid inventory count");
    std::sort(result.begin(),result.end(),[](const auto& a,const auto& b){return a.native_id<b.native_id;});
    for (std::size_t i=1;i<result.size();++i) need(result[i-1].native_id!=result[i].native_id,"duplicate native row ID");
    return result;
}
std::vector<Pair> required(const Catalogue& table) {
    auto pairs=table.pairs();
    for (const auto* target:{"OH","O","H2"}) pairs.push_back({table.species_id("e"),table.species_id(target)});
    std::sort(pairs.begin(),pairs.end()); pairs.erase(std::unique(pairs.begin(),pairs.end()),pairs.end());
    return pairs;
}
} // namespace

std::string sha256(std::span<const unsigned char> bytes) {
    unsigned char out[EVP_MAX_MD_SIZE]; unsigned n=0;
    need(EVP_Digest(bytes.data(),bytes.size(),out,&n,EVP_sha256(),nullptr)==1 && n==32,"digest failure");
    static constexpr char hex[]="0123456789abcdef"; std::string result;
    for (unsigned i=0;i<n;++i) {result.push_back(hex[out[i]>>4]);result.push_back(hex[out[i]&15]);}
    return result;
}
std::string sha256(const std::string& text) {
    return sha256(std::span(reinterpret_cast<const unsigned char*>(text.data()),text.size()));
}
std::string census_digest(std::span<const InventoryRow> rows) {
    const auto sorted=canonical(rows); std::string data="all-species-native-census-v1";
    append_u64(data,sorted.size());
    for (const auto& row:sorted) {
        append_u64(data,row.native_id);append_u64(data,row.collision_cell);append_u64(data,row.species);
        append_u64(data,std::bit_cast<std::uint64_t>(row.physical_count==0?0.0:row.physical_count));
    }
    return sha256(data);
}

Catalogue::Catalogue(std::span<const unsigned char> bytes,const std::string& expected) {
    digest_valid(expected);blob_sha_=sha256(bytes);need(blob_sha_==expected,"table identity differs");
    need(bytes.size()>=16 && std::string(reinterpret_cast<const char*>(bytes.data()),8)=="PFCHEMV3","table magic");
    const auto version=read_le(bytes,8,4);
    need(version==1 || version==2,"table schema");const auto ns=read_le(bytes,12,4);
    need(ns>0 && ns<64 && 16+40*ns<=bytes.size(),"table directory");
    std::map<std::string,std::span<const unsigned char>> sections;
    std::vector<std::pair<std::uint64_t,std::uint64_t>> intervals;
    for (std::uint64_t i=0;i<ns;++i) {
        const auto at=16+40*i;std::string name;
        for (std::size_t j=0;j<24 && bytes[at+j];++j) name.push_back(char(bytes[at+j]));
        const auto off=read_le(bytes,at+24,8),len=read_le(bytes,at+32,8);
        need(!name.empty() && off%64==0 && off>=16+40*ns && off<=bytes.size() && len<=bytes.size()-off,"table section bounds");
        need(sections.emplace(name,bytes.subspan(off,len)).second,"duplicate section");
        intervals.emplace_back(off,off+len);
    }
    std::sort(intervals.begin(),intervals.end());
    for (std::size_t i=1;i<intervals.size();++i) need(intervals[i-1].second<=intervals[i].first,"overlapping sections");
    auto section=[&](const char* name){const auto it=sections.find(name);need(it!=sections.end(),std::string("missing section ")+name);return it->second;};
    // Scope reads encounter coverage only. Native's qualified C3 loader has
    // already validated all v2 descriptor/source/model/event operands.
    if (version==2) {
        need(!section("cdf_desc").empty() && section("cdf_desc").size()%256==0 &&
             !section("cdf_quantile").empty() && section("cdf_quantile").size()%8==0 &&
             section("cdf_fragment_eV").size()==section("cdf_quantile").size(),"CDF sections");
    } else need(!sections.count("cdf_desc") && !sections.count("cdf_quantile") &&
                !sections.count("cdf_fragment_eV"),"legacy table contains CDF sections");
    const auto meta=section("meta");need(meta.size()==16,"meta size");
    const auto count=read_le(meta,0,4),groups=read_le(meta,4,4),procs=read_le(meta,8,4),regions=read_le(meta,12,4);
    need(count>0 && count<=64 && groups<=count*count,"species/group count");
    const auto sp=section("species"),gr=section("groups"),lookup=section("group_of");
    need(sp.size()==64*count && gr.size()==64*groups && lookup.size()==4*count*count,"registry section size");
    const auto strings=section("strings");std::vector<std::string> lines;std::string line;
    for (unsigned char c:strings) {if (c=='\n') {lines.push_back(line);line.clear();}else {need(c!=0,"NUL name");line.push_back(char(c));}}
    need(line.empty() && lines.size()==count+2*procs+regions,"registry names");
    names_.assign(lines.begin(),lines.begin()+count);
    std::set<std::string> unique;
    std::string registry="table-species-registry-v1";append_u64(registry,count);
    for (std::size_t i=0;i<count;++i) {
        need(!names_[i].empty() && names_[i].find_first_of(" \t\r\n")==std::string::npos && unique.insert(names_[i]).second,"species name");
        append_text(registry,names_[i]);registry.append(reinterpret_cast<const char*>(sp.data()+64*i),64);
    }
    registry_sha_=sha256(registry);
    const auto bo=section("bin_off"),br=section("bin_region"),bounds=section("bin_row");
    need(bo.size()==4*groups && br.size()==bounds.size() && br.size()%4==0,"coverage section size");
    struct GroupInfo {Pair pair; PairStatus status; bool resolved;};std::vector<GroupInfo> info;
    std::map<Pair,std::uint32_t> index;
    for (std::uint32_t g=0;g<groups;++g) {
        const Pair pair{std::uint32_t(read_le(gr,64*g,4)),std::uint32_t(read_le(gr,64*g+4,4))};
        need(pair.projectile<count && pair.target<count && index.emplace(pair,g).second,"group registry");
        const auto n=read_le(gr,64*g+12,4),off=read_le(bo,4*g,4);
        need(n>=2 && off<=br.size()/4 && n+1<=br.size()/4-off,"coverage range");
        bool all_unknown=true,any_resolved=false;
        for (std::uint64_t j=0;j<=n;++j) {
            const auto region=read_le(br,4*(off+j),4),bound=read_le(bounds,4*(off+j),4);
            need(region==none || region<regions,"coverage region ID");
            any_resolved=any_resolved || region==none;
            all_unknown=all_unknown && region!=none && bound==none;
        }
        info.push_back({pair,all_unknown?PairStatus::entirely_unbounded:PairStatus::supported,any_resolved});
    }
    for (std::uint32_t p=0;p<count;++p) for (std::uint32_t t=0;t<count;++t) {
        const auto it=index.find({p,t}); const auto expected_index=it==index.end()?none:it->second;
        need(read_le(lookup,4*(std::uint64_t(p)*count+t),4)==expected_index,"group lookup differs from registry");
    }
    std::sort(info.begin(),info.end(),[](const auto& a,const auto& b){return a.pair<b.pair;});
    for (const auto& entry:info) {pairs_.push_back(entry.pair);status_.push_back(entry.status);resolved_.push_back(entry.resolved);}
}
std::uint32_t Catalogue::species_id(const std::string& name) const {
    const auto it=std::find(names_.begin(),names_.end(),name);need(it!=names_.end(),"required species absent: "+name);
    return std::uint32_t(it-names_.begin());
}
PairStatus Catalogue::status(Pair pair) const {
    const auto it=std::lower_bound(pairs_.begin(),pairs_.end(),pair);
    return it==pairs_.end() || *it!=pair ? PairStatus::absent : status_[std::size_t(it-pairs_.begin())];
}
bool Catalogue::has_resolved_interval(Pair pair) const {
    const auto it=std::lower_bound(pairs_.begin(),pairs_.end(),pair);
    return it!=pairs_.end() && *it==pair && resolved_[std::size_t(it-pairs_.begin())];
}
std::string minimum_contract(const Catalogue& table,const std::string& native_registry,const std::string& configuration) {
    digest_valid(native_registry);digest_valid(configuration);
    std::ostringstream s;s<<"MINIMUM-WATER-ENCOUNTERS 1\nmodel_scope full_water_breakdown_candidate\n"
      <<"blob_sha256 "<<table.blob_sha256()<<"\nregistry_sha256 "<<table.registry_sha256()<<"\n"
      <<"native_registry_sha256 "<<native_registry<<"\nconfiguration_sha256 "<<configuration<<"\n"
      <<"census all_registered_species_before_target_filter\nzero_inventory inactive_unresolved\n"
      <<"existing_region_admission required\nnetwork_complete false\nphysical_model_validated false\n";
    const auto pairs=required(table);s<<"required_pairs "<<pairs.size()<<'\n';
    for (const auto pair:pairs) s<<"pair "<<table.species()[pair.projectile]<<' '<<table.species()[pair.target]<<'\n';
    return s.str();
}
Guard::Guard(Catalogue table,const std::string& contract,const std::string& expected):table_(std::move(table)) {
    digest_valid(expected);contract_sha_=sha256(contract);need(contract_sha_==expected,"scope contract identity differs");
    std::istringstream s(contract);std::string line;
    while(std::getline(s,line)) {
        const std::string n="native_registry_sha256 ",c="configuration_sha256 ";
        if(line.starts_with(n)) {need(native_registry_sha_.empty(),"duplicate native registry");native_registry_sha_=line.substr(n.size());}
        if(line.starts_with(c)) {need(configuration_sha_.empty(),"duplicate configuration");configuration_sha_=line.substr(c.size());}
    }
    need(contract==minimum_contract(table_,native_registry_sha_,configuration_sha_),"scope contract differs from the explicit minimum");
    required_=required(table_);
}
Admission Guard::check(std::span<const InventoryRow> rows,const StageIdentity& id) const {
    digest_valid(id.input_state_sha256);digest_valid(id.expected_census_sha256);
    need(id.native_registry_sha256==native_registry_sha_ && id.configuration_sha256==configuration_sha_,"native registry/configuration identity differs");
    need(id.native_row_count==rows.size() && id.all_registered_species_before_target_filter,"incomplete pre-filter census");
    need((id.stage==1 || id.stage==2) && std::isfinite(id.dt_s) && id.dt_s>0 &&
         std::isfinite(id.collision_cell_volume_m3) && id.collision_cell_volume_m3>0,"invalid stage interval/cell volume");
    Admission result;result.enforced=true;result.existing_region_admission_required=true;result.scope_sha256=contract_sha_;
    result.census_sha256=census_digest(rows);need(result.census_sha256==id.expected_census_sha256,"census identity differs");
    result.stage_binding_sha256=bind_stage(rows,id);
    digest_valid(id.expected_stage_binding_sha256);
    need(result.stage_binding_sha256==id.expected_stage_binding_sha256,"stage input binding differs");
    std::map<std::uint64_t,std::uint64_t> positive;
    for (const auto& row:rows) {
        need(row.species<table_.species().size(),"unregistered census species");
        if (row.physical_count>0) positive[row.collision_cell]|=std::uint64_t(1)<<row.species;
    }
    for (const auto& [cell,present]:positive) for (const auto pair:required_) {
        if ((present&(std::uint64_t(1)<<pair.projectile)) && (present&(std::uint64_t(1)<<pair.target))) {
            const auto status=table_.status(pair);
            if (status!=PairStatus::supported) result.gaps.push_back({cell,pair,status});
        }
    }
    result.admitted=result.gaps.empty();
    return result;
}
std::string Guard::bind_stage(std::span<const InventoryRow> rows,const StageIdentity& id) const {
    std::string binding="minimum-water-stage-admission-v1";
    for (const auto& value:{contract_sha_,table_.blob_sha256(),table_.registry_sha256(),id.native_registry_sha256,
                           id.configuration_sha256,id.input_state_sha256,census_digest(rows)}) append_text(binding,value);
    append_u64(binding,id.native_row_count);append_u64(binding,id.step);append_u64(binding,id.stage);
    append_u64(binding,std::bit_cast<std::uint64_t>(id.dt_s));append_u64(binding,std::bit_cast<std::uint64_t>(id.collision_cell_volume_m3));
    append_u64(binding,id.all_registered_species_before_target_filter?1:0);
    return sha256(binding);
}
Admission admit(Mode mode,const Guard* guard,std::span<const InventoryRow> rows,const StageIdentity& identity) {
    switch(mode) {
    case Mode::legacy_parent: case Mode::legacy_diagnostic: return Admission{};
    case Mode::full_water_breakdown_candidate: need(guard!=nullptr,"strict mode requires pinned scope guard");return guard->check(rows,identity);
    }
    throw std::runtime_error("water scope: unknown mode");
}
} // namespace water_scope
