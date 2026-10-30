#include "raikou/mfem_gap_pic.h"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <vector>
using namespace raikou::mfem_backend;

// Native bytes are valid only with the frozen header, compiler and ABI.
// This program constructs no GapPic, initializes no MPI and calls no solver.
static_assert(std::is_trivially_copyable_v<GapPicOptions>);
static_assert(std::is_trivially_copyable_v<GapPicDiagnostics>);
static_assert(std::is_trivially_copyable_v<GapParticle>);

int main(int argc, char** argv) {
 try {
    if(argc != 2) throw std::runtime_error("Expected one copied checkpoint directory");
    const std::filesystem::path directory(argv[1]);
    std::ifstream meta(directory/"restart.meta");
    std::string name; std::uint64_t step; int next_frame;
    meta >> name >> step >> next_frame;
    if(!meta || std::filesystem::path(name).filename() != name)
        throw std::runtime_error("Invalid copied restart metadata");
    std::ifstream input(directory/name, std::ios::binary);
    const auto get=[&](auto& value){input.read(reinterpret_cast<char*>(&value),sizeof(value));};
    std::uint64_t magic=0, next_id=0, count=0;
    GapPicOptions options; GapPicDiagnostics ledger;
    get(magic); get(options); get(ledger); get(next_id); get(count);
    if(!input || magic!=0x5241494b50494332ULL || count>2000000 || options.nx<1 ||
       options.ny<1 || options.nx>4096 || options.ny>4096 || options.depth_m<=0)
        throw std::runtime_error("Invalid native header or ABI mismatch");
    const auto expected=sizeof(magic)+sizeof(options)+sizeof(ledger)+sizeof(next_id)+sizeof(count)+count*sizeof(GapParticle);
    if(std::filesystem::file_size(directory/name)!=expected)
        throw std::runtime_error("Native checkpoint byte count does not match ABI");
    std::set<std::uint64_t> ids;
    long double weight=0, weight2=0, charge=0, kinetic=0, midweight=0;
    double minimum=std::numeric_limits<double>::infinity(),maximum=0;
    std::vector<double> cells(options.nx*options.ny,0), midcells(cells.size(),0), w2cells(cells.size(),0);
    std::vector<std::uint64_t> counts(cells.size(),0), midcounts(cells.size(),0);
    const double hx=options.width_m/options.nx, hy=options.gap_m/options.ny;
    for(std::uint64_t index=0;index<count;++index) {
        GapParticle p; get(p);
        if(!input || !ids.insert(p.id).second || p.id>=next_id || !std::isfinite(p.weight) || p.weight<=0 ||
           p.charge_C!=-pic_e || p.mass_kg!=pic_me)
            throw std::runtime_error("Invalid particle or unexpected species");
        for(double x:p.position_m) if(!std::isfinite(x)) throw std::runtime_error("Nonfinite position");
        for(double u:p.proper_velocity_m_s) if(!std::isfinite(u)) throw std::runtime_error("Nonfinite velocity");
        if(p.position_m[0]<0 || p.position_m[0]>options.width_m || p.position_m[1]<0 || p.position_m[1]>options.gap_m ||
           (options.finite_axial_extent && (p.position_m[2]<options.short_z_m || p.position_m[2]>options.feed_z_m)))
            throw std::runtime_error("Particle outside geometry");
        const int i=std::clamp(int(p.position_m[0]/hx),0,options.nx-1),j=std::clamp(int(p.position_m[1]/hy),0,options.ny-1);
        const auto cell=j*options.nx+i;
        weight+=p.weight;weight2+=static_cast<long double>(p.weight)*p.weight;charge+=static_cast<long double>(p.charge_C)*p.weight;
        minimum=std::min(minimum,p.weight);maximum=std::max(maximum,p.weight);
        cells[cell]+=p.weight;w2cells[cell]+=p.weight*p.weight;++counts[cell];
        if(std::abs(p.position_m[2])<.0005){midcells[cell]+=p.weight;++midcounts[cell];midweight+=p.weight;}
        const double u2=p.proper_velocity_m_s[0]*p.proper_velocity_m_s[0]+p.proper_velocity_m_s[1]*p.proper_velocity_m_s[1]+p.proper_velocity_m_s[2]*p.proper_velocity_m_s[2];
        kinetic+=static_cast<long double>(p.mass_kg*u2/(std::sqrt(1+u2/(pic_c*pic_c))+1))*p.weight;
    }
    if(!input || input.peek()!=std::char_traits<char>::eof()) throw std::runtime_error("Truncated or trailing native checkpoint bytes");
    const double volume=hx*hy*options.depth_m;
    auto maximum_count=*std::max_element(counts.begin(),counts.end());
    std::size_t occupied=0,midoccupied=0,singletons=0,midsingletons=0,underfive=0,midunderfive=0;
    std::uint64_t midmacro=0;
    for(std::size_t i=0;i<counts.size();++i){occupied+=counts[i]>0;midoccupied+=midcounts[i]>0;singletons+=counts[i]==1;midsingletons+=midcounts[i]==1;underfive+=counts[i]>0&&counts[i]<5;midunderfive+=midcounts[i]>0&&midcounts[i]<5;midmacro+=midcounts[i];}
    const long double defect=charge+ledger.lower_absorbed_charge_C+ledger.upper_absorbed_charge_C+ledger.short_absorbed_charge_C-
        ledger.lower_injected_charge_C-ledger.upper_injected_charge_C-ledger.authored_charge_C+ledger.escaped_charge_C;
    double pending=0;
    for(int f=0;f<4;++f){double span=0;meta>>span;std::string rest;std::getline(meta,rest);std::string integral,temperature;
        std::getline(meta,integral);std::getline(meta,temperature);
        std::istringstream row(integral);double value=0,total=0;std::size_t patches=0;
        while(row>>value){total+=value;++patches;}if(patches==0)throw std::runtime_error("Invalid reservoir metadata");
        pending+=2*total*(f<2?.005:.0002)/patches*options.depth_m/pic_e;
    }
    double elapsed=0, voltage=0;meta>>elapsed>>voltage;
    if(!meta || elapsed!=ledger.time_s || voltage!=options.upper_potential_V)
        throw std::runtime_error("Checkpoint and restart metadata disagree");
    std::cout<<std::setprecision(17)
        <<"{\n\"checkpoint\":\""<<name<<"\",\n\"step\":"<<step<<",\n\"elapsed_time_s\":"<<ledger.time_s
        <<",\n\"abi\":{\"options_bytes\":"<<sizeof(options)<<",\"diagnostics_bytes\":"<<sizeof(ledger)<<",\"particle_bytes\":"<<sizeof(GapParticle)<<"}"
        <<",\n\"macro_particles\":"<<count<<",\n\"active_physical_electrons\":"<<weight<<",\n\"active_charge_C\":"<<charge
        <<",\n\"weight_min\":"<<(count?minimum:0)<<",\n\"weight_mean\":"<<(count?weight/count:0)<<",\n\"weight_max\":"<<maximum
        <<",\n\"weight_below_one_possible\":true,\n\"weight_effective_sample_size\":"<<(weight2>0?weight*weight/weight2:0)
        <<",\n\"active_kinetic_energy_J\":"<<kinetic<<",\n\"ledger_charge_difference_C\":"<<charge-ledger.particle_charge_C
        <<",\n\"particle_charge_defect_C\":"<<defect<<",\n\"gauss_defect_C\":"<<ledger.gauss_defect_C
        <<",\n\"maximum_nodal_continuity_defect_C\":"<<ledger.maximum_nodal_continuity_defect_C
        <<",\n\"injected_physical_electrons\":"<<-(ledger.lower_injected_charge_C+ledger.upper_injected_charge_C)/pic_e
        <<",\n\"lower_absorbed_physical_electrons\":"<<-ledger.lower_absorbed_charge_C/pic_e
        <<",\n\"upper_absorbed_physical_electrons\":"<<-ledger.upper_absorbed_charge_C/pic_e
        <<",\n\"short_absorbed_physical_electrons\":"<<-ledger.short_absorbed_charge_C/pic_e
        <<",\n\"escaped_physical_electrons\":"<<-ledger.escaped_charge_C/pic_e
        <<",\n\"axial_escaped_physical_electrons_subset\":"<<-ledger.axial_escaped_charge_C/pic_e
        <<",\n\"source_reservoir_physical_electrons\":"<<pending
        <<",\n\"injected_kinetic_energy_J\":"<<ledger.injected_energy_J
        <<",\n\"short_absorbed_energy_J\":"<<ledger.short_absorbed_energy_J
        <<",\n\"field_energy_J\":"<<ledger.field_energy_J
        <<",\n\"all_xy_cells\":"<<cells.size()<<",\n\"occupied_xy_cells\":"<<occupied<<",\n\"singleton_xy_cells\":"<<singletons
        <<",\n\"occupied_xy_cells_below_five_macros\":"<<underfive
        <<",\n\"maximum_xy_macro_count\":"<<maximum_count
        <<",\n\"height_average_peak_density_cm3\":"<<*std::max_element(cells.begin(),cells.end())/volume/1e6
        <<",\n\"midplane_slab_macro_particles\":"<<midmacro<<",\n\"midplane_slab_physical_electrons\":"<<midweight
        <<",\n\"midplane_slab_occupied_xy_cells\":"<<midoccupied
        <<",\n\"midplane_slab_singleton_xy_cells\":"<<midsingletons
        <<",\n\"midplane_slab_occupied_xy_cells_below_five_macros\":"<<midunderfive
        <<",\n\"midplane_slab_maximum_xy_macro_count\":"<<*std::max_element(midcounts.begin(),midcounts.end())
        <<",\n\"midplane_slab_peak_density_cm3\":"<<*std::max_element(midcells.begin(),midcells.end())/(hx*hy*.001)/1e6
        <<",\n\"density_sampling_depth_m\":0.001,\n\"depth_m\":"<<options.depth_m<<"\n}\n";
    return 0;
 }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 2;}
}
