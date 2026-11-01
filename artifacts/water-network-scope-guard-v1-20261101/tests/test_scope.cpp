#include "water_network_scope.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

using namespace water_scope;
int main(int argc,char** argv) {
    if(argc!=4) return 2;
    try {
        std::ifstream f(argv[1],std::ios::binary);
        std::vector<unsigned char> bytes((std::istreambuf_iterator<char>(f)),{});
        Catalogue table(bytes,argv[2]);
        const auto native=sha256("synthetic-native-registry"), config=sha256("synthetic-test-configuration");
        const auto contract=minimum_contract(table,native,config);
        Guard guard(table,contract,sha256(contract));
        int checks=0,failed=0;
        auto check=[&](bool ok,const std::string& name){++checks;failed+=!ok;std::cout<<(ok?"PASS ":"FAIL ")<<name<<'\n';};
        auto refuse=[&](auto action,const std::string& name){bool threw=false;try {action();}catch(const std::exception&){threw=true;}check(threw,name);};
        auto rows_for=[&](const std::string& p,const std::string& t){return std::vector<InventoryRow>{{1,7,table.species_id(p),1},{2,7,table.species_id(t),2}};};
        auto identity=[&](const std::vector<InventoryRow>& rows){
            StageIdentity id{sha256("synthetic-accepted-input-state"),native,config,census_digest(rows),"",rows.size(),9,1,1e-14,1e-12,true};
            id.expected_stage_binding_sha256=guard.bind_stage(rows,id);return id;
        };
        for(const auto* target:{"OH","O","H2"}) {
            auto rows=rows_for("e",target);const auto id=identity(rows);const auto before=census_digest(rows);
            const auto d=admit(Mode::full_water_breakdown_candidate,&guard,rows,id);
            check(d.enforced && !d.admitted && d.gaps.size()==1 && d.gaps[0].status==PairStatus::absent,
                  std::string("positive missing e+")+target+" refuses without energy assumptions");
            check(before==census_digest(rows),std::string("refusal preserves inventory e+")+target);
        }
        auto rows=rows_for("e","OH");rows[1].physical_count=0;
        auto d=guard.check(rows,identity(rows));
        check(d.admitted && !d.network_complete && !d.physical_model_validated,"zero target count is inactive, not resolved");
        rows.pop_back();check(guard.check(rows,identity(rows)).admitted,"absent target inventory inactive");
        rows=rows_for("e","OH");rows[1].collision_cell=8;
        check(guard.check(rows,identity(rows)).admitted,"species in distinct cells do not create an encounter");
        rows=rows_for("e","OH");rows[1].physical_count=std::numeric_limits<double>::denorm_min();
        check(!guard.check(rows,identity(rows)).admitted,"positive subnormal target count never rounds into absent inventory");
        rows=rows_for("e","H");
        check(table.has_resolved_interval({table.species_id("e"),table.species_id("H")}),"actual supported e+H table has resolved threshold interval");
        d=guard.check(rows,identity(rows));
        check(d.admitted && d.existing_region_admission_required && !d.network_complete && !d.physical_model_validated,
              "present supported pair preserved; regional guard still mandatory");
        for(const auto pair:{std::pair{"e","OH+"},std::pair{"H3O+","OH-"}}) {
            rows=rows_for(pair.first,pair.second);d=guard.check(rows,identity(rows));
            check(!d.admitted && d.gaps.size()==1 && d.gaps[0].status==PairStatus::entirely_unbounded,
                  std::string("existing unknown recombination remains refused: ")+pair.first+"+"+pair.second);
        }
        rows=rows_for("e","H2O");auto id=identity(rows);
        check(guard.check(rows,id).admitted,"supported e+H2O minimum gate preserved");
        std::reverse(rows.begin(),rows.end());
        check(guard.check(rows,id).admitted && census_digest(rows)==id.expected_census_sha256,"canonical inventory ignores row ordering");
        rows[0].physical_count+=1;
        refuse([&]{guard.check(rows,id);},"changed target count with old identity refuses");
        rows=rows_for("e","H2O");id=identity(rows);rows[1].species=table.species_id("OH");
        refuse([&]{guard.check(rows,id);},"changed species with old identity refuses");
        rows=rows_for("e","H2O");id=identity(rows);rows.pop_back();
        refuse([&]{guard.check(rows,id);},"filtered-out target refuses original row census");
        rows=rows_for("e","H2O");id=identity(rows);id.all_registered_species_before_target_filter=false;
        refuse([&]{guard.check(rows,id);},"table-filtered census is not accepted");
        for(int field=0;field<5;++field) {
            id=identity(rows);
            if(field==0)id.step++;
            if(field==1)id.stage=2;
            if(field==2)id.dt_s*=2;
            if(field==3)id.input_state_sha256=sha256("different-state");
            if(field==4)id.collision_cell_volume_m3*=2;
            refuse([&]{guard.check(rows,id);},"changed stage input refuses field "+std::to_string(field));
        }
        id=identity(rows);id.native_registry_sha256=sha256("changed registry");
        refuse([&]{guard.check(rows,id);},"native registry identity mismatch refuses");
        id=identity(rows);id.configuration_sha256=sha256("changed config");
        refuse([&]{guard.check(rows,id);},"configuration identity mismatch refuses");
        id=identity(rows);id.dt_s=0;id.expected_stage_binding_sha256=guard.bind_stage(rows,id);
        refuse([&]{guard.check(rows,id);},"zero-duration stage is invalid, not an admission waiver");
        auto wrong=bytes;wrong.back()^=1;
        refuse([&]{Catalogue t(wrong,argv[2]);},"changed actual table bytes refuse old table pin");
        refuse([&]{Guard g(table,contract+"extra\n",sha256(contract));},"contract byte tampering refuses");
        auto weak=contract;const auto at=weak.find("pair e OH\n");check(at!=std::string::npos,"required missing pair exists in explicit contract");
        weak.erase(at,std::string("pair e OH\n").size());
        refuse([&]{Guard g(table,weak,sha256(weak));},"self-consistent weakened contract still refuses minimum scope");
        for(double weight:{-1.0,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}) {
            rows=rows_for("e","H2O");rows[1].physical_count=weight;
            refuse([&]{census_digest(rows);},"invalid physical count refuses");
        }
        rows=rows_for("e","H2O");rows[1].native_id=rows[0].native_id;
        refuse([&]{census_digest(rows);},"duplicate native ID refuses");
        rows=rows_for("e","H2O");rows[1].species=999;
        refuse([&]{guard.check(rows,identity(rows));},"unregistered species refuses even with recomputed census");
        StageIdentity ignored{};std::vector<InventoryRow> bad{{0,0,999,std::numeric_limits<double>::quiet_NaN()}};
        for(auto mode:{Mode::legacy_parent,Mode::legacy_diagnostic}) {
            const auto bypass=admit(mode,nullptr,bad,ignored);
            check(!bypass.enforced && !bypass.admitted && !bypass.network_complete,"legacy mode untouched and supplies no full-network admission");
        }
        refuse([&]{admit(static_cast<Mode>(99),nullptr,{},ignored);},"unknown mode refuses");
        refuse([&]{admit(Mode::full_water_breakdown_candidate,nullptr,{},ignored);},"strict opt-in requires pinned guard");
        std::ofstream out(argv[3]);out<<contract;out.close();if(!out)throw std::runtime_error("contract output failed");
        std::cout<<"CHECKS "<<checks<<" FAILED "<<failed<<"\nTABLE "<<table.blob_sha256()<<"\nREGISTRY "<<table.registry_sha256()
                 <<"\nSYNTHETIC_CONTRACT "<<sha256(contract)<<"\n";
        return failed?1:0;
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 2;}
}
