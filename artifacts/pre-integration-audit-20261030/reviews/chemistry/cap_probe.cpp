#include "generic_water_bridge.hpp"
#include <fstream>
#include <iostream>
namespace g=raikou::generic_water;namespace pc=raikou::pic_collisions;using J=g::Json;
std::vector<std::string> split(const std::string&s){std::vector<std::string> out;std::istringstream f(s);std::string x;while(std::getline(f,x,','))out.push_back(x);return out;}
int main(int argc,char**argv){try{
 if(argc!=5)return 2;
 auto data=std::make_shared<const raikou::kinetics::StrictCollisionDataset>(raikou::kinetics::StrictCollisionDataset::read(argv[1]));J binding;{std::ifstream f(argv[2]);f>>binding;}
 constexpr double amu=1.66053906660e-27,M=18.01056468403*amu,H=1.00782503223*amu;
 g::Bridge b(data,binding.at("dataset_sha256"),binding.at("kernel_manifest_sha256"),M,H,g::SitePolicy::PreserveReactantSites);
 J join=J::array();const auto channels=b.comparison_channels(g::Scope::CompareSupportedSubset);
 for(double E:{std::nextafter(1000.,0.),1000.,std::nextafter(1000.,INFINITY)}){J rates=J::array();double sum=0;for(const auto&ch:channels){const double value=ch.sigma_m2(E);sum+=value;rates.push_back({{"channel",ch.id},{"sigma_m2",value}});}join.push_back({{"energy_eV",E},{"subset_sum_m2",sum},{"rates",rates}});}
 J probes=J::array();for(std::size_t count:{257,1024,1025}){
  std::vector<g::PhysicalPacket> electrons;for(std::size_t i=0;i<count;++i){double E=20+.01*i,x=E*g::qe/(g::me*g::c*g::c),u=g::c*std::sqrt(x*(x+2));electrons.push_back({{"e/"+std::to_string(i),"electron",7,1,g::me,{u,0,0}},{.004,.004,.0125}});}
  const std::vector<g::PhysicalPacket> water={{{"water","H2O",7,1,M,{}},{.004,.004,.0125}}};
  J r={{"distinct_projectiles",count},{"coverage_cap",count==257?256:1024},{"counterfactual_static_fixture",true}};
  try{auto result=b.candidates(electrons,water,{{7,1e-12}},g::Scope::CompareSupportedSubset,2000,16000,5000,count==257?256:1024);r["admitted"]=true;r["coverage_records"]=result.coverage.size();r["candidate_rows"]=result.rows.size();}
  catch(const std::exception&e){r["admitted"]=false;r["reason"]=e.what();}probes.push_back(r);
 }
 std::ifstream f(argv[3]);std::string line;std::getline(f,line);const auto names=split(line);std::map<std::string,std::size_t> col;for(std::size_t i=0;i<names.size();++i)col[names[i]]=i;
 std::set<double> unique;std::size_t particles=0;long double weight=0;double maximum=0;
 pc::Packet electron,target;target.mass_kg=M;
 while(std::getline(f,line)){const auto row=split(line);if(row.size()!=names.size())throw std::runtime_error("CSV field count");electron.mass_kg=std::stod(row.at(col.at("mass_kg")));electron.proper_velocity_m_s={std::stod(row.at(col.at("ux_m_s"))),std::stod(row.at(col.at("uy_m_s"))),std::stod(row.at(col.at("uz_m_s")))};double E=pc::relative(electron,target).target_frame_energy_eV;unique.insert(E);maximum=std::max(maximum,E);weight+=std::stod(row.at(col.at("weight")));++particles;}
 J result={{"one_keV_join",join},{"coverage_cap_probes",probes},{"old_fixture_stationary_H2O_relative_energy_inventory",{{"particles",particles},{"distinct_binary64_energies",unique.size()},{"maximum_eV",maximum},{"physical_electrons",static_cast<double>(weight)},{"new_water_EEDF",false},{"velocity_columns","ux_m_s,uy_m_s,uz_m_s captured pre-kick operands"},{"purpose","representative broad-spectrum input only; no new history"}}},{"physical_time_advanced",false},{"native_or_field_solve_executed",false}};
 std::ofstream(argv[4])<<result.dump(2)<<'\n';std::cout<<result.dump(2)<<'\n';
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
