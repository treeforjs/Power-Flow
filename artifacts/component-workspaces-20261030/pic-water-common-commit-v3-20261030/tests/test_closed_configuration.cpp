#include "closed_parent_joint_step.hpp"
#include <iostream>
using namespace raikou::native_water;
int main(){
  std::size_t checks=0;
  auto require=[&](bool ok){++checks;if(!ok)throw std::runtime_error("closed configuration check failed");};
  const std::string parent(64,'1'),field(64,'2');
  const ClosedParentInputs inputs{std::string(64,'3'),std::string(64,'4'),4096,100000};
  const auto base=closed_parent_configuration_sha256(parent,field,inputs);
  const auto source=closed_parent_inactive_source_sha256(base);
  require(pc::digest(base)&&pc::digest(source)&&base!=source);
  require(base==closed_parent_configuration_sha256(parent,field,inputs));
  require(source==closed_parent_inactive_source_sha256(base));
  auto changed=[&](const std::string&p,const std::string&f,const ClosedParentInputs&i){
    const auto config=closed_parent_configuration_sha256(p,f,i);
    require(config!=base);require(closed_parent_inactive_source_sha256(config)!=source);
  };
  changed(std::string(64,'5'),field,inputs);changed(parent,std::string(64,'6'),inputs);
  auto in=inputs;in.species_sha256=std::string(64,'7');changed(parent,field,in);
  in=inputs;in.frame_sha256=std::string(64,'8');changed(parent,field,in);
  in=inputs;--in.maximum_native_records;changed(parent,field,in);
  in=inputs;--in.maximum_carriers;changed(parent,field,in);
  auto refuses=[&](auto fn){bool rejected=false;try{fn();}catch(const std::exception&){rejected=true;}require(rejected);};
  refuses([&]{(void)closed_parent_configuration_sha256("",field,inputs);});
  refuses([&]{(void)closed_parent_configuration_sha256(parent,std::string(64,'G'),inputs);});
  in=inputs;in.frame_sha256="";refuses([&]{(void)closed_parent_configuration_sha256(parent,field,in);});
  in=inputs;in.species_sha256="";refuses([&]{(void)closed_parent_configuration_sha256(parent,field,in);});
  in=inputs;in.maximum_native_records=0;refuses([&]{(void)closed_parent_configuration_sha256(parent,field,in);});
  in=inputs;in.maximum_carriers=0;refuses([&]{(void)closed_parent_configuration_sha256(parent,field,in);});
  refuses([&]{(void)closed_parent_inactive_source_sha256("unknown");});
  std::cout<<Json{{"status","passed"},{"checks",checks},
    {"scope","pure closed configuration identity and input refusal"},
    {"accepted_owners_constructed",0},{"full_joint_test",false},
    {"physical_time_advanced_s",0}}.dump()<<'\n';
}
