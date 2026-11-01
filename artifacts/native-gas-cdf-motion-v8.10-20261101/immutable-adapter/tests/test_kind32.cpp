#include "native_gas_stage.hpp"
#include <iostream>
int main(){unsigned checks=0;auto check=[&](bool x){++checks;if(!x)throw std::runtime_error("kind32 mapping");};
 check(native_gas_stage::internal_carrier(32,3)==pfg::kNone);
 check(native_gas_stage::internal_carrier(7,3)==2);
 check(native_gas_stage::internal_carrier(3,3)==2);
 for(unsigned c:{0u,1u,2u,4u,5u,6u}){bool refused=false;try{(void)native_gas_stage::internal_carrier(32,c);}catch(const std::invalid_argument&){refused=true;}check(refused);}
 bool refused=false;try{(void)native_gas_stage::internal_carrier(31,3);}catch(const std::invalid_argument&){refused=true;}check(refused);
 std::cout<<"kind32 mapping controls "<<checks<<"\n";
}
