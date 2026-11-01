#include <iostream>
#include <iomanip>
#include "chem_event.h"
#include "operands.hpp"
int main(){
    std::cout<<std::setprecision(17);
    for(double E:energies){
        double total=0;
        for(const auto& o:orbitals)if(E>o[0]){
            const auto s=pfc::sdcs_shape(o[0],o[1],o[2],E,5.29177211e-11,13.605,510998.95069,0);
            total+=pfc::sdcs_total(s);
        }
        std::cout<<E<<' '<<total<<'\n';
    }
}
