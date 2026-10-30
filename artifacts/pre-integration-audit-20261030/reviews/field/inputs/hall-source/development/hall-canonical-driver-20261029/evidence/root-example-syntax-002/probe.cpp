#include "native_parent_subflow.hpp"
#include "root_joint_prepare_example.hpp"
void probe(raikou::mfem_backend::PicArena&a,raikou::native_water::ParentSubflow&p,const axial::Mesh&g,const axial::Fields&f,axial::State&s,hall_driver::Thermal&t,const hall_driver::CurrentTrace&c,const std::string&h){hall_driver::root_two_exposure_example(a,p,g,f,s,t,c,1e-14,h,h,h,0,100000,4096);}
