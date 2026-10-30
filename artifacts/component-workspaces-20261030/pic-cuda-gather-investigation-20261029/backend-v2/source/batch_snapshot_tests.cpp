#include "batch_snapshot.hpp"
#include "batch_compare.hpp"
#include <iostream>

int main(int argc,char** argv) {
 try {
    using namespace pic_batch;namespace io=snapshot_io;namespace gp=gather_push;
    io::require(argc==2,"Expected synthetic fixture for CPU schema tests");
    auto s=HostSnapshot::load(argv[1],false);
    auto reject=[&](const HostSnapshot& broken,bool require_real=false){bool rejected=false;try{broken.validate(require_real);}catch(const std::exception&){rejected=true;}io::require(rejected,"Invalid snapshot passed validation");};
    reject(s,true);
    auto bad=s;bad.parameters.magic=0;reject(bad);
    bad=s;bad.parameters.dt=0;reject(bad);
    bad=s;bad.native_cells[0].dof[0]=-1;reject(bad);
    bad=s;bad.psi[0]=std::numeric_limits<double>::quiet_NaN();reject(bad);
    bad=s;bad.particles[0].weight=0;reject(bad);
    bad=s;bad.inputs[0].u.v[0]+=1;reject(bad);
    std::vector<gp::Outputs> computed;computed.reserve(s.inputs.size());const auto tables=s.tables();
    for(const auto& input:s.inputs)computed.push_back(gp::evaluate(tables,input));
    const auto c=batch_compare::compare(s.reference,computed,s.particles);io::require(c.failures==0,"Shared CPU reference gate changed");
    computed[0].u.v[0]=std::numeric_limits<double>::quiet_NaN();io::require(batch_compare::compare(s.reference,computed,s.particles).failures>0,"NaN bypassed shared comparison");
    std::cout<<"{\"passed\":true,\"GPU_executed\":false,\"synthetic_rejected_by_default\":true,\"invalid_snapshots_rejected\":6}\n";
 }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
