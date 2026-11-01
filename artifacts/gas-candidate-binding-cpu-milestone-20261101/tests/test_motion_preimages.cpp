#include "pf_gas_v1/gas_owner.hpp"
#include <array>
#include <iostream>
#include <limits>
using namespace pf_gas_v1;
int main(int argc,char **argv) {
  const bool baseline=argc==2 && std::string(argv[1])=="baseline";
  Configuration c;c.species={{"gas",1,0,true}};c.grid.low={0,0,0};c.grid.high={2,1,1};
  c.grid.n={2,1,1};c.grid.volume=1;c.transport=Transport::carrier_shares;
  CensusRow a{};a.native_id=1;a.weight=1;a.position_m[0]=a.position_m[1]=a.position_m[2]=0.5;
  auto b=a;b.native_id=2;const std::array<CensusRow,2> before={a,b};int n=0;
  auto test=[&](auto f){bool rejected=false;try{f();}catch(const Refusal&){rejected=true;}
    if(rejected==baseline)throw std::runtime_error("preimage outcome mismatch at "+std::to_string(n));
    ++n;};
  try {
    const std::array<CensusRow,2> duplicate={a,a};
    test([&]{census(c,duplicate);});
    auto bad=a;bad.weight=0;test([&]{census(c,{&bad,1});});
    const std::array<CensusRow,3> duplicate_after={a,b,a};
    test([&]{movers(c,before,duplicate_after,{});});
    auto changed=before;changed[0].weight=2;test([&]{movers(c,before,changed,{});});
    const std::uint64_t id=1;test([&]{movers(c,before,before,{&id,1});});
    const std::uint64_t unknown=3;test([&]{movers(c,before,before,{&unknown,1});});
    test([&]{Exact exact;exact.add_product(std::numeric_limits<double>::infinity(),1);});
    Owner owner(c);owner.begin_step(1);owner.take(1,1,{},{});
    MoverRecord m{};m.native_id=1;m.weight=1;m.cell_count=1;m.to_cell=1;
    std::array<MoverRecord,2> mm={m,m};
    test([&]{owner.transport(1,mm);});
    mm[1].native_id=2;test([&]{owner.transport(1,mm);});
    std::cout<<(baseline?"REPRODUCED ":"REFUSED ")<<n<<" invalid motion/exact inputs\n";
  }catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
}
