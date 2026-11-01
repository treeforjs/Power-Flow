#pragma once
#include "pf_gas_v1/gas_owner.hpp"
#include <algorithm>
#include <stdexcept>

// Pure evidence predicate, NOT a publication capability. The caller must obtain
// this receipt from its authentic live common-motion token, retain that binding,
// and recheck the common native/gas/Hall/source publication capabilities.
namespace powerflow::conductor_correction_guard_v1 {
namespace g=::pf_gas_v1;
struct Evidence {std::size_t conductor_exits=0,exterior_exits=0;};
inline void require(bool ok,const char*why){if(!ok)throw std::invalid_argument(why);}
inline bool digest(const std::string&s){return s.size()==64&&std::all_of(s.begin(),s.end(),[](char c){return (c>='0'&&c<='9')||(c>='a'&&c<='f');});}
inline Evidence require_no_unrouted_conductor_correction(const g::TransportReceipt&receipt,
    std::uint64_t expected_step,const g::Configuration&configuration){
  require(configuration.transport==g::Transport::census_conservative,"conductor correction: census transport required");
  require(configuration.execution==g::Execution::host,"conductor correction: qualified host receipt required");
  require(receipt.accepted_step==expected_step,"conductor correction: wrong transport step");
  require(digest(receipt.inputs_sha256)&&digest(receipt.candidate_before_sha256)&&digest(receipt.candidate_after_sha256),"conductor correction: missing transport identities");
  require(configuration.max_transport_carriers>0&&configuration.max_transport_exit_bytes>=sizeof(g::ExactExitShare)&&
    receipt.exits.size()<=configuration.max_transport_carriers&&
    receipt.exits.size()<=configuration.max_transport_exit_bytes/sizeof(g::ExactExitShare),"conductor correction: bounded exit payload required");
  Evidence proof;std::uint64_t previous=0;
  for(const auto&exit:receipt.exits){
    require(exit.native_id>previous,"conductor correction: noncanonical or duplicate native ID");previous=exit.native_id;
    require(exit.species<configuration.species.size()&&configuration.species[exit.species].charge==0,"conductor correction: unknown or charged correction owner");
    require(exit.destination==g::ExitDestination::conductor||exit.destination==g::ExitDestination::exterior,"conductor correction: unknown destination");
    const bool conductor=exit.destination==g::ExitDestination::conductor;
    for(const auto&quantity:exit.digits){
      // Gas transport publishes base-2^32 normalized digits. The highest digit
      // is signed; all lower digits must be unsigned 32-bit canonical values.
      for(std::size_t k=0;k+1<quantity.size();++k)
        require(quantity[k]>=0&&quantity[k]<=0xffffffffll,"conductor correction: noncanonical exact digits");
      // Never round or aggregate. A one-digit subnormal-scale share matters,
      // and opposite signed shares at distinct surfaces cannot cancel here.
      if(conductor)require(std::all_of(quantity.begin(),quantity.end(),[](std::int64_t d){return d==0;}),
        "conductor correction: nonzero per-exit correction needs a qualified local thermal/momentum destination");
    }
    if(conductor)++proof.conductor_exits;else ++proof.exterior_exits;
  }
  return proof;
}
} // namespace powerflow::conductor_correction_guard_v1
