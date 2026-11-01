#pragma once
#include "raikou/mfem_pic_full_time.h"
namespace native_gas_common {
// nlohmann's ordinary container destructor allocates a flattening stack. These
// proof trees have bounded structural depth, so drain recursively before normal
// destruction, including persistent-allocation-failure unwinding.
inline void drain_proof_json(raikou::mfem_backend::ArenaJson&j)noexcept {
 using J=raikou::mfem_backend::ArenaJson;
 if(j.is_object()){auto&o=j.get_ref<J::object_t&>();for(auto&x:o)drain_proof_json(x.second);o.clear();}
 else if(j.is_array()){auto&a=j.get_ref<J::array_t&>();for(auto&x:a)drain_proof_json(x);a.clear();}
}
inline void clone_proof_json(raikou::mfem_backend::ArenaJson&to,const raikou::mfem_backend::ArenaJson&from){
 using J=raikou::mfem_backend::ArenaJson;
 // Caller retains to under a drain guard. No nonempty JSON temporaries are
 // created between a potentially throwing allocation and their destruction.
 drain_proof_json(to);
 if(from.is_object()){to=J::object();for(const auto&[key,value]:from.items())clone_proof_json(to[key],value);}
 else if(from.is_array()){to=J::array();for(const auto&value:from){to.push_back(nullptr);clone_proof_json(to.back(),value);}}
 else to=from;
}
class GuardedProofJson {
 // Construct storage before insertion. Null-to-object operator[] promotion in
 // the frozen JSON changes the type before its allocation can fail.
 raikou::mfem_backend::ArenaJson value_=raikou::mfem_backend::ArenaJson::object();
public:
 GuardedProofJson()=default;
 explicit GuardedProofJson(raikou::mfem_backend::ArenaJson&&value)noexcept:value_(std::move(value)){}
 GuardedProofJson(const GuardedProofJson&)=delete;GuardedProofJson&operator=(const GuardedProofJson&)=delete;
 ~GuardedProofJson()noexcept{drain_proof_json(value_);}
 auto&get()noexcept{return value_;}const auto&get()const noexcept{return value_;}
 raikou::mfem_backend::ArenaJson release()noexcept{return std::move(value_);}
};
}
