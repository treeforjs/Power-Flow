#pragma once
#include <cudss.h>
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unistd.h>
namespace pic_cudss_tuning {
struct ConfigValue {
 cudssStatus_t status;size_t written;std::int32_t raw;
 bool verified() const{return status==CUDSS_STATUS_SUCCESS && written==sizeof(raw);}
 void json(std::ostream& out) const {if(verified())out<<raw;else out<<"null";}
};
class ConfigDiagnostics {
 inline static std::size_t next_id=0;
 std::ofstream events;
 static std::string hex(const void* value,size_t bytes) {
  const char* digits="0123456789abcdef";const auto* p=static_cast<const unsigned char*>(value);std::string out;
  for(size_t i=0;i<std::min(bytes,size_t(32));++i){out+=digits[p[i]>>4];out+=digits[p[i]&15];}return out;
 }
 void event(const char* operation,const char* parameter,const char* type,cudssStatus_t status,
            const void* buffer,size_t bytes,const size_t* written) {
  if(!enabled())return;
  events<<"{\"operation\":\""<<operation<<"\",\"parameter\":\""<<parameter<<"\",\"scalar_type\":\""<<type
        <<"\",\"return_code\":"<<int(status)<<",\"buffer_bytes\":"<<bytes<<",\"raw_size_written\":";
  if(written)events<<*written;else events<<"null";
  events<<",\"size_written_valid\":"<<(written && status==CUDSS_STATUS_SUCCESS?"true":"false")
        <<",\"raw_hex_bytes_recorded\":"<<std::min(bytes,size_t(32))<<",\"raw_bytes_hex\":\""<<hex(buffer,bytes)<<"\"}\n";
  events.flush();if(!events)throw std::runtime_error("Cannot write cuDSS config diagnostics");
 }
public:
 const std::size_t id=next_id++;
 std::string prefix;
 explicit ConfigDiagnostics(int rows) {
  const char* root=std::getenv("RAIKOU_PIC_CUDSS_AUDIT_PREFIX");if(!root)return;
  if(!*root || *root!='/')throw std::runtime_error("CUDSS_AUDIT_PREFIX must be absolute and nonempty");
  prefix=std::string(root)+"."+std::to_string(id);
  if(std::filesystem::exists(prefix+".config.jsonl") || std::filesystem::exists(prefix+".json"))throw std::runtime_error("Use a fresh CUDSS_AUDIT_PREFIX");
  events.open(prefix+".config.jsonl");
  events<<"{\"diagnostics_only\":true,\"solver_id\":"<<id<<",\"pid\":"<<::getpid()<<",\"matrix_rows\":"<<rows<<"}\n";
  events.flush();if(!events)throw std::runtime_error("Cannot initialize cuDSS config diagnostics");
 }
 bool enabled() const{return !prefix.empty();}
 template<class T>cudssStatus_t set(cudssConfig_t config,cudssConfigParam_t param,const char* name,T value) {
  static_assert(std::is_same_v<T,int> || std::is_same_v<T,cudssAlgType_t>);
  static_assert(sizeof(T)==4);
  const auto status=cudssConfigSet(config,param,&value,sizeof(value));
  event("ConfigSet",name,std::is_same_v<T,int>?"int":"cudssAlgType_t",status,&value,sizeof(value),nullptr);
  return status; // Caller retains the original mandatory setter error handling.
 }
 template<class T>ConfigValue get(cudssConfig_t config,cudssConfigParam_t param,const char* name) {
  static_assert(std::is_same_v<T,int> || std::is_same_v<T,cudssAlgType_t>);
  static_assert(sizeof(T)==4);
  T value;std::memset(&value,0xa5,sizeof(value));size_t written=0;
  const auto status=cudssConfigGet(config,param,&value,sizeof(value),&written);
  event("ConfigGet",name,std::is_same_v<T,int>?"int":"cudssAlgType_t",status,&value,sizeof(value),&written);
  std::int32_t raw;std::memcpy(&raw,&value,sizeof(raw));return {status,written,raw};
 }
 cudssStatus_t data_get(cudssHandle_t handle,cudssData_t data,cudssDataParam_t param,const char* name,
                       void* value,size_t bytes,size_t* written) {
  const auto status=cudssDataGet(handle,data,param,value,bytes,written);
  event("DataGet",name,"see_parameter_contract",status,value,bytes,written);return status;
 }
};
}
