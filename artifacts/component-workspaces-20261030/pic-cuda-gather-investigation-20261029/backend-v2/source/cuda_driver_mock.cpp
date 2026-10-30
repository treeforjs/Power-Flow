// CPU-only link-time CUDA Driver API mock. This target does NOT link libcuda.
// Pageable copies deliberately remain queued until an explicit context fence;
// destroying a stream deliberately does not drain it. These adversarial but
// API-permitted semantics catch the two reviewed ordering/lifetime defects.
#include "cuda_driver_mock.hpp"
#include "gather_push_core.hpp"
#include <cstring>
#include <cstdlib>
#include <deque>
#include <functional>
#include <set>
#include <memory>

namespace mock_cuda {
std::vector<std::string> trace;
std::map<std::string,int> calls;
namespace {
std::map<std::string,std::set<int>> faults;
std::vector<CUcontext> stack;
struct StreamState {CUcontext context;std::deque<std::function<void()>> work;};
std::map<CUstream,StreamState> streams;
std::map<CUcontext,std::deque<std::function<void()>>> defaults;
std::map<CUevent,CUstream> events;
std::set<void*> device_memory,host_memory;
std::size_t next_handle=100,unsafe=0;
int references=0;
bool inject(const char* name){trace.emplace_back(name);return faults[name].contains(++calls[name]);}
CUcontext current(){return stack.empty()?nullptr:stack.back();}
void drain(std::deque<std::function<void()>>& queue){while(!queue.empty()){auto job=std::move(queue.front());queue.pop_front();job();}}
void release_check(){if(queued())++unsafe;}
}
CUcontext primary(){return reinterpret_cast<CUcontext>(1);}
CUcontext foreign(){return reinterpret_cast<CUcontext>(2);}
void set_current(CUcontext context){if(stack.empty()){if(context)stack.push_back(context);}else if(context)stack.back()=context;else stack.pop_back();}
void fail(const std::string& name,int occurrence){faults[name].insert(occurrence);}
int primary_refs(){return references;}
std::size_t queued(){std::size_t n=0;for(auto& [_,s]:streams)n+=s.work.size();for(auto& [_,q]:defaults)n+=q.size();return n;}
std::size_t live_allocations(){return device_memory.size()+host_memory.size();}
std::size_t unsafe_releases(){return unsafe;}
void reset(bool with_primary) {
    // Test-owner recovery of intentionally quarantined mock resources. Pending
    // CPU jobs are drained before their fake allocations are released.
    for(auto& [_,q]:defaults)drain(q);
    for(auto& [_,s]:streams)drain(s.work);
    for(void* p:device_memory)std::free(p);
    for(void* p:host_memory)std::free(p);
    device_memory.clear();host_memory.clear();streams.clear();defaults.clear();events.clear();faults.clear();trace.clear();calls.clear();stack.clear();
    unsafe=0;references=with_primary?1:0;next_handle=100;if(with_primary)stack.push_back(primary());
}
}
using namespace mock_cuda;
extern "C" {
CUresult CUDAAPI cuInit(unsigned){return inject("init")?CUDA_ERROR_UNKNOWN:CUDA_SUCCESS;}
CUresult CUDAAPI cuGetErrorString(CUresult,const char** text){*text="mock injected failure";return CUDA_SUCCESS;}
CUresult CUDAAPI cuDeviceGet(CUdevice* d,int ordinal){*d=ordinal;return CUDA_SUCCESS;}
CUresult CUDAAPI cuCtxGetCurrent(CUcontext* c){if(inject("get_current"))return CUDA_ERROR_UNKNOWN;*c=current();return CUDA_SUCCESS;}
CUresult CUDAAPI cuDevicePrimaryCtxRetain(CUcontext* c,CUdevice){if(inject("retain_primary"))return CUDA_ERROR_UNKNOWN;++references;*c=primary();return CUDA_SUCCESS;}
CUresult CUDAAPI cuDevicePrimaryCtxRelease(CUdevice){if(inject("release_primary"))return CUDA_ERROR_UNKNOWN;--references;return CUDA_SUCCESS;}
CUresult CUDAAPI cuCtxPushCurrent(CUcontext c){if(inject("push_current"))return CUDA_ERROR_UNKNOWN;stack.push_back(c);return CUDA_SUCCESS;}
CUresult CUDAAPI cuCtxPopCurrent(CUcontext* c){if(inject("pop_current"))return CUDA_ERROR_UNKNOWN;if(stack.empty())return CUDA_ERROR_INVALID_CONTEXT;*c=stack.back();stack.pop_back();return CUDA_SUCCESS;}
CUresult CUDAAPI cuCtxSynchronize(){if(inject("context_sync"))return CUDA_ERROR_UNKNOWN;drain(defaults[current()]);for(auto& [_,s]:streams)if(s.context==current())drain(s.work);return CUDA_SUCCESS;}
CUresult CUDAAPI cuMemAlloc(CUdeviceptr* p,std::size_t n){if(inject("device_alloc"))return CUDA_ERROR_OUT_OF_MEMORY;void* m=std::calloc(1,n);if(!m)return CUDA_ERROR_OUT_OF_MEMORY;device_memory.insert(m);*p=reinterpret_cast<CUdeviceptr>(m);return CUDA_SUCCESS;}
CUresult CUDAAPI cuMemFree(CUdeviceptr p){release_check();if(inject("device_free"))return CUDA_ERROR_UNKNOWN;void* m=reinterpret_cast<void*>(p);device_memory.erase(m);std::free(m);return CUDA_SUCCESS;}
CUresult CUDAAPI cuMemHostAlloc(void** p,std::size_t n,unsigned){if(inject("host_alloc"))return CUDA_ERROR_OUT_OF_MEMORY;*p=std::calloc(1,n);if(!*p)return CUDA_ERROR_OUT_OF_MEMORY;host_memory.insert(*p);return CUDA_SUCCESS;}
CUresult CUDAAPI cuMemFreeHost(void* p){release_check();if(inject("host_free"))return CUDA_ERROR_UNKNOWN;host_memory.erase(p);std::free(p);return CUDA_SUCCESS;}
CUresult CUDAAPI cuMemcpyHtoD(CUdeviceptr to,const void* from,std::size_t n){
    auto staged=std::make_shared<std::vector<char>>(n);std::memcpy(staged->data(),from,n);
    defaults[current()].push_back([=]{std::memcpy(reinterpret_cast<void*>(to),staged->data(),n);});
    return inject("static_copy")?CUDA_ERROR_UNKNOWN:CUDA_SUCCESS;
}
CUresult CUDAAPI cuMemcpyHtoDAsync(CUdeviceptr to,const void* from,std::size_t n,CUstream stream){streams.at(stream).work.push_back([=]{std::memcpy(reinterpret_cast<void*>(to),from,n);});return inject("h2d")?CUDA_ERROR_UNKNOWN:CUDA_SUCCESS;}
CUresult CUDAAPI cuMemcpyDtoHAsync(void* to,CUdeviceptr from,std::size_t n,CUstream stream){streams.at(stream).work.push_back([=]{std::memcpy(to,reinterpret_cast<void*>(from),n);});return inject("d2h")?CUDA_ERROR_UNKNOWN:CUDA_SUCCESS;}
CUresult CUDAAPI cuModuleLoadDataEx(CUmodule* m,const void*,unsigned,CUjit_option*,void**){if(inject("module_load"))return CUDA_ERROR_UNKNOWN;*m=reinterpret_cast<CUmodule>(++next_handle);return CUDA_SUCCESS;}
CUresult CUDAAPI cuModuleGetFunction(CUfunction* f,CUmodule,const char*){if(inject("get_function"))return CUDA_ERROR_UNKNOWN;*f=reinterpret_cast<CUfunction>(99);return CUDA_SUCCESS;}
CUresult CUDAAPI cuModuleUnload(CUmodule){release_check();return inject("module_unload")?CUDA_ERROR_UNKNOWN:CUDA_SUCCESS;}
CUresult CUDAAPI cuStreamCreate(CUstream* s,unsigned){if(inject("stream_create"))return CUDA_ERROR_UNKNOWN;*s=reinterpret_cast<CUstream>(++next_handle);streams.emplace(*s,StreamState{current(),{}});return CUDA_SUCCESS;}
CUresult CUDAAPI cuStreamSynchronize(CUstream s){if(inject("stream_sync"))return CUDA_ERROR_UNKNOWN;drain(streams.at(s).work);return CUDA_SUCCESS;}
CUresult CUDAAPI cuStreamDestroy(CUstream s){release_check();if(inject("stream_destroy"))return CUDA_ERROR_UNKNOWN;streams.erase(s);return CUDA_SUCCESS;}
CUresult CUDAAPI cuEventCreate(CUevent* e,unsigned){if(inject("event_create"))return CUDA_ERROR_UNKNOWN;*e=reinterpret_cast<CUevent>(++next_handle);events[*e]=nullptr;return CUDA_SUCCESS;}
CUresult CUDAAPI cuEventRecord(CUevent e,CUstream s){events[e]=s;return inject("event_record")?CUDA_ERROR_UNKNOWN:CUDA_SUCCESS;}
CUresult CUDAAPI cuEventSynchronize(CUevent e){if(inject("event_sync"))return CUDA_ERROR_UNKNOWN;drain(streams.at(events.at(e)).work);return CUDA_SUCCESS;}
CUresult CUDAAPI cuEventElapsedTime(float* ms,CUevent,CUevent){if(inject("event_elapsed"))return CUDA_ERROR_UNKNOWN;*ms=.1f;return CUDA_SUCCESS;}
CUresult CUDAAPI cuEventDestroy(CUevent e){release_check();if(inject("event_destroy"))return CUDA_ERROR_UNKNOWN;events.erase(e);return CUDA_SUCCESS;}
CUresult CUDAAPI cuLaunchKernel(CUfunction,unsigned,unsigned,unsigned,unsigned,unsigned,unsigned,unsigned,CUstream s,void** args,void**) {
    const auto table=*static_cast<gather_push::Tables*>(args[0]);
    const auto inputs=*static_cast<CUdeviceptr*>(args[1]),outputs=*static_cast<CUdeviceptr*>(args[2]);const auto n=*static_cast<unsigned*>(args[3]);
    streams.at(s).work.push_back([=]{auto* in=reinterpret_cast<const gather_push::Inputs*>(inputs);auto* out=reinterpret_cast<gather_push::Outputs*>(outputs);for(unsigned i=0;i<n;++i)out[i]=gather_push::evaluate(table,in[i]);if(inject("kernel_status"))out[0].status=gather_push::invalid_input;});
    return inject("launch")?CUDA_ERROR_UNKNOWN:CUDA_SUCCESS;
}
// These APIs must never be used by the backend, even during injected failures.
CUresult CUDAAPI cuCtxDestroy(CUcontext){inject("FORBIDDEN_context_destroy");return CUDA_ERROR_UNKNOWN;}
CUresult CUDAAPI cuDevicePrimaryCtxReset(CUdevice){inject("FORBIDDEN_primary_reset");return CUDA_ERROR_UNKNOWN;}
}
