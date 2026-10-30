#include "cuda_force_runtime.hpp"
#include "runtime_guards.hpp"
#include <cuda.h>
#include <array>
#include <chrono>
#include <thread>
#include <cstdio>

namespace pic_batch {
namespace gp=gather_push;
namespace io=snapshot_io;
using namespace raikou::mfem_backend;
using Clock=std::chrono::steady_clock;
namespace {
void cuda_ok(CUresult result,const char* operation) {
    if(result!=CUDA_SUCCESS) {
        const char* detail=nullptr;cuGetErrorString(result,&detail);
        throw std::runtime_error(std::string(operation)+": "+(detail?detail:"CUDA error"));
    }
}
double elapsed(Clock::time_point begin,Clock::time_point end) {return std::chrono::duration<double,std::milli>(end-begin).count();}
void cleanup_issue(CleanupReport& report,const char* operation,CUresult result) noexcept {
    if(result==CUDA_SUCCESS)return;
    if(report.issue_count<report.issues.size())report.issues[report.issue_count]={operation,int(result)};
    ++report.issue_count;
}
struct Context {
    CUcontext value=nullptr;CUdevice device=0;
    bool retained=false,initial_push=false,cleanup_push=false,closed=false;
    bool current_was_replaced=false;
    ContextPolicy policy;
    CleanupReport& report;CleanupReport* sink;
    std::size_t logged_issues=0;bool logged_retained=false;
    std::thread::id thread=std::this_thread::get_id();
    Context(BackendOptions options,CleanupReport& r):policy(options.context),report(r),sink(options.cleanup_report) {
        try {
            cuda_ok(cuInit(0),"cuInit");CUcontext current=nullptr;
            cuda_ok(cuCtxGetCurrent(&current),"cuCtxGetCurrent");
            if(policy==ContextPolicy::RequireCurrentSharedPrimary)
                io::require(current!=nullptr,"cuDSS-first shared-primary current context is required; no private fallback");
            else io::require(current==nullptr,"StandalonePrimary requires no current context; integration must borrow cuDSS first");
            cuda_ok(cuDeviceGet(&device,0),"cuDeviceGet");
            cuda_ok(cuDevicePrimaryCtxRetain(&value,device),"cuDevicePrimaryCtxRetain");retained=true;
            if(policy==ContextPolicy::RequireCurrentSharedPrimary) {
                io::require(current==value,"Current context is not device0 shared primary; no context switch permitted");
            } else {
                cuda_ok(cuCtxPushCurrent(value),"cuCtxPushCurrent standalone primary");initial_push=true;
            }
        }catch(...){close();throw;}
    }
    void publish() noexcept {
        if(report.resources_retained)report.resources_released=false;
        if(sink){*sink=report;return;}
        for(std::size_t i=logged_issues;i<report.issue_count && i<report.issues.size();++i)
            std::fprintf(stderr,"PIC CUDA cleanup: %s (status%d)\n",report.issues[i].operation,report.issues[i].code);
        logged_issues=report.issue_count;
        if(report.resources_retained && !logged_retained) {
            std::fprintf(stderr,"PIC CUDA cleanup retained resources because safe release was not established; abort this candidate process.\n");logged_retained=true;
        }
    }
    void current() const {CUcontext now=nullptr;cuda_ok(cuCtxGetCurrent(&now),"cuCtxGetCurrent");io::require(now==value,"Batched-force CUDA context changed");}
    bool enter_cleanup() noexcept {
        if(closed || !value)return false;
        if(std::this_thread::get_id()!=thread) {cleanup_issue(report,"cleanup owner thread",CUDA_ERROR_INVALID_CONTEXT);return false;}
        CUcontext now=nullptr;const auto got=cuCtxGetCurrent(&now);cleanup_issue(report,"cleanup cuCtxGetCurrent",got);
        if(got==CUDA_SUCCESS && now==value){report.context_established=true;return true;}
        current_was_replaced=true;
        const auto pushed=cuCtxPushCurrent(value);cleanup_issue(report,"cleanup cuCtxPushCurrent",pushed);
        if(pushed!=CUDA_SUCCESS)return false;
        cleanup_push=true;report.context_established=true;return true;
    }
    void close() noexcept {
        if(closed)return;
        closed=true;
        if(std::this_thread::get_id()!=thread) {
            cleanup_issue(report,"close owner thread",CUDA_ERROR_INVALID_CONTEXT);
            report.resources_retained=true;publish();return;
        }
        bool restored=true;
        if(cleanup_push) {
            CUcontext popped=nullptr;const auto result=cuCtxPopCurrent(&popped);
            cleanup_issue(report,"cleanup cuCtxPopCurrent",result);restored=result==CUDA_SUCCESS && popped==value;
            if(!restored)report.resources_retained=true;
        }
        if(initial_push) {
            if(current_was_replaced) {
                // Do not pop a context another component placed on the stack.
                cleanup_issue(report,"standalone caller changed context",CUDA_ERROR_INVALID_CONTEXT);
                report.resources_retained=true;
            } else {
                CUcontext popped=nullptr;const auto result=cuCtxPopCurrent(&popped);
                cleanup_issue(report,"standalone cuCtxPopCurrent",result);restored=result==CUDA_SUCCESS && popped==value;
                if(!restored)report.resources_retained=true;
            }
        }
        report.caller_context_restored=restored;
        if(retained && !report.resources_retained) {
            const auto result=cuDevicePrimaryCtxRelease(device);
            cleanup_issue(report,"cuDevicePrimaryCtxRelease own reference",result);
            report.primary_reference_released=result==CUDA_SUCCESS;
            if(result!=CUDA_SUCCESS)report.resources_retained=true;
        }
        publish();
    }
    ~Context(){close();}
};
// One resource owner centralizes completion and release order. Its destructor
// also runs when Impl construction throws, unlike Impl::~Impl itself.
struct Resources {
    struct Device {CUdeviceptr value=0;std::size_t bytes=0;};
    struct Host {void* value=nullptr;};
    std::array<Device,12> device;
    Host inputs,phi,outputs;
    CUmodule module=nullptr;CUstream stream=nullptr;std::array<CUevent,4> events{};
    Context& context;CleanupReport& report;
    bool closed=false,static_inflight=false,stream_inflight=false;
    Resources(Context& c,CleanupReport& r):context(c),report(r){}
    void allocate_device(std::size_t slot,std::size_t bytes) {
        cuda_ok(cuMemAlloc(&device[slot].value,bytes),"cuMemAlloc");device[slot].bytes=bytes;
    }
    void allocate_host(Host& host,std::size_t bytes) {cuda_ok(cuMemHostAlloc(&host.value,bytes,0),"cuMemHostAlloc");}
    template<class Handle,class Release> void release(Handle& handle,Release function,const char* operation) noexcept {
        if(!handle)return;
        const auto result=function(handle);cleanup_issue(report,operation,result);
        if(result!=CUDA_SUCCESS)report.resources_retained=true;
        handle={};
    }
    void close() noexcept {
        if(closed)return;
        closed=true;report.attempted=true;
        if(!context.enter_cleanup()){report.resources_retained=true;context.publish();return;}
        bool drained=true;
        // Stream destruction is not a completion boundary. Always drain once
        // at close, before releasing any module, pinned memory or device data.
        if(stream) {
            const auto result=cuStreamSynchronize(stream);cleanup_issue(report,"close cuStreamSynchronize",result);
            drained=result==CUDA_SUCCESS;
        }
        if(!drained || static_inflight) {
            report.used_context_fallback=true;
            const auto result=cuCtxSynchronize();cleanup_issue(report,"close cuCtxSynchronize",result);
            drained=result==CUDA_SUCCESS;
        }
        if(!drained) {
            // Context state cannot establish completion. Quarantine allocations
            // and the retained reference rather than risk use-after-free. The
            // isolated caller must abort; no borrowed context is reset/destroyed.
            report.resources_retained=true;context.publish();return;
        }
        report.stream_drained=true;static_inflight=false;stream_inflight=false;
        for(auto& event:events)release(event,cuEventDestroy,"cuEventDestroy");
        release(stream,cuStreamDestroy,"cuStreamDestroy after drain");
        release(module,cuModuleUnload,"cuModuleUnload after drain");
        release(outputs.value,cuMemFreeHost,"cuMemFreeHost outputs after drain");
        release(phi.value,cuMemFreeHost,"cuMemFreeHost phi after drain");
        release(inputs.value,cuMemFreeHost,"cuMemFreeHost inputs after drain");
        for(auto& buffer:device)release(buffer.value,cuMemFree,"cuMemFree after drain");
        report.resources_released=!report.resources_retained;context.publish();
    }
    ~Resources(){close();}
};
enum Slot {SelfCells,NativeCells,Psi,Az0,Az1,Rate0,Rate1,Xs,Ys,Phi,Inputs,Outputs,Slots};
std::array<std::size_t,Slots> allocation_sizes(const RuntimeSnapshot& s,std::size_t capacity) {
    return {s.self_cells.size()*sizeof(gp::Cell),s.native_cells.size()*sizeof(gp::Cell),
        s.psi.size()*8,s.az0.size()*8,s.az1.size()*8,s.rate0.size()*8,s.rate1.size()*8,s.xs.size()*8,s.ys.size()*8,
        s.phi.size()*8,capacity*sizeof(gp::Inputs),capacity*sizeof(gp::Outputs)};
}
}
struct CudaForceRuntime::Impl {
    CleanupReport cleanup;
    Context context;
    Resources resources;
    CUfunction kernel=nullptr;
    std::thread::id thread=std::this_thread::get_id();
    Parameters static_parameters;
    std::size_t capacity,self_nodes,native_nodes,total_bytes=0,pending_count=0;
    std::vector<GapParticle> pending_inputs;
    bool pending=false,poisoned=false;
    BatchTiming last;
    Impl(const std::filesystem::path& path,const RuntimeSnapshot& s,std::size_t cap,BackendOptions options)
        :context(options,cleanup),resources(context,cleanup),static_parameters(s.parameters),capacity(cap),self_nodes(s.phi.size()),native_nodes(s.az0.size()),pending_inputs(cap) {
        const auto sizes=allocation_sizes(s,capacity);
        for(auto bytes:sizes)total_bytes+=bytes;
        io::require(total_bytes<=128*1024*1024,"Batched-force128MiB explicit device allocation guard exceeded");
        for(std::size_t i=0;i<resources.device.size();++i)resources.allocate_device(i,sizes[i]);
        resources.allocate_host(resources.inputs,capacity*sizeof(gp::Inputs));resources.allocate_host(resources.phi,self_nodes*8);resources.allocate_host(resources.outputs,capacity*sizeof(gp::Outputs));
        auto ptx=io::read<char>(path);ptx.push_back('\0');cuda_ok(cuModuleLoadDataEx(&resources.module,ptx.data(),0,nullptr,nullptr),"cuModuleLoadDataEx");
        cuda_ok(cuModuleGetFunction(&kernel,resources.module,"gather_push_kernel"),"cuModuleGetFunction");
        cuda_ok(cuStreamCreate(&resources.stream,CU_STREAM_NON_BLOCKING),"cuStreamCreate");
        for(auto& event:resources.events)cuda_ok(cuEventCreate(&event,CU_EVENT_BLOCKING_SYNC),"cuEventCreate");
    }
    ~Impl(){close();}
    CleanupReport close() noexcept {resources.close();context.close();return cleanup;}
    void ready() const {io::require(!resources.closed,"Batched-force backend is closed");io::require(!poisoned,"Batched-force backend is poisoned; reconstruct it");io::require(std::this_thread::get_id()==thread,"Batched-force owner thread changed");context.current();}
    void upload(Slot slot,const void* source,std::size_t bytes) {io::require(bytes==resources.device[slot].bytes,"Static snapshot size changed; reconstruct backend");cuda_ok(cuMemcpyHtoD(resources.device[slot].value,source,bytes),"cuMemcpyHtoD static snapshot");}
};
CudaForceRuntime::CudaForceRuntime(const std::filesystem::path& path,const RuntimeSnapshot& s,std::size_t capacity,BackendOptions options) {
    // Validate before the constructor above can initialize CUDA.
    s.validate();io::require(capacity>=1 && capacity<=2000000,"Invalid persistent particle capacity");
    std::size_t total=0;for(auto n:allocation_sizes(s,capacity))total+=n;
    io::require(total<=128*1024*1024,"Batched-force128MiB explicit allocation guard exceeded before CUDA initialization");
    p_=std::make_unique<Impl>(path,s,capacity,options);set_snapshot(s);
}
CudaForceRuntime::~CudaForceRuntime()=default;
void CudaForceRuntime::set_snapshot(const RuntimeSnapshot& s) {
    auto& b=*p_;b.ready();io::require(!b.pending,"Discard or commit pending batch before snapshot refresh");
    try {
        s.validate();io::require(s.phi.size()==b.self_nodes && s.az0.size()==b.native_nodes,"Snapshot nodal counts changed");
        b.resources.static_inflight=true;
        b.upload(SelfCells,s.self_cells.data(),s.self_cells.size()*sizeof(gp::Cell));b.upload(NativeCells,s.native_cells.data(),s.native_cells.size()*sizeof(gp::Cell));
        const std::array<std::pair<Slot,const std::vector<double>*>,7> fields={{{Psi,&s.psi},{Az0,&s.az0},{Az1,&s.az1},{Rate0,&s.rate0},{Rate1,&s.rate1},{Xs,&s.xs},{Ys,&s.ys}}};
        for(const auto& [slot,field]:fields)b.upload(slot,field->data(),field->size()*8);
        cuda_ok(cuCtxSynchronize(),"set_snapshot static publication fence");
        b.resources.static_inflight=false;
        b.static_parameters=s.parameters;
    }catch(...){b.poisoned=true;b.close();throw;}
}
std::span<const gp::Outputs> CudaForceRuntime::evaluate(std::span<const GapParticle> particles,
        std::span<const double> values,const Parameters& p) {
    const auto start=Clock::now();
    auto& b=*p_;b.ready();io::require(!b.pending,"Discard or commit prior batch first");
    try {
        const auto n=particles.size();
        io::require(n>0 && n<=b.capacity && values.size()==b.self_nodes,"Batched-force shape/capacity mismatch");
        const auto& fixed=b.static_parameters;
        io::require(p.magic==snapshot_magic && p.flags==1 && p.self_nx==fixed.self_nx && p.self_ny==fixed.self_ny &&
            p.native_nx==fixed.native_nx && p.native_ny==fixed.native_ny && p.width==fixed.width && p.height==fixed.height &&
            p.depth==fixed.depth && p.short_z==fixed.short_z,"Refresh static snapshot after geometry/layout changes");
        io::require(std::isfinite(p.fraction) && p.fraction>=0 && p.fraction<=1 && std::isfinite(p.voltage) && std::isfinite(p.dt) && p.dt>0,"Invalid midpoint parameters");
        auto* packed=static_cast<gp::Inputs*>(b.resources.inputs.value);
        for(std::size_t i=0;i<n;++i) {
            const auto& q=particles[i];
            io::require(gp::finite(io::vec(q.position_m)) && gp::finite(io::vec(q.proper_velocity_m_s)) && std::isfinite(q.charge_C) &&
                std::isfinite(q.mass_kg) && q.mass_kg>0 && std::isfinite(q.weight) && q.weight>0 &&
                std::isfinite(q.charge_C*q.weight) && std::isfinite(q.charge_C/q.mass_kg) && gp::finite(gp::norm2(io::vec(q.proper_velocity_m_s))),"Nonfinite/invalid input particle");
            io::require(q.position_m[0]>=0 && q.position_m[0]<=p.width && q.position_m[1]>=0 && q.position_m[1]<=p.height,"Particle outside supported domain");
            b.pending_inputs[i]=q;packed[i]={q.position_m[0],q.position_m[1],q.position_m[2],io::vec(q.proper_velocity_m_s),q.charge_C/q.mass_kg};
        }
        for(double v:values)io::require(std::isfinite(v),"Nonfinite recovered phi");
        std::memcpy(b.resources.phi.value,values.data(),values.size_bytes());const auto packed_at=Clock::now();
        gp::Tables table{{b.resources.device[SelfCells].value?reinterpret_cast<const gp::Cell*>(b.resources.device[SelfCells].value):nullptr,p.self_nx,p.self_ny},
            {reinterpret_cast<const gp::Cell*>(b.resources.device[NativeCells].value),p.native_nx,p.native_ny},
            reinterpret_cast<const double*>(b.resources.device[Phi].value),reinterpret_cast<const double*>(b.resources.device[Psi].value),
            reinterpret_cast<const double*>(b.resources.device[Az0].value),reinterpret_cast<const double*>(b.resources.device[Az1].value),
            reinterpret_cast<const double*>(b.resources.device[Rate0].value),reinterpret_cast<const double*>(b.resources.device[Rate1].value),
            reinterpret_cast<const double*>(b.resources.device[Xs].value),reinterpret_cast<const double*>(b.resources.device[Ys].value),
            p.width,p.height,p.depth,p.short_z,p.fraction,p.voltage,p.dt};
        auto din=b.resources.device[Inputs].value,dout=b.resources.device[Outputs].value;auto count=static_cast<unsigned>(n);
        void* arguments[]={&table,&din,&dout,&count};
        b.resources.stream_inflight=true;
        cuda_ok(cuEventRecord(b.resources.events[0],b.resources.stream),"event start");
        cuda_ok(cuMemcpyHtoDAsync(din,b.resources.inputs.value,n*sizeof(gp::Inputs),b.resources.stream),"H2D particles");
        cuda_ok(cuMemcpyHtoDAsync(b.resources.device[Phi].value,b.resources.phi.value,values.size_bytes(),b.resources.stream),"H2D phi");
        cuda_ok(cuEventRecord(b.resources.events[1],b.resources.stream),"event after H2D");
        cuda_ok(cuLaunchKernel(b.kernel,(count+127)/128,1,1,128,1,1,0,b.resources.stream,arguments,nullptr),"gather_push_kernel");
        cuda_ok(cuEventRecord(b.resources.events[2],b.resources.stream),"event after kernel");
        cuda_ok(cuMemcpyDtoHAsync(b.resources.outputs.value,dout,n*sizeof(gp::Outputs),b.resources.stream),"D2H force results");
        cuda_ok(cuEventRecord(b.resources.events[3],b.resources.stream),"event after D2H");
        cuda_ok(cuEventSynchronize(b.resources.events[3]),"final batch synchronization");b.resources.stream_inflight=false;const auto synchronized=Clock::now();
        float stages[3];for(int i=0;i<3;++i)cuda_ok(cuEventElapsedTime(&stages[i],b.resources.events[i],b.resources.events[i+1]),"event elapsed");
        const std::span<const gp::Outputs> result{static_cast<const gp::Outputs*>(b.resources.outputs.value),n};
        validate_runtime_outputs(result,particles);
        const auto validated=Clock::now();b.last={elapsed(start,packed_at),stages[0],stages[1],stages[2],elapsed(start,synchronized),elapsed(synchronized,validated),elapsed(start,validated)};
        b.pending_count=n;b.pending=true;return result;
    }catch(...){b.pending=false;b.poisoned=true;b.close();throw;}
}
void CudaForceRuntime::commit_velocities(std::span<GapParticle> particles) {
    auto& b=*p_;b.ready();io::require(b.pending && particles.size()==b.pending_count,"No validated batch to commit");
    for(std::size_t i=0;i<particles.size();++i)io::require(std::memcmp(&particles[i],&b.pending_inputs[i],sizeof(GapParticle))==0,"Particle state changed before transactional commit");
    const auto* result=static_cast<const gp::Outputs*>(b.resources.outputs.value);
    for(std::size_t i=0;i<particles.size();++i)particles[i].proper_velocity_m_s=io::vec(result[i].u);
    b.pending=false;
}
void CudaForceRuntime::discard(){p_->ready();p_->pending=false;}
CleanupReport CudaForceRuntime::close() noexcept{return p_->close();}
const BatchTiming& CudaForceRuntime::timing() const{return p_->last;}
std::size_t CudaForceRuntime::explicit_device_bytes() const{return p_->total_bytes;}
}
