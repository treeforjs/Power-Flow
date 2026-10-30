#include "cuda_batched_force.hpp"
#include <cuda.h>
#include <array>
#include <chrono>
#include <thread>

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
struct Context {
    CUcontext value=nullptr;bool owned=false,pushed=false;
    Context() {
        cuda_ok(cuInit(0),"cuInit");cuda_ok(cuCtxGetCurrent(&value),"cuCtxGetCurrent");
        if(!value) {CUdevice device;cuda_ok(cuDeviceGet(&device,0),"cuDeviceGet");cuda_ok(cuCtxCreate(&value,CU_CTX_SCHED_BLOCKING_SYNC,device),"cuCtxCreate");owned=true;}
    }
    void current() const {CUcontext now;cuda_ok(cuCtxGetCurrent(&now),"cuCtxGetCurrent");io::require(now==value,"Batched-force CUDA context changed");}
    void cleanup_context() noexcept {CUcontext now=nullptr;if(cuCtxGetCurrent(&now)==CUDA_SUCCESS && now!=value && cuCtxPushCurrent(value)==CUDA_SUCCESS)pushed=true;}
    ~Context() {if(pushed){CUcontext discarded;cuCtxPopCurrent(&discarded);}if(owned && value)cuCtxDestroy(value);}
};
struct DeviceBuffer {
    CUdeviceptr value=0;std::size_t bytes=0;
    void allocate(std::size_t n) {io::require(!value && n>0,"Invalid persistent device allocation");cuda_ok(cuMemAlloc(&value,n),"cuMemAlloc");bytes=n;}
    ~DeviceBuffer() {if(value)cuMemFree(value);}
};
struct HostBuffer {
    void* value=nullptr;
    void allocate(std::size_t n) {cuda_ok(cuMemHostAlloc(&value,n,0),"cuMemHostAlloc");}
    ~HostBuffer() {if(value)cuMemFreeHost(value);}
};
struct Module {CUmodule value=nullptr;~Module(){if(value)cuModuleUnload(value);}};
struct Stream {CUstream value=nullptr;~Stream(){if(value)cuStreamDestroy(value);}};
struct Event {CUevent value=nullptr;~Event(){if(value)cuEventDestroy(value);}};
enum Slot {SelfCells,NativeCells,Psi,Az0,Az1,Rate0,Rate1,Xs,Ys,Phi,Inputs,Outputs,Slots};
std::array<std::size_t,Slots> allocation_sizes(const HostSnapshot& s,std::size_t capacity) {
    return {s.self_cells.size()*sizeof(gp::Cell),s.native_cells.size()*sizeof(gp::Cell),
        s.psi.size()*8,s.az0.size()*8,s.az1.size()*8,s.rate0.size()*8,s.rate1.size()*8,s.xs.size()*8,s.ys.size()*8,
        s.phi.size()*8,capacity*sizeof(gp::Inputs),capacity*sizeof(gp::Outputs)};
}
}
struct CudaBatchedForce::Impl {
    Context context; // first constructed, last destroyed
    std::array<DeviceBuffer,Slots> device;
    HostBuffer inputs,phi,outputs;
    Module module;Stream stream;std::array<Event,4> events;
    CUfunction kernel=nullptr;
    std::thread::id thread=std::this_thread::get_id();
    Parameters static_parameters;
    std::size_t capacity,self_nodes,native_nodes,total_bytes=0,pending_count=0;
    std::vector<GapParticle> pending_inputs;
    bool pending=false,poisoned=false;
    BatchTiming last;
    batch_compare::Comparison comparison;
    Impl(const std::filesystem::path& path,const HostSnapshot& s,std::size_t cap)
        :static_parameters(s.parameters),capacity(cap),self_nodes(s.phi.size()),native_nodes(s.az0.size()),pending_inputs(cap) {
        const auto sizes=allocation_sizes(s,capacity);
        for(auto bytes:sizes)total_bytes+=bytes;
        io::require(total_bytes<=128*1024*1024,"Batched-force128MiB explicit device allocation guard exceeded");
        for(std::size_t i=0;i<device.size();++i)device[i].allocate(sizes[i]);
        inputs.allocate(capacity*sizeof(gp::Inputs));phi.allocate(self_nodes*8);outputs.allocate(capacity*sizeof(gp::Outputs));
        auto ptx=io::read<char>(path);ptx.push_back('\0');cuda_ok(cuModuleLoadDataEx(&module.value,ptx.data(),0,nullptr,nullptr),"cuModuleLoadDataEx");
        cuda_ok(cuModuleGetFunction(&kernel,module.value,"gather_push_kernel"),"cuModuleGetFunction");
        cuda_ok(cuStreamCreate(&stream.value,CU_STREAM_NON_BLOCKING),"cuStreamCreate");
        for(auto& e:events)cuda_ok(cuEventCreate(&e.value,CU_EVENT_BLOCKING_SYNC),"cuEventCreate");
    }
    ~Impl(){context.cleanup_context();}
    void ready() const {io::require(!poisoned,"Batched-force backend is poisoned; reconstruct it");io::require(std::this_thread::get_id()==thread,"Batched-force owner thread changed");context.current();}
    void upload(Slot slot,const void* source,std::size_t bytes) {io::require(bytes==device[slot].bytes,"Static snapshot size changed; reconstruct backend");cuda_ok(cuMemcpyHtoD(device[slot].value,source,bytes),"cuMemcpyHtoD static snapshot");}
};
CudaBatchedForce::CudaBatchedForce(const std::filesystem::path& path,const HostSnapshot& s,std::size_t capacity) {
    // Validate before the constructor above can initialize CUDA.
    s.validate(true);io::require(capacity>=s.particles.size() && capacity<=2000000,"Invalid persistent particle capacity");
    std::size_t total=0;for(auto n:allocation_sizes(s,capacity))total+=n;
    io::require(total<=128*1024*1024,"Batched-force128MiB explicit allocation guard exceeded before CUDA initialization");
    p_=std::make_unique<Impl>(path,s,capacity);set_snapshot(s);
}
CudaBatchedForce::~CudaBatchedForce()=default;
void CudaBatchedForce::set_snapshot(const HostSnapshot& s) {
    auto& b=*p_;b.ready();io::require(!b.pending,"Discard or commit pending batch before snapshot refresh");
    try {
        s.validate(true);io::require(s.phi.size()==b.self_nodes && s.az0.size()==b.native_nodes,"Snapshot nodal counts changed");
        b.upload(SelfCells,s.self_cells.data(),s.self_cells.size()*sizeof(gp::Cell));b.upload(NativeCells,s.native_cells.data(),s.native_cells.size()*sizeof(gp::Cell));
        const std::array<std::pair<Slot,const std::vector<double>*>,7> fields={{{Psi,&s.psi},{Az0,&s.az0},{Az1,&s.az1},{Rate0,&s.rate0},{Rate1,&s.rate1},{Xs,&s.xs},{Ys,&s.ys}}};
        for(const auto& [slot,field]:fields)b.upload(slot,field->data(),field->size()*8);
        b.static_parameters=s.parameters;
    }catch(...){b.poisoned=true;throw;}
}
std::span<const gp::Outputs> CudaBatchedForce::evaluate(std::span<const GapParticle> particles,
        std::span<const double> values,const Parameters& p,std::span<const gp::Outputs> reference) {
    const auto start=Clock::now();
    auto& b=*p_;b.ready();io::require(!b.pending,"Discard or commit prior batch first");
    try {
        const auto n=particles.size();
        io::require(n>0 && n<=b.capacity && reference.size()==n && values.size()==b.self_nodes,"Batched-force shape/capacity mismatch");
        const auto& fixed=b.static_parameters;
        io::require(p.magic==snapshot_magic && p.flags==1 && p.self_nx==fixed.self_nx && p.self_ny==fixed.self_ny &&
            p.native_nx==fixed.native_nx && p.native_ny==fixed.native_ny && p.width==fixed.width && p.height==fixed.height &&
            p.depth==fixed.depth && p.short_z==fixed.short_z,"Refresh static snapshot after geometry/layout changes");
        io::require(std::isfinite(p.fraction) && p.fraction>=0 && p.fraction<=1 && std::isfinite(p.voltage) && std::isfinite(p.dt) && p.dt>0,"Invalid midpoint parameters");
        auto* packed=static_cast<gp::Inputs*>(b.inputs.value);
        for(std::size_t i=0;i<n;++i) {
            const auto& q=particles[i];const auto& r=reference[i];
            io::require(gp::finite(io::vec(q.position_m)) && gp::finite(io::vec(q.proper_velocity_m_s)) && std::isfinite(q.charge_C) &&
                std::isfinite(q.mass_kg) && q.mass_kg>0 && std::isfinite(q.weight) && q.weight>0 &&
                std::isfinite(q.charge_C*q.weight) && std::isfinite(q.charge_C/q.mass_kg) && gp::finite(gp::norm2(io::vec(q.proper_velocity_m_s))),"Nonfinite/invalid input particle");
            io::require(q.position_m[0]>=0 && q.position_m[0]<=p.width && q.position_m[1]>=0 && q.position_m[1]<=p.height,"Particle outside supported domain");
            io::require(r.status==0 && gp::finite(r.E) && gp::finite(r.B) && gp::finite(r.u) && gp::finite(gp::norm2(r.u)),"Invalid independent CPU reference");
            b.pending_inputs[i]=q;packed[i]={q.position_m[0],q.position_m[1],q.position_m[2],io::vec(q.proper_velocity_m_s),q.charge_C/q.mass_kg};
        }
        for(double v:values)io::require(std::isfinite(v),"Nonfinite recovered phi");
        std::memcpy(b.phi.value,values.data(),values.size_bytes());const auto packed_at=Clock::now();
        gp::Tables table{{b.device[SelfCells].value?reinterpret_cast<const gp::Cell*>(b.device[SelfCells].value):nullptr,p.self_nx,p.self_ny},
            {reinterpret_cast<const gp::Cell*>(b.device[NativeCells].value),p.native_nx,p.native_ny},
            reinterpret_cast<const double*>(b.device[Phi].value),reinterpret_cast<const double*>(b.device[Psi].value),
            reinterpret_cast<const double*>(b.device[Az0].value),reinterpret_cast<const double*>(b.device[Az1].value),
            reinterpret_cast<const double*>(b.device[Rate0].value),reinterpret_cast<const double*>(b.device[Rate1].value),
            reinterpret_cast<const double*>(b.device[Xs].value),reinterpret_cast<const double*>(b.device[Ys].value),
            p.width,p.height,p.depth,p.short_z,p.fraction,p.voltage,p.dt};
        auto din=b.device[Inputs].value,dout=b.device[Outputs].value;auto count=static_cast<unsigned>(n);
        void* arguments[]={&table,&din,&dout,&count};
        cuda_ok(cuEventRecord(b.events[0].value,b.stream.value),"event start");
        cuda_ok(cuMemcpyHtoDAsync(din,b.inputs.value,n*sizeof(gp::Inputs),b.stream.value),"H2D particles");
        cuda_ok(cuMemcpyHtoDAsync(b.device[Phi].value,b.phi.value,values.size_bytes(),b.stream.value),"H2D phi");
        cuda_ok(cuEventRecord(b.events[1].value,b.stream.value),"event after H2D");
        cuda_ok(cuLaunchKernel(b.kernel,(count+127)/128,1,1,128,1,1,0,b.stream.value,arguments,nullptr),"gather_push_kernel");
        cuda_ok(cuEventRecord(b.events[2].value,b.stream.value),"event after kernel");
        cuda_ok(cuMemcpyDtoHAsync(b.outputs.value,dout,n*sizeof(gp::Outputs),b.stream.value),"D2H force results");
        cuda_ok(cuEventRecord(b.events[3].value,b.stream.value),"event after D2H");
        cuda_ok(cuEventSynchronize(b.events[3].value),"final batch synchronization");const auto synchronized=Clock::now();
        float stages[3];for(int i=0;i<3;++i)cuda_ok(cuEventElapsedTime(&stages[i],b.events[i].value,b.events[i+1].value),"event elapsed");
        const std::span<const gp::Outputs> result{static_cast<const gp::Outputs*>(b.outputs.value),n};
        b.comparison=batch_compare::compare(reference,result,particles);
        io::require(b.comparison.failures==0,"GPU batch failed independent finite/componentwise/weighted CPU reference gates");
        const auto validated=Clock::now();b.last={elapsed(start,packed_at),stages[0],stages[1],stages[2],elapsed(start,synchronized),elapsed(synchronized,validated),elapsed(start,validated)};
        b.pending_count=n;b.pending=true;return result;
    }catch(...){b.pending=false;b.poisoned=true;throw;}
}
void CudaBatchedForce::commit_velocities(std::span<GapParticle> particles) {
    auto& b=*p_;b.ready();io::require(b.pending && particles.size()==b.pending_count,"No validated batch to commit");
    for(std::size_t i=0;i<particles.size();++i)io::require(std::memcmp(&particles[i],&b.pending_inputs[i],sizeof(GapParticle))==0,"Particle state changed before transactional commit");
    const auto* result=static_cast<const gp::Outputs*>(b.outputs.value);
    for(std::size_t i=0;i<particles.size();++i)particles[i].proper_velocity_m_s=io::vec(result[i].u);
    b.pending=false;
}
void CudaBatchedForce::discard(){p_->ready();p_->pending=false;}
const BatchTiming& CudaBatchedForce::timing() const{return p_->last;}
const batch_compare::Comparison& CudaBatchedForce::comparison() const{return p_->comparison;}
std::size_t CudaBatchedForce::explicit_device_bytes() const{return p_->total_bytes;}
}
