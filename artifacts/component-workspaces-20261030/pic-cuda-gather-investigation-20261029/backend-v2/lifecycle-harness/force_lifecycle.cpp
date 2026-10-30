// Qualification only. No application hook, production action, or device fault injection.
#include "native_reference.hpp"
#include "cuda_batched_force.hpp"
#include <cuda.h>
#include "cudss_poisson.hpp" // original frozen solver, unchanged
#include <chrono>
#include <iomanip>
#include <iostream>
#include <sched.h>
#include <sys/resource.h>

using namespace pic_batch;
namespace io=snapshot_io;
namespace gp=gather_push;
using Clock=std::chrono::steady_clock;

// Link-time observation of this executable/static backend's Driver API calls.
// Dynamic cuDSS/runtime internals are not interposed or counted. No CUDA call
// is replaced with a mock; only attempts to reset/destroy the borrowed primary
// are rejected. The harness destroys its own deliberately different context.
struct DriverAudit {
    CUcontext borrowed=nullptr,last_retained=nullptr;
    int retains=0,releases=0,launches=0,h2d=0,d2h=0;
    int reset_attempts=0,borrowed_destroy_attempts=0,owned_context_destroys=0;
} audit;
extern "C" {
CUresult CUDAAPI __real_cuDevicePrimaryCtxRetain(CUcontext*,CUdevice);
CUresult CUDAAPI __real_cuDevicePrimaryCtxRelease_v2(CUdevice);
CUresult CUDAAPI __real_cuCtxDestroy_v2(CUcontext);
CUresult CUDAAPI __real_cuLaunchKernel(CUfunction,unsigned,unsigned,unsigned,unsigned,unsigned,unsigned,unsigned,CUstream,void**,void**);
CUresult CUDAAPI __real_cuMemcpyHtoDAsync_v2(CUdeviceptr,const void*,size_t,CUstream);
CUresult CUDAAPI __real_cuMemcpyDtoHAsync_v2(void*,CUdeviceptr,size_t,CUstream);
CUresult CUDAAPI __wrap_cuDevicePrimaryCtxRetain(CUcontext* value,CUdevice device) {
    const auto r=__real_cuDevicePrimaryCtxRetain(value,device);
    if(r==CUDA_SUCCESS){++audit.retains;audit.last_retained=*value;}return r;
}
CUresult CUDAAPI __wrap_cuDevicePrimaryCtxRelease_v2(CUdevice device) {
    const auto r=__real_cuDevicePrimaryCtxRelease_v2(device);
    if(r==CUDA_SUCCESS)++audit.releases;
    return r;
}
CUresult CUDAAPI __wrap_cuDevicePrimaryCtxReset_v2(CUdevice) {
    ++audit.reset_attempts;return CUDA_ERROR_INVALID_CONTEXT;
}
CUresult CUDAAPI __wrap_cuCtxDestroy_v2(CUcontext value) {
    if(value==audit.borrowed && value){++audit.borrowed_destroy_attempts;return CUDA_ERROR_INVALID_CONTEXT;}
    const auto r=__real_cuCtxDestroy_v2(value);if(r==CUDA_SUCCESS)++audit.owned_context_destroys;return r;
}
CUresult CUDAAPI __wrap_cuLaunchKernel(CUfunction f,unsigned gx,unsigned gy,unsigned gz,unsigned bx,unsigned by,unsigned bz,unsigned shared,CUstream stream,void** args,void** extra) {
    ++audit.launches;return __real_cuLaunchKernel(f,gx,gy,gz,bx,by,bz,shared,stream,args,extra);
}
CUresult CUDAAPI __wrap_cuMemcpyHtoDAsync_v2(CUdeviceptr target,const void* source,size_t n,CUstream stream) {
    ++audit.h2d;return __real_cuMemcpyHtoDAsync_v2(target,source,n,stream);
}
CUresult CUDAAPI __wrap_cuMemcpyDtoHAsync_v2(void* target,CUdeviceptr source,size_t n,CUstream stream) {
    ++audit.d2h;return __real_cuMemcpyDtoHAsync_v2(target,source,n,stream);
}
}
namespace {
void cuda_ok(CUresult code,const char* operation) {
    if(code==CUDA_SUCCESS)return;
    const char* detail=nullptr;cuGetErrorString(code,&detail);
    throw std::runtime_error(std::string(operation)+": "+(detail?detail:"CUDA error"));
}
CUcontext current() {CUcontext value=nullptr;cuda_ok(cuCtxGetCurrent(&value),"cuCtxGetCurrent");return value;}
int primary_active() {unsigned flags=0;int active=-1;cuda_ok(cuDevicePrimaryCtxGetState(0,&flags,&active),"cuDevicePrimaryCtxGetState");return active;}
std::string json_string(const std::string& value) {
    std::ostringstream out;out<<'"';
    for(unsigned char c:value) {
        if(c=='"' || c=='\\')out<<'\\'<<char(c);
        else if(c<32)out<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<int(c)<<std::dec;
        else out<<char(c);
    }
    out<<'"';return out.str();
}
std::string pointer(CUcontext value) {std::ostringstream out;out<<"0x"<<std::hex<<reinterpret_cast<std::uintptr_t>(value);return out.str();}
struct Event {std::string stage,context;int active=0,retains=0,releases=0,launches=0,h2d=0,d2h=0;};
struct Solve {std::string stage;double norm=0,residual=0,limit=0,relative=0;};
struct Receipt {
    CleanupReport shared,wrong,lost,standalone;
    std::vector<Event> events;
    std::vector<Solve> solves;
    bool wrong_rejected=false,lost_rejected=false,passed=false;
    int standalone_active_before=-1,standalone_active_after=-1;
    std::size_t particles=0,rows=0,forces_completed=0;
    std::string mode,error;
    void event(const std::string& stage) {
        events.push_back({stage,pointer(current()),primary_active(),audit.retains,audit.releases,audit.launches,audit.h2d,audit.d2h});
    }
    void write(const std::filesystem::path& output) const {
        std::filesystem::create_directories(output);
        std::ofstream file(output/"lifecycle.json");file<<std::setprecision(17)<<std::boolalpha;
        file<<"{\"passed\":"<<passed<<",\"mode\":"<<json_string(mode)<<",\"original_error\":"<<json_string(error)
            <<",\"particles\":"<<particles<<",\"manufactured_Poisson_rows\":"<<rows
            <<",\"actual_midpoint_fixture_required\":true,\"forces_completed\":"<<forces_completed<<",\"GPU_kernel_dispatched\":"<<(audit.launches>0)
            <<",\"full_PIC_replay\":false,\"production_actions\":false,\"device_fault_injection\":false"
            <<",\"borrowed_primary\":"<<json_string(pointer(audit.borrowed))<<",\"successful_own_retains\":"<<audit.retains
            <<",\"successful_own_releases\":"<<audit.releases<<",\"reset_attempts\":"<<audit.reset_attempts
            <<",\"borrowed_destroy_attempts\":"<<audit.borrowed_destroy_attempts<<",\"owned_context_destroys\":"<<audit.owned_context_destroys
            <<",\"wrong_current_rejected_before_dispatch\":"<<wrong_rejected<<",\"lost_current_rejected_before_dispatch\":"<<lost_rejected
            <<",\"standalone_primary_active_before\":"<<standalone_active_before<<",\"standalone_primary_active_after\":"<<standalone_active_after
            <<",\"cleanup\":{";
        const std::array<std::pair<const char*,const CleanupReport*>,4> reports={{{"shared",&shared},{"wrong",&wrong},{"lost",&lost},{"standalone",&standalone}}};
        for(std::size_t k=0;k<reports.size();++k) {
            const auto& [name,r]=reports[k];file<<(k?",":"")<<json_string(name)<<":{\"attempted\":"<<r->attempted
                <<",\"context_established\":"<<r->context_established<<",\"stream_drained\":"<<r->stream_drained
                <<",\"used_context_fallback\":"<<r->used_context_fallback<<",\"resources_released\":"<<r->resources_released
                <<",\"resources_retained\":"<<r->resources_retained<<",\"caller_context_restored\":"<<r->caller_context_restored
                <<",\"primary_reference_released\":"<<r->primary_reference_released<<",\"issue_count\":"<<r->issue_count<<",\"issues\":[";
            for(std::size_t i=0;i<r->issue_count && i<r->issues.size();++i)
                file<<(i?",":"")<<"{\"operation\":"<<json_string(r->issues[i].operation?r->issues[i].operation:"unknown")<<",\"code\":"<<r->issues[i].code<<'}';
            file<<"]}";
        }
        file<<"},\"events\":[";
        for(std::size_t i=0;i<events.size();++i) {
            const auto& e=events[i];file<<(i?",":"")<<"{\"stage\":"<<json_string(e.stage)<<",\"current\":"<<json_string(e.context)
                <<",\"primary_active\":"<<e.active<<",\"own_retains\":"<<e.retains<<",\"own_releases\":"<<e.releases
                <<",\"launches\":"<<e.launches<<",\"H2D_dispatches\":"<<e.h2d<<",\"D2H_dispatches\":"<<e.d2h<<'}';
        }
        file<<"],\"solves\":[";
        for(std::size_t i=0;i<solves.size();++i) {
            const auto& s=solves[i];file<<(i?",":"")<<"{\"stage\":"<<json_string(s.stage)<<",\"rhs_norm\":"<<s.norm
                <<",\"residual_norm\":"<<s.residual<<",\"guard_limit\":"<<s.limit<<",\"relative_residual\":"<<s.relative<<'}';
        }
        file<<"],\"scope\":\"Original frozen cuDSS wrapper on manufactured Q1 Poisson at actual mesh resolution; exact original residual threshold. Actual captured midpoint force plus fresh original CPU references and unchanged gates. Not a deposited-charge field solve or evolving application replay. Driver retain/release audit covers statically linked backend calls, not dynamic runtime/cuDSS internals.\"}\n";
        file.close();io::require(bool(file),"Cannot persist lifecycle/cleanup receipt");
    }
};
void clean(const CleanupReport& r) {
    io::require(r.attempted && r.context_established && r.stream_drained && r.resources_released && !r.resources_retained &&
        r.caller_context_restored && r.primary_reference_released && r.issue_count==0,"Lifecycle cleanup disqualified run; see lifecycle.json");
}
void balanced(int before_retains,int before_releases) {
    io::require(audit.retains==before_retains+1 && audit.releases==before_releases+1,"Backend own primary reference not retained/released exactly once");
}
void shared_current() {io::require(current()==audit.borrowed && primary_active()==1,"cuDSS shared primary identity/lifetime changed");}

class ManufacturedPoisson {
    mfem::Mesh mesh;
    mfem::H1_FECollection collection;
    mfem::FiniteElementSpace space;
    mfem::ConstantCoefficient epsilon;
    mfem::BilinearForm form;
    mfem::GridFunction trial;
    mfem::LinearForm zero;
    mfem::OperatorPtr matrix;
    mfem::Vector x,b,known,rhs;
    mfem::SparseMatrix* sparse=nullptr;
    std::unique_ptr<raikou_pic_parallel::CudssPoisson> solver;
public:
    explicit ManufacturedPoisson(const HostSnapshot& s)
        :mesh(mfem::Mesh::MakeCartesian2D(s.parameters.self_nx,s.parameters.self_ny,mfem::Element::QUADRILATERAL,true,s.parameters.width,s.parameters.height,false)),
         collection(1,2),space(&mesh,&collection),epsilon(8.8541878188e-12),form(&space),trial(&space),zero(&space) {
        form.AddDomainIntegrator(new mfem::DiffusionIntegrator(epsilon));form.Assemble();form.Finalize();
        mfem::Array<int> marker(mesh.bdr_attributes.Max());marker=0;
        mfem::Array<int> vertices;
        for(int boundary=0;boundary<mesh.GetNBE();++boundary) {
            mesh.GetBdrElementVertices(boundary,vertices);
            const double y0=mesh.GetVertex(vertices[0])[1],y1=mesh.GetVertex(vertices[1])[1];
            if(y0==y1 && (y0==0 || y0==s.parameters.height))marker[mesh.GetBdrAttribute(boundary)-1]=1;
        }
        mfem::Array<int> essential;space.GetEssentialTrueDofs(marker,essential);trial=0.;zero=0.;
        form.FormLinearSystem(essential,trial,zero,matrix,x,b);
        sparse=dynamic_cast<mfem::SparseMatrix*>(matrix.Ptr());io::require(sparse!=nullptr,"Expected manufactured sparse Q1 operator");
        io::require(s.phi.size()==std::size_t(sparse->Height()),"Manufactured system resolution differs from actual captured mesh");
        known.SetSize(sparse->Width());std::copy(s.phi.begin(),s.phi.end(),known.GetData());
        // Manufactured RHS, intentionally NOT particle charge deposition.
        // Eliminate boundary values to avoid mixing unit diagonal with eps0.
        for(int d:essential)known[d]=0.;
        rhs.SetSize(sparse->Height());sparse->Mult(known,rhs);
        io::require(std::isfinite(rhs.Norml2()) && rhs.Norml2()>0,"Invalid manufactured Poisson RHS");
        solver=std::make_unique<raikou_pic_parallel::CudssPoisson>(*sparse);
    }
    void solve(const std::string& stage,Receipt& receipt) {
        solver->Solve(rhs,x);
        for(int i=0;i<x.Size();++i)io::require(std::isfinite(x[i]),"Nonfinite cuDSS solution");
        mfem::Vector error(rhs.Size());sparse->Mult(x,error);error-=rhs;
        const double residual=error.Norml2(),norm=rhs.Norml2(),relative=norm>0?residual/norm:residual;
        const double limit=std::max(1e-30,2e-13*norm); // exact original LinearSystem guard
        receipt.solves.push_back({stage,norm,residual,limit,relative});
        io::require(residual<=limit && std::isfinite(relative),"Original CUDA Poisson residual/finite guard failed");
        receipt.rows=rhs.Size();receipt.event(stage);
    }
};
void force(CudaBatchedForce& backend,const HostSnapshot& snapshot,native_reference::NativeReference& native,
        const std::filesystem::path& output,const std::string& stage,Receipt& receipt) {
    std::vector<gp::Outputs> reference(snapshot.particles.size());
    native.generate(snapshot.particles,reference);
    io::require(std::memcmp(reference.data(),snapshot.reference.data(),reference.size()*sizeof(gp::Outputs))==0,"Fresh original CPU E/B/u differs from true midpoint capture");
    auto live=snapshot.particles;
    const auto result=backend.evaluate(snapshot.particles,snapshot.phi,snapshot.parameters,reference);
    io::require(backend.comparison().failures==0,"Original independent-reference gates failed");
    backend.commit_velocities(live);
    for(std::size_t i=0;i<live.size();++i) {
        auto restored=live[i];restored.proper_velocity_m_s=snapshot.particles[i].proper_velocity_m_s;
        io::require(std::memcmp(&restored,&snapshot.particles[i],sizeof(restored))==0,"Force commit changed particle metadata");
        io::require(std::memcmp(live[i].proper_velocity_m_s.data(),result[i].u.v,24)==0,"Force commit differs from qualified result");
    }
    io::write(output/(stage+"-device.bin"),std::vector<gp::Outputs>(result.begin(),result.end()));
    std::ofstream comparison(output/(stage+"-comparison.json"));batch_compare::print_comparison(backend.comparison(),comparison);
    comparison.close();io::require(bool(comparison),"Cannot save force qualification comparison");
    ++receipt.forces_completed;receipt.event(stage);
}
bool reject_without_dispatch(CudaBatchedForce& backend,const HostSnapshot& snapshot) {
    const int kernels=audit.launches,h2d=audit.h2d,d2h=audit.d2h;bool rejected=false;
    try{backend.evaluate(snapshot.particles,snapshot.phi,snapshot.parameters,snapshot.reference);}
    catch(const std::runtime_error& error){rejected=std::string(error.what())=="Batched-force CUDA context changed";}
    io::require(rejected && audit.launches==kernels && audit.h2d==h2d && audit.d2h==d2h,"Wrong/lost context did not reject before asynchronous dispatch");return true;
}
class OwnedWrongContext {
    CUcontext value=nullptr,previous=nullptr;
public:
    OwnedWrongContext():previous(current()){cuda_ok(cuCtxCreate(&value,0,0),"cuCtxCreate owned wrong context");io::require(value!=audit.borrowed,"Wrong context unexpectedly equals shared primary");}
    void close() {
        if(!value)return;
        io::require(current()==value,"Owned wrong context was not restored by backend close");
        CUcontext popped=nullptr;cuda_ok(cuCtxPopCurrent(&popped),"cuCtxPopCurrent owned wrong context");
        io::require(popped==value && current()==previous,"Owned wrong context stack was not balanced");
        cuda_ok(cuCtxDestroy(value),"cuCtxDestroy owned wrong context only");value=nullptr;
    }
    ~OwnedWrongContext(){if(value){cuCtxSetCurrent(previous);cuCtxDestroy(value);}}
};
}
int main(int argc,char** argv) {
    Receipt receipt;std::filesystem::path output;
    // All cleanup sinks outlive their backends AND the outer catch.
    try {
        io::require(argc==6 && std::string(argv[4])=="--parent-authorized" &&
            (std::string(argv[5])=="--shared-lifecycle" || std::string(argv[5])=="--standalone-primary"),
            "Expected ACTUAL_FIXTURE PTX OUTPUT --parent-authorized --shared-lifecycle|--standalone-primary");
        output=argv[3];receipt.mode=argv[5];std::filesystem::create_directories(output);
        cpu_set_t mask;CPU_ZERO(&mask);
        io::require(sched_getaffinity(0,sizeof(mask),&mask)==0 && CPU_COUNT(&mask)==1 && CPU_ISSET(7,&mask) && getpriority(PRIO_PROCESS,0)>=15,
            "Lifecycle qualification requires parent-owned CPU7/nice15 and external GPU/core7 leases");
        const auto start=Clock::now();
        const auto snapshot=HostSnapshot::load(argv[1],true);receipt.particles=snapshot.particles.size();
        native_reference::NativeReference native(snapshot,argv[1],"/home/treefor/sim-workspace/powerflow/outputs/gap-pulse-midplane110kV.csv");
        auto bounded=[&]{io::require(std::chrono::duration<double>(Clock::now()-start).count()<30,"Lifecycle internal30-second guard exceeded");};
        if(receipt.mode=="--shared-lifecycle") {
            ManufacturedPoisson poisson(snapshot);poisson.solve("solve1_before_force",receipt);
            audit.borrowed=current();io::require(audit.borrowed!=nullptr,"cuDSS did not establish current shared primary");shared_current();
            const int r=audit.retains,d=audit.releases;
            CudaBatchedForce backend(argv[2],snapshot,snapshot.particles.size(),{ContextPolicy::RequireCurrentSharedPrimary,&receipt.shared});
            io::require(audit.last_retained==audit.borrowed,"Backend did not retain cuDSS shared primary");
            shared_current();force(backend,snapshot,native,output,"force1",receipt);shared_current();bounded();
            poisson.solve("solve2_between_forces",receipt);shared_current();
            force(backend,snapshot,native,output,"force2",receipt);shared_current();
            receipt.shared=backend.close();clean(receipt.shared);balanced(r,d);shared_current();receipt.event("shared_explicit_close");
            poisson.solve("solve3_after_force_close",receipt);shared_current();bounded();
            {
                const int wr=audit.retains,wd=audit.releases;
                CudaBatchedForce wrong(argv[2],snapshot,snapshot.particles.size(),{ContextPolicy::RequireCurrentSharedPrimary,&receipt.wrong});
                OwnedWrongContext owned;receipt.event("wrong_current_selected");
                receipt.wrong_rejected=reject_without_dispatch(wrong,snapshot);
                receipt.wrong=wrong.close();clean(receipt.wrong);balanced(wr,wd);receipt.event("wrong_close_restored_owned_context");
                owned.close();shared_current();receipt.event("wrong_owned_context_destroyed_only");
            }
            poisson.solve("solve4_after_wrong_context_probe",receipt);shared_current();bounded();
            {
                const int lr=audit.retains,ld=audit.releases;
                CudaBatchedForce lost(argv[2],snapshot,snapshot.particles.size(),{ContextPolicy::RequireCurrentSharedPrimary,&receipt.lost});
                cuda_ok(cuCtxSetCurrent(nullptr),"cuCtxSetCurrent null for lost-context rejection");receipt.event("lost_current_selected");
                receipt.lost_rejected=reject_without_dispatch(lost,snapshot);
                receipt.lost=lost.close();clean(receipt.lost);balanced(lr,ld);
                io::require(current()==nullptr,"Cleanup did not restore null caller context");receipt.event("lost_close_restored_null");
                cuda_ok(cuCtxSetCurrent(audit.borrowed),"Restore caller shared primary after lost-context probe");shared_current();
            }
            poisson.solve("solve5_after_lost_context_probe",receipt);shared_current();bounded();
        } else {
            cuda_ok(cuInit(0),"cuInit standalone probe");io::require(current()==nullptr,"Standalone process already has a current context");
            receipt.standalone_active_before=primary_active();receipt.event("standalone_before");const int r=audit.retains,d=audit.releases;
            CudaBatchedForce backend(argv[2],snapshot,snapshot.particles.size(),{ContextPolicy::StandalonePrimary,&receipt.standalone});
            io::require(current()==audit.last_retained && primary_active()==1,"Standalone primary was not selected");
            force(backend,snapshot,native,output,"standalone_force1",receipt);bounded();
            force(backend,snapshot,native,output,"standalone_force2",receipt);
            receipt.standalone=backend.close();clean(receipt.standalone);balanced(r,d);
            io::require(current()==nullptr,"Standalone close did not restore original null current context");
            receipt.standalone_active_after=primary_active();
            io::require(receipt.standalone_active_after==receipt.standalone_active_before,"Standalone own-reference lifetime changed primary activity");receipt.event("standalone_after_close");bounded();
        }
        io::require(audit.reset_attempts==0 && audit.borrowed_destroy_attempts==0 && audit.retains==audit.releases,"Forbidden shared-context lifetime operation or unbalanced backend references");
        receipt.passed=true;receipt.write(output);return 0;
    }catch(const std::exception& error) {
        receipt.passed=false;receipt.error=error.what();
        if(!output.empty())try{receipt.write(output);}catch(const std::exception& report){std::cerr<<"Receipt write failed: "<<report.what()<<'\n';}
        std::cerr<<error.what()<<'\n';return 1;
    }
}
