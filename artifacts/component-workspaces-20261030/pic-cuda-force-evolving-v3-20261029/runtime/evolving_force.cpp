#include "evolving_force.hpp"
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <numeric>
#include <sched.h>
#include <sys/resource.h>
#include <thread>

namespace pic_evolving {
namespace io=snapshot_io;namespace gp=gather_push;
using namespace pic_batch;using Clock=std::chrono::steady_clock;
namespace {
constexpr std::size_t capacity=65536,trace_budget=2ULL*1024*1024*1024;
constexpr std::uint64_t step_budget=256;
std::uint64_t initial_step=0;
bool restart_bounded=false;
struct Session {
    GapPic* owner=nullptr;Mode mode=Mode::Off;Provider provider;
    std::thread::id thread;
    std::filesystem::path ptx,output;
    RuntimeSnapshot snapshot;
    CleanupReport cleanup; // outlives backend and caller's outer catch
    std::unique_ptr<CudaForceRuntime> backend;
    const mfem::FiniteElementSpace* self_space=nullptr;
    std::uint64_t layout_epoch=0,unit_epoch=0,lower_frame=0,upper_frame=0,steps=0,refreshes=0;
    double lower_time=0,upper_time=0,relative=0,dt=0,reference_ms=0,api_ms=0,commit_ms=0;
    std::size_t trace_bytes=0;
    batch_compare::Comparison comparison;
    Clock::time_point force_start;
    std::ofstream trace,history;
    bool closed=false,failed=false;
    void bytes(const void* data,std::size_t n) {
        io::require(n<=trace_budget-trace_bytes,"Two-GiB isolated full-state trace budget exceeded");
        trace.write(static_cast<const char*>(data),n);io::require(bool(trace),"Full-state trace write failed");trace_bytes+=n;
    }
    template<class T>void scalar(const T& value){bytes(&value,sizeof(value));}
    template<class T>void array(std::span<const T> values){bytes(values.data(),values.size_bytes());}
} session;
double milliseconds(Clock::time_point a,Clock::time_point b){return std::chrono::duration<double,std::milli>(b-a).count();}
bool owns(const GapPic* owner){
    if(session.owner!=owner || session.closed)return false;
    io::require(std::this_thread::get_id()==session.thread,"Evolving force owner thread changed");
    io::require(!session.failed,"Evolving force session failed; no retry or native fallback");return true;
}
void clean(const CleanupReport& r) {
    io::require(r.attempted && r.context_established && r.stream_drained && r.resources_released && !r.resources_retained &&
                r.caller_context_restored && r.primary_reference_released && !r.issue_count,"Evolving candidate cleanup disqualified; see force-cleanup.json");
}
void static_fields(const NativeSelection& n,const mfem::GridFunction& phi,const GapPicOptions& o) {
    auto& s=session.snapshot;
    io::require(n.psi && n.az0 && n.az1 && n.rate0 && n.rate1 && n.elements && n.xs && n.ys,
                "Missing typed live native/unit state");
    io::require(n.layout_epoch>0 && n.unit_epoch>0 && n.upper_frame==n.lower_frame+1 &&
                std::isfinite(n.lower_time) && std::isfinite(n.upper_time) && n.lower_time<n.upper_time,
                "Invalid explicit native topology/unit/frame epochs");
    if(session.self_space) {
        io::require(phi.FESpace()==session.self_space && session.layout_epoch==n.layout_epoch && session.unit_epoch==n.unit_epoch,
                    "One-ps candidate requires unchanged self/unit topology and epoch; reconstruct explicitly");
    } else {
        session.self_space=phi.FESpace();
        std::vector<int> order(phi.FESpace()->GetMesh()->GetNE());std::iota(order.begin(),order.end(),0);
        s.self_cells=io::export_grid(phi,order);
        const auto unit_cells=io::export_grid(*n.psi,order);
        io::require(unit_cells.size()==s.self_cells.size() && std::memcmp(unit_cells.data(),s.self_cells.data(),unit_cells.size()*sizeof(gp::Cell))==0,
                    "Actual self/unit layouts differ");
        s.psi=io::values(*n.psi);session.layout_epoch=n.layout_epoch;session.unit_epoch=n.unit_epoch;
    }
    s.native_cells=io::export_grid(*n.az0,*n.elements);
    for(const auto* field:{n.az1,n.rate0,n.rate1}) {
        const auto cells=io::export_grid(*field,*n.elements);
        io::require(cells.size()==s.native_cells.size() && std::memcmp(cells.data(),s.native_cells.data(),cells.size()*sizeof(gp::Cell))==0,
                    "Native Az/rate lower/upper layouts differ");
    }
    s.az0=io::values(*n.az0);s.az1=io::values(*n.az1);s.rate0=io::values(*n.rate0);s.rate1=io::values(*n.rate1);
    s.xs=*n.xs;s.ys=*n.ys;
    s.parameters.self_nx=o.nx;s.parameters.self_ny=o.ny;
    s.parameters.native_nx=int(s.xs.size())-1;s.parameters.native_ny=int(s.ys.size())-1;
    s.parameters.width=o.width_m;s.parameters.height=o.gap_m;s.parameters.depth=o.depth_m;s.parameters.short_z=o.short_z_m;
    session.lower_frame=n.lower_frame;session.upper_frame=n.upper_frame;session.lower_time=n.lower_time;session.upper_time=n.upper_time;
}
}
Mode requested_mode() {
    const std::string value=std::getenv("RAIKOU_PIC_FORCE_MODE")?std::getenv("RAIKOU_PIC_FORCE_MODE"):"off";
    if(value=="off")return Mode::Off;
    if(value=="shadow")return Mode::Shadow;
    if(value=="optimized")return Mode::OptimizedCandidate;
    throw std::invalid_argument("Expected explicit off/shadow/optimized force mode");
}
void candidate_environment() {
    cpu_set_t mask;CPU_ZERO(&mask);
    io::require(sched_getaffinity(0,sizeof(mask),&mask)==0 && CPU_COUNT(&mask)==1 && CPU_ISSET(7,&mask) && getpriority(PRIO_PROCESS,0)>=15,
                "Evolving GPU candidate requires parent CPU7/nice15 and GPU-then-core7 leases");
    const char* token=std::getenv("RAIKOU_PIC_BOUNDED_REPLAY");
    io::require(token && std::string(token)=="parent-approved-one-ps-v3","Deployment hard disabled; explicit parent one-ps replay token required");
    for(const char* name:{"OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"})
        io::require(std::getenv(name) && std::string(std::getenv(name))=="1","Evolving candidate requires one thread");
    static_cast<void>(requested_mode());
}
void validate_output(const std::filesystem::path& output) {
    const auto root=std::filesystem::weakly_canonical("/home/treefor/projects/runtimes/pic-cuda-force-evolving-v3-20261029/qualification-runs");
    const auto path=std::filesystem::weakly_canonical(output);
    const auto relative=path.lexically_relative(root);
    io::require(!relative.empty() && relative!="." && *relative.begin()!=".." && path!=root,
                "Deployment hard disabled: output must stay inside owned candidate qualification-runs");
}
void bound_restart(double origin,double relative_start,double end,std::uint64_t step) {
    const double start=origin+relative_start;
    io::require(step==437411 && std::isfinite(start) && std::abs(start-2.6478441131125727e-08)<=1e-22 &&
                std::isfinite(end) && end>start && end-start<=1.00000001e-12,
                "Deployment hard disabled: exact late seed and at most one-ps replay required");
    initial_step=step;restart_bounded=true;
}
void configure(GapPic& owner,Mode mode,const std::filesystem::path& ptx,const std::filesystem::path& output,Provider provider) {
    io::require(restart_bounded && !session.owner && provider,"Unbounded or duplicate evolving force configuration");
    validate_output(output);
    session.owner=&owner;session.thread=std::this_thread::get_id();session.mode=mode;session.ptx=ptx;session.output=output;session.provider=std::move(provider);
    session.trace.open(output/"evolving-state.bin",std::ios::binary);session.history.open(output/"force-stages.csv");
    io::require(bool(session.trace) && bool(session.history),"Cannot open bounded evolving qualification traces");
    std::ofstream schema(output/"trace-schema.json");
    schema<<"{\"version\":3,\"endian\":\"little\",\"particle_bytes\":"<<sizeof(GapParticle)<<",\"diagnostics_bytes\":"<<sizeof(GapPicDiagnostics)
        <<",\"record_order\":\"uint64 magic0x5049435354455033, uint64 step, double dt, raw GapPicDiagnostics, uint64 particle_count/phi_count/horizontal_count/vertical_count, raw GapParticle[], double phi[], horizontal[], vertical[], four double reservoir spans, four repetitions of uint64 size/double integral[]/double weighted_temperature[]\"}\n";
    schema.close();io::require(bool(schema),"Cannot write exact native ABI trace schema");
    session.history<<"relative_midpoint_s,dt_s,mode,reference_ms,GPU_API_ms,commit_ms,whole_force_ms,lower_frame,upper_frame,voltage_V,fraction,comparison_failures,max_E_error,max_B_error,max_u_error,weighted_u_error,weighted_energy_error\n";
}
bool active(const GapPic* owner){return owns(owner) && session.mode!=Mode::Off;}
bool shadow(const GapPic* owner){return owns(owner) && session.mode==Mode::Shadow;}
void begin_force(const GapPic* owner,double relative_midpoint,double dt) {
    if(!owns(owner))return;
    session.relative=relative_midpoint;session.dt=dt;session.reference_ms=0;session.api_ms=0;session.commit_ms=0;
    session.comparison={};session.force_start=Clock::now();
}
void reference_ready(const GapPic* owner){if(shadow(owner))session.reference_ms=milliseconds(session.force_start,Clock::now());}
void force(GapPic* owner,std::span<GapParticle> particles,const mfem::GridFunction& phi,const GapPicOptions& o,
           double relative_midpoint,double dt,std::span<const gp::Outputs> reference) {
    try {
    io::require(active(owner),"No explicit force candidate binding");
    io::require(particles.size()<=capacity,"One-ps fixed particle capacity exceeded; no fallback/substepping");
    if(particles.empty()){io::require(reference.empty(),"Empty force/reference mismatch");return;}
    io::require((session.mode==Mode::Shadow && reference.size()==particles.size()) ||
                (session.mode==Mode::OptimizedCandidate && reference.empty()),"Shadow reference policy/mode mismatch");
    const auto n=session.provider(relative_midpoint); // exact original selector, explicitly typed
    io::require(n.absolute_midpoint==n.origin+relative_midpoint && std::isfinite(n.absolute_midpoint) &&
                std::isfinite(n.fraction) && n.fraction>=0 && n.fraction<=1 && std::isfinite(n.voltage) &&
                n.lower_time<n.upper_time && n.fraction==std::clamp((n.absolute_midpoint-n.lower_time)/(n.upper_time-n.lower_time),0.,1.),
                "Invalid typed midpoint clock/parameters");
    auto& s=session.snapshot;s.parameters.flags=1;s.parameters.dt=dt;s.parameters.fraction=n.fraction;s.parameters.voltage=n.voltage;
    s.phi=io::values(phi); // exact current recovered phi, never a saved/cached midpoint
    const bool refresh=!session.backend || session.lower_frame!=n.lower_frame || session.upper_frame!=n.upper_frame ||
                       session.lower_time!=n.lower_time || session.upper_time!=n.upper_time;
    io::require(!session.backend || (session.layout_epoch==n.layout_epoch && session.unit_epoch==n.unit_epoch &&
                session.self_space==phi.FESpace()),"Static/unit epoch changed; isolated candidate disqualified");
    if(refresh) {
        static_fields(n,phi,o);s.validate();
        if(session.backend)session.backend->set_snapshot(s);
        else session.backend=std::make_unique<CudaForceRuntime>(session.ptx,s,capacity,BackendOptions{ContextPolicy::RequireCurrentSharedPrimary,&session.cleanup});
        ++session.refreshes;
    }
    const auto begin=Clock::now();const auto result=session.backend->evaluate(particles,s.phi,s.parameters);
    session.api_ms=milliseconds(begin,Clock::now());
    if(session.mode==Mode::Shadow) {
        session.comparison=batch_compare::compare(reference,result,particles);
        if(session.comparison.failures){session.backend->discard();throw std::runtime_error("Evolving shadow independent finite/componentwise/weighted reference gates failed");}
    }
    const auto commit=Clock::now();session.backend->commit_velocities(particles);session.commit_ms=milliseconds(commit,Clock::now());
    }catch(...){session.failed=true;throw;}
}
void end_force(const GapPic* owner) {
    if(!owns(owner))return;
    const auto& p=session.snapshot.parameters;const auto& c=session.comparison;
    session.history<<std::setprecision(17)<<session.relative<<','<<session.dt<<','<<int(session.mode)<<','<<session.reference_ms<<','<<session.api_ms<<','<<session.commit_ms
        <<','<<milliseconds(session.force_start,Clock::now())<<','<<session.lower_frame<<','<<session.upper_frame<<','<<p.voltage<<','<<p.fraction
        <<','<<c.failures<<','<<c.max_E<<','<<c.max_B<<','<<c.max_u<<','<<c.weighted_u<<','<<c.weighted_energy<<'\n';
    io::require(bool(session.history),"Force stage trace failed");
}
void observe_step(const GapPic& pic,std::uint64_t step,double dt) {
    io::require(owns(&pic) && step>initial_step && step-initial_step<=step_budget,"One-ps replay step/owner bound exceeded");
    const auto& particles=pic.particles();const auto& phi=pic.potential();
    io::require(particles.size()<=capacity,"One-ps trace capacity exceeded");
    const std::uint64_t magic=0x5049435354455033ULL,n=particles.size(),nodes=phi.Size(),h=pic.horizontal_transfers_C().size(),v=pic.vertical_transfers_C().size();
    session.scalar(magic);session.scalar(step);session.scalar(dt);session.scalar(pic.diagnostics());
    session.scalar(n);session.scalar(nodes);session.scalar(h);session.scalar(v);
    session.array<GapParticle>(particles);session.array<double>({phi.HostRead(),std::size_t(phi.Size())});
    session.array<double>(pic.horizontal_transfers_C());session.array<double>(pic.vertical_transfers_C());++session.steps;
}
void observe_sources(const std::array<double,4>& spans,const std::array<std::span<const double>,4>& integral,
                     const std::array<std::span<const double>,4>& temperature) {
    io::require(session.owner && !session.closed && !session.failed && std::this_thread::get_id()==session.thread,"Invalid source observer owner/state");
    session.array<double>(spans);
    for(int i=0;i<4;++i){const std::uint64_t size=integral[i].size();io::require(size==temperature[i].size() && size<=4096,"Source trace size guard");session.scalar(size);session.array<double>(integral[i]);session.array<double>(temperature[i]);}
}
void shutdown(const char* original_error) {
    if(!session.owner || session.closed)return;
    session.closed=true;
    if(session.backend)session.cleanup=session.backend->close();
    session.trace.close();session.history.close();
    const auto& r=session.cleanup;const bool okay=!original_error && !session.failed && bool(session.trace) && bool(session.history) && (!session.backend ||
        (r.attempted && r.context_established && r.stream_drained && r.resources_released && !r.resources_retained && r.caller_context_restored && r.primary_reference_released && !r.issue_count));
    std::ofstream file(session.output/"force-cleanup.json");file<<std::boolalpha
        <<"{\"passed\":"<<okay<<",\"mode\":"<<int(session.mode)<<",\"optimized_mode_is_unqualified_candidate\":true,\"deployment_authorized\":false,\"production_resume_authorized\":false"
        <<",\"steps_recorded\":"<<session.steps<<",\"trace_bytes\":"<<session.trace_bytes<<",\"static_uploads\":"<<session.refreshes
        <<",\"attempted\":"<<r.attempted<<",\"context_established\":"<<r.context_established<<",\"stream_drained\":"<<r.stream_drained
        <<",\"resources_released\":"<<r.resources_released<<",\"resources_retained\":"<<r.resources_retained<<",\"caller_context_restored\":"<<r.caller_context_restored
        <<",\"primary_reference_released\":"<<r.primary_reference_released<<",\"issue_count\":"<<r.issue_count<<",\"issues\":[";
    for(std::size_t i=0;i<r.issue_count && i<r.issues.size();++i)file<<(i?",":"")<<"{\"operation\":"<<std::quoted(r.issues[i].operation)<<",\"code\":"<<r.issues[i].code<<'}';
    file<<"],\"original_error\":"<<std::quoted(original_error?original_error:"")<<"}\n";file.close();
    io::require(bool(file) && bool(session.trace) && bool(session.history),"Evolving qualification report/trace close failed");
    if(!original_error)io::require(okay,"Evolving candidate state/cleanup/trace disqualified");
    if(!original_error && session.backend)clean(r);
}
}
