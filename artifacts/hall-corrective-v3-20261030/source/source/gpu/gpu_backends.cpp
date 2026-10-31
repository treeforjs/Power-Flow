// GPU backends of Hall corrective v3 (see gpu_backends.hpp). v3: the motion provider
// runs the v3 engine (hall_gpu_push_v3.cpp: v1 engine with first-hit truncation from
// sink_core.hpp) and assembles survivors and sink records; the field solver is the v2 one.
// v2 description:
//
// Field solver: algorithm, kernels, stopping rule and launch geometry of
// hall-gpu-resident-v1-20261030/field (src/gpu_cg.hpp, src/kernels.hpp,
// factored operator), qualified there as bitwise equal to the sealed apply and
// bitwise reproducible. The driver wrapper of field/ (its own context, stream
// and allocation code) is replaced by the common hall_gpu::Session: one device
// session, one non-blocking ordered stream, pinned staging, drained teardown.
// Device-to-device copies use a copy kernel on the same stream; graphs are not
// used (gated kernels make the result independent of batching and graphs).
//
// Motion provider: hall-gpu-resident-v1-20261030/push MotionEngine, unchanged,
// on the same session; the host resolves refusals in the sealed order and
// assembles the v2 receipt (MotionAssembly) with the sealed record and seal.
#include "gpu_backends.hpp"
#include "common/hall_gpu_common.hpp"
#include "push/hall_adapter.hpp"
#include "hall_gpu_push_v3.hpp"
#include <chrono>
#include <cstring>
#include <map>

namespace hall_corrective_v3::gpu {
namespace {
using Clock=std::chrono::steady_clock;
double since(Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();}
constexpr int TPB=256,NBLK=512,BATCH=4;
struct Ctl {
    double rz,rr,rnorm,alpha,beta,pap,target,bnorm,tnorm,spare;
    std::int32_t it,done,error,maxit,corr,pad0,pad1,pad2;
};
static_assert(sizeof(Ctl)==112,"control block layout");
// Kernels of hall-gpu-resident-v1-20261030/field/src/kernels.hpp (factored form),
// unchanged, plus k_copy replacing the driver-level device-to-device copy.
const char* kernel_source=R"CUDA(
#define TPB 256
#define NBLK 512
struct Ctl {
    double rz, rr, rnorm, alpha, beta, pap, target, bnorm, tnorm, spare;
    int it, done, error, maxit, corr, pad0, pad1, pad2;
};
__device__ __forceinline__ bool stopped(const Ctl* c) { return c->done != 0 || c->error != 0; }
__device__ __forceinline__ double block_tree(double* sh, double v) {
    const int t = threadIdx.x;
    sh[t] = v;
    __syncthreads();
    for (int s = TPB / 2; s > 0; s >>= 1) {
        if (t < s) sh[t] = sh[t] + sh[t + s];
        __syncthreads();
    }
    return sh[0];
}
__device__ __forceinline__ double final_tree(double* sh, const double* part) {
    const int t = threadIdx.x;
    sh[t] = part[t];
    __syncthreads();
    for (int s = NBLK / 2; s > 0; s >>= 1) {
        if (t < s) sh[t] = sh[t] + sh[t + s];
        __syncthreads();
    }
    return sh[0];
}
#define ROWS(i, n) for (int i = blockIdx.x * TPB + threadIdx.x; i < (n); i += NBLK * TPB)

extern "C" __global__ void k_copy(int n, const double* __restrict__ src, double* __restrict__ dst) {
    ROWS(i, n) dst[i] = src[i];
}
extern "C" __global__ void k_dot_self(int n, const double* __restrict__ b, double* part) {
    __shared__ double sh[TPB];
    double s = 0.0;
    ROWS(i, n) s += b[i] * b[i];
    const double v = block_tree(sh, s);
    if (threadIdx.x == 0) part[blockIdx.x] = v;
}
extern "C" __global__ void k_fin_target(const double* part, Ctl* c, int maxit) {
    __shared__ double sh[NBLK];
    const double v = final_tree(sh, part);
    if (threadIdx.x == 0) {
        c->bnorm = sqrt(v);
        c->target = fmax(1e-30, 2e-13 * c->bnorm);
        c->done = 0; c->error = 0; c->it = 0; c->maxit = maxit; c->corr = 0;
    }
}
extern "C" __global__ void k_init(int n, const double* __restrict__ r, const double* __restrict__ diag,
                                  double* __restrict__ z, double* __restrict__ p, double* part_rz, double* part_rr) {
    __shared__ double sh[TPB];
    double s1 = 0.0, s2 = 0.0;
    ROWS(i, n) {
        const double ri = r[i];
        const double zi = ri / diag[i];
        z[i] = zi; p[i] = zi;
        s1 += ri * zi; s2 += ri * ri;
    }
    double v = block_tree(sh, s1);
    if (threadIdx.x == 0) part_rz[blockIdx.x] = v;
    __syncthreads();
    v = block_tree(sh, s2);
    if (threadIdx.x == 0) part_rr[blockIdx.x] = v;
}
extern "C" __global__ void k_fin_init(const double* part_rz, const double* part_rr, Ctl* c, int corr, double* hist_rnorm) {
    __shared__ double sh[NBLK];
    const double rz = final_tree(sh, part_rz);
    __syncthreads();
    const double rr = final_tree(sh, part_rr);
    if (threadIdx.x == 0) {
        c->rz = rz; c->rr = rr; c->rnorm = sqrt(rr); c->it = 0; c->corr = corr; c->error = 0;
        c->done = (0 < c->maxit && c->rnorm > c->target * .2) ? 0 : 1;
        hist_rnorm[corr * (c->maxit + 1)] = c->rnorm;
    }
}
extern "C" __global__ void k_face(int nb, const int* __restrict__ fe, const double* __restrict__ fc, const double* __restrict__ H,
                                  const double* __restrict__ sqrtM, const unsigned char* __restrict__ fixed,
                                  const double* __restrict__ v, double* __restrict__ Bh, const Ctl* c, int gated) {
    if (gated && stopped(c)) return;
    const int b = blockIdx.x * TPB + threadIdx.x;
    if (b >= nb) return;
    double s = 0.0;
    for (int k = 0; k < 4; ++k) {
        const int e = fe[4 * b + k];
        const double E = fixed[e] ? 0.0 : v[e] / sqrtM[e];
        s += fc[4 * b + k] * E;
    }
    Bh[b] = s * H[b];
}
extern "C" __global__ void k_edge(int n, const int* __restrict__ ctp, const int* __restrict__ ctf, const double* __restrict__ ctc,
                                  const double* __restrict__ Bh, const double* __restrict__ sqrtM, const double* __restrict__ a, double q,
                                  const unsigned char* __restrict__ fixed, const double* __restrict__ v, double* __restrict__ Ap,
                                  double* part, const Ctl* c, int gated, int with_dot) {
    if (gated && stopped(c)) return;
    __shared__ double sh[TPB];
    double s = 0.0;
    ROWS(i, n) {
        double y;
        if (fixed[i]) y = v[i];
        else {
            const double E = v[i] / sqrtM[i];
            double K = 0.0;
            const int end = ctp[i + 1];
            for (int j = ctp[i]; j < end; ++j) K += ctc[j] * Bh[ctf[j]];
            y = v[i] + (a[i] * E + q * K) / sqrtM[i];
        }
        Ap[i] = y;
        s += v[i] * y;
    }
    if (with_dot) {
        const double r = block_tree(sh, s);
        if (threadIdx.x == 0) part[blockIdx.x] = r;
    }
}
extern "C" __global__ void k_fin_alpha(const double* part, Ctl* c, double* hist_alpha) {
    if (stopped(c)) return;
    __shared__ double sh[NBLK];
    const double pap = final_tree(sh, part);
    if (threadIdx.x == 0) {
        c->pap = pap;
        if (!(pap > 0.0 && isfinite(pap))) { c->error = 1; return; }
        c->alpha = c->rz / pap;
        hist_alpha[c->corr * c->maxit + c->it] = c->alpha;
    }
}
extern "C" __global__ void k_update(int n, double* __restrict__ x, double* __restrict__ r, const double* __restrict__ p,
                                    const double* __restrict__ Ap, const double* __restrict__ diag, double* __restrict__ z,
                                    double* part_rz, double* part_rr, const Ctl* c) {
    if (stopped(c)) return;
    __shared__ double sh[TPB];
    const double alpha = c->alpha;
    double s1 = 0.0, s2 = 0.0;
    ROWS(i, n) {
        x[i] = x[i] + alpha * p[i];
        const double ri = r[i] - alpha * Ap[i];
        r[i] = ri;
        const double zi = ri / diag[i];
        z[i] = zi;
        s1 += ri * zi; s2 += ri * ri;
    }
    double v = block_tree(sh, s1);
    if (threadIdx.x == 0) part_rz[blockIdx.x] = v;
    __syncthreads();
    v = block_tree(sh, s2);
    if (threadIdx.x == 0) part_rr[blockIdx.x] = v;
}
extern "C" __global__ void k_fin_beta(const double* part_rz, const double* part_rr, Ctl* c, double* hist_beta, double* hist_rnorm) {
    if (stopped(c)) return;
    __shared__ double sh[NBLK];
    const double nrz = final_tree(sh, part_rz);
    __syncthreads();
    const double rr = final_tree(sh, part_rr);
    if (threadIdx.x == 0) {
        const double beta = nrz / c->rz;
        c->beta = beta; c->rz = nrz; c->rr = rr; c->rnorm = sqrt(rr);
        hist_beta[c->corr * c->maxit + c->it] = beta;
        const int it = c->it + 1;
        c->it = it;
        hist_rnorm[c->corr * (c->maxit + 1) + it] = c->rnorm;
        c->done = (it < c->maxit && c->rnorm > c->target * .2) ? 0 : 1;
    }
}
extern "C" __global__ void k_direction(int n, const double* __restrict__ z, double* __restrict__ p, const Ctl* c) {
    if (stopped(c)) return;
    const double beta = c->beta;
    const int i = blockIdx.x * TPB + threadIdx.x;
    if (i < n) p[i] = z[i] + beta * p[i];
}
extern "C" __global__ void k_residual(int n, const double* __restrict__ b, const double* __restrict__ Ax, double* __restrict__ r, double* part) {
    __shared__ double sh[TPB];
    double s = 0.0;
    ROWS(i, n) {
        const double ri = b[i] - Ax[i];
        r[i] = ri;
        s += ri * ri;
    }
    const double v = block_tree(sh, s);
    if (threadIdx.x == 0) part[blockIdx.x] = v;
}
extern "C" __global__ void k_fin_true(const double* part, Ctl* c) {
    __shared__ double sh[NBLK];
    const double v = final_tree(sh, part);
    if (threadIdx.x == 0) c->tnorm = sqrt(v);
}
)CUDA";
} // namespace

// ------------------------------------------------------------------ device
struct Device::Impl {
    hall_gpu::Session session;
    bool closed=false;std::string teardown;
    explicit Impl(const hall_gpu::SessionOptions&o):session(o){}
};
static hall_gpu::SessionOptions session_options(const DeviceOptions&o){
    need(!o.uuid.empty()&&o.nvrtc_library.is_absolute(),"GPU device requires a UUID and an absolute NVRTC library path");
    hall_gpu::SessionOptions s;s.require_uuid=o.uuid;s.nvrtc_library=o.nvrtc_library;s.blocking_sync=true;return s;
}
Device::Device(const DeviceOptions&o):p_(std::make_unique<Impl>(session_options(o))){}
Device::~Device(){(void)close_json();}
hall_gpu::Session& Device::session(){need(!p_->closed,"GPU device already closed");return p_->session;}
std::string Device::device_json()const{return hall_gpu::device_json(p_->session.device());}
std::string Device::memory_json(){return hall_gpu::memory_json(p_->session.memory());}
std::string Device::close_json()noexcept{
    if(!p_)return "null";
    if(!p_->closed){p_->closed=true;try{p_->teardown=hall_gpu::teardown_json(p_->session.close());}catch(...){p_->teardown="{\"error\":\"teardown report\"}";}}
    return p_->teardown;
}

// ------------------------------------------------------------------ field solver
namespace {
class GpuFieldSolver final:public FieldSolveBackend {
    std::shared_ptr<Device> device_;hall_gpu::Session&s_;hall_gpu::Module module_;
    std::map<std::string,hall_gpu::Kernel> k_;
    std::string name_="gpu1-factored-cg",descriptor_;
    // operator structure (per mesh)
    const Mesh* mesh_=nullptr;int n_=0,nb_=0,maxit_=0;
    hall_gpu::DeviceBuffer fe_,fc_,H_,ctp_,ctf_,ctc_,sqrtM_,a_,fixed_,diag_,b_,x_,r_,rt_,z_,p_,Ap_,Bh_,part1_,part2_,ctl_,ha_,hb_,hr_;
    hall_gpu::HostBuffer stage_,ctl_host_,x_host_;
    Vec a_cache_,diag_cache_;
public:
    GpuSolveStats stats;
    explicit GpuFieldSolver(std::shared_ptr<Device>d):device_(std::move(d)),s_(device_->session()){
        hall_gpu::CompileRequest req;req.main={"hall_v2_field_cg.cu",kernel_source};
        module_=s_.compile_and_load(req);
        for(const char*k:{"k_copy","k_dot_self","k_fin_target","k_init","k_fin_init","k_face","k_edge","k_fin_alpha","k_update",
                          "k_fin_beta","k_direction","k_residual","k_fin_true"})k_[k]=s_.kernel(module_,k);
        const auto&rec=module_.record;
        Digest dg("hall-v2-field-solve-backend-v1");dg.text(name_);
        dg.text("deterministic FP64 Jacobi CG on the factored operator (curl rows, transposed curl rows ascending face order, "
                "sqrt(M), a=(.5 dt)(S+Y), q=(.25 dt) dt, PEC mask); row-per-thread products; fixed reduction tree "
                "NBLK=512 x TPB=256 grid-stride, shared-memory halving, one finishing block; no atomics; recursive stop at "
                "0.2*target, device true residual, at most three correction cycles");
        dg.text(bytes_sha256(kernel_source));
        dg.u64(rec.options.size());for(const auto&o:rec.options)dg.text(o);
        dg.i64(rec.target_cc);dg.text(rec.image_kind);dg.i64(rec.nvrtc_major);dg.i64(rec.nvrtc_minor);
        descriptor_=dg.hex();
    }
    const std::string& name()const noexcept override{return name_;}
    const std::string& descriptor_sha256()const noexcept override{return descriptor_;}
    const hall_gpu::CompileRecord& compile_record()const noexcept{return module_.record;}
    void launch(const char*k,unsigned grid,unsigned block,std::initializer_list<void*> args){
        std::vector<void*>a(args);s_.launch(k_.at(k),{grid},{block},0,a.data());}
    template<class T>void upload_static(hall_gpu::DeviceBuffer&dst,const std::vector<T>&v,const char*tag){
        const std::size_t bytes=std::max<std::size_t>(v.size()*sizeof(T),8);dst=s_.alloc(bytes,tag);
        if(v.empty())return;need(v.size()*sizeof(T)<=stage_.bytes,"GPU staging buffer too small");
        std::memcpy(stage_.pointer,v.data(),v.size()*sizeof(T));s_.upload(dst,0,stage_,0,v.size()*sizeof(T));s_.synchronize("static operator upload");
    }
    void build(const Mesh&g,int maxit){
        need(mesh_==nullptr,"GPU field solver is bound to one mesh");mesh_=&g;n_=g.ne;nb_=g.nb;maxit_=maxit;
        std::vector<std::int32_t>fe(4*std::size_t(g.nb)),ctp(g.ne+1,0),ctf,count(g.ne+1,0);Vec fc(4*std::size_t(g.nb)),H(g.nb),ctc,sqrtM(g.ne);
        std::vector<std::uint8_t>fixed(g.ne);
        for(int b=0;b<g.nb;++b){need(g.C[b].count==4,"axial curl rows have four entries");H[b]=g.H[b];
            for(int k=0;k<4;++k){fe[4*b+k]=g.C[b].index[k];fc[4*b+k]=g.C[b].coefficient[k];++count[g.C[b].index[k]+1];}}
        for(int e=0;e<g.ne;++e){ctp[e+1]=ctp[e]+count[e+1];sqrtM[e]=std::sqrt(g.M[e]);fixed[e]=g.fixed_E[e]?1:0;}
        ctf.resize(ctp[g.ne]);ctc.resize(ctp[g.ne]);std::vector<std::int32_t>fill(ctp.begin(),ctp.end()-1);
        for(int b=0;b<g.nb;++b)for(int k=0;k<4;++k){const int e=g.C[b].index[k];ctf[fill[e]]=b;ctc[fill[e]]=g.C[b].coefficient[k];++fill[e];}
        std::size_t biggest=std::max({fe.size()*4,fc.size()*8,ctf.size()*4,ctc.size()*8,std::size_t(g.ne)*8,std::size_t(g.nb)*8});
        stage_=s_.alloc_pinned(biggest,"field staging");ctl_host_=s_.alloc_pinned(sizeof(Ctl),"field control");x_host_=s_.alloc_pinned(std::size_t(g.ne)*8,"field solution");
        upload_static(fe_,fe,"face edges");upload_static(fc_,fc,"face coefficients");upload_static(H_,H,"H");
        upload_static(ctp_,ctp,"transposed curl ptr");upload_static(ctf_,ctf,"transposed curl faces");upload_static(ctc_,ctc,"transposed curl coefficients");
        upload_static(sqrtM_,sqrtM,"sqrt M");upload_static(fixed_,fixed,"PEC mask");
        const std::size_t vb=std::size_t(g.ne)*8;
        for(auto*b:{&a_,&diag_,&b_,&x_,&r_,&rt_,&z_,&p_,&Ap_})*b=s_.alloc(vb,"field vector");
        Bh_=s_.alloc(std::size_t(g.nb)*8,"face work");part1_=s_.alloc(NBLK*8,"partials");part2_=s_.alloc(NBLK*8,"partials");
        ctl_=s_.alloc(sizeof(Ctl),"control");
        ha_=s_.alloc(std::size_t(3)*maxit*8,"alpha history");hb_=s_.alloc(std::size_t(3)*maxit*8,"beta history");hr_=s_.alloc(std::size_t(3)*(maxit+1)*8,"residual history");
        s_.zero(ha_,0,ha_.bytes);s_.zero(hb_,0,hb_.bytes);s_.zero(hr_,0,hr_.bytes);s_.synchronize("field solver ready");
    }
    void upload_vector(hall_gpu::DeviceBuffer&dst,const Vec&v){std::memcpy(stage_.pointer,v.data(),v.size()*8);s_.upload(dst,0,stage_,0,v.size()*8);s_.synchronize("field vector upload");}
    void apply(std::uint64_t v,std::uint64_t out,int gated,int with_dot){
        int n=n_,nb=nb_;double q=q_;
        auto fe=fe_.address,fc=fc_.address,H=H_.address,sq=sqrtM_.address,fx=fixed_.address,Bh=Bh_.address,ctl=ctl_.address,
             ctp=ctp_.address,ctf=ctf_.address,ctc=ctc_.address,a=a_.address,part=part1_.address;
        launch("k_face",unsigned((nb_+TPB-1)/TPB),TPB,{&nb,&fe,&fc,&H,&sq,&fx,&v,&Bh,&ctl,&gated});
        launch("k_edge",NBLK,TPB,{&n,&ctp,&ctf,&ctc,&Bh,&sq,&a,&q,&fx,&v,&out,&part,&ctl,&gated,&with_dot});
    }
    double q_=0;
    void read(){s_.download(ctl_host_,0,ctl_,0,sizeof(Ctl));s_.synchronize("field control readback");
        if(ctl_host_.as<Ctl>()->error)throw std::runtime_error("non-SPD field operator (device CG)");}
    void iteration(){
        int n=n_;auto x=x_.address,r=r_.address,p=p_.address,Ap=Ap_.address,dg=diag_.address,z=z_.address,p1=part1_.address,p2=part2_.address,
            ctl=ctl_.address,ha=ha_.address,hb=hb_.address,hr=hr_.address;
        apply(p,Ap,1,1);
        launch("k_fin_alpha",1,NBLK,{&p1,&ctl,&ha});
        launch("k_update",NBLK,TPB,{&n,&x,&r,&p,&Ap,&dg,&z,&p1,&p2,&ctl});
        launch("k_fin_beta",1,NBLK,{&p1,&p2,&ctl,&hb,&hr});
        launch("k_direction",unsigned((n_+TPB-1)/TPB),TPB,{&n,&z,&p,&ctl});
    }
    void tail(){int n=n_;auto x=x_.address,Ap=Ap_.address,b=b_.address,rt=rt_.address,p1=part1_.address,ctl=ctl_.address;
        apply(x,Ap,0,0);launch("k_residual",NBLK,TPB,{&n,&b,&Ap,&rt,&p1});launch("k_fin_true",1,NBLK,{&p1,&ctl});}
    void init_from(std::uint64_t src,int corr){int n=n_;auto r=r_.address,dg=diag_.address,z=z_.address,p=p_.address,p1=part1_.address,p2=part2_.address,
        ctl=ctl_.address,hr=hr_.address;
        launch("k_copy",NBLK,TPB,{&n,&src,&r});
        launch("k_init",NBLK,TPB,{&n,&r,&dg,&z,&p,&p1,&p2});
        launch("k_fin_init",1,NBLK,{&p1,&p2,&ctl,&corr,&hr});}
    Vec solve(const Apply&,const Vec&rhs,int maxit,const SolveContext&c,SolveRecord&rec)override{
        const auto t0=Clock::now();
        if(!mesh_)build(c.mesh,maxit);
        need(&c.mesh==mesh_&&maxit==maxit_&&rhs.size()==std::size_t(n_)&&c.diagonal.size()==rhs.size()&&c.conduction.size()==rhs.size(),
             "GPU field solver operator binding");
        Vec a(n_);for(int e=0;e<n_;++e)a[e]=.5*c.dt*(c.conduction[e]+c.mesh.Y[e]);
        q_=.25*c.dt*c.dt;
        if(a!=a_cache_){upload_vector(a_,a);a_cache_=a;}
        if(c.diagonal!=diag_cache_){upload_vector(diag_,c.diagonal);diag_cache_=c.diagonal;}
        upload_vector(b_,rhs);
        int n=n_,mi=maxit_,corr0=0;auto b=b_.address,p1=part1_.address,ctl=ctl_.address;
        launch("k_dot_self",NBLK,TPB,{&n,&b,&p1});
        launch("k_fin_target",1,NBLK,{&p1,&ctl,&mi});
        s_.zero(x_,0,x_.bytes);
        init_from(b_.address,corr0);
        for(int k=0;k<BATCH;++k)iteration();
        tail();read();
        rec.backend=name_;rec.iterations=0;rec.corrections=0;
        std::vector<int>its;
        for(int corr=0;corr<3;++corr){
            rec.corrections=corr+1;bool needs_tail=false;
            if(corr>0){init_from(rt_.address,corr);ctl_host_.as<Ctl>()->done=0;needs_tail=true;}
            for(int guard=0;!ctl_host_.as<Ctl>()->done;++guard){
                need(guard<=maxit_/BATCH+2,"device CG did not stop within maxit");
                for(int k=0;k<BATCH;++k)iteration();read();needs_tail=true;}
            if(needs_tail){tail();read();}
            const Ctl&C=*ctl_host_.as<Ctl>();its.push_back(C.it);rec.iterations+=C.it;rec.backend_residual=C.tnorm;
            if(C.tnorm<=C.target){
                s_.download(x_host_,0,x_,0,std::size_t(n_)*8);s_.synchronize("field solution download");
                Vec x(x_host_.as<double>(),x_host_.as<double>()+n_);
                ++stats.solves;stats.total_s+=since(t0);stats.last_iterations=its;stats.last_device_residual=C.tnorm;
                return x;}
        }
        throw std::runtime_error("field solve independent residual failed (device CG)");
    }
};
} // namespace
std::shared_ptr<FieldSolveBackend> make_gpu_field_solver(std::shared_ptr<Device>d){return std::make_shared<GpuFieldSolver>(std::move(d));}
std::string gpu_field_solver_compile_json(const FieldSolveBackend&b){
    const auto*g=dynamic_cast<const GpuFieldSolver*>(&b);need(g!=nullptr,"not the GPU field solver");return hall_gpu::compile_json(g->compile_record());}
GpuSolveStats gpu_field_solver_stats(const FieldSolveBackend&b){
    const auto*g=dynamic_cast<const GpuFieldSolver*>(&b);need(g!=nullptr,"not the GPU field solver");return g->stats;}

// ------------------------------------------------------------------ motion provider
namespace hp=hall_gpu::push;namespace ha=hall_gpu::push::hall_adapter;namespace h3=hall_gpu::push_v3;
struct GpuMotionProvider::Impl {
    std::shared_ptr<Device> device;const Mesh&mesh;std::string mesh_sha256;std::unique_ptr<h3::MotionEngineV3> engine;std::string name="gpu-push-v3-sinks";
    const MotionContext* context=nullptr;std::span<const native::ArenaRateEntry> census;ha::OrderedCensus ordered;GpuMotionStats stats;
    bool over_cap=false;
    Impl(std::shared_ptr<Device>d,const MotionContext&c):device(std::move(d)),mesh(*c.mesh),mesh_sha256(c.mesh_sha256){}
};
GpuMotionProvider::GpuMotionProvider(std::shared_ptr<Device>d,const MotionContext&c,std::size_t capacity):p_(std::make_unique<Impl>(std::move(d),c)){
    need(c.mesh!=nullptr,"GPU motion provider requires a motion context");
    hp::EngineOptions eo;eo.capacity=std::max<std::size_t>(capacity,1);
    h3::SinkSupport sup;for(int k=0;k<5;++k){for(int d=0;d<3;++d){sup.lo[k][d]=c.sink_lo[k][d];sup.hi[k][d]=c.sink_hi[k][d];}sup.supported[k]=c.sink_supported[k];}
    sup.aligned=c.sink_support_aligned;
    p_->engine=std::make_unique<h3::MotionEngineV3>(p_->device->session(),ha::mesh_tables(p_->mesh),ha::geometry(c.options),sup,eo);
}
GpuMotionProvider::~GpuMotionProvider()=default;
const std::string& GpuMotionProvider::name()const noexcept{return p_->name;}
const GpuMotionStats& GpuMotionProvider::stats()const noexcept{return p_->stats;}
std::string GpuMotionProvider::compile_json()const{return hall_gpu::compile_json(p_->engine->compile_record());}
void GpuMotionProvider::begin(const MotionContext&c,std::span<const native::ArenaRateEntry> census){
    const auto t=Clock::now();
    need(c.mesh==&p_->mesh&&c.mesh_sha256==p_->mesh_sha256,"GPU motion provider bound to another mesh or geometry");
    p_->context=&c;p_->census=census;
    // Census order/identity and the Count2 packet comparison on the host, as the
    // sealed provider does; carriers then stay on the device for every iterate.
    // Above the configured Hall cap nothing is uploaded: motion() then refuses in
    // MotionAssembly::admit with the sealed text at the sealed position.
    p_->over_cap=census.size()>c.maximum_carriers;
    if(p_->over_cap){p_->stats.upload_carriers_s+=since(t);return;}
    p_->ordered=ha::order_census(census,false);
    p_->engine->set_carriers(p_->ordered.view());
    p_->stats.upload_carriers_s+=since(t);
}
void GpuMotionProvider::end()noexcept{p_->context=nullptr;p_->census={};}
hall_driver::MotionReceipt GpuMotionProvider::motion(const hall_driver::MotionBinding&binding,const Vec&E,const Vec&B,FreeCharge&after){
    auto&I=*p_;need(I.context!=nullptr,"motion provider used outside a proposal");
    const auto t0=Clock::now();
    auto admitted=MotionAssembly::admit(*I.context,I.census,binding,{},E,B);
    need(!I.over_cap,"motion census exceeds the configured Hall carrier cap");
    const auto t1=Clock::now();I.stats.admit_s+=std::chrono::duration<double>(t1-t0).count();
    const auto fl=MotionAssembly::flights(*I.context,admitted,binding,newborns_);
    if(fl.count)I.engine->set_flights(fl.flight,fl.face_body,fl.newborn);else I.engine->set_flights({},{},{});
    I.engine->set_field(E,B);auto v=I.engine->run(binding.dt_s,false);
    const auto t2=Clock::now();I.stats.device_run_s+=std::chrono::duration<double>(t2-t1).count();
    for(const auto&ph:v.gpu_phases)if(ph.name=="gpu_total_run")I.stats.last_gpu_total_ms=ph.ms;
    I.stats.last_host_fold_ms=v.host_fold_ms;I.stats.last_edge_records=v.edge_records;
    // Refusal order of the sealed provider: the host Count2 packet comparison sits
    // inside the first per-carrier loop after the interior check and before the
    // plasma operand and charge deposit of the same carrier.
    const long rp=I.ordered.packet_failure_rank;
    if(rp>=0){bool packet_first=!v.refusal.refused||v.refusal.stage!="census"||v.refusal.rank>rp||
        (v.refusal.rank==rp&&(v.refusal.code==hgp::unresolved_plasma_operand||v.refusal.code==hgp::charge_outside_mesh));
        if(packet_first)throw std::runtime_error("motion divergent Count2 packet carriers");}
    if(v.refusal.refused)throw std::runtime_error(v.refusal.message);
    const std::size_t n=v.n;need(n==admitted.rows.size(),"device motion carrier count");
    MotionArrays m;m.endpoints.reserve(n);m.paths.reserve(n);m.birth_scope=birth_scope_;
    const native::PicVector shift{-I.context->options.width_m/2,-I.context->options.gap_m/2,0};
    for(std::size_t i=0;i<n;++i){const auto&row=*admitted.rows[i];auto e=row.native;
        if(fl.newborn[i]){NewbornRecord nb;nb.native_id=e.id;nb.birth_time_s=fl.tau[i];nb.flight_s=fl.flight[i];nb.origin=fl.face_body[i];nb.axis=fl.face_axis[i];
            nb.start_m=row.native.position_m;nb.packet_id=row.packet_id;nb.state=row.state;nb.charge_C=e.charge_C;nb.weight=e.weight;
            nb.birth_weighted_kinetic_J=e.weight*native::pic_kinetic_energy_J(row.native);m.births.push_back(std::move(nb));}
        e.position_m={v.end_x[i],v.end_y[i],v.end_z[i]};e.proper_velocity_m_s={v.end_ux[i],v.end_uy[i],v.end_uz[i]};
        if(v.sink_kind[i]>=0){
            // Sink record from the device outputs; host values with the CPU expressions.
            SinkRecord k;k.native_id=e.id;k.packet_id=row.packet_id;k.state=row.state;k.boundary=v.sink_kind[i];
            k.collected=k.boundary==0||k.boundary==1||k.boundary==3;k.kind=k.collected?1:2;k.conductor=k.boundary==3?2:(k.collected?k.boundary:-1);
            k.axis=v.sink_axis[i];k.plane_m=v.sink_plane[i];k.fraction=v.sink_fraction[i];
            k.segment_m={v.sink_segment_x[i],v.sink_segment_y[i],v.sink_segment_z[i]};k.hit_point_m=e.position_m;
            for(int d=0;d<3;++d)k.field_point_m[d]=e.position_m[d]+shift[d];
            k.incident_proper_velocity_m_s=e.proper_velocity_m_s;k.flight_s=v.sink_flight[i];k.discarded_flight_s=v.sink_remaining[i];
            const double K1=native::pic_kinetic_energy_J(e),K0=native::pic_kinetic_energy_J(row.native);
            k.charge_C=e.charge_C;k.mass_kg=e.mass_kg;k.weight=e.weight;k.kinetic_J=K1;k.weighted_kinetic_J=e.weight*K1;k.initial_weighted_kinetic_J=e.weight*K0;
            k.gyro=v.gyro[i];k.displacement_cells=v.dcells[i];
            m.sinks.push_back(std::move(k));
            m.paths.push_back({e.id,{v.start_x[i],v.start_y[i],v.start_z[i]},{v.disp_x[i],v.disp_y[i],v.disp_z[i]},v.sink_flight[i]});
            continue;}
        m.endpoints.push_back(e);
        m.paths.push_back({e.id,{v.start_x[i],v.start_y[i],v.start_z[i]},{v.disp_x[i],v.disp_y[i],v.disp_z[i]},fl.flight[i]});}
    m.current_A_m.assign(v.current_A_m.begin(),v.current_A_m.end());
    m.initial_charge_C.assign(v.initial_charge_C.begin(),v.initial_charge_C.end());
    m.endpoint_charge_C.assign(v.endpoint_charge_C.begin(),v.endpoint_charge_C.end());
    const auto&b=v.budget;m.budget={b.kinetic_change_J,b.field_work_J,b.work_defect_J,b.continuity_C,b.plasma_dt,b.maximum_gyro,b.maximum_displacement_cells};
    // charge_of_endpoints is evaluated after the motion succeeded; same order here.
    if(v.charge_refusal.refused)throw std::runtime_error(v.charge_refusal.message);
    need(v.free_net_hi.size()==std::size_t(I.mesh.nn)&&v.free_net_lo.size()==std::size_t(I.mesh.nn)&&v.free_gross_C.size()==std::size_t(I.mesh.nn),
         "device endpoint free-charge dimensions");
    after.net_C.resize(I.mesh.nn);after.gross_C.assign(v.free_gross_C.begin(),v.free_gross_C.end());
    for(std::size_t k=0;k<after.net_C.size();++k){after.net_C[k].hi=v.free_net_hi[k];after.net_C[k].lo=v.free_net_lo[k];}
    auto r=MotionAssembly::assemble(*I.context,admitted,binding,std::move(m));
    I.stats.assemble_s+=since(t2);I.stats.total_s+=since(t0);++I.stats.runs;
    return r;
}
} // namespace hall_corrective_v3::gpu
