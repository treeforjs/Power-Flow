#pragma once
// Stage 1 test fixture for coordinator v5. Defines the test-only access
// structures that the v5 sources declare but never define.
#include "genesis_v5.hpp"
#include "gpu/gpu_backends.hpp"
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <exception>
#include <execinfo.h>
#include <sys/wait.h>
#include <unistd.h>

extern "C" {
extern int stage1_injection_report_fd;
long long stage1_allocation_count() noexcept;
void stage1_set_allocation_countdown(long long) noexcept;
long long stage1_allocation_countdown() noexcept;
long long stage1_injected_allocation_failures() noexcept;
}

namespace raikou::native_water::v5 {
struct TrialTestAccess {
  using Trial = ClosedParentTrialV5;
  static void set_can_publish_fault(Trial& t,int k){t.can_publish_fault_=k;}
  static int refusal(const Trial& t){return t.refusal_point_;}
  static void set_publish_probe(Trial& t,void (*p)(int)){t.publish_probe_=p;}
  static void retain(Trial& t,std::unique_ptr<m::PreparedArena>* n,
                     std::unique_ptr<hc::JointProposal>* f,PreparedParentState* r){
    t.retain_native_=n;t.retain_field_heat_=f;t.retain_random_=r;
  }
  static std::unique_ptr<m::PreparedArena>& native(Trial& t){return t.native_;}
  static std::unique_ptr<hc::JointProposal>& field_heat(Trial& t){return t.field_heat_;}
  static const hc::JointProposal& proposal(const Trial& t){return *t.field_heat_;}
  static PreparedParentState& random(Trial& t){return t.random_;}
  static ClosedParentInputsV5& inputs(Trial& t){return t.inputs_;}
  static std::unique_ptr<m::PrivateFullTimeStep>& private_step(Trial& t){return t.private_step_;}
  static const std::vector<m::ArenaRateEntry>& census(const Trial& t){return t.post_first_census_;}
  static const common_commit::ExpectedStep& expected(const Trial& t){return t.expected_;}
  static bool consumed(const Trial& t){return t.consumed_;}
  static Trial shell(){return Trial();}
  static void copy_contract(Trial& to,const Trial& from){
    to.inputs_=from.inputs_;to.expected_=from.expected_;to.post_first_census_=from.post_first_census_;
  }
};
struct GenesisTestAccess {
  static hc::ChargeLedger declared_ledger(const hc::Model& model,const axial::State& seed){
    return JointGenesisV5::declared_ledger(model,seed,"stage1 genesis: DC conductor/PEC induced charge of the C19 first-knot seed",
                                           "stage1 genesis: DC exterior termination charge of the C19 first-knot seed");}
  using PacketHook=JointGenesisV5::PacketHook;
  using LedgerHook=JointGenesisV5::LedgerHook;
  static JointOwnersV5 admit(const GenesisRequestV5& rq,const w::ParentBridge& bridge,
                             PacketHook packets,LedgerHook ledger,int* reached){
    return JointGenesisV5::admit_with_hooks(rq,bridge,packets,ledger,reached);
  }
};
} // namespace raikou::native_water::v5

namespace stage1 {
namespace nw = raikou::native_water;
namespace v5 = raikou::native_water::v5;
namespace hg = ::hall_corrective_v2::gpu;
namespace m = raikou::mfem_backend;
namespace h = ::hall_driver;
namespace hc = ::hall_corrective_v2;
namespace pc = raikou::pic_collisions;
namespace w = raikou::pic_water;
namespace fs = std::filesystem;
using J = nw::Json;
using TA = v5::TrialTestAccess;
using steady = std::chrono::steady_clock;

inline constexpr double mn = 2.9915e-26;
inline constexpr double dt = 1e-14;
extern int forbidden_solve_calls, forbidden_advance_calls;

inline double since(steady::time_point a){return std::chrono::duration<double>(steady::now()-a).count();}

// Exactly hall_corrective test_owner.cpp geometry(): the native frame is the
// Hall field box translated by (+w/2, +g/2, 0), z in (-depth/2, +depth/2).
inline m::GapPicOptions hall_geometry(const axial::Mesh& g){
  m::GapPicOptions o;o.nx=o.ny=g.nx;o.width_m=o.gap_m=g.side;o.depth_m=g.depth;
  o.finite_foils=o.finite_axial_extent=true;o.short_z_m=-g.depth/2;o.feed_z_m=g.depth/2;
  o.lower_foil={g.side/2-.0025,g.side/2+.0025,g.side/2-.0027,g.side/2-.0025};
  o.upper_foil={g.side/2-.0025,g.side/2+.0025,g.side/2+.0025,g.side/2+.0027};
  return o;
}
// States and channel exactly as V3 test_native_parent_subflow.cpp.
inline std::vector<m::ArenaStateSpec> states(){
  return {{"electron",0,0,-1,w::me,0,0},{"H2O",2,1,0,mn,0,0},
          {"H2O+",2,1,1,mn-w::me+w::threshold_eV*pc::qe/(pc::c*pc::c),0,w::threshold_eV*pc::qe}};
}
inline std::vector<m::ArenaChannel> channels(const nw::Binding& b){
  return {{w::parent_channel,"electron","H2O",{"electron","electron","H2O+"},
           b.dataset_sha256,b.kernel_sha256,b.kernel_binding_file_sha256,w::parent_model}};
}
inline m::PicVector electron_u(double eV){
  const double x=eV*pc::qe/(w::me*pc::c*pc::c);return {pc::c*std::sqrt(x*(x+2)),0,0};
}
// Plan position (native frame): Hall field origin, on the k=3 Hall node plane
// z=0. Hall corrective v2 admits charge-conserving motion there (deposit-
// rounding term of the endpoint guard). centre_stage1 is the Stage 1 position
// (1 mm above the plane) used only for the bit-for-bit comparison with Stage 1.
inline constexpr m::PicVector centre_plan{.0045,.0045,0};
inline constexpr m::PicVector centre_stage1{.0045,.0045,.001};
inline constexpr m::PicVector centre=centre_plan;
// P1: one electron and one H2O+ packet co-located, counts {1,0}: every node
// starts exactly neutral. 100 eV electron, thermal (about 300 K) ion.
inline std::vector<m::FullTimeAuthoredPacket> packets_p1(m::PicVector x=centre){
  return {{"e0","electron","stage1/authored-electron-100eV",{1,0},x,electron_u(100)},
          {"i0","H2O+","stage1/authored-thermal-ion",{1,0},x,{0,372,0}}};
}
// P2: P1 plus a high-weight H2O packet in the same collision cell.
inline std::vector<m::FullTimeAuthoredPacket> packets_p2(double water_weight,m::PicVector x=centre){
  auto p=packets_p1(x);
  p.push_back({"w0","H2O","stage1/authored-water-reservoir",{water_weight,0},x,{0,0,0}});
  return p;
}

// Backend selection (configuration): field-solve backend cpu-pcg | gpu1 (bound into
// the Hall configuration identity) and motion provider cpu | gpu (identical
// records; execution envelope only). One GPU device session for both.
struct Backends {
  std::string field="cpu",push="cpu",uuid,nvrtc;
  std::shared_ptr<hg::Device> device;
  std::shared_ptr<hc::FieldSolveBackend> solver;
  std::shared_ptr<hg::Device> gpu(){if(!device)device=std::make_shared<hg::Device>(hg::DeviceOptions{uuid,nvrtc});return device;}
};
inline Backends backends;
// Carrier caps (configuration): Hall CarrierLimits and native v5 FullTimeConfiguration.
struct Caps {std::size_t carriers=100000;std::uint64_t hall_budget=1ull<<30,native_budget=2ull<<30;};
inline Caps caps;
struct Inputs {
  std::shared_ptr<const hc::Model> model;
  std::shared_ptr<const w::Data> data;
  std::shared_ptr<const w::ParentBridge> bridge;
  nw::Binding binding;
};
inline Inputs load_inputs(const std::string& trace,const std::string& data,const std::string& binding){
  Inputs in;
  if(backends.field=="gpu"){if(!backends.solver)backends.solver=hg::make_gpu_field_solver(backends.gpu());}
  else if(backends.field=="cpu")backends.solver=hc::make_cpu_pcg();
  else throw std::runtime_error("unknown field backend");
  in.model=std::make_shared<const hc::Model>(hc::MeshSpec{},trace,backends.solver,hc::CarrierLimits{caps.carriers,caps.hall_budget});
  std::ifstream f(binding);J b;f>>b;
  if(!f)throw std::runtime_error("binding read");
  in.binding={b.at("dataset_sha256"),b.at("kernel_sha256"),b.at("kernel_binding_file_sha256")};
  in.data=std::make_shared<const w::Data>(w::Data::read(data));
  in.bridge=std::make_shared<const w::ParentBridge>(in.data,w::parent_model,in.binding.dataset_sha256,in.binding.kernel_sha256,mn);
  return in;
}
inline v5::GenesisRequestV5 request(const Inputs& in,std::vector<m::FullTimeAuthoredPacket> packets){
  v5::GenesisRequestV5 rq;
  rq.model=in.model;rq.options=hall_geometry(in.model->mesh());rq.seed=in.model->dc_seed();
  rq.states=states();rq.channels=channels(in.binding);rq.packets=std::move(packets);
  rq.grid_counts={45,45,5};rq.binding=in.binding;rq.limits={};rq.clock_seed=123;rq.product_seed=456;
  rq.origin_provenance="stage1/authored-neutral-resident-population-at-first-C19-knot";
  rq.conductor_declaration="stage1 genesis: DC conductor/PEC induced charge of the C19 first-knot seed";
  rq.exterior_declaration="stage1 genesis: DC exterior termination charge of the C19 first-knot seed";
  rq.maximum_native_records=4096;rq.maximum_carriers=caps.carriers;rq.native_memory_budget_bytes=caps.native_budget;
  // Above 1e5 carriers the parent chemistry carrier limit follows the configured cap
  // (it is part of the chemistry configuration identity); default runs are unchanged.
  if(caps.carriers>rq.limits.carriers)rq.limits.carriers=caps.carriers;
  return rq;
}
inline void attach_motion(v5::JointOwnersV5& f,std::size_t capacity){
  if(backends.push=="gpu")f.motion=std::make_shared<hg::GpuMotionProvider>(backends.gpu(),f.owner->motion_context(),capacity);
  else if(backends.push!="cpu")throw std::runtime_error("unknown motion provider");
}
inline v5::JointOwnersV5 genesis(const Inputs& in,std::vector<m::FullTimeAuthoredPacket> packets,std::size_t capacity=1024){
  auto f=v5::JointGenesisV5::admit(request(in,std::move(packets)),*in.bridge);
  attach_motion(f,capacity);return f;
}
inline v5::ClosedParentTrialV5 prepare(v5::JointOwnersV5& f,int fault=0,double step=dt){
  return v5::prepare_closed_parent_joint_v5(*f.arena,*f.chemistry,*f.owner,*f.motion,f.inputs,step,fault);
}
inline bool can_publish(const v5::JointOwnersV5& f,const v5::ClosedParentTrialV5& t){
  return v5::can_publish_closed_parent_joint_v5(*f.arena,*f.chemistry,*f.owner,t);
}
inline bool publish(v5::JointOwnersV5& f,v5::ClosedParentTrialV5& t){
  return v5::try_publish_closed_parent_joint_v5(*f.arena,*f.chemistry,*f.owner,t);
}

inline std::string read_bytes(const fs::path& p){
  std::ifstream f(p,std::ios::binary);
  if(!f)throw std::runtime_error("read "+p.string());
  return {std::istreambuf_iterator<char>(f),{}};
}
inline void write_text(const fs::path& p,const std::string& s){
  std::ofstream f(p,std::ios::binary);f<<s;f.close();
  if(!f)throw std::runtime_error("write "+p.string());
}
inline J boundary_json(const m::FullTimeAcceptedBoundary& b){
  return {{"accepted_step",b.accepted_step},{"accepted_arena_revision",b.accepted_arena_revision},
    {"time_s",b.time_s},{"accepted_arena_sha256",b.accepted_arena_sha256},{"field_sha256",b.field_sha256},
    {"heat_sha256",b.heat_sha256},{"source_sha256",b.source_sha256},{"rng_sha256",b.rng_sha256},
    {"mesh_sha256",b.mesh_sha256},{"material_sha256",b.material_sha256},{"frame_sha256",b.frame_sha256},
    {"common_clock_sha256",b.common_clock_sha256},{"authored_full_time_genesis",b.authored_full_time_genesis}};
}

// "Byte-identical" evidence: native save() bytes, state_prefix() and boundary
// JSON; Hall save() bytes and field/heat hashes; chemistry checkpoint dump.
inline long snapshot_counter=0;
inline J snapshot(const v5::JointOwnersV5& f,const fs::path& scratch,const std::string& tag,
                  const fs::path& keep=fs::path()){
  const auto unique=tag+"-"+std::to_string(::getpid())+"-"+std::to_string(++snapshot_counter);
  J s;
  const auto native=scratch/(unique+"-native.json");
  try{
    f.arena->save(native.string());const auto bytes=read_bytes(native);
    s["native_save_sha256"]=nw::sha(bytes);s["native_save_bytes"]=bytes.size();
    if(!keep.empty())fs::copy_file(native,keep/"native-save.json",fs::copy_options::overwrite_existing);
  }catch(const std::exception& e){s["native_save_error"]=e.what();}
  std::error_code ignored;fs::remove(native,ignored);
  s["native_state_prefix"]=f.arena->state_prefix();
  s["native_boundary"]=boundary_json(*f.arena->full_time_boundary());
  const auto hall=scratch/(unique+"-hall");
  try{
    f.owner->save(hall);J files=J::object();
    std::vector<fs::path> paths;
    for(const auto& e:fs::recursive_directory_iterator(hall))if(e.is_regular_file())paths.push_back(e.path());
    std::sort(paths.begin(),paths.end());
    for(const auto& p:paths)files[fs::relative(p,hall).string()]=nw::sha(read_bytes(p));
    s["hall_save_sha256"]=files;
    if(!keep.empty())fs::copy(hall,keep/"hall-save",fs::copy_options::recursive|fs::copy_options::overwrite_existing);
  }catch(const std::exception& e){s["hall_save_error"]=e.what();}
  fs::remove_all(hall,ignored);
  s["hall_field_sha256"]=f.owner->field_sha256();s["hall_heat_sha256"]=f.owner->heat_sha256();
  s["hall_time_s"]=f.owner->time_s();s["hall_revision"]=f.owner->revision();
  s["chemistry_checkpoint"]=f.chemistry->checkpoint().dump();
  const auto digest=nw::sha(s.dump());s["digest"]=digest;
  return s;
}

struct Results {
  J cases=J::array();
  void record(const std::string& test,const std::string& name,bool passed,J evidence=J::object()){
    cases.push_back({{"test",test},{"name",name},{"passed",passed},{"evidence",std::move(evidence)}});
    std::cout<<test<<" | "<<name<<": "<<(passed?"pass":"FAIL")<<std::endl;
  }
  J summary()const{
    J by=J::object();
    for(const auto& c:cases){
      auto& t=by[c.at("test").get<std::string>()];
      if(t.is_null())t={{"checks",0},{"passed",0},{"failed",0}};
      t["checks"]=t["checks"].get<int>()+1;
      if(c.at("passed").get<bool>())t["passed"]=t["passed"].get<int>()+1;
      else t["failed"]=t["failed"].get<int>()+1;
    }
    return by;
  }
};
// Runs f expecting an exception; returns the message (empty if none thrown).
template<class F>std::string refusal_message(F&& f){
  try{f();}catch(const std::exception& e){return std::string("refused: ")+e.what();}
  catch(...){return "refused: non-standard exception";}
  return std::string();
}

// In a child, std::terminate writes a marker and the symbolized backtrace to
// the report pipe and exits with status 86 (backtrace() is preloaded so the
// handler itself does not need operator new).
inline int child_report_fd=-1;
[[noreturn]] inline void report_terminate() noexcept{
  static const char marker[]="\n@@terminate-backtrace@@\n";
  void* frames[96];const int count=::backtrace(frames,96);
  if(child_report_fd>=0){
    if(::write(child_report_fd,marker,sizeof marker-1)<0){}
    ::backtrace_symbols_fd(frames,count,child_report_fd);
  }
  ::_exit(86);
}
struct Child {int code=-1,signal=0;std::string out;};
inline Child in_child(const std::function<std::string()>& body){
  std::cout.flush();std::cerr.flush();std::fflush(nullptr);
  int fd[2];
  if(::pipe(fd)!=0)throw std::runtime_error("pipe");
  const pid_t pid=::fork();
  if(pid<0)throw std::runtime_error("fork");
  if(pid==0){
    ::close(fd[0]);int code=0;std::string out;
    {void* warm[4];(void)::backtrace(warm,4);}
    child_report_fd=fd[1];std::set_terminate(report_terminate);
    try{out=body();}catch(const std::exception& e){out=std::string("exception: ")+e.what();code=3;}
    catch(...){out="non-standard exception";code=4;}
    const char* p=out.data();std::size_t n=out.size();
    while(n){const auto k=::write(fd[1],p,n);if(k<=0)break;p+=k;n-=static_cast<std::size_t>(k);}
    ::close(fd[1]);std::fflush(nullptr);::_exit(code);
  }
  ::close(fd[1]);Child c;char buffer[65536];
  for(;;){const auto k=::read(fd[0],buffer,sizeof buffer);if(k<=0)break;c.out.append(buffer,static_cast<std::size_t>(k));}
  ::close(fd[0]);int status=0;::waitpid(pid,&status,0);
  if(WIFEXITED(status))c.code=WEXITSTATUS(status);
  else if(WIFSIGNALED(status))c.signal=WTERMSIG(status);
  return c;
}
inline J child_json(const Child& c){
  J j={{"exit_code",c.code},{"signal",c.signal}};
  std::string text=c.out;
  const auto end=text.rfind("@@end-injection@@\n");
  if(end!=std::string::npos){j["injection_backtrace"]=text.substr(0,end);text=text.substr(end+18);}
  try{j["output"]=J::parse(text);}catch(...){j["output_text"]=text.substr(0,12000);}
  return j;
}
} // namespace stage1
