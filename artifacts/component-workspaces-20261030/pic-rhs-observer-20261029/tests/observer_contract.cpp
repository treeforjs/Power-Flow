#include "../source/src/rhs_observer.hpp"
#include <cassert>
#include <cmath>
#include <iostream>

namespace {
constexpr std::array<int,5> rows{0,2,5,7,8};
constexpr std::array<int,8> cols{0,1,0,1,2,1,2,3};
constexpr std::array<double,8> matrix{2,0,0,4,-1,-1,4,3};
constexpr std::array<int,2> essential{0,3};
void emit(rhs_observer::Recorder& observer,std::uint64_t step,bool populated=true) {
 const double time=step*1e-14;
 observer.driver_step(step,time,1e-14);
 observer.begin_advance(populated,populated?4:0,time-1e-10,1e-14,0,10);
 if(populated)for(unsigned half=1;half<=2;++half) {
  observer.set_half(half);
  std::array<double,4> trial{10,1000.+step+half,2000.+step+half,0};
  if(half==2)trial={10,1,2,0};
  std::array<double,4> charge{0,2,7,0},rhs{20,2,7,0},answer{10,1,2,0},residual{0,0,0,0};
  const auto token=observer.before_form(trial.data(),charge.data(),4,essential.data(),2,0x111);
  // Synthetic contract stand-in for the native FormLinearSystem overwrite.
  trial[1]=trial[2]=0;
  observer.after_form(token,rows.data(),cols.data(),matrix.data(),4,8,trial.data(),rhs.data(),0x222);
  observer.accepted(rhs.data(),answer.data(),residual.data(),4,0,std::sqrt(453.),77,0x333);
  trial.fill(-999);charge.fill(-999);rhs.fill(-999);answer.fill(-999);residual.fill(-999);
 }
 observer.end_advance();
}
template<class T>std::vector<T> read(const std::filesystem::path& path) {
 std::ifstream file(path,std::ios::binary);
 assert(file);const auto bytes=std::filesystem::file_size(path);assert(bytes%sizeof(T)==0);
 std::vector<T> data(bytes/sizeof(T));file.read(reinterpret_cast<char*>(data.data()),bytes);assert(file);
 return data;
}
std::uint64_t disk_bytes(const std::filesystem::path& dir) {
 std::uint64_t bytes=0;
 for(const auto& file:std::filesystem::recursive_directory_iterator(dir)) {
  if(file.is_regular_file())bytes+=file.file_size();
 }
 return bytes;
}
}
int main(int argc,char** argv) {
 assert(argc==2);const std::filesystem::path root(argv[1]);
 std::filesystem::create_directories(root);
 {
  rhs_observer::Recorder observer(root/"normal",8,256*1024,192*1024,true);
  // Startup and unit-bias-like calls without an advance are excluded.
  assert(!observer.before_form(nullptr,nullptr,0,nullptr,0,0).selected);
  observer.accepted(nullptr,nullptr,nullptr,0,0,0,1,1);
  for(unsigned i=0;i<12;++i)emit(observer,100+i,i!=2 && i!=7);
  observer.finish();observer.finish();
  auto projected=read<double>(root/"normal/call_0001_projected_previous_phi.f64");
  auto incoming=read<double>(root/"normal/call_0001_incoming_zeroed_X.f64");
  auto accepted=read<double>(root/"normal/call_0001_accepted_solution.f64");
  assert((projected==std::vector<double>{10,1101,2101,0}));
  assert((incoming==std::vector<double>{10,0,0,0}));
  assert((accepted==std::vector<double>{10,1,2,0}));
  assert((read<int>(root/"normal/matrix.rows.i32")==std::vector<int>(rows.begin(),rows.end())));
  assert((read<double>(root/"normal/matrix.values.f64")==std::vector<double>(matrix.begin(),matrix.end())));
  assert((read<std::uint8_t>(root/"normal/essential_mask.u8")==std::vector<std::uint8_t>{1,0,0,1}));
  assert(std::filesystem::exists(root/"normal/call_0016_guard_residual.f64"));
  assert(!std::filesystem::exists(root/"normal/call_0017_guard_residual.f64"));
  assert(disk_bytes(root/"normal")<=256*1024);
 }
 {
  rhs_observer::Recorder observer(root/"capacity",8,32768,32000,true);
  for(unsigned i=0;i<10;++i)emit(observer,200+i);
  observer.finish();assert(disk_bytes(root/"capacity")<=32768);
  assert(std::filesystem::exists(root/"capacity/call_0002_guard_residual.f64"));
  assert(!std::filesystem::exists(root/"capacity/call_0003_guard_residual.f64"));
 }
 {
  rhs_observer::Recorder observer(root/"matrix_capacity",8,32768,32368,true);
  emit(observer,300);observer.finish();
  assert(disk_bytes(root/"matrix_capacity")<=32768);
  assert(!std::filesystem::exists(root/"matrix_capacity/matrix.values.f64"));
  assert(!std::filesystem::exists(root/"matrix_capacity/call_0001_projected_previous_phi.f64"));
 }
 {
  rhs_observer::Recorder observer(root/"huge_dimensions",8,32768,32000,true);
  observer.driver_step(400,0,1e-14);observer.begin_advance(true,1,0,1e-14,0,10);
  observer.set_half(1);
  assert(!observer.before_form(nullptr,nullptr,1ULL<<40,nullptr,0,0).selected);
  observer.end_advance();observer.finish();assert(disk_bytes(root/"huge_dimensions")<=32768);
 }
 {
  std::filesystem::create_directories(root/"nonempty");
  {std::ofstream sentinel(root/"nonempty/sentinel");sentinel<<"untouched";}
  rhs_observer::Recorder observer(root/"nonempty",8,32768,32000,true);
  bool rejected=false;try{observer.finish();}catch(const std::runtime_error&){rejected=true;}
  assert(rejected && std::filesystem::file_size(root/"nonempty/sentinel")==9);
 }
 {
  rhs_observer::Recorder observer(root/"mutated_matrix",8,256*1024,192*1024,true);
  emit(observer,500);
  observer.driver_step(501,0,1e-14);observer.begin_advance(true,1,0,1e-14,0,10);
  observer.set_half(1);
  std::array<double,4> values{10,20,30,0},rhs{20,2,7,0};
  auto changed=matrix;changed[3]+=1;
  auto token=observer.before_form(values.data(),rhs.data(),4,essential.data(),2,0x111);
  bool rejected=false;
  try{observer.after_form(token,rows.data(),cols.data(),changed.data(),4,8,values.data(),rhs.data(),0x222);}
  catch(const std::runtime_error&){rejected=true;}
  assert(rejected);observer.end_advance();observer.finish();
  assert(!std::filesystem::exists(root/"mutated_matrix/call_0003_projected_previous_phi.f64"));
 }
 std::cout<<"CPU-only observer contracts passed: first8/16, startup/empty exclusion, true pre-form copy, deep copies, exact CSR/zeros/mask, byte caps, huge-dimension rejection, no overwrite, matrix mutation detection.\n";
}
