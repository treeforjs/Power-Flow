#pragma once
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace rhs_observer {
static_assert(sizeof(double)==8 && std::numeric_limits<double>::is_iec559);
static_assert(sizeof(int)==4 && std::endian::native==std::endian::little);
inline constexpr std::uint64_t output_limit=128ULL*1024*1024;
inline constexpr std::uint64_t metadata_reserve=2ULL*1024*1024;
inline constexpr unsigned requested_steps=8;
struct Driver {
 bool valid=false;std::uint64_t step=0;double time_s=0,dt_s=0;
};
struct Token {bool selected=false;std::size_t index=0;};
struct Call {
 std::uint64_t sequence=0,step_ordinal=0,driver_step=0,particles=0;
 int half=0,linear_system_id=-1;
 std::uintptr_t field_state=0,matrix_address=0,linear_system_address=0;
 double driver_time_s=0,dt_s=0,ledger_time_s=0,lower_V=0,upper_V=0;
 double residual_norm=0,rhs_norm=0;
 bool formed=false,accepted=false;
 // Actual pre-form trial, original full charge RHS, actual constrained RHS,
 // actual post-form X, checked solution, and the existing guard residual.
 std::array<std::vector<double>,6> vectors;
};
class Recorder {
 std::filesystem::path directory;
 unsigned limit_steps;
 std::uint64_t limit_bytes,reserve_bytes,reserved_payload=0;
 Driver driver;
 bool active=false,selected=false,limited=false,incomplete=false,finished=false;
 std::uint64_t seen_steps=0,complete_steps=0,particles=0;
 double ledger_time=0,lower_V=0,upper_V=0;
 unsigned half=0,started=0,accepted_calls=0;
 std::size_t step_first_call=0;
 std::uint64_t step_vector_reservation=0;
 std::vector<Call> calls;
 std::vector<std::int32_t> rows,columns,essential;
 std::vector<double> coefficients;
 std::vector<std::uint8_t> mask;
 std::uintptr_t captured_matrix=0;
 bool matrix_present=false;
 bool synthetic_fixture=false;
 static constexpr std::array<const char*,6> vector_names{{
  "projected_previous_phi","full_charge_rhs","constrained_rhs",
  "incoming_zeroed_X","accepted_solution","guard_residual"}};
 static void require(bool condition,const char* message) {
  if(!condition)throw std::runtime_error(std::string("RHS observer: ")+message);
 }
 template<class T> static std::vector<T> copy(const T* data,std::size_t size) {
  if(size==0)return {};
  require(data!=nullptr,"null input buffer");
  return std::vector<T>(data,data+size);
 }
 template<class T> static bool same(const std::vector<T>& saved,const T* data,std::size_t size) {
  return saved.size()==size && (size==0 || std::memcmp(saved.data(),data,size*sizeof(T))==0);
 }
 bool reserve(std::uint64_t bytes) {
  if(limit_bytes<reserve_bytes || bytes>limit_bytes-reserve_bytes-reserved_payload) {
   limited=true;return false;
  }
  reserved_payload+=bytes;return true;
 }
 void abandon_step() noexcept {
  calls.resize(step_first_call);
  reserved_payload-=step_vector_reservation;
  step_vector_reservation=0;selected=false;
 }
 static std::string call_stem(std::uint64_t index) {
  std::ostringstream name;name<<"call_"<<std::setfill('0')<<std::setw(4)<<index;return name.str();
 }
 static std::string address(std::uintptr_t value) {
  std::ostringstream text;text<<"0x"<<std::hex<<value;return text.str();
 }
 std::uint64_t binary_bytes() const {
  std::uint64_t size=matrix_present?(rows.size()+columns.size()+essential.size())*4ULL+
                                    coefficients.size()*8ULL+mask.size():0;
  for(const auto& call:calls)for(const auto& vector:call.vectors)size+=vector.size()*8ULL;
  return size;
 }
 std::string manifest() const {
  std::ostringstream out;out<<std::setprecision(17);
  const bool complete=complete_steps==limit_steps && calls.size()==2ULL*limit_steps && !limited && !incomplete;
  out<<"{\"schema\":\"pic_rhs_observation_v1\",\"status\":"
     <<std::quoted(complete?"complete":limited?"capacity_limited":"incomplete")
     <<",\"requested_populated_steps\":"<<limit_steps<<",\"captured_populated_steps\":"<<complete_steps
     <<",\"populated_advances_seen\":"<<seen_steps<<",\"captured_solves\":"<<calls.size()
     <<",\"output_limit_bytes\":"<<limit_bytes<<",\"metadata_reserved_bytes\":"<<reserve_bytes
     <<",\"payload_bytes\":"<<binary_bytes()
     <<",\"synthetic_contract_fixture\":"<<(synthetic_fixture?"true":"false")
     <<",\"startup_and_unit_bias_excluded\":true,\"selection\":\"First populated main-driver advance entries; both original half-step solves retained\","
     <<"\"producer_header\":"<<std::quoted(std::string(__FILE__))
     <<",\"original_source\":\"/home/treefor/projects/runtimes/pic-cpu-candidate-v2-20261029/source\","
     <<"\"solver_contract\":{\"backend\":\"cudss\",\"matrix_and_vector_precision\":\"FP64\","
     <<"\"ir_n_steps\":2,\"settings_source\":\"Unchanged source configuration; no extra solver getter or solve was called\","
     <<"\"guard\":\"norm(Ax-b)<=max(1e-30,2e-13*norm(b)); unchanged CPU guard\"},"
     <<"\"matrix\":{\"present\":"<<(matrix_present?"true":"false")<<",\"rows\":"<<mask.size()
     <<",\"stored_entries_including_zeros\":"<<coefficients.size()
     <<",\"storage\":\"Exact native constrained MFEM CSR, including original row/column order and stored zeros\","
     <<"\"row_file\":\"matrix.rows.i32\",\"column_file\":\"matrix.cols.i32\",\"value_file\":\"matrix.values.f64\","
     <<"\"essential_file\":\"essential.i32\",\"essential_mask_file\":\"essential_mask.u8\","
     <<"\"essential_list_entries\":"<<essential.size()<<",\"matrix_address\":"<<std::quoted(address(captured_matrix))
     <<",\"matrix_and_mask_bitwise_checked_each_call\":true},"
     <<"\"projected_prior_is_actual_full_trial\":true,\"raw_unprojected_potential_exported\":false,"
     <<"\"vector_definitions\":{\"projected_previous_phi\":\"Actual trial after original ProjectBdrCoefficient and finite-foil assignments, copied immediately before FormLinearSystem. Original projection may also touch DOFs outside the final essential mask. Never reconstructed from zeroed X.\","
     <<"\"full_charge_rhs\":\"Existing full_rhs copy of the assembled particle source before boundary elimination.\","
     <<"\"constrained_rhs\":\"Actual b after native FormLinearSystem and before the backend call.\","
     <<"\"incoming_zeroed_X\":\"Actual X after default copy_interior=0, before the backend modifies it.\","
     <<"\"accepted_solution\":\"Actual answer after existing CPU residual/finite guards passed.\","
     <<"\"guard_residual\":\"Existing CPU guard error vector A*answer-rhs; no extra matrix-vector evaluation.\"},"
     <<"\"binary_format\":{\"floating\":\"little-endian IEEE754 binary64\",\"indices\":\"little-endian int32\",\"mask\":\"uint8,1=essential\"},\"calls\":[";
  for(std::size_t i=0;i<calls.size();++i) {
   if(i)out<<',';
   const auto& call=calls[i];
   out<<"{\"sequence\":"<<call.sequence<<",\"captured_step_ordinal\":"<<call.step_ordinal
      <<",\"driver_step_index\":"<<call.driver_step<<",\"half\":"<<call.half
      <<",\"phase\":"<<std::quoted(call.half==1?"after_first_half_drift":"after_second_half_drift")
      <<",\"system_role\":\"main_system_within_populated_advance\",\"linear_system_id\":"<<call.linear_system_id
      <<",\"linear_system_address\":"<<std::quoted(address(call.linear_system_address))
      <<",\"field_state_address\":"<<std::quoted(address(call.field_state))
      <<",\"matrix_address\":"<<std::quoted(address(call.matrix_address))
      <<",\"particles_at_advance_entry\":"<<call.particles
      <<",\"driver_start_time_s\":"<<call.driver_time_s<<",\"dt_s\":"<<call.dt_s
      <<",\"ledger_time_s_at_advance_entry\":"<<call.ledger_time_s
      <<",\"lower_potential_V\":"<<call.lower_V<<",\"upper_potential_V\":"<<call.upper_V
      <<",\"residual_norm\":"<<call.residual_norm<<",\"rhs_norm\":"<<call.rhs_norm
      <<",\"unchanged_cpu_guard_passed\":"<<(call.accepted?"true":"false")<<",\"vectors\":{";
   for(std::size_t v=0;v<vector_names.size();++v) {
    if(v)out<<',';
    out<<std::quoted(vector_names[v])<<":{\"file\":"
       <<std::quoted(call_stem(call.sequence)+"_"+vector_names[v]+".f64")
       <<",\"entries\":"<<call.vectors[v].size()<<",\"bytes\":"<<call.vectors[v].size()*8ULL<<'}';
   }
   out<<"}}";
  }
  out<<"],\"no_extra_solver_or_residual_evaluation\":true,\"no_physics_or_checkpoint_data_changed\":true}\n";
  return out.str();
 }
 void write_file(const std::string& name,const void* data,std::size_t bytes) const {
  std::ofstream output(directory/name,std::ios::binary);
  require(bool(output),"cannot create capture file");
  if(bytes)output.write(static_cast<const char*>(data),static_cast<std::streamsize>(bytes));
  output.close();require(bool(output),"capture file write failed");
 }
public:
 explicit Recorder(std::filesystem::path path={},unsigned steps=requested_steps,
                   std::uint64_t bytes=output_limit,std::uint64_t metadata=metadata_reserve,
                   bool synthetic=false)
  :directory(std::move(path)),limit_steps(steps),limit_bytes(bytes),reserve_bytes(metadata),
   synthetic_fixture(synthetic) {
  require(directory.empty() || directory.is_absolute(),"capture directory must be absolute");
  require(limit_bytes>=reserve_bytes,"metadata reservation exceeds output cap");
 }
 bool enabled() const {return !directory.empty();}
 void driver_step(std::uint64_t index,double time,double dt) {
  if(enabled())driver={true,index,time,dt};
 }
 void begin_advance(bool populated,std::uint64_t particle_count,double ledger,double dt,
                    double lower,double upper) {
  if(!enabled() || !populated)return;
  require(!active,"nested populated advance");
  require(driver.valid,"missing main-driver step metadata");
  require(driver.dt_s==dt,"driver/advance timestep mismatch");
  ++seen_steps;active=true;selected=seen_steps<=limit_steps && !limited;
  particles=particle_count;ledger_time=ledger;lower_V=lower;upper_V=upper;
  half=started=accepted_calls=0;step_first_call=calls.size();step_vector_reservation=0;
 }
 void set_half(unsigned value) {
  if(!active)return;
  require(value==half+1 && value<=2,"invalid advance half ordering");half=value;
 }
 Token before_form(const double* projected,const double* charge,std::size_t n,
                   const int* essential_dofs,std::size_t essential_count,std::uintptr_t field_state) {
  if(!active || !selected)return {};
  require(half==started+1 && started<2,"unexpected or repeated solve in advance half");
  if(half==1) {
   require(n<=std::numeric_limits<std::uint64_t>::max()/96ULL,"vector size overflow");
   const std::uint64_t reservation=n*96ULL; // two calls, six FP64 vectors each
   if(!reserve(reservation)){selected=false;return {};}
   step_vector_reservation=reservation;
   if(essential_count>(limit_bytes-reserve_bytes-reserved_payload)/4ULL) {
    limited=true;abandon_step();return {};
   }
  }
  Call call;
  call.sequence=calls.size()+1;call.step_ordinal=seen_steps;call.driver_step=driver.step;
  call.particles=particles;call.half=half;call.field_state=field_state;
  call.driver_time_s=driver.time_s;call.dt_s=driver.dt_s;call.ledger_time_s=ledger_time;
  call.lower_V=lower_V;call.upper_V=upper_V;
  call.vectors[0]=copy(projected,n);call.vectors[1]=copy(charge,n);
  if(matrix_present)require(same(essential,essential_dofs,essential_count),"essential list changed");
  else essential=copy(essential_dofs,essential_count);
  calls.push_back(std::move(call));++started;
  return {true,calls.size()-1};
 }
 void after_form(Token token,const int* matrix_rows,const int* matrix_columns,
                 const double* values,std::size_t n,std::size_t nnz,
                 const double* guess,const double* rhs,std::uintptr_t matrix_address) {
  if(!token.selected || !selected)return;
  require(token.index<calls.size() && token.index+1==calls.size(),"invalid capture token");
  auto& call=calls[token.index];
  require(call.vectors[0].size()==n,"pre/post-form vector size changed");
  if(!matrix_present) {
   if(nnz>(limit_bytes-reserve_bytes)/12ULL || n>(limit_bytes-reserve_bytes)/4ULL) {
    limited=true;abandon_step();return;
   }
   const std::uint64_t matrix_bytes=(n+1)*4ULL+nnz*12ULL+essential.size()*4ULL+n;
   if(!reserve(matrix_bytes)){abandon_step();return;}
   rows=copy(matrix_rows,n+1);columns=copy(matrix_columns,nnz);coefficients=copy(values,nnz);
   mask.assign(n,0);
   for(auto dof:essential){require(dof>=0 && std::size_t(dof)<n,"essential DOF out of range");mask[dof]=1;}
   captured_matrix=matrix_address;matrix_present=true;
  } else {
   require(same(rows,matrix_rows,n+1) && same(columns,matrix_columns,nnz) &&
           same(coefficients,values,nnz),"native matrix changed during capture");
   require(captured_matrix==matrix_address,"native matrix identity changed");
  }
  call.matrix_address=matrix_address;
  call.vectors[2]=copy(rhs,n);call.vectors[3]=copy(guess,n);call.formed=true;
 }
 void accepted(const double* rhs,const double* answer,const double* residual,std::size_t n,
               double residual_norm,double rhs_norm,int system_id,std::uintptr_t system_address) {
  if(!active || !selected)return;
  require(!calls.empty(),"accepted solve has no pre-form observation");
  auto& call=calls.back();
  require(call.formed && !call.accepted && call.half==int(half),"accepted solve ordering mismatch");
  require(same(call.vectors[2],rhs,n),"backend RHS differs from observed constrained RHS");
  call.vectors[4]=copy(answer,n);call.vectors[5]=copy(residual,n);
  call.residual_norm=residual_norm;call.rhs_norm=rhs_norm;call.linear_system_id=system_id;
  call.linear_system_address=system_address;call.accepted=true;++accepted_calls;
 }
 void end_advance() noexcept {
  if(!active)return;
  if(selected) {
   if(started==2 && accepted_calls==2 && half==2)++complete_steps;
   else {incomplete=true;abandon_step();}
  }
  active=false;selected=false;
 }
 void finish() {
  if(!enabled() || finished)return;
  require(!active,"finalized during an advance");
  for(const auto& call:calls)require(call.accepted,"incomplete call in output");
  const std::string text=manifest();
  require(text.size()<=reserve_bytes,"manifest exceeds reserved metadata cap");
  require(binary_bytes()+text.size()<=limit_bytes,"total output exceeds byte cap");
  if(std::filesystem::exists(directory))
   require(std::filesystem::is_directory(directory) && std::filesystem::is_empty(directory),
           "capture output directory is not empty");
  else std::filesystem::create_directories(directory);
  if(matrix_present) {
   write_file("matrix.rows.i32",rows.data(),rows.size()*4);
   write_file("matrix.cols.i32",columns.data(),columns.size()*4);
   write_file("matrix.values.f64",coefficients.data(),coefficients.size()*8);
   write_file("essential.i32",essential.data(),essential.size()*4);
   write_file("essential_mask.u8",mask.data(),mask.size());
  }
  for(const auto& call:calls)for(std::size_t v=0;v<vector_names.size();++v)
   write_file(call_stem(call.sequence)+"_"+vector_names[v]+".f64",
              call.vectors[v].data(),call.vectors[v].size()*8);
  write_file("capture.json",text.data(),text.size());finished=true;
 }
};
inline Recorder& recorder() {
 static Recorder instance([]{
  const char* path=std::getenv("RAIKOU_PIC_RHS_OBSERVE_DIR");
  return path?std::filesystem::path(path):std::filesystem::path{};
 }());
 return instance;
}
inline void driver_step(std::uint64_t step,double time,double dt) {recorder().driver_step(step,time,dt);}
class AdvanceScope {
public:
 AdvanceScope(bool populated,std::uint64_t particles,double ledger,double dt,double lower,double upper) {
  recorder().begin_advance(populated,particles,ledger,dt,lower,upper);
 }
 ~AdvanceScope(){recorder().end_advance();}
 AdvanceScope(const AdvanceScope&)=delete;AdvanceScope& operator=(const AdvanceScope&)=delete;
};
inline void set_half(unsigned half){recorder().set_half(half);}
inline Token before_form(const double* projected,const double* charge,std::size_t n,
                         const int* essential,std::size_t count,std::uintptr_t field) {
 return recorder().before_form(projected,charge,n,essential,count,field);
}
inline void after_form(Token token,const int* rows,const int* columns,const double* values,
                       std::size_t n,std::size_t nnz,const double* guess,const double* rhs,
                       std::uintptr_t matrix) {
 recorder().after_form(token,rows,columns,values,n,nnz,guess,rhs,matrix);
}
inline void accepted(const double* rhs,const double* answer,const double* residual,std::size_t n,
                     double residual_norm,double rhs_norm,int system,std::uintptr_t address) {
 recorder().accepted(rhs,answer,residual,n,residual_norm,rhs_norm,system,address);
}
inline void finish(){recorder().finish();}
}
