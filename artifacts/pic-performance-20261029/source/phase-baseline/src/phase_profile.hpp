#pragma once
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

// Timing-only state. No field, particle, solver, or checkpoint data is stored here.
namespace pic_phase {
using Clock=std::chrono::steady_clock;
enum class Phase : std::size_t {
 DriverStep,DriverSelection,ThermalSource,Diagnostics,PeriodicOutput,
 Advance,DriftFirst,SolveFieldsFirst,ForceLoop,DriftSecond,SolveFieldsSecond,
 AdvanceEmpty,CudaFirst,CudaSecond,CudaOther,ResidualFirst,ResidualSecond,
 ResidualOther,FactorSetup,Count
};
constexpr std::size_t count=static_cast<std::size_t>(Phase::Count);
struct Definition {const char* name;const char* meaning;};
inline constexpr std::array<Definition,count> definitions{{
 {"driver_step","Complete successful while-condition and driver loop body; exclusive time is remaining driver glue."},
 {"driver_selection","Field history/boundary updates and initial timestep limits; nested diagnostics excluded from exclusive time."},
 {"thermal_source_reservoirs","Thermal interpolation, source timestep limits, reservoir accumulation, and injection."},
 {"diagnostics","Every diagnostics() call inside a step, including the loop condition, selection, guard, and progress/save calls."},
 {"periodic_output","Conditional frame/checkpoint/progress work inside a step; nested diagnostics excluded from exclusive time."},
 {"advance","advance() validation, transfer reset, and time update; nested drift/solve/force excluded from exclusive time."},
 {"drift_first","First half drift, including charge-conserving transfer bookkeeping and absorption."},
 {"solve_fields_first","First solve_fields(); exclusive time excludes nested CUDA solve, linear residual guard, and factor setup."},
 {"force_gather_boris","Whole particle force loop: self/external gathers, imposed fields, and Boris update."},
 {"drift_second","Second half drift, including charge-conserving transfer bookkeeping and absorption."},
 {"solve_fields_second","Second solve_fields(); exclusive time excludes nested CUDA solve, linear residual guard, and factor setup."},
 {"advance_empty","Whole advance_empty() call."},
 {"cuda_solve_first","CUDA RHS transfer, unchanged cuDSS SOLVE, and answer transfer inside first advance solve_fields."},
 {"cuda_solve_second","CUDA RHS transfer, unchanged cuDSS SOLVE, and answer transfer inside second advance solve_fields."},
 {"cuda_solve_other_step","CUDA transfer plus SOLVE occurring in a step outside the two advance solve_fields scopes."},
 {"cpu_residual_first","LinearSystem CPU residual vector/norm, solver bookkeeping, and unchanged guards in first solve_fields."},
 {"cpu_residual_second","LinearSystem CPU residual vector/norm, solver bookkeeping, and unchanged guards in second solve_fields."},
 {"cpu_residual_other_step","LinearSystem CPU residual/guards in a step outside the two advance solve_fields scopes."},
 {"cuda_factor_setup_in_step","Any cuDSS constructor/factorization that actually occurs during a driver step; startup constructors excluded."}
}};
struct Entry {double inclusive=0,exclusive=0;std::uint64_t calls=0;};
struct StepRecord {
 std::uint64_t step=0;
 double physical_start_s=0,dt_s=0,wall_s=0,top_level_child_s=0;
 std::array<Entry,count> phases{};
};
struct Scope;
inline thread_local StepRecord* current_step=nullptr;
inline thread_local Scope* current_scope=nullptr;
inline StepRecord working,terminal;
inline Clock::time_point step_start;
inline std::vector<StepRecord> steps;
inline std::array<Entry,count> totals{};
inline double total_wall_s=0;
inline bool terminal_present=false;

struct Scope {
 Phase phase;StepRecord* record=nullptr;Scope* parent=nullptr;
 Clock::time_point start;double child_s=0;bool running=false;
 explicit Scope(Phase value):phase(value),record(current_step) {
  if(record){parent=current_scope;current_scope=this;start=Clock::now();running=true;}
 }
 void stop() {
  if(!running)return;
  const double elapsed=std::chrono::duration<double>(Clock::now()-start).count();
  auto& entry=record->phases[static_cast<std::size_t>(phase)];
  entry.inclusive+=elapsed;entry.exclusive+=elapsed-child_s;++entry.calls;
  if(parent)parent->child_s+=elapsed;else record->top_level_child_s+=elapsed;
  current_scope=parent;running=false;
 }
 ~Scope(){stop();}
 Scope(const Scope&)=delete;Scope& operator=(const Scope&)=delete;
};
inline Phase solve_phase(Phase first,Phase second,Phase other) {
 for(auto* scope=current_scope;scope;scope=scope->parent) {
  if(scope->phase==Phase::SolveFieldsFirst)return first;
  if(scope->phase==Phase::SolveFieldsSecond)return second;
 }
 return other;
}
inline void begin_step(std::uint64_t index) {
 working=StepRecord{};working.step=index;
 current_step=&working;current_scope=nullptr;step_start=Clock::now();
}
inline void interval(double start,double dt) {
 working.physical_start_s=start;working.dt_s=dt;
}
inline void finish_record(StepRecord& record) {
 record.wall_s=std::chrono::duration<double>(Clock::now()-step_start).count();
 record.phases[static_cast<std::size_t>(Phase::DriverStep)]={record.wall_s,record.wall_s-record.top_level_child_s,1};
 current_step=nullptr;current_scope=nullptr;
}
inline void end_step() {
 finish_record(working);
 total_wall_s+=working.wall_s;
 for(std::size_t i=0;i<count;++i) {
  totals[i].inclusive+=working.phases[i].inclusive;
  totals[i].exclusive+=working.phases[i].exclusive;
  totals[i].calls+=working.phases[i].calls;
 }
 steps.push_back(working);
}
inline void end_terminal_condition() {
 finish_record(working);terminal=working;terminal_present=true;
}
inline void write_entries(std::ostream& out,const std::array<Entry,count>& entries) {
 out<<'{';
 for(std::size_t i=0;i<count;++i) {
  if(i)out<<',';
  out<<std::quoted(definitions[i].name)<<":{\"seconds_inclusive\":"<<entries[i].inclusive
     <<",\"seconds_exclusive\":"<<entries[i].exclusive<<",\"calls\":"<<entries[i].calls<<'}';
 }
 out<<'}';
}
inline double exclusive_sum(const std::array<Entry,count>& entries) {
 double result=0;for(const auto& entry:entries)result+=entry.exclusive;return result;
}
inline void write(const std::string& path) {
 std::ofstream out(path);out<<std::setprecision(17);
 out<<"{\"schema\":\"pic_additive_step_phases_v1\",\"clock\":\"steady_clock_wall_seconds\","
    <<"\"window\":\"Each successful driver while-condition plus its loop body. Startup/load-checkpoint/unit-bias solves and final output are excluded.\","
    <<"\"additive_rule\":\"Sum seconds_exclusive only. Inclusive times overlap through nesting and must not be added.\","
    <<"\"step_count\":"<<steps.size()<<",\"step_wall_seconds\":"<<total_wall_s
    <<",\"additive_error_seconds\":"<<exclusive_sum(totals)-total_wall_s<<",\"definitions\":{";
 for(std::size_t i=0;i<count;++i){if(i)out<<',';out<<std::quoted(definitions[i].name)<<':'<<std::quoted(definitions[i].meaning);}
 out<<"},\"totals\":";write_entries(out,totals);out<<",\"steps\":[";
 for(std::size_t i=0;i<steps.size();++i) {
  if(i)out<<',';const auto& step=steps[i];
  out<<"{\"step\":"<<step.step<<",\"physical_start_s\":"<<step.physical_start_s
     <<",\"dt_s\":"<<step.dt_s<<",\"wall_seconds\":"<<step.wall_s
     <<",\"additive_error_seconds\":"<<exclusive_sum(step.phases)-step.wall_s<<",\"phases\":";
  write_entries(out,step.phases);out<<'}';
 }
 out<<"],\"terminal_loop_condition\":";
 if(terminal_present){out<<"{\"wall_seconds\":"<<terminal.wall_s<<",\"phases\":";write_entries(out,terminal.phases);out<<'}';}
 else out<<"null";
 out<<"}\n";
}
}
