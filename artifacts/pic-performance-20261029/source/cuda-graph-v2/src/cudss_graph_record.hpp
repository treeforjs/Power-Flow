#pragma once
// Process-local diagnostics only. No CUDA calls and no per-solve file I/O.
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <unistd.h>

namespace pic_graph {
class RunRecord {
 struct Solver {
  int rows;
  std::string requested_mode,status="initializing",reason;
  bool closed=false,first_rhs_bitwise_qualified=false;
  std::size_t nodes=0,kernels=0,copies=0,direct_solves=0,graph_replays=0;
 };
 bool started=false;
 std::string path,run_id,context,requested_mode,run_status="running",reason;
 std::vector<Solver> solvers;
 static std::string quote(const std::string& value) {
  const char* hex="0123456789abcdef";std::string result="\"";
  for(unsigned char c:value) {
   if(c=='\"' || c=='\\'){result+='\\';result+=char(c);}
   else if(c<32){result+="\\u00";result+=hex[c>>4];result+=hex[c&15];}
   else result+=char(c);
  }
  return result+'\"';
 }
 bool all_qualified() const {
  if(run_status=="failed" || solvers.empty())return false;
  for(const auto& s:solvers)if(s.status!="qualified")return false;
  return true;
 }
 std::string graph_status() const {
  if(run_status=="failed")return "failed";
  bool fallback=false,pending=solvers.empty(),off=false,qualified=false;
  for(const auto& s:solvers) {
   if(s.status=="failed" || s.status=="unsupported")return "failed";
   fallback|=s.status=="fallback";
   off|=s.status=="off";qualified|=s.status=="qualified";
   pending|=s.status=="initializing" || s.status=="pending";
  }
  if(fallback)return "fallback";
  if(pending)return "pending";
  if(off && qualified)return "mixed";
  return off?"off":"qualified";
 }
 void write() const {
  if(path.empty())return;
  const std::string temporary=path+".tmp."+std::to_string(::getpid());
  std::ofstream output(temporary);
  if(!output)throw std::runtime_error("Cannot open graph run record temporary file");
  bool required=!solvers.empty();
  for(const auto& s:solvers)required&=s.requested_mode=="required";
  output<<"{\"schema_version\":2,\"run_id\":"<<quote(run_id)
        <<",\"pid\":"<<::getpid()<<",\"context\":"<<quote(context)
        <<",\"requested_mode\":"<<quote(requested_mode)
        <<",\"run_status\":"<<quote(run_status)<<",\"reason\":"<<quote(reason)
        <<",\"graph_status\":"<<quote(graph_status())
        <<",\"all_solvers_graph_qualified\":"<<(all_qualified()?"true":"false")
        <<",\"all_solvers_requested_required\":"<<(required?"true":"false")
        <<",\"outer_process_receipt_required\":true,\"solvers\":[";
  for(std::size_t i=0;i<solvers.size();++i) {
   const auto& s=solvers[i];if(i)output<<',';
   output<<"{\"solver_id\":"<<i<<",\"matrix_rows\":"<<s.rows
         <<",\"requested_mode\":"<<quote(s.requested_mode)<<",\"status\":"<<quote(s.status)
         <<",\"reason\":"<<quote(s.reason)<<",\"closed\":"<<(s.closed?"true":"false")
         <<",\"ir_steps\":2,\"graph_nodes\":"<<s.nodes<<",\"kernel_nodes\":"<<s.kernels
         <<",\"memcpy_nodes\":"<<s.copies<<",\"direct_solve_calls\":"<<s.direct_solves
         <<",\"graph_replays\":"<<s.graph_replays
         <<",\"first_rhs_bitwise_qualified\":"<<(s.first_rhs_bitwise_qualified?"true":"false")<<'}';
  }
  output<<"]}\n";output.close();
  if(!output){std::remove(temporary.c_str());throw std::runtime_error("Cannot write graph run record");}
  if(std::rename(temporary.c_str(),path.c_str())!=0) {
   std::remove(temporary.c_str());throw std::runtime_error("Cannot replace graph run record");
  }
 }
public:
 static RunRecord& instance(){static RunRecord record;return record;}
 void begin(const std::string& run_context) {
  if(started)throw std::runtime_error("Graph run record already started");
  started=true;context=run_context;
  const char* record_path=std::getenv("RAIKOU_PIC_CUDSS_GRAPH_RECORD");
  const std::string requested_path=record_path?record_path:"";
  if(!requested_path.empty() && requested_path.front()!='/')throw std::runtime_error("Graph run record path must be absolute");
  path=requested_path;
  const char* mode=std::getenv("RAIKOU_PIC_CUDSS_GRAPH");requested_mode=mode?mode:"off";
  const auto now=std::chrono::system_clock::now().time_since_epoch();
  run_id=std::to_string(::getpid())+":"+std::to_string(std::chrono::duration_cast<std::chrono::nanoseconds>(now).count());
  write(); // Replace a previous run immediately, including off and zero-solver runs.
 }
 std::size_t add(int rows,const std::string& mode) {
  if(!started)begin("solver_only");
  if(run_status!="running")throw std::runtime_error("Cannot add solver after graph run ended");
  solvers.push_back(Solver{rows,mode});write();return solvers.size()-1;
 }
 void update(std::size_t id,const char* status,const char* detail,
             std::size_t nodes,std::size_t kernels,std::size_t copies,bool first_rhs) {
  auto& s=solvers.at(id);
  if(s.status=="failed" || s.status=="unsupported")return;
  s.status=status;s.reason=detail;s.nodes=nodes;s.kernels=kernels;s.copies=copies;
  s.first_rhs_bitwise_qualified=first_rhs;write();
 }
 void direct_solve(std::size_t id){++solvers.at(id).direct_solves;}
 void graph_replay(std::size_t id){++solvers.at(id).graph_replays;}
 void fail(std::size_t id,const char* detail) noexcept {
  try {
   if(run_status!="failed")reason=detail;
   run_status="failed";
   if(id<solvers.size()) {
    auto& s=solvers[id];if(s.status!="failed"){s.status="failed";s.reason=detail;}
   }
   write();
  } catch(...) {std::cerr<<"Failed to persist graph run failure record\n";}
 }
 void close(std::size_t id) noexcept {
  try {solvers.at(id).closed=true;write();}
  catch(...) {fail(id,"Cannot persist solver close record");}
 }
 void complete() {
  for(const auto& s:solvers)if(!s.closed) {
   fail(std::size_t(-1),"Graph run completed with an unclosed solver");
   throw std::runtime_error("Graph run completed with an unclosed solver");
  }
  if(run_status=="failed") {write();throw std::runtime_error("Graph run has a recorded failure");}
  run_status="completed";
  try {write();}catch(...) {fail(std::size_t(-1),"Cannot persist graph run completion");throw;}
 }
};
}
