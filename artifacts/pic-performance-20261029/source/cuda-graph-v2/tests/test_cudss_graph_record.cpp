#include "cudss_graph_record.hpp"
#include <filesystem>

int main(int argc,char** argv) {
 try {
  if(argc!=3)throw std::runtime_error("Expected scenario and absolute output path");
  const std::string scenario=argv[1],path=argv[2];
  ::setenv("RAIKOU_PIC_CUDSS_GRAPH_RECORD",path.c_str(),1);
  ::setenv("RAIKOU_PIC_CUDSS_GRAPH",scenario=="off"?"off":"required",1);
  auto& run=pic_graph::RunRecord::instance();
  if(scenario=="invalid_path") {
   ::setenv("RAIKOU_PIC_CUDSS_GRAPH_RECORD","relative-record-must-not-exist.json",1);
   bool rejected=false;try{run.begin(scenario);}catch(const std::exception&){rejected=true;}
   run.fail(std::size_t(-1),"Rejected relative path");
   if(!rejected || std::filesystem::exists("relative-record-must-not-exist.json"))
    throw std::runtime_error("Relative record path was not rejected safely");
   return 0;
  }
  {std::ofstream old(path);old<<"{\"run_status\":\"completed\",\"graph_status\":\"qualified\",\"run_id\":\"old-run\"}\n";}
  run.begin(scenario);
  std::filesystem::copy_file(path,path+".initial",std::filesystem::copy_options::overwrite_existing);
  auto a=run.add(130321,scenario=="off"?"off":scenario=="fallback"?"auto":"required");
  if(scenario=="off") {
   run.update(a,"off","Graph disabled",0,0,0,false);run.direct_solve(a);run.close(a);run.complete();
  } else if(scenario=="unclosed") {
   bool rejected=false;try{run.complete();}catch(const std::exception&){rejected=true;}
   if(!rejected)throw std::runtime_error("Unclosed solver accepted");run.close(a);
  } else {
   auto b=run.add(130321,"required");
   run.update(a,scenario=="fallback"?"fallback":"qualified","first solver",40,37,2,scenario!="fallback");
   run.direct_solve(a);
   if(scenario!="fallback")run.graph_replay(a);
   run.update(b,"qualified","second solver",40,37,2,true);run.graph_replay(b);
   if(scenario=="late_failure") {
    run.fail(a,"Late launch failure: \"quoted\"\nline");
    run.update(a,"qualified","must not erase failure",40,37,2,true);
   }
   run.close(a);run.close(b);
   if(scenario=="outer_failure")run.fail(std::size_t(-1),"Outer scientific guard failed");
   bool rejected=false;try{run.complete();}catch(const std::exception&){rejected=true;}
   if(rejected!=(scenario=="late_failure" || scenario=="outer_failure"))
    throw std::runtime_error("Unexpected completion outcome");
  }
  std::cout<<scenario<<": passed\n";return 0;
 } catch(const std::exception& error) {std::cerr<<error.what()<<'\n';return 1;}
}
