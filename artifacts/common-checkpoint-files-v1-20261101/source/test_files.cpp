#include "checkpoint_files.hpp"
#include <iostream>
namespace c=powerflow::common_checkpoint_files_v1;using J=c::Json;namespace fs=std::filesystem;
int interrupt_point=0;void inject(int point){if(point==interrupt_point)throw std::runtime_error("injected interrupted checkpoint");}
int main(int argc,char**argv){
 if(argc!=2){return 2;}
 const fs::path out=argv[1];fs::create_directories(out);J tests=J::array();unsigned failed=0;
 auto check=[&](const char*name,bool ok){tests.push_back({{"name",name},{"passed",ok}});failed+=!ok;std::cout<<(ok?"PASS ":"FAIL ")<<name<<'\n';};
 auto refuses=[&](const char*name,auto f){bool threw=false;try{f();}catch(const std::exception&){threw=true;}check(name,threw);};
 auto write=[](const fs::path&p,const std::string&s){std::ofstream(p,std::ios::binary)<<s;};
 auto seed=[&](const fs::path&p){fs::create_directories(p/"hall");c::write_text(p/"native.json","native checkpoint\n");c::write_text(p/"hall/manifest.json","nested owner manifest\n");};
 try{
  check("known SHA256",c::hash_text("abc")=="ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
  const auto root=out/"history";fs::create_directory(root);auto tmp=root/(c::step_name(1)+".tmp");seed(tmp);
  const std::string large(3u<<20,'x');c::write_text(tmp/"large.bin",large);
  check("streamed multi-block digest matches exact text digest",c::hash_file(tmp/"large.bin",large.size())==c::hash_text(large));
  refuses("stream file byte cap",[&]{c::hash_file(tmp/"large.bin",large.size()-1);});
  auto limits=c::Limits{};limits.files=2;refuses("file count cap",[&]{c::file_table(tmp,limits);});
  limits=c::Limits{};limits.total_bytes=large.size();refuses("total byte cap",[&]{c::file_table(tmp,limits);});
  const auto table=c::file_table(tmp);check("nested owner manifest is covered",table.contains("hall/manifest.json")&&table.size()==3);
  const J state={{"step",1},{"common_root",std::string(64,'a')}};const auto committed=c::commit(tmp,root/c::step_name(1),state);
  check("coherent manifest verifies exact state",c::verify(committed)==state);
  check("newest committed checkpoint",c::newest_verified(root)==committed);
  refuses("existing accepted checkpoint never overwritten",[&]{c::commit(root/"fresh.tmp",committed,state);});
  write(committed/"hall/manifest.json","different nested owner manifest\n");refuses("corrupted nested manifest refused",[&]{c::verify(committed);});
  write(committed/"hall/manifest.json","nested owner manifest\n");check("original files restore verifier parity",c::verify(committed)==state);
  write(committed/"extra","extra");refuses("added file refused",[&]{c::verify(committed);});fs::remove(committed/"extra");
  fs::rename(committed/"native.json",committed/"native.moved");refuses("missing or renamed file refused",[&]{c::verify(committed);});fs::rename(committed/"native.moved",committed/"native.json");
  fs::create_symlink(committed/"native.json",committed/"alias");refuses("symlink refused",[&]{c::verify(committed);});fs::remove(committed/"alias");
  auto manifest=c::read_text(committed/"manifest.json",8u<<20);write(committed/"manifest.json","{}");refuses("malformed manifest refused",[&]{c::verify(committed);});write(committed/"manifest.json",manifest);
  limits=c::Limits{};limits.manifest_bytes=1;refuses("manifest read cap",[&]{c::verify(committed,limits);});
  for(int point=1;point<=4;++point){const auto step=std::uint64_t(point+1);tmp=root/(c::step_name(step)+".tmp");seed(tmp);interrupt_point=point;
   refuses("injected commit interruption observed",[&]{c::commit(tmp,root/c::step_name(step),{{"step",step}},c::Limits{},inject);});
   check("interrupted checkpoint only visible after commit point",fs::exists(root/c::step_name(step))==(point==4));
   check("incomplete directories cannot supersede accepted history",c::newest_verified(root)==(point==4?root/c::step_name(step):committed));
  }
  const auto newest=root/c::step_name(5);write(newest/"native.json","corrupt newest");refuses("corrupt newest is refused without silently falling back",[&]{c::newest_verified(root);});
  refuses("step name overflow refused",[]{c::step_name(1000000000000ULL);});
  check("step names sort chronologically",c::step_name(9)<c::step_name(10));
  tmp=out/"sealed.tmp";seed(tmp);c::write_text(tmp/"manifest.json","already sealed");refuses("presealed temporary directory refused",[&]{c::commit(tmp,out/"sealed",{});});
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';check("unexpected exception",false);}
 const J result={{"passed",failed==0},{"failed",failed},{"cases",tests},{"scope","durable file-set component only; does not authorize restored owners or continuation"},{"hash_buffer_bytes",65536}};
 std::ofstream(out/"result.json")<<result.dump(2)<<'\n';return failed?1:0;
}
