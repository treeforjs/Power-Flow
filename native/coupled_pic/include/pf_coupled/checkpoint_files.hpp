#pragma once
#include "raikou/json.h"
#include <openssl/evp.h>
#include <array>
#include <algorithm>
#include <cerrno>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <fcntl.h>
#include <memory>
#include <set>
#include <unistd.h>
namespace powerflow::common_checkpoint_files_v1 {
namespace fs=std::filesystem;using Json=raikou::Json;
inline constexpr char schema[]="common-checkpoint-files-v1";
struct Limits {std::uintmax_t file_bytes=std::uintmax_t(8)<<30,total_bytes=std::uintmax_t(32)<<30;std::size_t files=4096,manifest_bytes=8u<<20;};
inline void need(bool yes,const std::string&why){if(!yes)throw std::runtime_error("common checkpoint: "+why);}
inline std::string hex(const unsigned char*p,std::size_t size){static constexpr char h[]="0123456789abcdef";std::string s(2*size,'0');for(std::size_t k=0;k<size;++k){s[2*k]=h[p[k]>>4];s[2*k+1]=h[p[k]&15];}return s;}
class Hash {
 std::unique_ptr<EVP_MD_CTX,decltype(&EVP_MD_CTX_free)> context_{EVP_MD_CTX_new(),EVP_MD_CTX_free};
public:
 Hash(){need(bool(context_)&&EVP_DigestInit_ex(context_.get(),EVP_sha256(),nullptr)==1,"SHA256 initialization failed");}
 void add(const void*p,std::size_t n){need(EVP_DigestUpdate(context_.get(),p,n)==1,"SHA256 update failed");}
 std::string finish(){std::array<unsigned char,EVP_MAX_MD_SIZE>bytes{};unsigned n=0;need(EVP_DigestFinal_ex(context_.get(),bytes.data(),&n)==1&&n==32,"SHA256 final failed");return hex(bytes.data(),n);}
};
inline std::string hash_text(const std::string&s){Hash h;h.add(s.data(),s.size());return h.finish();}
inline void ordinary_file(const fs::path&p){need(fs::is_regular_file(fs::symlink_status(p)),"missing/nonregular file "+p.string());}
inline std::string hash_file(const fs::path&p,std::uintmax_t cap){
 ordinary_file(p);const auto size=fs::file_size(p);need(size<=cap,"file cap "+p.string());
 std::ifstream in(p,std::ios::binary);need(bool(in),"open "+p.string());Hash hash;std::array<char,65536>buffer{};std::uintmax_t total=0;
 while(in){in.read(buffer.data(),buffer.size());const auto n=in.gcount();need(n>=0&&std::uintmax_t(n)<=cap-total,"file grew above cap");total+=n;hash.add(buffer.data(),n);}
 need(in.eof()&&!in.bad()&&total==size&&fs::file_size(p)==size,"file changed or read failed "+p.string());return hash.finish();
}
inline std::string read_text(const fs::path&p,std::size_t cap){ordinary_file(p);const auto size=fs::file_size(p);need(size<=cap,"text cap");std::string s(size,'\0');std::ifstream in(p,std::ios::binary);in.read(s.data(),s.size());need(bool(in)&&in.peek()==EOF,"short or growing text");return s;}
inline void sync(const fs::path&p,bool directory){const int fd=::open(p.c_str(),O_RDONLY|O_CLOEXEC|O_NOFOLLOW|(directory?O_DIRECTORY:0));need(fd>=0,"fsync open "+p.string());const int rc=::fsync(fd),saved=errno;::close(fd);need(rc==0,"fsync "+p.string()+": "+std::strerror(saved));}
inline void write_text(const fs::path&p,const std::string&bytes){need(!fs::exists(p),"refusing overwrite "+p.string());{std::ofstream out(p,std::ios::binary);out.write(bytes.data(),bytes.size());out.close();need(bool(out),"write "+p.string());}sync(p,false);}
inline Json file_table(const fs::path&dir,const Limits&limits={}){
 need(fs::is_directory(fs::symlink_status(dir)),"ordinary checkpoint directory required");
 std::vector<fs::path>files;std::uintmax_t total=0;
 for(const auto&e:fs::recursive_directory_iterator(dir)){
  const auto status=e.symlink_status();need(fs::is_regular_file(status)||fs::is_directory(status),"nonregular checkpoint entry "+e.path().string());
  if(!fs::is_regular_file(status)||e.path()==dir/"manifest.json")continue;
  need(files.size()<limits.files,"file count cap");const auto size=e.file_size();need(size<=limits.file_bytes&&size<=limits.total_bytes-total,"checkpoint byte cap");total+=size;files.push_back(e.path());
 }
 std::sort(files.begin(),files.end());Json table=Json::object();
 for(const auto&p:files)table[p.lexically_relative(dir).generic_string()]={{"bytes",fs::file_size(p)},{"sha256",hash_file(p,limits.file_bytes)}};
 return table;
}
inline void sync_tree(const fs::path&dir){std::vector<fs::path>dirs{dir};for(const auto&e:fs::recursive_directory_iterator(dir)){const auto st=e.symlink_status();if(fs::is_directory(st))dirs.push_back(e.path());else{ordinary_file(e.path());sync(e.path(),false);}}for(auto p=dirs.rbegin();p!=dirs.rend();++p)sync(*p,true);}
// The caller first serializes a coherent immutable accepted owner set into tmp.
// This layer verifies files/durability only; it cannot authorize owner restore.
inline fs::path commit(const fs::path&tmp,const fs::path&final,const Json&state,const Limits&limits={},void(*hook)(int)=nullptr){
 need(tmp.parent_path()==final.parent_path()&&tmp!=final&&!fs::exists(final),"fresh sibling checkpoint destination required");
 need(!fs::exists(tmp/"manifest.json")&&!fs::exists(tmp/"manifest.json.tmp"),"fresh unsealed temporary directory required");
 need(state.is_object(),"common state must be an object");
 Json body={{"schema",schema},{"state",state},{"files",file_table(tmp,limits)}};
 Json envelope={{"body",body},{"body_sha256",hash_text(body.dump())}};const auto encoded=envelope.dump(1)+"\n";
 need(encoded.size()<=limits.manifest_bytes,"manifest cap");if(hook)hook(1);
 write_text(tmp/"manifest.json.tmp",encoded);if(hook)hook(2);
 fs::rename(tmp/"manifest.json.tmp",tmp/"manifest.json");sync_tree(tmp);if(hook)hook(3);
 need(!fs::exists(final),"checkpoint destination appeared");fs::rename(tmp,final);sync(final.parent_path(),true);if(hook)hook(4);return final;
}
inline Json verify(const fs::path&dir,const Limits&limits={}){
 const auto envelope=Json::parse(read_text(dir/"manifest.json",limits.manifest_bytes));
 need(envelope.is_object()&&envelope.size()==2&&envelope.contains("body")&&envelope.contains("body_sha256"),"manifest envelope");
 const auto&body=envelope.at("body");need(hash_text(body.dump())==envelope.at("body_sha256").get<std::string>(),"manifest seal");
 need(body.is_object()&&body.size()==3&&body.at("schema")==schema&&body.at("state").is_object()&&body.at("files").is_object(),"manifest schema");
 need(body.at("files")==file_table(dir,limits),"checkpoint file set/digest differs");return body.at("state");
}
inline std::string step_name(std::uint64_t step){need(step<=999999999999ULL,"checkpoint step name range");auto value=std::to_string(step);return "ckpt-"+std::string(12-value.size(),'0')+value;}
inline fs::path newest_verified(const fs::path&root,const Limits&limits={}){
 need(fs::is_directory(fs::symlink_status(root)),"checkpoint root required");std::vector<fs::path>done;
 for(const auto&e:fs::directory_iterator(root)){const auto name=e.path().filename().string();if(name.rfind("ckpt-",0)!=0)continue;
  if(name.size()==21&&name.substr(17)==".tmp")continue;
  need(name.size()==17&&std::all_of(name.begin()+5,name.end(),[](char c){return c>='0'&&c<='9';})&&fs::is_directory(e.symlink_status()),"unexpected checkpoint entry "+name);done.push_back(e.path());}
 need(!done.empty(),"no committed checkpoint");std::sort(done.begin(),done.end());verify(done.back(),limits);return done.back();
}
} // namespace powerflow::common_checkpoint_files_v1
