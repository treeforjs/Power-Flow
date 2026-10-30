#include "collision_scheduler.hpp"
#include <chrono>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sys/resource.h>
namespace pc=raikou::pic_collisions;
std::uint64_t bits(double x){std::uint64_t v;std::memcpy(&v,&x,8);return v;}
void need(bool x,const char*m){if(!x)throw std::runtime_error(m);}
template<class T>void read(std::ifstream&in,T&x){in.read(reinterpret_cast<char*>(&x),sizeof x);need(bool(in),"truncated binary input");}
std::vector<pc::Cell> cells(const std::string&path){std::ifstream in(path,std::ios::binary);need(bool(in),"cell file missing");std::vector<pc::Cell>out;while(in.peek()!=EOF){need(out.size()<1000000,"cell input cap");pc::Cell c;read(in,c.id);read(in,c.volume_m3);out.push_back(c);}return out;}
std::vector<pc::Packet> packets(const std::string&path){std::ifstream in(path,std::ios::binary);need(bool(in),"particle file missing");std::vector<pc::Packet>out;while(in.peek()!=EOF){need(out.size()<100000,"particle input cap");std::uint64_t id;pc::Packet p;read(in,id);read(in,p.cell);read(in,p.weight);read(in,p.mass_kg);for(auto&u:p.proper_velocity_m_s)read(in,u);p.id="actual-like/e/"+std::to_string(id);p.state="electron";out.push_back(p);}return out;}
bool same(const pc::Candidate&a,const pc::Candidate&b){
 if(std::tie(a.projectile_id,a.target_id,a.channel_id,a.data_sha256,a.kernel_sha256,a.cell)!=std::tie(b.projectile_id,b.target_id,b.channel_id,b.data_sha256,b.kernel_sha256,b.cell))return false;
 const double av[]={a.event_weight,a.projectile_before,a.target_before,a.projectile_remainder.high,a.projectile_remainder.low,a.target_remainder.high,a.target_remainder.low,a.rate_s,a.energy_eV,a.flux_m_s,a.sigma_m2,a.volume_m3};
 const double bv[]={b.event_weight,b.projectile_before,b.target_before,b.projectile_remainder.high,b.projectile_remainder.low,b.target_remainder.high,b.target_remainder.low,b.rate_s,b.energy_eV,b.flux_m_s,b.sigma_m2,b.volume_m3};
 for(int i=0;i<12;++i)if(bits(av[i])!=bits(bv[i]))return false;
 return true;
}
void run(const char*name,const std::vector<pc::Packet>&e,const std::vector<pc::Packet>&n,const std::vector<pc::Cell>&full,const std::vector<pc::Cell>&sparse){
 std::vector<std::uint64_t> queries;pc::Channel ch{"manufactured-constant","electron","manufactured-target",std::string(64,'a'),std::string(64,'b'),pc::Observable::ExclusiveEvent,[&](double E){queries.push_back(bits(E));return 1e-20;}};
 pc::CandidateWork wf,ws;auto start=std::chrono::steady_clock::now();auto first=pc::candidates(e,n,full,{ch},1000000,1000000,1000000,{},&wf);const auto full_s=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();auto expected=queries;queries.clear();start=std::chrono::steady_clock::now();auto second=pc::candidates(e,n,sparse,{ch},1000000,1000000,1000000,{},&ws);const auto sparse_s=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();need(first.size()==second.size()&&queries==expected,"candidate count/query order differs");
 double r1=0,r2=0,w1=0,w2=0;for(std::size_t i=0;i<first.size();++i){need(same(first[i],second[i]),"candidate field bit mismatch");r1+=first[i].rate_s;r2+=second[i].rate_s;w1+=first[i].rate_s*first[i].event_weight;w2+=second[i].rate_s*second[i].event_weight;}need(bits(r1)==bits(r2)&&bits(w1)==bits(w2),"total or weighted rate differs");need(wf.pair_visits==ws.pair_visits&&wf.channel_queries==ws.channel_queries,"sparse view changed pair/query work");
 auto missing=sparse;need(!missing.empty(),"missing occupied control empty");missing.erase(missing.begin());queries.clear();bool refused=false;try{pc::candidates(e,n,missing,{ch});}catch(const pc::Error&){refused=true;}need(refused&&queries.empty(),"unknown occupied bin not refused before query");
 std::cout<<"{\"case\":\""<<name<<"\",\"status\":\"PASS\",\"projectiles\":"<<e.size()<<",\"manufactured_targets\":"<<n.size()<<",\"exhaustive_cells\":"<<full.size()<<",\"occupied_cells\":"<<sparse.size()<<",\"candidate_rows\":"<<first.size()<<",\"pair_visits\":"<<wf.pair_visits<<",\"query_calls\":"<<expected.size()<<",\"exhaustive_wall_s\":"<<full_s<<",\"sparse_wall_s\":"<<sparse_s<<",\"candidate_fields_and_query_order_bitwise_equal\":true,\"total_and_weighted_rate_bitwise_equal\":true,\"unknown_occupied_bin_refused_before_query\":true,\"scheduler_event_selected_or_committed\":false}"<<'\n';
}
int main(int argc,char**argv){try{need(argc==2,"expected actual-like input directory");std::cout<<std::setprecision(17);
 std::vector<pc::Cell>full,sparse;for(std::uint64_t i=0;i<128;++i)full.push_back({i,1e-12*(1+i%7)});for(auto i:{5,17,29,61,127})sparse.push_back(full.at(i));std::vector<pc::Packet>e,n;
 for(std::size_t i=0;i<64;++i)e.push_back({"manufactured/e/"+std::to_string(i),"electron",sparse[i%sparse.size()].id,(i%3==0)?.25:1.,9.1093837139e-31,{1e6+i*1000.,0,0}});
 for(std::size_t i=0;i<11;++i)n.push_back({"manufactured/n/"+std::to_string(i),"manufactured-target",sparse[i%sparse.size()].id,8.,2.9915e-26,{0,i*100.,0}});
 run("small-variable-volume",e,n,full,sparse);
 const std::string dir=argv[1];full=cells(dir+"/atlas-cells.bin");sparse=cells(dir+"/occupied-cells.bin");e=packets(dir+"/projectiles.bin");n.clear();for(const auto&cell:sparse)n.push_back({"manufactured/n/"+std::to_string(cell.id),"manufactured-target",cell.id,8.,2.9915e-26,{0,0,0}});run("preserved-electron-layout-manufactured-targets",e,n,full,sparse);
 struct rusage usage{};getrusage(RUSAGE_SELF,&usage);std::cout<<"{\"process_peak_RSS_KiB\":"<<usage.ru_maxrss<<",\"sizeof_Cell\":"<<sizeof(pc::Cell)<<",\"physical_evolution\":false,\"GPU\":false}"<<'\n';
 }catch(const std::exception&error){std::cerr<<error.what()<<'\n';return 1;}}
