#include "collision_scheduler.hpp"
#include "legacy_collision_scheduler.hpp"
#include <bit>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sys/resource.h>
#include <sched.h>
namespace now=raikou::pic_collisions;
namespace old=raikou::pic_collisions_legacy;
using Clock=std::chrono::steady_clock;
constexpr double me=9.1093837139e-31,mn=2.9915e-26;
const std::string hash(64,'a'),kernel(64,'b');
int checks=0,comparisons=0;
void check(bool ok,const char* msg){++checks;if(!ok)throw std::runtime_error(msg);}
struct Def {std::string id,p="electron",n="H2O";now::Observable kind=now::Observable::ExclusiveEvent;std::function<double(double)> sigma;};
struct Input {std::vector<now::Packet> p,n;std::vector<now::Cell> cells;std::vector<Def> channels;};
using Trace=std::vector<std::pair<std::string,std::uint64_t>>;
std::uint64_t bits(double x){return std::bit_cast<std::uint64_t>(x);}
template<class C>std::string bytes(const C& r){
    std::string result;
    auto str=[&](const std::string& s){const std::uint64_t n=s.size();result.append(reinterpret_cast<const char*>(&n),sizeof n);result+=s;};
    for(const auto* s:{&r.projectile_id,&r.target_id,&r.channel_id,&r.data_sha256,&r.kernel_sha256})str(*s);
    result.append(reinterpret_cast<const char*>(&r.cell),sizeof r.cell);
    for(double x:{r.event_weight,r.projectile_before,r.target_before,r.projectile_remainder.high,r.projectile_remainder.low,r.target_remainder.high,r.target_remainder.low,r.rate_s,r.energy_eV,r.flux_m_s,r.sigma_m2,r.volume_m3}){const auto b=bits(x);result.append(reinterpret_cast<const char*>(&b),sizeof b);}
    return result;
}
template<class T>std::uint64_t fingerprint(const std::vector<T>& rows){std::uint64_t h=1469598103934665603ULL;for(const auto& r:rows)for(unsigned char c:bytes(r)){h^=c;h*=1099511628211ULL;}return h;}
std::vector<now::Channel> channels_now(const Input& x,Trace* trace=nullptr){
 std::vector<now::Channel> result;for(const auto& d:x.channels){auto fn=d.sigma;auto id=d.id;result.push_back({d.id,d.p,d.n,hash,kernel,d.kind,[fn,id,trace](double e){if(trace)trace->emplace_back(id,bits(e));return fn(e);}});}return result;
}
struct LegacyInput {std::vector<old::Packet> p,n;std::vector<old::Cell> cells;std::vector<old::Channel> channels;};
LegacyInput legacy(const Input& x,Trace* trace=nullptr){
 LegacyInput result;for(const auto* list:{&x.p,&x.n})for(const auto& r:*list){auto& dest=list==&x.p?result.p:result.n;dest.push_back({r.id,r.state,r.cell,r.weight,r.mass_kg,r.proper_velocity_m_s});}
 for(const auto& c:x.cells)result.cells.push_back({c.id,c.volume_m3});
 for(const auto& d:x.channels){auto fn=d.sigma;auto id=d.id;result.channels.push_back({d.id,d.p,d.n,hash,kernel,static_cast<old::Observable>(d.kind),[fn,id,trace](double e){if(trace)trace->emplace_back(id,bits(e));return fn(e);}});}return result;
}
struct Answer {std::vector<now::Candidate> rows;std::string error;Trace trace;now::CandidateWork work;};
Answer run(const Input& x,std::size_t pair=1000000,std::size_t query=1000000,std::size_t row=1000000,const now::CandidateInputLimits& caps={}){
 Answer a;try{a.rows=now::candidates(x.p,x.n,x.cells,channels_now(x,&a.trace),pair,query,row,caps,&a.work);}catch(const std::exception& e){a.error=e.what();}return a;
}
Answer compare(const Input& x,std::size_t pair=1000000,std::size_t query=1000000,std::size_t row=1000000){
 auto current=run(x,pair,query,row);Trace trace;auto previous=legacy(x,&trace);std::vector<old::Candidate> rows;std::string error;
 try{rows=old::candidates(previous.p,previous.n,previous.cells,previous.channels,pair,query,row);}catch(const std::exception& e){error=e.what();}
 check(error==current.error,"legacy/indexed exception mismatch");check(trace==current.trace,"legacy/indexed support query sequence mismatch");check(rows.size()==current.rows.size(),"legacy/indexed candidate count mismatch");
 for(std::size_t i=0;i<rows.size();++i)check(bytes(rows[i])==bytes(current.rows[i]),"legacy/indexed exact candidate byte mismatch");
 ++comparisons;return current;
}
double proper(double eV){const double x=eV*now::qe/(me*now::c*now::c);return now::c*std::sqrt(x*(x+2));}
Input basic(){Input x;x.cells={{7,1e-12}};for(int i=0;i<2;++i)x.p.push_back({"e"+std::to_string(i),"electron",7,.031334406939830223*(i+1),me,{proper(100.+10*i),0,0}});for(int i=0;i<3;++i)x.n.push_back({"n"+std::to_string(i),"H2O",7,1e20,mn,{0,0,0}});x.channels={{"z","electron","H2O",now::Observable::ExclusiveEvent,[](double){return 1e-20;}},{"a","electron","H2O",now::Observable::ExclusiveEvent,[](double){return 2e-20;}}};return x;}
std::size_t input_bytes(const Input& x){std::size_t n=0;for(const auto* list:{&x.p,&x.n})for(const auto& r:*list)n+=r.id.size()+r.state.size();for(const auto& c:x.channels)n+=c.id.size()+c.p.size()+c.n.size()+hash.size()+kernel.size();return n;}
void correctness(){
 std::mt19937_64 rng(0x873218982ULL);
 for(int test=0;test<300;++test){
  Input x;const int nc=1+rng()%8,np=rng()%13,nn=rng()%14,nch=rng()%9;
  for(int i=0;i<nc;++i)x.cells.push_back({static_cast<std::uint64_t>(101+3*i),1e-12*(i+1)});
  for(int i=0;i<np;++i)x.p.push_back({"e"+std::to_string(i),i%3?"electron":"beam",x.cells[rng()%nc].id,std::array<double,4>{.25,.031334406939830223,64.,1e20}[rng()%4],me,{proper(14.+rng()%900),double(rng()%1000)-500,double(rng()%1000)-500}});
  for(int i=0;i<nn;++i)x.n.push_back({"n"+std::to_string(i),i%3?"H2O":"H",x.cells[rng()%nc].id,std::array<double,4>{1.,.031334406939830223,64.,1e20}[rng()%4],mn,{double(rng()%1000)-500,double(rng()%1000)-500,0}});
  for(int i=0;i<nch;++i)x.channels.push_back({"ch"+std::to_string(i),i%3?"electron":"beam",i%2?"H2O":"H",now::Observable::ExclusiveEvent,[i](double e){return i%3==0?0.:(e<30?0.:(1+i)*1e-20);}});
  const auto a=compare(x);std::shuffle(x.p.begin(),x.p.end(),rng);std::shuffle(x.n.begin(),x.n.end(),rng);std::shuffle(x.cells.begin(),x.cells.end(),rng);std::shuffle(x.channels.begin(),x.channels.end(),rng);const auto b=compare(x);
  check(a.error.empty()&&b.error.empty(),"valid seeded input rejected");check(fingerprint(a.rows)==fingerprint(b.rows),"input permutation changes canonical candidates");
 }
 const auto x=basic();const auto base=compare(x,6,12,12);check(base.rows.size()==12,"dense independent row count");check(base.work.pair_visits==6&&base.work.channel_queries==12&&base.work.candidate_rows==12,"independent dense work count");
 Input unusual=x;unusual.p[0].id=std::string("e\0inside",8);unusual.n[0].id=std::string(1300,'q');unusual.channels[0].id=std::string("ch\0inside",9);check(compare(unusual).rows.size()==12,"embedded NUL/long identities change candidates");
 check(compare(x,5,12,12).error=="collision pair work limit","pair exact boundary");check(compare(x,6,11,12).error=="collision query work limit","query exact boundary");check(compare(x,6,12,11).error=="collision candidate allocation limit","row exact boundary");
 compare(x,SIZE_MAX,SIZE_MAX,SIZE_MAX);
 Input same=x;same.p.resize(1);same.n.resize(1);same.n[0].proper_velocity_m_s=same.p[0].proper_velocity_m_s;same.channels.resize(1);same.channels[0].sigma=[](double e)->double{if(e==0)throw std::runtime_error("manufactured support gap");return 1e-20;};
 auto zero=compare(same);check(zero.error=="manufactured support gap"&&zero.trace.size()==1,"zero flux support query retained");
 same.channels[0].sigma=[](double){return 0.;};check(compare(same).rows.empty(),"zero cross section preserved");
 same=x;same.channels.clear();same.p[0].proper_velocity_m_s={1e308,0,0};same.n[0].proper_velocity_m_s={-1e308,0,0};check(compare(same).error=="relative kinematics overflow","unmatched pair kinematic rejection retained");
 auto invalid=[&](auto mutation){Input z=x;mutation(z);const auto a=compare(z);check(!a.error.empty(),"invalid legacy input silently admitted");};
 invalid([](Input& z){z.cells.push_back(z.cells[0]);});invalid([](Input& z){z.p[0].id=z.n[0].id;});invalid([](Input& z){z.p[0].cell=999;});invalid([](Input& z){z.cells[0].volume_m3=0;});invalid([](Input& z){z.p[0].mass_kg=0;});invalid([](Input& z){z.p[0].proper_velocity_m_s[0]=NAN;});invalid([](Input& z){z.channels[0].kind=now::Observable::InclusiveProductYield;});invalid([](Input& z){z.channels[0].sigma=[](double){return -1.;};});invalid([](Input& z){z.channels[0].sigma=[](double){return NAN;};});
 now::CandidateInputLimits caps;caps.projectiles=x.p.size();caps.targets=x.n.size();caps.cells=x.cells.size();caps.channels=x.channels.size();caps.input_string_bytes=input_bytes(x);
 std::size_t payload=0;for(const auto& r:base.rows)payload+=5+r.projectile_id.size()+r.target_id.size()+r.channel_id.size()+r.data_sha256.size()+r.kernel_sha256.size();caps.candidate_string_bytes=payload;
 check(run(x,6,12,12,caps).rows.size()==12,"exact input/payload caps rejected");
 auto early=[&](const now::CandidateInputLimits& c,const char* expected){auto a=run(x,6,12,12,c);check(a.error==expected,"input cap rejection mismatch");check(a.work.cell_lookups==0&&a.work.target_index_entries==0&&a.work.channel_queries==0,"input cap did work before rejection");};
 auto c=caps;c.projectiles--;early(c,"collision projectile input limit");c=caps;c.targets--;early(c,"collision target input limit");c=caps;c.cells--;early(c,"collision cell input limit");c=caps;c.channels--;early(c,"collision channel input limit");c=caps;c.input_string_bytes--;early(c,"collision input string byte limit");
 c=caps;c.candidate_string_bytes--;auto limited=run(x,6,12,12,c);check(limited.error=="collision candidate string byte limit"&&limited.work.candidate_rows==11,"payload budget silently truncates or overallocates");
 c=caps;c.candidate_string_bytes=4;check(run(x,6,12,12,c).error=="collision candidate string byte limit","sub-terminator payload cap underflow");
 Input empty;now::CandidateInputLimits none{0,0,0,0,0,0};check(run(empty,0,0,0,none).error.empty(),"zero-cap empty input rejected");
}
struct Bench {std::string label;std::size_t np,nn,cells,channels,rows;std::uint64_t cartesian;double indexed_ms,legacy_ms;now::CandidateWork work;std::uint64_t checksum;};
Bench bench(std::string label,const Input& x,bool run_legacy){
 auto previous=legacy(x);auto current_channels=channels_now(x);now::CandidateWork work;const auto counted=now::candidates(x.p,x.n,x.cells,current_channels,SIZE_MAX,SIZE_MAX,1000000,{},&work);const auto want=fingerprint(counted);const auto rows=counted.size();
 auto measure=[&](bool historical){std::vector<double> times;for(int repeat=0;repeat<3;++repeat){const auto start=Clock::now();std::uint64_t got=0;if(historical){auto r=old::candidates(previous.p,previous.n,previous.cells,previous.channels,SIZE_MAX,SIZE_MAX,1000000);const auto stop=Clock::now();times.push_back(std::chrono::duration<double,std::milli>(stop-start).count());got=fingerprint(r);}else{auto r=now::candidates(x.p,x.n,x.cells,current_channels,SIZE_MAX,SIZE_MAX,1000000);const auto stop=Clock::now();times.push_back(std::chrono::duration<double,std::milli>(stop-start).count());got=fingerprint(r);}check(got==want,"benchmark exact candidate checksum mismatch");}std::sort(times.begin(),times.end());return times[1];};
 const double old_ms=run_legacy?measure(true):-1;const double new_ms=measure(false);
 return {std::move(label),x.p.size(),x.n.size(),x.cells.size(),x.channels.size(),rows,static_cast<std::uint64_t>(x.p.size())*x.n.size(),new_ms,old_ms,work,want};
}
Input disjoint(int n){Input x;x.cells={{1,1e-12},{2,1e-12}};for(int i=0;i<n;++i){x.p.push_back({"e"+std::to_string(i),"electron",1,1,me,{proper(100),0,0}});x.n.push_back({"n"+std::to_string(i),"H2O",2,64,mn,{0,0,0}});}x.channels={{"parent","electron","H2O",now::Observable::ExclusiveEvent,[](double){return 1e-20;}}};return x;}
Input matched(int n){auto x=disjoint(n);x.cells.clear();for(int i=0;i<n;++i){const auto cell=static_cast<std::uint64_t>(i);x.cells.push_back({cell,1e-12*(1+i%7)});x.p[i].cell=cell;x.n[i].cell=cell;}return x;}
Input dense(int n,int channels){auto x=disjoint(n);x.cells.resize(1);for(auto& target:x.n)target.cell=1;for(int i=1;i<channels;++i)x.channels.push_back({"ch"+std::to_string(i),"electron","H2O",now::Observable::ExclusiveEvent,[i](double){return (i+1)*1e-20;}});return x;}
void emit(std::ostream& o,const Bench& r){const auto& w=r.work;o<<"{\"label\":\""<<r.label<<"\",\"projectiles\":"<<r.np<<",\"targets\":"<<r.nn<<",\"cells\":"<<r.cells<<",\"channels\":"<<r.channels<<",\"rows\":"<<r.rows<<",\"cartesian_comparisons_legacy\":"<<r.cartesian<<",\"indexed_median_ms\":"<<r.indexed_ms<<",\"legacy_median_ms\":"<<r.legacy_ms<<",\"checksum\":"<<r.checksum<<",\"work\":{\"cell_lookups\":"<<w.cell_lookups<<",\"target_index_entries\":"<<w.target_index_entries<<",\"channel_index_entries\":"<<w.channel_index_entries<<",\"pair_visits\":"<<w.pair_visits<<",\"kinematics_calls\":"<<w.kinematics_calls<<",\"channel_queries\":"<<w.channel_queries<<",\"candidate_rows\":"<<w.candidate_rows<<",\"candidate_string_bytes\":"<<w.candidate_string_bytes<<"}}";}
int main(int argc,char** argv){try{
 check(argc==2,"result path required");cpu_set_t affinity;CPU_ZERO(&affinity);check(sched_getaffinity(0,sizeof affinity,&affinity)==0&&CPU_COUNT(&affinity)==1&&CPU_ISSET(7,&affinity),"CPU7 required");check(getpriority(PRIO_PROCESS,0)==15,"nice15 required");
 correctness();std::vector<Bench> records;
 for(int n:{1000,5000,10000})records.push_back(bench("disjoint_"+std::to_string(n),disjoint(n),true));
 auto large=disjoint(100000);auto sparse=run(large,0,0,0);check(sparse.error.empty()&&sparse.rows.empty()&&sparse.work.pair_visits==0&&sparse.work.kinematics_calls==0&&sparse.work.channel_queries==0,"large disjoint case exceeds zero local work budget");records.push_back(bench("disjoint_100000",large,false));
 for(int n:{1000,5000})records.push_back(bench("matched_"+std::to_string(n),matched(n),true));
 records.push_back(bench("dense_128_two_channels",dense(128,2),true));
 records.push_back(bench("dense_256_two_channels",dense(256,2),true));
 auto many=dense(16,1);for(int i=0;i<4096;++i)many.channels.push_back({"unused"+std::to_string(i),"unused"+std::to_string(i),"H2O",now::Observable::ExclusiveEvent,[](double){return 1e-20;}});records.push_back(bench("many_unmatched_channels",many,true));
 auto many_large=dense(64,1);many_large.channels=many.channels;records.push_back(bench("many_unmatched_channels_64",many_large,true));
 struct rusage usage{};check(getrusage(RUSAGE_SELF,&usage)==0,"resource measurement");std::ofstream o(argv[1]);o<<std::setprecision(17)<<"{\"passed\":true,\"checks\":"<<checks<<",\"exact_legacy_comparisons\":"<<comparisons<<",\"randomized_cases\":300,\"CPU\":7,\"nice\":15,\"threads\":1,\"max_rss_KiB\":"<<usage.ru_maxrss<<",\"user_cpu_s\":"<<usage.ru_utime.tv_sec+1e-6*usage.ru_utime.tv_usec<<",\"system_cpu_s\":"<<usage.ru_stime.tv_sec+1e-6*usage.ru_stime.tv_usec<<",\"benchmarks\":[";for(std::size_t i=0;i<records.size();++i){if(i)o<<',';emit(o,records[i]);}o<<"],\"physical_kernel_change\":false,\"native_mutation\":false}\n";if(!o)throw std::runtime_error("result write");std::cout<<checks<<" independent indexed checks passed, "<<comparisons<<" exact legacy comparisons\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
