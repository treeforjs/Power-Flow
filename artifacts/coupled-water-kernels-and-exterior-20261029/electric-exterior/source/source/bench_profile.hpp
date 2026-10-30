#pragma once
#include <chrono>
#include <fstream>
#include <iomanip>
#include <map>
#include <string>
#include <sys/resource.h>
#include <sched.h>
namespace pic_bench {
using Clock=std::chrono::steady_clock;
struct Entry { double seconds=0; unsigned long long calls=0; };
inline std::map<std::string,Entry> entries;
struct Timer {
 std::string name; Clock::time_point start=Clock::now();
 explicit Timer(std::string n):name(std::move(n)){}
 ~Timer(){auto& e=entries[name];e.seconds+=std::chrono::duration<double>(Clock::now()-start).count();++e.calls;}
};
inline void resources(const std::string& path) {
 struct rusage r{};getrusage(RUSAGE_SELF,&r);cpu_set_t cpus;CPU_ZERO(&cpus);sched_getaffinity(0,sizeof(cpus),&cpus);
 std::ofstream f(path);f<<std::setprecision(17)<<"{\"user_cpu_s\":"<<r.ru_utime.tv_sec+r.ru_utime.tv_usec*1e-6<<",\"system_cpu_s\":"<<r.ru_stime.tv_sec+r.ru_stime.tv_usec*1e-6<<",\"max_rss_KiB\":"<<r.ru_maxrss<<",\"allowed_logical_cpus\":[";
 bool first=true;for(int i=0;i<CPU_SETSIZE;++i)if(CPU_ISSET(i,&cpus)){if(!first)f<<',';first=false;f<<i;}f<<"]}\n";
}
inline void write(const std::string& path) {
 std::ofstream f(path);f<<std::setprecision(17)<<"{\n";bool first=true;
 for(const auto& [name,e]:entries){if(!first)f<<",\n";first=false;f<<'"'<<name<<"\":{\"seconds\":"<<e.seconds<<",\"calls\":"<<e.calls<<'}';}
 f<<"\n}\n";
}
}
