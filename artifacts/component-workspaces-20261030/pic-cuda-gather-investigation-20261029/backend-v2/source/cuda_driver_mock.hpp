#pragma once
#include <cuda.h>
#include <map>
#include <string>
#include <vector>
namespace mock_cuda {
extern std::vector<std::string> trace;
extern std::map<std::string,int> calls;
void reset(bool current_primary=true);
void fail(const std::string&,int occurrence);
void set_current(CUcontext);
CUcontext primary();
CUcontext foreign();
int primary_refs();
std::size_t queued();
std::size_t live_allocations();
std::size_t unsafe_releases();
}
