#pragma once
#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
namespace pic_cudss_tuning {
struct Structure {
 std::size_t isolated=0,components=0,largest_component=0,max_degree=0,bandwidth=0;
 Structure(const std::vector<int>& rows,const std::vector<int>& cols) {
  const std::size_t n=rows.size()-1;std::vector<std::vector<int>> adjacent(n);
  for(std::size_t i=0;i<n;++i)for(int k=rows[i];k<rows[i+1];++k) {
   const int j=cols[k];if(j==int(i))continue;
   adjacent[i].push_back(j);adjacent[j].push_back(int(i));bandwidth=std::max(bandwidth,std::size_t(j-int(i)));
  }
  std::vector<unsigned char> seen(n);std::vector<int> queue;queue.reserve(n);
  for(std::size_t i=0;i<n;++i) {
   max_degree=std::max(max_degree,adjacent[i].size());isolated+=adjacent[i].empty();if(seen[i])continue;
   ++components;queue.clear();queue.push_back(int(i));seen[i]=1;
   for(std::size_t p=0;p<queue.size();++p)for(int j:adjacent[queue[p]])if(!seen[j]){seen[j]=1;queue.push_back(j);}
   largest_component=std::max(largest_component,queue.size());
  }
 }
 void json(std::ostream& out) const {
  out<<"\"isolated_diagonal_rows\":"<<isolated<<",\"connected_components\":"<<components
     <<",\"largest_component_rows\":"<<largest_component<<",\"max_offdiagonal_degree\":"<<max_degree<<",\"original_index_bandwidth\":"<<bandwidth;
 }
};
template<class T>void binary(const std::string& path,const std::vector<T>& values) {
 if(std::filesystem::exists(path))throw std::runtime_error("Audit output already exists: "+path);
 std::ofstream out(path,std::ios::binary);out.write(reinterpret_cast<const char*>(values.data()),values.size()*sizeof(T));
 out.close();if(!out)throw std::runtime_error("Cannot write audit array: "+path);
}
inline void matrix_arrays(const std::string& prefix,const std::vector<int>& rows,const std::vector<int>& cols,const std::vector<double>& values) {
 binary(prefix+".rows.i32",rows);binary(prefix+".cols.i32",cols);binary(prefix+".values.f64",values);
}
}
