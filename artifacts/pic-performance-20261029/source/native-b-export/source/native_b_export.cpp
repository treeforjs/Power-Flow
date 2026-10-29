// Compile the exact sealed samplers, retaining no live reference to its driver.
#define main sealed_pic_driver_main_unused
#include "/home/treefor/projects/runtimes/mfem-gpu-pic-20261029/source/examples/mfem_heating_history_pic.cpp"
#undef main
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <stdexcept>

namespace {
struct Frame {
 int index;std::filesystem::path source;
};
double frame_time(const std::filesystem::path& source) {
 std::ifstream input(source);std::string line;std::getline(input,line);
 if(!std::getline(input,line))throw std::runtime_error("Frame has no samples");
 return std::stod(split(line).at(0));
}
void verify(Fields& fields,const std::filesystem::path& directory,double reference_time,
            const std::filesystem::path& output) {
 fields.set(reference_time);const int interval=fields.current;
 const double t0=fields.times[interval],t1=fields.times[interval+1],mid=t0+.5*(t1-t0);
 NativeField lower(directory,t0),upper(directory,t1);
 const std::array<std::array<double,2>,7> points{{
  {{.0000125,.0000125}},{{.0089875,.0000125}},{{.0000125,.0089875}},
  {{.0089875,.0089875}},{{.0045125,.0045125}},{{.0025125,.0035125}},
  {{.0065125,.0055125}}}};
 double max_error=0,max_value=0;std::size_t comparisons=0;
 std::ofstream report(output);report<<std::setprecision(17)
  <<"{\"method\":\"Original Fields::value compared with direct original NativeField::magnetic at both endpoints and the interval midpoint\","
  <<"\"lower_time_s\":"<<t0<<",\"upper_time_s\":"<<t1<<",\"midpoint_time_s\":"<<mid
  <<",\"samples\":[";
 bool first=true;
 for(double time:{t0,mid,t1}) {
  const double q=(time-t0)/(t1-t0);
  for(const auto& point:points) {
   auto a=lower.magnetic(point[0],point[1]),b=upper.magnetic(point[0],point[1]);
   auto actual=fields.value(point[0],point[1],time,true);
   PicVector expected{};
   for(int component=0;component<3;++component) {
    expected[component]=(1-q)*a[component]+q*b[component];
    if(!std::isfinite(actual[component]) || !std::isfinite(expected[component]))
     throw std::runtime_error("Nonfinite verification sample");
    const double error=std::abs(actual[component]-expected[component]);
    max_error=std::max(max_error,error);max_value=std::max(max_value,std::abs(expected[component]));
    if(error>1e-13*std::max(1.,std::abs(expected[component])))
     throw std::runtime_error("Original Fields and direct NativeField endpoint/midpoint mismatch");
    ++comparisons;
   }
   if(!first)report<<',';first=false;
   report<<"{\"time_s\":"<<time<<",\"x_m\":"<<point[0]<<",\"y_m\":"<<point[1]<<",\"actual_T\":["
         <<actual[0]<<','<<actual[1]<<','<<actual[2]<<"],\"expected_T\":["
         <<expected[0]<<','<<expected[1]<<','<<expected[2]<<"]}";
  }
 }
 report<<"],\"component_comparisons\":"<<comparisons<<",\"max_absolute_error_T\":"<<max_error
       <<",\"max_absolute_reference_T\":"<<max_value<<",\"passed\":true}\n";
 report.close();if(!report)throw std::runtime_error("Verification report write failed");
}
void export_frame(Fields& fields,const Frame& frame,const std::filesystem::path& directory) {
 const std::string stem="frame_"+std::to_string(frame.index)+"_B";
 const auto temporary=directory/(stem+".csv.tmp"),output=directory/(stem+".csv");
 std::ifstream input(frame.source);std::ofstream csv(temporary);std::string line;
 if(!input || !csv)throw std::runtime_error("Cannot open frame input or magnetic export");
 std::getline(input,line);
 const auto columns=split(line);
 if(columns.size()<3 || columns[0]!="time_s" || columns[1]!="x_m" || columns[2]!="y_m")
  throw std::runtime_error("Unexpected frame columns");
 csv<<std::setprecision(17)<<"time_s,time_ns,x_m,y_m,Bx_T,By_T,Bz_T,Bmag_T\n";
 std::array<double,4> minimum,maximum;
 minimum.fill(std::numeric_limits<double>::infinity());
 maximum.fill(-std::numeric_limits<double>::infinity());
 std::size_t rows=0;double time=0;
 while(std::getline(input,line)) {
  if(line.empty())throw std::runtime_error("Empty row in native frame");
  const auto values=split(line);
  const double current_time=std::stod(values.at(0)),x=std::stod(values.at(1)),y=std::stod(values.at(2));
  if(rows==0)time=current_time;
  if(current_time!=time || !std::isfinite(time) || !std::isfinite(x) || !std::isfinite(y))
   throw std::runtime_error("Inconsistent or nonfinite frame time/coordinates");
  if(rows>=360ULL*360 || std::abs(x-(double(rows%360)+.5)*(.009/360))>1e-15 ||
     std::abs(y-(double(rows/360)+.5)*(.009/360))>1e-15)
   throw std::runtime_error("Frame is not the original 360x360 snapshot cell-centre grid");
  auto magnetic=fields.value(x,y,time,true);
  const double norm=std::hypot(magnetic[0],magnetic[1],magnetic[2]);
  const std::array<double,4> field{magnetic[0],magnetic[1],magnetic[2],norm};
  for(int component=0;component<4;++component) {
   if(!std::isfinite(field[component]))throw std::runtime_error("Nonfinite magnetic export");
   minimum[component]=std::min(minimum[component],field[component]);
   maximum[component]=std::max(maximum[component],field[component]);
  }
  csv<<time<<','<<time*1e9<<','<<x<<','<<y<<','<<field[0]<<','<<field[1]<<','<<field[2]<<','<<field[3]<<'\n';
  ++rows;
 }
 if(rows!=360ULL*360)throw std::runtime_error("Native frame has an incomplete grid");
 csv.close();if(!csv)throw std::runtime_error("Magnetic CSV write failed");
 std::filesystem::rename(temporary,output);
 std::ofstream report(directory/(stem+".json"));
 report<<std::setprecision(17)<<"{\"frame_index\":"<<frame.index<<",\"source_frame\":"
       <<std::quoted(frame.source.string())<<",\"time_s\":"<<time<<",\"time_ns\":"<<time*1e9
       <<",\"rows\":"<<rows<<",\"nx\":360,\"ny\":360,\"unit\":\"Tesla\",\"all_finite\":true,"
       <<"\"native_interval_index\":"<<fields.current<<",\"native_lower_time_s\":"<<fields.times[fields.current]
       <<",\"native_upper_time_s\":"<<fields.times[fields.current+1]<<",\"native_fraction\":"<<fields.fraction
       <<",\"components\":{";
 const std::array<const char*,4> names{"Bx_T","By_T","Bz_T","Bmag_T"};
 for(int i=0;i<4;++i){if(i)report<<',';report<<std::quoted(names[i])<<":{\"min\":"<<minimum[i]<<",\"max\":"<<maximum[i]<<'}';}
 report<<"},\"csv\":"<<std::quoted(output.string())<<"}\n";
 report.close();if(!report)throw std::runtime_error("Magnetic stats write failed");
 std::cout<<std::setprecision(17)<<"frame "<<frame.index<<" time_ns="<<time*1e9
          <<" rows="<<rows<<" |B|_T=["<<minimum[3]<<','<<maximum[3]<<"] "<<output<<'\n'<<std::flush;
}
}
int main(int argc,char** argv) {
 try {
  if(argc!=4)throw std::runtime_error("Usage: native_b_export FIELD_DIR JOBS_TSV OUTPUT_DIR");
  const std::filesystem::path directory(argv[1]),output(argv[3]);
  std::ifstream jobs(argv[2]);std::string line;std::vector<Frame> frames;
  while(std::getline(jobs,line)) {
   const auto tab=line.find('\t');if(tab==std::string::npos)throw std::runtime_error("Invalid export job");
   frames.push_back({std::stoi(line.substr(0,tab)),line.substr(tab+1)});
  }
  if(frames.empty())throw std::runtime_error("No native magnetic export jobs");
  Fields fields(directory);
  verify(fields,directory,frame_time(frames.front().source),output/"endpoint_midpoint_verification.json");
  for(const auto& frame:frames)export_frame(fields,frame,output);
  return 0;
 } catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
