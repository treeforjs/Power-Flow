#include "fixture.hpp"
int main(int argc,char** argv){try{
 if(argc!=3)return 2;
 const std::string input=argv[1];J result=J::array();
 for(const auto& pair:{std::pair<std::string,std::string>{"case-1.json","test/initial-joint-state"},{"neutral.json","neutral/initial"}}){m::GapPic pic(options());auto restored=m::PicArena::restore(pic,input+"/"+pair.first,states(),channels(),pair.second);result.push_back({{"file",pair.first},{"packets",restored->snapshot()->packets.size()},{"sinks",restored->snapshot()->sinks.size()},{"accepted",true}});}
 std::ofstream out(argv[2]);out<<result.dump(2)<<'\n';if(!out)return 3;return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
