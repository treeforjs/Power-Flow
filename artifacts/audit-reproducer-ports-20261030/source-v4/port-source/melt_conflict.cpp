// Cross-component check: can the exact source stock produced by owner v4 for the
// audit A07/C1 input be consumed by the SS304 melt/ablation proposal v1?
// The melt proposal's own fixture helpers (scenario, view, action, request) are
// reused unchanged by including its test file with its entry point renamed.
// Built with the melt proposal's own compiler flags.
//
// usage: melt_conflict <v4-accepted-stock-carriers.txt> <output-directory>
#define main melt_proposal_test_main
#include "test_proposal.cpp"
#undef main
#include <cstdlib>

namespace {
std::vector<double> read_carriers(const char* path){
 std::ifstream in(path);a::need(bool(in),"carrier file");std::vector<double> w;std::string line;
 while(std::getline(in,line))if(!line.empty())w.push_back(std::strtod(line.c_str(),nullptr));
 return w;
}
J attempt_stocks(const a::Evaluator& e,const std::vector<a::Stock>& stocks){
 J r;try{(void)e.energy_at(stocks,400.);r["admitted"]=true;}catch(const std::exception& x){r["admitted"]=false;r["message"]=x.what();}return r;
}
}

int main(int argc,char** argv){
 try{
  a::need(argc==3,"usage: melt_conflict <carriers> <outdir>");std::filesystem::create_directories(argv[2]);
  const auto carriers=read_carriers(argv[1]);const auto v4stock=Count::from_positive_carriers(carriers);
  auto s=scenario();a::Evaluator e(s);J out;
  out["v4_stock_carriers"]=J::array();for(double w:carriers){std::ostringstream h;h<<std::hexfloat<<w;out["v4_stock_carriers"].push_back(h.str());}
  out["v4_stock_two_limb_export"]=[&]{J r;try{(void)v4stock.require_two_limb_values();r["exported"]=true;}catch(const std::exception& x){r["exported"]=false;r["message"]=x.what();}return r;}();
  // Control: the melt proposal's own wet fixture (water count 4, two-limb) is admitted.
  auto control=view(e);out["control_fixture_water_count_4"]=attempt_stocks(e,control.stocks);
  // The audit Count2 input {1e15,-0.01} is itself two-limb and is admitted.
  auto audit_input=control.stocks;audit_input[2].count=Count(1e15,-.01);out["audit_count2_input_stock"]=attempt_stocks(e,audit_input);
  // The v4 accepted stock after the exact 1e-30 debit, substituted into the same fixture.
  auto substituted=control.stocks;substituted[2].count=v4stock;out["v4_accepted_stock_in_melt_fixture"]=attempt_stocks(e,substituted);
  // The melt proposal's own release of the same audit amount from the same audit input.
  J own;try{
   auto v=control;v.stocks[2].count=Count(1e15,-.01);v.thermal_energy_J=e.energy_at(v.stocks,400.);
   auto r=request(e,v,{action(s,"audit-c1-release","water",std::strtod("1.0000000000000024e-30",nullptr))});
   auto p=e.evaluate(v,r);own["admitted"]=true;own["stock_after_carriers"]=J::array();for(double w:p.stocks()[2].count.positive_carriers()){std::ostringstream h;h<<std::hexfloat<<w;own["stock_after_carriers"].push_back(h.str());}
  }catch(const std::exception& x){own["admitted"]=false;own["message"]=x.what();}
  out["melt_own_release_of_audit_amount"]=own;
  // A multi-word stock that IS two-limb representable is admitted (the refusal is about scale, not ExactParticleCount itself).
  auto two_limb=control.stocks;two_limb[2].count=Count(4.).debit(.5);out["two_limb_debited_stock"]=attempt_stocks(e,two_limb);
  const bool conflict=!out["v4_stock_two_limb_export"]["exported"].get<bool>()&&!out["v4_accepted_stock_in_melt_fixture"]["admitted"].get<bool>()&&out["control_fixture_water_count_4"]["admitted"].get<bool>()&&out["audit_count2_input_stock"]["admitted"].get<bool>();
  out["conflict_confirmed"]=conflict;
  std::ofstream f(std::filesystem::path(argv[2])/"melt-conflict.json");f<<out.dump(1)<<"\n";
  std::cout<<"conflict_confirmed="<<(conflict?"true":"false")<<"\n";return 0;
 }catch(const std::exception& x){std::cerr<<"harness error: "<<x.what()<<"\n";return 1;}
}
