#include "../snapshot/source/native_parent_subflow.hpp"
#include <fstream>
#include <iostream>
namespace n=raikou::native_water;
static_assert(!std::is_copy_constructible_v<n::ParentSubflow>);
static_assert(!std::is_copy_constructible_v<n::PreparedParentStage>);
static_assert(std::is_nothrow_move_constructible_v<n::PreparedParentStage>);
static_assert(!std::is_copy_constructible_v<n::PreparedParentState>);
static_assert(std::is_nothrow_move_constructible_v<n::PreparedParentState>);
int main(int argc,char**argv){
  if(argc!=2)return 2;
  n::m::GapPicOptions o;o.width_m=o.gap_m=.009;o.depth_m=.025;
  o.finite_axial_extent=true;o.short_z_m=0;o.feed_z_m=.025;
  n::Binding b{std::string(64,'a'),std::string(64,'b'),std::string(64,'c')};n::Limits l;
  const auto initial=n::ParentSubflow::initial_full_time_rng_sha256(o,{45,45,5},b,l,17,29);
  std::size_t checks=0;
  auto require=[&](bool yes,const char* why){++checks;if(!yes)throw std::runtime_error(why);};
  require(n::pc::digest(initial),"initial hash format");
  require(initial==n::ParentSubflow::initial_full_time_rng_sha256(o,{45,45,5},b,l,17,29),"initial hash repeatability");
  for(std::uint64_t seed=18;seed<26;++seed)require(initial!=n::ParentSubflow::initial_full_time_rng_sha256(o,{45,45,5},b,l,seed,29),"clock seed unbound");
  for(std::uint64_t seed=30;seed<38;++seed)require(initial!=n::ParentSubflow::initial_full_time_rng_sha256(o,{45,45,5},b,l,17,seed),"product seed unbound");
  require(initial!=n::ParentSubflow::initial_full_time_rng_sha256(o,{45,45,6},b,l,17,29),"grid unbound");
  auto q=l;++q.reactions;require(initial!=n::ParentSubflow::initial_full_time_rng_sha256(o,{45,45,5},b,q,17,29),"reaction cap unbound");
  q=l;++q.carriers;require(initial!=n::ParentSubflow::initial_full_time_rng_sha256(o,{45,45,5},b,q,17,29),"carrier cap unbound");
  q=l;++q.cells;require(initial!=n::ParentSubflow::initial_full_time_rng_sha256(o,{45,45,5},b,q,17,29),"cell cap unbound");
  q=l;++q.pairs;require(initial!=n::ParentSubflow::initial_full_time_rng_sha256(o,{45,45,5},b,q,17,29),"pair cap unbound");
  q=l;++q.queries;require(initial!=n::ParentSubflow::initial_full_time_rng_sha256(o,{45,45,5},b,q,17,29),"query cap unbound");
  q=l;++q.candidates;require(initial!=n::ParentSubflow::initial_full_time_rng_sha256(o,{45,45,5},b,q,17,29),"candidate cap unbound");
  auto binding=b;binding.dataset_sha256[0]='d';require(initial!=n::ParentSubflow::initial_full_time_rng_sha256(o,{45,45,5},binding,l,17,29),"dataset unbound");
  binding=b;binding.kernel_sha256[0]='d';require(initial!=n::ParentSubflow::initial_full_time_rng_sha256(o,{45,45,5},binding,l,17,29),"kernel unbound");
  binding=b;binding.kernel_binding_file_sha256[0]='d';require(initial!=n::ParentSubflow::initial_full_time_rng_sha256(o,{45,45,5},binding,l,17,29),"binding file unbound");
  std::ofstream out(argv[1]);out<<n::Json{{"status","passed_declared_header_and_initial_hash_only"},{"runtime_checks",checks},{"compile_time_token_traits",5},{"initial_rng_sha256",initial},{"native_full_time_implementation_linked",false},{"native_transactions_executed",0}}.dump(2)<<'\n';
  std::cout<<checks<<" initial hash assertions and five token traits pass\n";
}
