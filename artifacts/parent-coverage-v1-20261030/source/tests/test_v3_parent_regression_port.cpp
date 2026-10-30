#include "parent_coverage_subflow.hpp"
#include <iostream>
#include <fstream>
namespace nw=raikou::native_water;
namespace m=raikou::mfem_backend;
namespace pc=raikou::pic_collisions;
namespace w=raikou::pic_water;namespace cv=raikou::native_water::coverage_v1;
using J=nw::Json;
J checks=J::array(),receipts=J::array();
void check(bool ok,const std::string& name){checks.push_back({{"name",name},{"passed",ok}});if(!ok)throw std::runtime_error(name);}
template<class F>void reject(F f,const std::string& s){bool threw=false;try{f();}catch(const std::exception&){threw=true;}check(threw,s);}
constexpr double mn=2.9915e-26;
m::GapPicOptions options(){
  m::GapPicOptions o;o.nx=o.ny=45;o.width_m=o.gap_m=.009;o.depth_m=.025;
  o.finite_foils=true;o.finite_axial_extent=true;o.short_z_m=0;o.feed_z_m=.025;
  o.lower_foil={.002,.007,.0018,.002};o.upper_foil={.002,.007,.007,.0072};return o;
}
m::PicVector electron_u(double E){double x=E*pc::qe/(w::me*pc::c*pc::c);return {pc::c*std::sqrt(x*(x+2)),0,0};}
std::vector<m::ArenaStateSpec> states(){return {
  {"electron",0,0,-1,w::me,0,0},{"H2O",2,1,0,mn,0,0},
  {"H2O+",2,1,1,mn-w::me+w::threshold_eV*pc::qe/(pc::c*pc::c),0,w::threshold_eV*pc::qe}};}
std::vector<m::ArenaChannel> channels(const nw::Binding& b){return {{w::parent_channel,"electron","H2O",{"electron","electron","H2O+"},b.dataset_sha256,b.kernel_sha256,b.kernel_binding_file_sha256,w::parent_model}};}
struct Case {
  m::GapPic pic{options()};std::unique_ptr<m::PicArena> arena;
  Case(const nw::Binding& b,double ew=1,double nwgt=16,double E=100,bool separate=false){
    m::GapParticle e,n;e.weight=ew;e.mass_kg=w::me;e.position_m={.004,.004,.0125};e.proper_velocity_m_s=electron_u(E);
    n.weight=nwgt;n.mass_kg=mn;n.charge_C=0;n.position_m={separate?.005:.00403,.004,.0125};
    auto ids=pic.add_authored_batch({e,n});pic.solve_fields();
    arena=std::make_unique<m::PicArena>(pic,states(),channels(b),std::vector<m::ArenaPacket>{
      {"e","electron","authored_fixture","",{ew,0},{ids[0]}},
      {"n","H2O","authored_fixture","",{nwgt,0},{ids[1]}}},"native-water-subflow-fixture/v1");
  }
};
nw::Interval interval(m::PicArena& arena,const cv::ParentSubflow& flow,double duration,std::uint64_t sequence=0){
  nw::Interval i;i.sequence=sequence;i.end_s=duration;i.common_clock_sha256=nw::sha("frozen-native-parent-fixture/"+std::to_string(sequence));
  i.phase.state_prefix=arena.state_prefix();i.phase.mesh_prefix=flow.grid_identity();i.phase.clock_prefix=i.common_clock_sha256;return i;
}
long double weight(const m::PicArena& a,const std::string& state){
  long double total=0;for(const auto& [_,p]:a.snapshot()->packets)if(p.state==state)for(double part:m::PicArena::positive_carriers(p.count))total+=static_cast<long double>(part);return total;
}
int main(int argc,char**argv){
  try{
    if(argc!=5)throw std::runtime_error("dataset binding output-directory native-source-identity");
    std::ifstream f(argv[2]);J binding;f>>binding;
    nw::Binding b{binding.at("dataset_sha256"),binding.at("kernel_sha256"),binding.at("kernel_binding_file_sha256")};
    const std::string dir=argv[3];
    auto data=std::make_shared<const w::Data>(w::Data::read(argv[1]));w::ParentBridge bridge(data,w::parent_model,b.dataset_sha256,b.kernel_sha256,mn);const cv::CoverageScope scope(data,b.dataset_sha256,bridge);
    nw::Limits limits;
    // Exact geometric volume and axial grouping independently determine rate.
    nw::Grid grid(options(),{45,45,5},limits.cells);
    check(grid.locate({.004,.004,.0125})==grid.locate({.00403,.004,.0125}),"physical nearby pair same3D cell");
    check(grid.locate({.004,.004,.0125})!=grid.locate({.004,.004,.018}),"axially separated pairs are distinct");
    check(std::abs(grid.cells[0].volume_m3-2e-10)<1e-25,"physical3D volume includes depth once");
    check(grid.cells.size()==(45*45-2*25)*5,"fully solid foil cells excluded from collision volume");
    reject([&]{nw::Grid cut(options(),{44,45,5},limits.cells);},"unresolved partial solid collision cells rejected");
    reject([&]{grid.locate({.004,.0019,.0125});},"resident inside a solid collision cell rejected");
    reject([&]{grid.locate({.004,.002,.0125});},"material boundary needs an admitted vacuum-side transport stage");
    reject([&]{nw::Grid bad(options(),{UINT32_MAX,UINT32_MAX,UINT32_MAX},100);},"grid allocation cap before allocation");
    {
      Case c(b);cv::ParentSubflow flow(*c.arena,grid.counts,b,limits,bridge,scope,123,456);
      auto cp=flow.checkpoint();auto base=c.arena->state_prefix();auto native=m::PicArena::fingerprint(c.pic);
      auto request=interval(*c.arena,flow,1e8);auto trial=flow.prepare(request);
      check(flow.checkpoint()==cp&&c.arena->state_prefix()==base&&m::PicArena::fingerprint(c.pic)==native,"multi-event prepare leaves all accepted owners unchanged");
      check(trial.audit().at("events").size()>1,"finite-rate rebuilt populations produce multiple native events");
      check(flow.can_publish(trial)&&flow.try_publish(trial),"single common native+RNG publication");
      const auto completion=c.arena->stochastic_completion(cv::ParentSubflow::integration_owner());
      check(completion&&completion->common_clock_sha256==request.common_clock_sha256,
            "accepted reaction exposure retains exact common clock identity");
      check(!flow.try_publish(trial),"duplicate common publication refused");
      check(c.pic.diagnostics().time_s==0,"SSA waiting times never advance physical trajectory time");
      const auto ions=weight(*c.arena,"H2O+");
      check(ions>=2&&weight(*c.arena,"H2O")+ions==16&&weight(*c.arena,"electron")==1+ions,"native H2O depletion electron multiplication andcharge balance");
      double last=0;for(const auto& event:trial.audit().at("events")){
        double t=event.at("collision_elapsed_s");check(t>last&&t<1e8,"SSA eventclock strict order");last=t;
        check(event.at("physical_stage_s")==0,"each event bound to unchangedphysical stage");
      }
      receipts.push_back(trial.audit());
      reject([&]{flow.prepare(request);},"same interval cannot be applied twice");
      c.arena->save(dir+"/native");auto saved=flow.checkpoint();std::ofstream(dir+"/subflow.json")<<saved.dump(2);
      m::GapPic restored{options()};auto a=m::PicArena::restore(restored,dir+"/native",states(),channels(b),"native-water-subflow-fixture/v1");
      cv::ParentSubflow after(*a,grid.counts,b,limits,bridge,scope,9,10);after.restore(saved);
      check(after.checkpoint()==saved,"native andboth RNG streams restartexactly");
      auto bad=saved;bad["grid"]=nw::sha("foreign");reject([&]{after.restore(bad);},"foreign collisiongeometry restart refused");
      auto left=flow.prepare(interval(*c.arena,flow,1e8,1));auto right=after.prepare(interval(*a,after,1e8,1));
      check(left.audit()==right.audit(),"next native subflow receipt byte-equivalent afterrestart");
      check(flow.try_publish(left)&&after.try_publish(right)&&flow.checkpoint()==after.checkpoint(),"next accepted native+RNG state matches restart");
    }
    {
      Case c(b);auto small=limits;small.reactions=1;
      cv::ParentSubflow flow(*c.arena,grid.counts,b,small,bridge,scope,123,456);auto cp=flow.checkpoint();auto base=c.arena->state_prefix();
      reject([&]{flow.prepare(interval(*c.arena,flow,1e8));},"reaction cap aborts whole subflow ratherthan truncating");
      check(cp==flow.checkpoint()&&base==c.arena->state_prefix(),"late reactioncap restores nativeand RNGstate");
    }
    for(double E:{13.,8000.}){
      Case c(b,1,16,E);cv::ParentSubflow flow(*c.arena,grid.counts,b,limits,bridge,scope,1,2);auto old=flow.checkpoint();
      reject([&]{flow.prepare(interval(*c.arena,flow,1e8));},"unsupported incident energy cannot be silentlyignored");
      check(flow.checkpoint()==old,"support error leaves all owners unchanged");
    }
    {
      // The first event is supported; its new electrons land in the unresolved
      // onset gap. Rebuilding rates must reject the whole trial, not hide it.
      Case c(b,1,16,38.3);cv::ParentSubflow flow(*c.arena,grid.counts,b,limits,bridge,scope,123,456);
      auto old=flow.checkpoint();auto native=m::PicArena::fingerprint(c.pic);
      reject([&]{flow.prepare(interval(*c.arena,flow,1e8));},"post-reaction unsupported energy aborts whole subflow");
      check(flow.checkpoint()==old&&m::PicArena::fingerprint(c.pic)==native,
            "late support failure restores native population fields and both RNG streams");
    }
    {
      Case c(b,.25,1e20);cv::ParentSubflow flow(*c.arena,grid.counts,b,limits,bridge,scope,123,456);
      auto trial=flow.prepare(interval(*c.arena,flow,1e-5));
      const auto events=trial.audit().at("events").size();
      check(events>=2,"finite-rate compensated reservoir has multiple real events");
      check(flow.try_publish(trial),"compensated native reservoir subflow commits");
      const auto& count=c.arena->snapshot()->packets.at("n").count;
      const auto limbs=count.require_two_limb_values();
      check(limbs[0]==1e20&&limbs[1]==-.25*events,
            "sub-ULP neutral depletion retained in canonical two-limb inventory");
      check(c.arena->rate_carriers("n").size()>1,"subsequent propensities use positive exact native carriers");
      receipts.push_back(trial.audit());
    }
    {
      Case c(b,1,16,100,true);cv::ParentSubflow flow(*c.arena,grid.counts,b,limits,bridge,scope,1,2);
      auto old=flow.checkpoint();auto before=c.arena->state_prefix();
      auto trial=flow.prepare(interval(*c.arena,flow,1e8));
      check(trial.audit().at("events").empty()&&trial.audit().at("pairs")==0,"disjoint cells have no candidatepairs");
      check(flow.try_publish(trial),"empty collision subflow acceptedwithout changing particles");
      check(before!=c.arena->state_prefix(),"zero-event exposure changes native accepted history");
      auto latest=flow.checkpoint();
      reject([&]{flow.restore(old);},"unmodified mixed-age no-event checkpoint rejected");
      auto relabeled=old;relabeled["arena"]=c.arena->state_prefix();
      reject([&]{flow.restore(relabeled);},"old RNG clock cannot be relabeled to newer no-event arena");
      check(flow.checkpoint()==latest,"failed mixed-age restores leave all accepted state unchanged");
      cv::ParentSubflow fresh(*c.arena,grid.counts,b,limits,bridge,scope,1,2);
      reject([&]{fresh.checkpoint();},"new un-restored scheduler cannot masquerade as accepted arena owner");
      fresh.restore(latest);check(fresh.checkpoint()==latest,"matching no-event native/RNG checkpoint restores");
      c.arena->save(dir+"/native-no-event");
      m::GapPic recovered{options()};auto a=m::PicArena::restore(recovered,dir+"/native-no-event",states(),channels(b),"native-water-subflow-fixture/v1");
      cv::ParentSubflow from_disk(*a,grid.counts,b,limits,bridge,scope,999,999);from_disk.restore(latest);
      check(from_disk.checkpoint()==latest,"zero-event exposure and RNG completion persist through native disk restart");
      auto wrong=b;wrong.kernel_binding_file_sha256=nw::sha("wrong binding");
      reject([&]{cv::ParentSubflow bad(*a,grid.counts,wrong,limits,bridge,scope,1,2);},"channel binding enforced even when no collision is possible");
      auto exhausted=latest;exhausted["next_interval"]=std::uint64_t(UINT64_MAX);
      reject([&]{from_disk.restore(exhausted);},"exhausted interval counter refused during restore");
    }
    // One equal-weight neutral packet can react at most once: exact Bernoulli
    // probability 1-exp(-lambda*tau). No post-event chemistry approximation.
    int occurred=0;constexpr int samples=128;
    pc::Packet e{"e","electron",0,1,w::me,electron_u(100)},n{"n","H2O",0,1,mn,{}};
    const double rate=pc::candidates({e},{n},{{0,grid.cells[0].volume_m3}},{bridge.channel()})[0].rate_s;
    for(int seed=0;seed<samples;++seed){
      Case c(b,1,1);cv::ParentSubflow flow(*c.arena,grid.counts,b,limits,bridge,scope,seed,42);
      auto trial=flow.prepare(interval(*c.arena,flow,std::log(2.)/rate));
      const auto count=trial.audit().at("events").size();check(count<=1,"single target packet cannot react twice");
      occurred+=int(count);check(flow.try_publish(trial),"actual finite-rate one-target trial commits");
      check(weight(*c.arena,"H2O")==1-count&&weight(*c.arena,"electron")==1+count,"trial native species inventory matches sampled event");
    }
    check(occurred>=44&&occurred<=84,"128 seeded native first-event samples agree with analytic half-probability (3.54sigma)");
    J result{{"status","passed"},{"checks",checks},{"check_count",checks.size()},
      {"native_source_identity",argv[4]},{"first_event_sample_count",samples},{"first_events",occurred},
      {"expected_probability",.5},{"receipts",receipts},{"conditional_parent_only",true},
      {"physical_trajectory_advanced",false},{"full_coupled_run",false},{"coverage_scope_sha256",scope.identity()},{"coverage_scope",scope.description()}};
    std::ofstream(dir+"/result.json")<<result.dump(2)<<'\n';std::cout<<"Passed "<<checks.size()<<" checks; "<<occurred<<"/"<<samples<<" first events\n";
  }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}
