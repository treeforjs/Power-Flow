#include "motion_receipt.hpp"
#include "../hall-finite-foil-em-20261029/axial_particles.hpp"
#include <iostream>

namespace h=hall_driver;namespace n=raikou::mfem_backend;
using Rows=std::vector<n::ArenaRateEntry>;
Json checks=Json::array();int forbidden_calls=0;
extern "C" void __wrap__ZN6raikou12mfem_backend6GapPic12solve_fieldsEv(void*){++forbidden_calls;throw std::runtime_error("native Poisson forbidden in motion qualification");}
extern "C" void __wrap__ZN6raikou12mfem_backend6GapPic7advanceEd(void*,double){++forbidden_calls;throw std::runtime_error("native D-K-D forbidden in motion qualification");}
void check(bool ok,const std::string&name,Json detail=Json::object()){
    checks.push_back({{"name",name},{"passed",ok},{"detail",detail}});if(!ok)throw std::runtime_error(name);}
template<class F>void reject(F f,const std::string&name){bool caught=false;try{f();}catch(const std::exception&){caught=true;}check(caught,name);}
n::GapPicOptions options(const axial::Mesh&g){n::GapPicOptions o;o.nx=o.ny=g.nx;o.width_m=o.gap_m=g.side;o.depth_m=g.depth;
    o.finite_foils=o.finite_axial_extent=true;o.short_z_m=-g.depth/2;o.feed_z_m=g.depth/2;
    o.lower_foil={g.side/2-.0025,g.side/2+.0025,g.side/2-.0027,g.side/2-.0025};
    o.upper_foil={g.side/2-.0025,g.side/2+.0025,g.side/2+.0025,g.side/2+.0027};return o;}
Rows authored(){Rows rows;std::uint64_t id=1;
    // Fresh manufactured full-time initial data using the actual v6.2 positive
    // Count2 carrier factory. No accepted legacy phase/checkpoint is relabeled.
    for(int species=0;species<3;++species){n::ArenaCount count{1e5,std::ldexp(1.,-40)};
        for(double weight:n::PicArena::positive_carriers(count)){n::GapParticle p;
            p.id=id++;p.position_m={.00421,.00433,.001};p.proper_velocity_m_s={1e5,-2e5,3e5};
            p.charge_C=species==0?-pic_e:species==1?pic_e:0;p.mass_kg=species==0?pic_me:species==1?hall::proton_mass:2.9915e-26;p.weight=weight;
            rows.push_back({"manufactured/packet/"+std::to_string(species),"manufactured/state/"+std::to_string(species),p});}}
    return rows;
}
h::MotionBinding binding(const axial::Mesh&g,const n::GapPicOptions&o,const Rows&rows,const axial::Vec&E,const axial::Vec&B,double t,double dt){
    h::MotionBinding b;b.phase=h::ParticlePhase::full_time_position_and_proper_velocity;b.begin_s=t;b.dt_s=dt;b.end_s=t+dt;
    b.staged_arena_sha256=h::canonical_census_sha256(rows);b.trial_field_sha256=h::total_field_sha256(E,B);b.mesh_sha256=h::motion_mesh_sha256(g,o);
    b.material_sha256=std::string(64,'a');b.species_sha256=std::string(64,'b');b.common_clock_sha256=std::string(64,'c');
    b.kinetic_policy_sha256=h::kinetic_policy_sha256(h::declared_kinetic_policy(o));return b;
}
void apply_endpoints(Rows&rows,const h::MotionReceipt&r){std::map<std::uint64_t,n::GapParticle>map;for(auto p:r.endpoints())map.emplace(p.id,p);
    for(auto&row:rows)row.native=map.at(row.native.id);}
void rotate(Rows&rows,double dt){const double a=3e12*dt,c=std::cos(a),s=std::sin(a);for(auto&row:rows){auto&u=row.native.proper_velocity_m_s;
    const auto old=u;u[0]=c*old[0]-s*old[1];u[1]=s*old[0]+c*old[1];}}
Rows composition(const axial::Mesh&g,const n::GapPicOptions&o,const axial::Vec&E,const axial::Vec&B,unsigned refinement){
    auto rows=authored();double t=7.4978805549199997e-10;
    for(unsigned big=0;big<8;++big)for(unsigned sub=0;sub<refinement;++sub){double dt=(big%2?1e-14:5e-15)/refinement;
        rotate(rows,dt/2);auto b=binding(g,o,rows,E,B,t,dt);auto p=h::prepare_canonical_motion(g,o,rows,b,{},E,B);apply_endpoints(rows,p);t=b.end_s;rotate(rows,dt/2);}
    return rows;
}
double udifference(const Rows&a,const Rows&b){long double sum=0;for(std::size_t i=0;i<a.size();++i)for(int d=0;d<3;++d){
    long double x=a[i].native.proper_velocity_m_s[d]-b[i].native.proper_velocity_m_s[d];sum+=x*x;}return std::sqrt(sum);}
int main(int argc,char**argv){MPI_Init(&argc,&argv);int rc=0;try{
    hall::require(argc==2,"output path required");std::filesystem::path out=argv[1];axial::Mesh g;auto o=options(g);auto rows=authored();
    {n::GapPic sentinel(o);reject([&]{sentinel.solve_fields();},"native Poisson trap is active");
        reject([&]{sentinel.advance(1e-14);},"native D-K-D trap is active");check(forbidden_calls==2,"both forbidden entrypoint probes intercepted");forbidden_calls=0;}
    axial::Vec E(g.ne),B(g.nb);n::PicVector e{1e4,-2e4,3e4},b{.4,.3,.2};
    for(int i=0;i<g.ne;++i)E[i]=e[g.edges[i].axis];
    for(int i=0;i<g.nb;++i)B[i]=b[i<g.nbx?0:i<g.nbx+g.nby?1:2];
    const double t=7.4978805549199997e-10,dt=1e-14;auto bind=binding(g,o,rows,E,B,t,dt);auto before=h::canonical_census_sha256(rows);
    auto r=h::prepare_canonical_motion(g,o,rows,bind,{},E,B);
    check(rows.size()==6&&r.endpoints().size()==6&&r.paths().size()==6,"all six positive Count2 carriers get one full-time proposal");
    check(h::canonical_census_sha256(rows)==before,"motion prepare leaves every input carrier unchanged");
    bool identity=true,momentum=true,neutral=true;std::map<std::string,n::GapParticle> packets;
    double max_relative_u=0;
    for(std::size_t i=0;i<rows.size();++i){const auto&p=rows[i].native;const auto&q=r.endpoints()[i];
        identity=identity&&p.id==q.id&&p.weight==q.weight&&p.mass_kg==q.mass_kg&&p.charge_C==q.charge_C;
        auto expected=n::pic_boris(p.proper_velocity_m_s,p.charge_C/p.mass_kg,e,b,dt);
        for(int d=0;d<3;++d)max_relative_u=std::max(max_relative_u,std::abs(expected[d]-q.proper_velocity_m_s[d])/3e5);
        auto [it,added]=packets.emplace(rows[i].packet_id,q);if(!added)identity=identity&&q.position_m==it->second.position_m&&q.proper_velocity_m_s==it->second.proper_velocity_m_s;
        if(p.charge_C==0){neutral=neutral&&p.proper_velocity_m_s==q.proper_velocity_m_s;
            for(int d=0;d<3;++d)neutral=neutral&&std::abs(q.position_m[d]-(p.position_m[d]+dt*hall::velocity(p.proper_velocity_m_s)[d]))<=1e-18;}}
    momentum=max_relative_u<=2e-14;
    check(identity,"native IDs species exact weights and Count2 co-kinematics preserved");
    check(momentum,"momentum equals unchanged native Boris in original field axes",{{"max_relative_u",max_relative_u}});
    check(neutral,"neutral q=0 shares full-time drift with unchanged momentum");
    check(r.budget().continuity_C<=1e-25,"canonical translated endpoints satisfy compatible charge continuity",{{"defect_C",r.budget().continuity_C}});
    check(std::abs(r.budget().work_defect_J)<=1e-20,"same compatible current closes particle-field work",{{"defect_J",r.budget().work_defect_J}});
    auto reordered=rows;std::reverse(reordered.begin(),reordered.end());check(h::canonical_census_sha256(reordered)==before,"census hash independent of caller ordering");
    auto rr=h::prepare_canonical_motion(g,o,reordered,bind,{},E,B);check(rr.seal_sha256()==r.seal_sha256(),"receipt seal and endpoint order deterministic");
    auto changed=rows;changed.front().native.weight=std::nextafter(changed.front().native.weight,INFINITY);
    check(h::canonical_census_sha256(changed)!=before,"census seal detects one-ulp carrier weight change");
    auto poison=o;poison.lower_potential_V=1e9;poison.upper_potential_V=-1e9;poison.magnetic_field_T={1e9,1e9,1e9};
    auto rp=h::prepare_canonical_motion(g,poison,rows,bind,{},E,B);check(rp.seal_sha256()==r.seal_sha256(),"native prescribed E/B excluded from sole total-field force");
    auto invalid=bind;invalid.phase=h::ParticlePhase::unbound;reject([&]{h::prepare_canonical_motion(g,o,rows,invalid,{},E,B);},"unbound phase cannot masquerade as full-time");
    invalid=bind;invalid.dt_s=2e-14;invalid.end_s=invalid.begin_s+invalid.dt_s;reject([&]{h::prepare_canonical_motion(g,o,rows,invalid,{},E,B);},"10 fs transport limit unchanged");
    invalid=bind;invalid.end_s=std::nextafter(invalid.end_s,INFINITY);reject([&]{h::prepare_canonical_motion(g,o,rows,invalid,{},E,B);},"unbound endpoint clock rejected");
    invalid=bind;invalid.common_clock_sha256.clear();reject([&]{h::prepare_canonical_motion(g,o,rows,invalid,{},E,B);},"missing common clock seal rejected");
    auto wrongE=E;wrongE.front()+=1;reject([&]{h::prepare_canonical_motion(g,o,rows,bind,{},wrongE,B);},"field changed after binding rejected");
    auto wrongGeometry=o;wrongGeometry.lower_foil.left_m+=1e-6;reject([&]{h::prepare_canonical_motion(g,wrongGeometry,rows,bind,{},E,B);},"foil/translation mismatch rejected");
    std::vector<std::uint64_t> active{rows[0].native.id};reject([&]{h::prepare_canonical_motion(g,o,rows,bind,active,E,B);},"active layer cannot be silently frozen or double-pushed");
    changed=rows;changed.back().native.id=changed.front().native.id;reject([&]{h::prepare_canonical_motion(g,o,changed,bind,{},E,B);},"duplicate carrier ID rejected");
    changed=rows;changed[1].native.position_m[0]+=1e-6;reject([&]{h::prepare_canonical_motion(g,o,changed,bind,{},E,B);},"divergent Count2 carrier trajectories rejected");
    changed=rows;for(auto&row:changed){row.native.position_m[1]=o.upper_foil.bottom_m-1e-12;row.native.proper_velocity_m_s={0,1e5,0};}
    reject([&]{h::prepare_canonical_motion(g,o,changed,bind,{},E,B);},"foil hit requires native source/boundary owner");
    changed=rows;for(auto&row:changed){row.native.position_m[2]=o.feed_z_m-1e-12;row.native.proper_velocity_m_s={0,0,1e5};}
    reject([&]{h::prepare_canonical_motion(g,o,changed,bind,{},E,B);},"native feed-plane escape rejected inside larger field box");
    changed=rows;for(auto&row:changed)row.native.weight=1e50;
    reject([&]{h::prepare_canonical_motion(g,o,changed,bind,{},E,B);},"full census plasma-frequency guard retained");
    auto highB=B;for(auto&v:highB)v*=1e3;invalid=binding(g,o,rows,E,highB,t,dt);
    reject([&]{h::prepare_canonical_motion(g,o,rows,invalid,{},E,highB);},"full-interval gyro guard retained");
    check(h::canonical_census_sha256(rows)==before,"all rejected motion attempts preserve original input");
    const auto policy=h::declared_kinetic_policy(o);auto neutral_particle=rows.back().native;
    auto stream_case=[&](const std::string&name,n::PicVector x,n::PicVector u,h::KineticSink kind,int axis,double plane){
        auto p=neutral_particle;p.position_m=x;p.proper_velocity_m_s=u;
        auto s=h::prepare_neutral_stream(policy,p,t,t+dt,dt);
        double expected=(plane-x[axis])/hall::velocity(u)[axis];
        check(s.hit.sink==kind&&s.endpoint.position_m[axis]==plane&&std::abs(s.flight_s-expected)<=8*std::numeric_limits<double>::epsilon()*expected&&
            s.flight_s>0&&s.flight_s<dt&&s.remaining_s==dt-s.flight_s,"neutral first-hit "+name);
        check(s.endpoint.id==p.id&&s.endpoint.weight==p.weight&&s.endpoint.mass_kg==p.mass_kg&&s.endpoint.charge_C==0&&
            s.endpoint.proper_velocity_m_s==p.proper_velocity_m_s&&s.incident_kinetic_J==p.weight*n::pic_kinetic_energy_J(p),"neutral sink identity and kinetic operand "+name);
    };
    stream_case("feed escape",{.0045,.0045,o.feed_z_m-1e-10},{0,0,1e5},h::KineticSink::axial_escape,2,o.feed_z_m);
    stream_case("short bridge",{.0045,.0045,o.short_z_m+1e-10},{0,0,-1e5},h::KineticSink::short_bridge,2,o.short_z_m);
    stream_case("nonbridge short-end escape",{.001,.0045,o.short_z_m+1e-10},{0,0,-1e5},h::KineticSink::axial_escape,2,o.short_z_m);
    stream_case("outer x escape",{o.width_m-1e-10,.0045,0},{1e5,0,0},h::KineticSink::outer_xy_escape,0,o.width_m);
    stream_case("upper foil",{.0045,o.upper_foil.bottom_m-1e-10,0},{0,1e5,0},h::KineticSink::upper_foil,1,o.upper_foil.bottom_m);
    stream_case("lower foil",{.0045,o.lower_foil.top_m+1e-10,0},{0,-1e5,0},h::KineticSink::lower_foil,1,o.lower_foil.top_m);
    auto free=h::prepare_neutral_stream(policy,neutral_particle,t,t+dt,dt);
    check(free.hit.sink==h::KineticSink::none&&free.flight_s==dt&&free.remaining_s==0&&free.end_s==t+dt,"neutral full-time duration uses authoritative dt at C19 epoch");
    n::PicVector corner{o.width_m-1e-4,.0045,o.feed_z_m-1e-4},D{o.width_m-corner[0],0,o.feed_z_m-corner[2]};
    auto tie=h::first_kinetic_hit(policy,corner,D);check(tie.fraction==1&&tie.sink==h::KineticSink::axial_escape&&tie.axis==2,"exact endpoint hit and axial tie order match native kinetic policy");
    reject([&]{h::prepare_neutral_stream(policy,rows.front().native,t,t+dt,dt);},"neutral helper refuses charged carrier");
    invalid=bind;invalid.kinetic_policy_sha256=std::string(64,'d');reject([&]{h::prepare_canonical_motion(g,o,rows,invalid,{},E,B);},"independent kinetic policy seal mismatch rejected");
    axial::Mesh wide(85,4,.017,3);axial::Vec Ew(wide.ne),Bw(wide.nb);
    for(int i=0;i<wide.ne;++i)Ew[i]=e[wide.edges[i].axis];for(int i=0;i<wide.nb;++i)Bw[i]=b[i<wide.nbx?0:i<wide.nbx+wide.nby?1:2];
    auto wb=binding(wide,o,rows,Ew,Bw,t,dt);auto wr=h::prepare_canonical_motion(wide,o,rows,wb,{},Ew,Bw);
    bool same=true;for(std::size_t i=0;i<rows.size();++i)for(int d=0;d<3;++d)same=same&&
        std::abs(wr.endpoints()[i].position_m[d]-r.endpoints()[i].position_m[d])<=1e-18&&
        std::abs(wr.endpoints()[i].proper_velocity_m_s[d]-r.endpoints()[i].proper_velocity_m_s[d])<=1e-9;
    check(same&&wb.mesh_sha256!=bind.mesh_sha256&&wb.kinetic_policy_sha256==bind.kinetic_policy_sha256,
        "wider field xy and axial padding preserve native policy and no-hit trajectory");
    check(wr.frame().native_to_field_m==r.frame().native_to_field_m&&
        wr.paths()[0].field_start_m[0]==rows[0].native.position_m[0]-o.width_m/2&&
        wr.paths()[0].field_start_m[1]==rows[0].native.position_m[1]-o.gap_m/2,
        "kinetic-centered xy transform applied exactly once in wider field mesh");
    for(bool is_neutral:{false,true})for(bool axial:{false,true}){
        Rows cross;for(const auto&row:rows)if((row.native.charge_C==0)==is_neutral&&(is_neutral||row.native.charge_C<0))cross.push_back(row);
        for(auto&row:cross){row.native.position_m={.0045,.0045,0};row.native.proper_velocity_m_s={0,0,0};
            row.native.position_m[axial?2:0]=(axial?o.feed_z_m:o.width_m)-1e-10;row.native.proper_velocity_m_s[axial?2:0]=1e5;}
        auto cb=binding(wide,o,cross,Ew,Bw,t,dt);
        reject([&]{h::prepare_canonical_motion(wide,o,cross,cb,{},Ew,Bw);},std::string(is_neutral?"neutral":"charged")+
            (axial?" feed":" xy")+" first hit rejected before endpoint admission despite wider field box");
    }
    auto c1=composition(g,o,E,B,1),c2=composition(g,o,E,B,2),c4=composition(g,o,E,B,4);
    double d1=udifference(c1,c2),d2=udifference(c2,c4),ratio=d1/d2;
    check(ratio>3.8&&ratio<4.2,"variable-step manufactured C/2 transport C/2 composition converges at second order",
        {{"successive_u_difference_m_s",{d1,d2}},{"ratio",ratio},{"actual_SSA",false}});
    check(forbidden_calls==0,"native Poisson and D-K-D entrypoints never called (link-time traps)");
    std::ofstream(out/"motion-receipt.json")<<h::motion_record(r).dump(2)<<'\n';
    std::ofstream(out/"motion-validation.json")<<Json{{"status","passed"},{"count",checks.size()},{"checks",checks},
        {"accepted_native_arena_commit",false},{"actual_SSA",false},{"fixture","fresh authored full-time Count2 carriers"}}.dump(2)<<'\n';
    std::cout<<checks.size()<<" motion checks passed\n";
    }catch(const std::exception&e){std::cerr<<e.what()<<'\n';rc=1;}MPI_Finalize();return rc;}
