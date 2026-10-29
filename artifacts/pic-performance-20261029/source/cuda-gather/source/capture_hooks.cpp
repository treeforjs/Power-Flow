#include "capture_hooks.hpp"
#include "snapshot_io.hpp"
#include <numeric>
#include <iomanip>

namespace pic_capture {
using namespace raikou::mfem_backend;
namespace io=snapshot_io;
namespace gp=gather_push;
namespace {
struct State {
    bool configured=false,active=false,done=false;
    std::filesystem::path directory;
    const mfem::GridFunction* psi=nullptr;
    double origin=0,midpoint=0,dt=0;
    const void* provider=nullptr;
    NativeWriter native=nullptr;
    std::size_t next=0;
    std::vector<GapParticle> before;
    std::vector<gp::Outputs> reference;
};
State state;
}
void configure(const std::filesystem::path& path,const mfem::GridFunction& psi,
               double origin,const void* provider,NativeWriter native) {
    io::require(!state.configured && !path.empty() && !std::filesystem::exists(path),"Capture requires a fresh output path");
    io::require(std::isfinite(origin) && provider && native,"Invalid capture configuration");
    state.configured=true;state.directory=path;state.psi=&psi;state.origin=origin;
    state.provider=provider;state.native=native;
}
void begin(const GapPic& pic,double dt) {
    auto& s=state;
    if(!s.configured || s.done || pic.particles().empty())return;
    io::require(!s.active,"Nested midpoint capture");
    s.before=pic.particles();s.reference.clear();s.reference.reserve(s.before.size());s.next=0;
    const auto diagnostics=pic.diagnostics();
    s.dt=dt;s.midpoint=diagnostics.time_s+.5*dt;
    std::filesystem::create_directories(s.directory);
    const auto& o=pic.options();std::vector<int> elements(o.nx*o.ny);std::iota(elements.begin(),elements.end(),0);
    const auto cells=io::export_grid(pic.potential(),elements),psi_cells=io::export_grid(*s.psi,elements);
    io::require(cells.size()==psi_cells.size() && std::memcmp(cells.data(),psi_cells.data(),cells.size()*sizeof(gp::Cell))==0,"Unit/self snapshot topology differs");
    io::write(s.directory/"self_cells.bin",cells);
    io::write(s.directory/"phi.bin",io::values(pic.potential()));
    io::write(s.directory/"psi.bin",io::values(*s.psi));
    io::write(s.directory/"particles.bin",s.before);
    io::write(s.directory/"diagnostics_before.bin",std::vector<GapPicDiagnostics>{diagnostics});
    io::write(s.directory/"options.bin",std::vector<GapPicOptions>{o});
    std::vector<gp::Inputs> inputs;inputs.reserve(s.before.size());
    for(const auto& p:s.before)inputs.push_back({p.position_m[0],p.position_m[1],p.position_m[2],io::vec(p.proper_velocity_m_s),p.charge_C/p.mass_kg});
    io::write(s.directory/"inputs.bin",inputs);
    s.active=true;
}
void record(const GapParticle& p,const PicVector& E,const PicVector& B) {
    auto& s=state;if(!s.active)return;
    io::require(s.next<s.before.size(),"Capture reference overflow");
    const auto& original=s.before[s.next++];
    io::require(p.id==original.id && p.position_m==original.position_m && p.charge_C==original.charge_C &&
                p.mass_kg==original.mass_kg && p.weight==original.weight,"Force loop changed identity/position/weight");
    s.reference.push_back({io::vec(E),io::vec(B),io::vec(p.proper_velocity_m_s),0});
}
void finish(const GapPic& pic,double dt) {
    auto& s=state;if(!s.active)return;
    io::require(s.next==s.before.size() && dt==s.dt && pic.particles().size()==s.before.size(),"Incomplete midpoint capture");
    io::write(s.directory/"reference.bin",s.reference);
    io::write(s.directory/"particles_after_kick.bin",pic.particles());
    const auto& o=pic.options();const auto d=pic.diagnostics();
    std::ofstream manifest(s.directory/"manifest.json.tmp");manifest<<std::setprecision(17)
        <<"{\n\"schema\":1,\"particles\":"<<s.before.size()<<",\"self_nx\":"<<o.nx<<",\"self_ny\":"<<o.ny
        <<",\"width\":"<<o.width_m<<",\"height\":"<<o.gap_m<<",\"depth\":"<<o.depth_m<<",\"short_z\":"<<o.short_z_m
        <<",\"dt\":"<<dt<<",\"ledger_time_s\":"<<d.time_s<<",\"absolute_native_time_s\":"<<s.origin+s.midpoint
        <<",\"real_midpoint_capture\":true,\"synthetic_self_fields\":false,\"capture_stage\":\"post_half_drift_post_first_guarded_solve_pre_kick\""
        <<",\"field_callbacks_added\":0,\"reference\":\"E/B and post-u copied from original force loop; no re-evaluation\""
        <<",\"potential_upper_boundary_V\":"<<o.upper_potential_V
        <<",\"poisson_residual_C_per_m\":"<<d.poisson_residual_C_per_m
        <<",\"maximum_nodal_continuity_defect_C\":"<<d.maximum_nodal_continuity_defect_C
        <<",\"particle_charge_C\":"<<d.particle_charge_C<<",\"deposited_charge_C\":"<<d.deposited_charge_C
        <<",\"not_a_restart_checkpoint\":true,\"GapParticle_bytes\":"<<sizeof(GapParticle)
        <<",\"GapPicOptions_bytes\":"<<sizeof(GapPicOptions)<<",\"GapPicDiagnostics_bytes\":"<<sizeof(GapPicDiagnostics);
    // The writer only reads Fields/Voltage state established by the original
    // callbacks. No set(), at(), electric(), magnetic() or extra solve occurs.
    s.native(s.provider,s.directory,s.origin+s.midpoint,manifest);
    manifest<<",\"gates\":{\"E_atol\":1e-6,\"B_atol\":1e-13,\"u_atol\":1e-5,\"field_rtol\":5e-12,\"u_rtol\":1e-11,\"weighted_relative_limit\":1e-12}\n}\n";
    manifest.close();io::require(bool(manifest),"Midpoint manifest write failed");
    std::filesystem::rename(s.directory/"manifest.json.tmp",s.directory/"manifest.json");
    s.active=false;s.done=true;
}
}
