#pragma once
// Hall corrective v3.5: ResidentOwner::restore_as_genesis -- genesis at a vacuum pre-run checkpoint (seal checked under
// its recorded configuration, heat-only source-patch declarations added, carrier limits and native geometry transferred
// for a vacuum state, revision rebased to 0 with the clock kept, provenance bound into the state and heat digests).
// Hall corrective v3.2: the v3.1 device-resident owner with the admitted Gauss ledger in the device
// state (bootstrap with a zero ledger) and checkpoint save/restore of the resident owner
// (ResidentOwner::save / ResidentOwner::restore; restore also into a Model whose configuration
// differs only in the field-step upper bound and/or the solve backend record, verified by
// recomputing the recorded configuration identity) and/or the motion forms (field state unaffected).
// Hall corrective v3.1: device-resident owner (ResidentOwner) with the protocol of the v3
// Owner (prepare_vacuum / preflight_vacuum_final / publish_vacuum_prevalidated, the closed
// coupled proposal prepare_closed_proposal / preflight_final / publish_prevalidated).
//
// The accepted and the candidate field states live on the device (gpu::ResidentField). The
// host keeps the thermal owner (ThermalV3, unchanged), the accepted free charge (host mirror
// of the motion result, for the census check), the kinetic ledgers, the digests and the
// slot metadata (time, revision, declarations). The physical arithmetic of a step equals
// the Hall v3 Model::prepare arithmetic bit for bit for the same increment; the solver is
// the gpu1 CG on the device right-hand side (or the sealed CPU PCG in the parity mode); the
// independent FP64 residual check runs on the host with the sealed arithmetic.
//
// Identities: the field state digest is the block-Merkle record hall-v3.1-field-state-merkle-v1
// and the midpoint role seal hall-v3.1-authoritative-midpoint-merkle-v1 (different values
// from the v3 streaming records; one-to-one with the same words). Heat, kinetic-ledger,
// conductivity and every motion/native-facing digest are the v3 records unchanged.
#include "root_api.hpp"
#include "gpu/resident_field.hpp"
#include "gpu/device_step_v34.hpp"
#include <ctime>
#include <optional>

namespace hall_corrective_v3 {
class ResidentOwner;
class ResidentJointProposal;
ResidentJointProposal prepare_closed_proposal_device(ResidentOwner&, const native::GapPicOptions&, std::span<const native::ArenaRateEntry>,
                                                     hall_driver::MotionBinding, const std::string&, const StepInputs&);
// Motion context of the resident owner: the binary v3 record digest only for the accepted iterate.
inline MotionContext deferred_v3(MotionContext c) { c.defer_record_v3 = true; return c; }

// Field-solve backend record of the resident step: bound into the Model configuration
// identity (name + descriptor digest). It is not callable through Model::prepare.
class ResidentBackend final : public FieldSolveBackend {
    std::string name_ = "gpu1-resident-v1", descriptor_;
public:
    explicit ResidentBackend(const std::string& field_descriptor_sha256) {
        Digest d("hall-v2-field-solve-backend-v1");
        d.text(name_); d.text(gpu::ResidentField::descriptor_text()); d.text(field_descriptor_sha256);
        descriptor_ = d.hex();
    }
    const std::string& name() const noexcept override { return name_; }
    const std::string& descriptor_sha256() const noexcept override { return descriptor_; }
    Vec solve(const Apply&, const Vec&, int, const SolveContext&, SolveRecord&) override {
        throw std::runtime_error("the resident field backend is driven by ResidentOwner, not by Model::prepare");
    }
    const char* residual_rule() const noexcept override { return gpu::ResidentField::residual_rule(); }   // v3.4 exact norms
};
// Descriptor digest of the resident kernels as compiled for this device (compiles once).
inline std::string resident_descriptor_sha256(std::shared_ptr<gpu::Device> device, const Model& probe_model) {
    gpu::ResidentField f(std::move(device), probe_model, gpu::ResidentOptions{});
    return f.descriptor_sha256();
}

class ResidentFieldHeatProposal {
    friend class ResidentOwner; friend class ResidentJointProposal;
    friend ResidentJointProposal prepare_closed_proposal(ResidentOwner&, const native::GapPicOptions&, std::span<const native::ArenaRateEntry>,
                                                         hall_driver::MotionBinding, const std::string&, MotionProvider&, const StepInputs&);
    friend ResidentJointProposal prepare_closed_proposal_device(ResidentOwner&, const native::GapPicOptions&, std::span<const native::ArenaRateEntry>,
                                                                hall_driver::MotionBinding, const std::string&, const StepInputs&);
    gpu::ResidentStep step_; gpu::SlotMeta meta_; ThermalTrialV3 heat_;
    std::weak_ptr<const OwnerLifetime> owner_; std::uint64_t generation_ = 0, token_ = 0;
    std::string before_, after_, configuration_, conductivity_, midpoint_seal_, before_heat_, after_heat_;
    FreeCharge after_free_; Vec current_; double dt_ = 0;
    bool device_free_ = false;   // v3.4: after_free_ not held on the host (device path)
    mutable bool final_approved_ = false; bool live_ = true;
    KineticLedgers ledgers_; std::string ledgers_sha256_;
    ResidentFieldHeatProposal() = default;
    void check_energy() const {
        const auto& r = step_.receipt;
        double transfers = std::abs(r.field_delta_J) + std::abs(r.plasma_work_J) + std::abs(heat_.thermal_delta_J) + std::abs(r.outer_loss_J) +
                           std::abs(r.outer_supply_J) + std::abs(r.circuit_work_J);
        if (heat_.loads) transfers += std::abs(heat_.wall_heat_J) + heat_.external_heating_J + heat_.external_cooling_J;
        need(std::isfinite(field_heat_balance_J()) && std::abs(field_heat_balance_J()) <= 2e-20 + 1e-10 * transfers + heat_.roundoff_bound_J,
             "unchanged field/heat component energy bound");
    }
public:
    ResidentFieldHeatProposal(const ResidentFieldHeatProposal&) = delete;
    ResidentFieldHeatProposal& operator=(const ResidentFieldHeatProposal&) = delete;
    ResidentFieldHeatProposal(ResidentFieldHeatProposal&&) noexcept = default;
    ResidentFieldHeatProposal& operator=(ResidentFieldHeatProposal&&) noexcept = default;
    const std::string& candidate_field_sha256() const noexcept { return after_; }
    const std::string& candidate_heat_sha256() const noexcept { return after_heat_; }
    double end_s() const noexcept { return meta_.time_s; }
    std::uint64_t endpoint_revision() const noexcept { return meta_.revision; }
    const std::string& midpoint_seal() const noexcept { return midpoint_seal_; }
    const std::string& configuration_sha256() const noexcept { return configuration_; }
    const axial::Receipt& field_receipt() const noexcept { return step_.receipt; }
    const Constraints& candidate_constraints() const noexcept { return step_.constraints; }
    const SolveRecord& solve_record() const noexcept { return step_.solve; }
    const gpu::ControlReceipt& control() const noexcept { return step_.control; }
    const gpu::ResidentStep& step() const noexcept { return step_; }
    const ThermalTrialV3& heat() const noexcept { return heat_; }
    const KineticLedgers& candidate_ledgers() const noexcept { return ledgers_; }
    const std::string& candidate_ledgers_sha256() const noexcept { return ledgers_sha256_; }
    const FreeCharge& candidate_free_charge() const noexcept { return after_free_; }
    double dt_s() const noexcept { return dt_; }
    double field_heat_balance_J() const noexcept {
        if (!heat_.loads) return step_.receipt.balance_J + heat_.thermal_delta_J - step_.receipt.joule_J;
        return step_.receipt.balance_J + heat_.thermal_delta_J - step_.receipt.joule_J - heat_.wall_heat_J - heat_.external_J;
    }
};

class ResidentJointProposal {
    friend class ResidentOwner;
    friend ResidentJointProposal prepare_closed_proposal(ResidentOwner&, const native::GapPicOptions&, std::span<const native::ArenaRateEntry>,
                                                         hall_driver::MotionBinding, const std::string&, MotionProvider&, const StepInputs&);
    friend ResidentJointProposal prepare_closed_proposal_device(ResidentOwner&, const native::GapPicOptions&, std::span<const native::ArenaRateEntry>,
                                                                hall_driver::MotionBinding, const std::string&, const StepInputs&);
    ResidentFieldHeatProposal field_heat_; hall_driver::MotionReceipt motion_; std::string source_;
    Vec midE_, midB_;   // host copy of the accepted fixed-point midpoint (equal to the device midpoint); v3.4 device path:
                        // filled when the owner keeps host record arrays (native v7 host admission), else empty
    bool device_ = false;          // v3.4: built by the device path (record arrays and census on the device)
    std::uint64_t seal_generation_ = 0;   // v3.4: DeviceStep::record_view().seal_generation of this proposal's seal
    std::string trial_digest_;     // v3.4: total-field digest of the accepted trial midpoint (computed once)
    std::vector<std::string> census_states_;   // v3.2: census v2 state table of the motion (native v7 forms)
    double combined_balance_ = 0; unsigned iterations_ = 0; ProposalTimings timings_; StepEnergy energy_;
    ResidentJointProposal(ResidentFieldHeatProposal&& p, hall_driver::MotionReceipt&& motion, std::string source, unsigned iterations, ProposalTimings timings,
                          Vec midE, Vec midB)
        : field_heat_(std::move(p)), motion_(std::move(motion)), source_(std::move(source)), midE_(std::move(midE)), midB_(std::move(midB)),
          iterations_(iterations), timings_(std::move(timings)) {
        // Same checks as JointProposal (root_api.hpp).
        const auto& r = field_heat_.step_.receipt; const auto& h = field_heat_.heat_;
        const bool births = !motion_.births().empty();
        const double work = motion_.budget().kinetic_change_J - r.plasma_work_J;
        // v3.2: the work tolerance of the motion (motion_work_rule), as JointProposal
        const double work_tolerance = motion_work_tolerance(motion_.work_scales());
        need(std::isfinite(work) && std::abs(work) <= work_tolerance, "unchanged particle/Maxwell work bound");
        combined_balance_ = field_heat_.field_heat_balance_J() + work;
        double transfers = std::abs(r.field_delta_J) + std::abs(motion_.budget().kinetic_change_J) + std::abs(h.thermal_delta_J) + std::abs(r.outer_loss_J) +
                           std::abs(r.outer_supply_J) + std::abs(r.circuit_work_J);
        if (births || h.loads) transfers += std::abs(h.wall_heat_J) + h.external_heating_J + h.external_cooling_J;
        need(std::isfinite(combined_balance_) && std::abs(combined_balance_) <= 2e-20 + 1e-10 * transfers + h.roundoff_bound_J + work_tolerance,
             "unchanged complete field/motion/heat energy bound");
        energy_ = StepEnergy::of(r, motion_, h, combined_balance_, transfers);
        const auto sinks = SinkLedger::of(motion_.sinks());
        const double foil = (sinks.boundary[0].kinetic_J + sinks.boundary[1].kinetic_J).rounded();
        need(std::abs(foil - h.wall_heat_J) <= 1e-12 * std::abs(foil) + 1e-300, "collected foil kinetic energy differs from the applied wall heat");
    }
public:
    ResidentJointProposal(const ResidentJointProposal&) = delete;
    ResidentJointProposal& operator=(const ResidentJointProposal&) = delete;
    ResidentJointProposal(ResidentJointProposal&&) noexcept = default;
    ResidentJointProposal& operator=(ResidentJointProposal&&) noexcept = default;
    const std::string& candidate_field_sha256() const noexcept { return field_heat_.candidate_field_sha256(); }
    const std::string& candidate_heat_sha256() const noexcept { return field_heat_.candidate_heat_sha256(); }
    double end_s() const noexcept { return field_heat_.end_s(); }
    std::uint64_t endpoint_revision() const noexcept { return field_heat_.endpoint_revision(); }
    const hall_driver::MotionReceipt& motion() const noexcept { return motion_; }
    double combined_balance_J() const noexcept { return combined_balance_; }
    const std::string& configuration_sha256() const noexcept { return field_heat_.configuration_sha256(); }
    unsigned iterations() const noexcept { return iterations_; }
    const ResidentFieldHeatProposal& field_heat() const noexcept { return field_heat_; }
    const ProposalTimings& timings() const noexcept { return timings_; }
    const StepEnergy& energy() const noexcept { return energy_; }
    const Vec& midpoint_E() const noexcept { return midE_; }
    const Vec& midpoint_B() const noexcept { return midB_; }
    // v3.2 native v7 admission inputs (API-007 amendment A1, HALL-LINK 3e): the trial midpoint E and B of the accepted
    // motion iterate (total_field_sha256 = binding().trial_field_sha256) and its plasma current.
    const Vec& trial_E() const noexcept { return midE_; }
    const Vec& trial_B() const noexcept { return midB_; }
    const Vec& current_A_m() const noexcept { return motion_.current_A_m(); }
    bool device_path() const noexcept { return device_; }
    // v3.4: the seal generation of the device step's record view that belongs to this proposal (native v7.1 lane)
    std::uint64_t device_seal_generation() const noexcept { return seal_generation_; }
    // v3.4: the accessors of JointProposal the coordinator uses (Owner-compatible names): the accepted trial midpoint
    // (host copies; present with host record arrays), its role seal (block-Merkle midpoint seal) and total-field digest.
    struct MidpointView {
        const ResidentJointProposal& p;
        const Vec& E() const noexcept { return p.midE_; }
        const Vec& B() const noexcept { return p.midB_; }
        const std::string& seal() const noexcept { return p.field_heat_.midpoint_seal(); }
        const std::string& total_field_sha256() const noexcept { return p.trial_digest_.empty() ? p.motion_.binding().trial_field_sha256 : p.trial_digest_; }
    };
    MidpointView midpoint() const noexcept { return MidpointView{*this}; }
    std::span<const SinkRecord> sinks() const noexcept { return motion_.sinks(); }
    std::span<const NewbornRecord> births() const noexcept { return motion_.births(); }
};

class ResidentOwner {
    friend ResidentJointProposal prepare_closed_proposal(ResidentOwner&, const native::GapPicOptions&, std::span<const native::ArenaRateEntry>,
                                                         hall_driver::MotionBinding, const std::string&, MotionProvider&, const StepInputs&);
    friend ResidentJointProposal prepare_closed_proposal_device(ResidentOwner&, const native::GapPicOptions&, std::span<const native::ArenaRateEntry>,
                                                                hall_driver::MotionBinding, const std::string&, const StepInputs&);
    std::shared_ptr<const Model> model_; std::shared_ptr<gpu::Device> device_; const native::GapPicOptions geometry_; const MotionContext motion_;
    std::unique_ptr<gpu::ResidentField> dev_;
    gpu::SlotMeta meta_; FreeCharge free_;
    std::unique_ptr<ThermalV3> thermal_; KineticLedgers ledgers_; std::string ledgers_sha256_;
    const std::shared_ptr<const OwnerLifetime> lifetime_ = std::make_shared<const OwnerLifetime>();
    std::uint64_t generation_ = 0, token_ = 0; bool vacuum_ = false;
    std::string field_sha256_, heat_sha256_;
    std::string heat_identity_;                 // v3.5: material binding digest, bound to the genesis provenance when rebased
    GenesisProvenance provenance_;              // v3.5: empty for a bootstrap or a restore of a non-genesis checkpoint
    std::unique_ptr<gpu::DeviceStep> step_;   // v3.4 device path (created at the first coupled proposal in device mode)
    mutable bool free_stale_ = false;          // v3.4: the accepted free charge lives on the device (host mirror refreshed on demand)
    mutable axial::State field_cache_; mutable ChargeLedger charge_cache_; mutable std::uint64_t field_cache_gen_ = 0;   // v3.4 field()/charge()
    void refresh_free() const {
        if (!free_stale_) return;
        axial::State st; ChargeLedger q; dev_->download(0, st, q, meta_);
        const_cast<FreeCharge&>(free_) = std::move(q.particles); free_stale_ = false;
    }
    static std::shared_ptr<const Model> checked_model(std::shared_ptr<const Model> m) {
        need(bool(m), "missing immutable model");
        need(m->solver().name() == "gpu1-resident-v1", "ResidentOwner requires a Model configured with the resident field backend");
        return m;
    }
    void clock_guard() const { need(meta_.time_s == thermal_->time_s() && meta_.revision == thermal_->revision(), "accepted field/heat clock mismatch"); }
    std::string conductivity_sha256() const { return Model::conductivity_sha256(thermal_->conductivity()); }
    void admit_census(std::span<const native::ArenaRateEntry> census) const {
        refresh_free();
        auto actual = charge_of_census(mesh(), geometry_, census, model_->limits().maximum_carriers);
        need(same_net_charge(free_, actual), "accepted charge does not match actual canonical census");
        (void)dev_->admit(0, meta_, &actual);
    }
    void interval_guard(double dt) const {
        need(std::isfinite(dt) && dt > 0 && dt <= model_->dt_max() && meta_.time_s + dt > meta_.time_s && meta_.time_s + dt <= model_->drive().end() &&
             meta_.revision < UINT64_MAX - 1, "authoritative field interval");
    }
    gpu::ResidentInputs inputs(double dt, const Vec& sigma) const {
        gpu::ResidentInputs in;
        in.dt = dt; in.sigma = &sigma;
        const auto& d = model_->drive();
        in.current = d.current(meta_.time_s + .5 * dt);
        const double i0 = d.current(meta_.time_s), i1 = d.current(meta_.time_s + dt);
        in.drive_rel_change = std::abs(i1 - i0) / std::max(std::abs(in.current), 1.);
        in.drive_knot_interval_s = d.knot_interval(meta_.time_s + .5 * dt);
        return in;
    }
    // Candidate identities and the thermal trial (v3 heat_for with the device candidate).
    ResidentFieldHeatProposal finish(gpu::ResidentStep&& step, double dt, const std::string& sigma_sha, std::span<const SurfaceLoad> loads,
                                     const hall_driver::MotionReceipt* motion, FreeCharge&& after, Vec&& current) {
        ResidentFieldHeatProposal p;
        p.meta_ = meta_; p.meta_.time_s = meta_.time_s + dt; ++p.meta_.revision;
        p.meta_.residue_present = step.residue_present;   // bank residue array present in the candidate ledger
        p.dt_ = dt; p.token_ = token_;
        p.before_ = field_sha256_; p.configuration_ = configuration_sha256(); p.conductivity_ = sigma_sha;
        p.step_ = std::move(step);
        p.after_ = dev_->state_digest(1, api, configuration_sha256(), p.meta_);
        p.midpoint_seal_ = dev_->midpoint_seal();
        need(p.before_ == field_sha256_ && p.configuration_ == configuration_sha256() && p.conductivity_ == conductivity_sha256(),
             "field preparation is stale or belongs to another configuration");
        p.heat_ = thermal_->prepare(p.step_.receipt, dt, p.after_, loads);
        KineticLedgers next = ledgers_;
        if (motion && !motion->v2_identical()) {
            const auto st = SinkLedger::of(motion->sinks()); next.sinks.add(st); next.births.add(BirthLedger::of(motion->births()));
            next.short_bridge_J = next.short_bridge_J + st.boundary[3].kinetic_J; next.escaped_J = next.escaped_J + st.escaped_kinetic_J();
        }
        if (p.heat_.loads) {
            next.loads += p.heat_.loads;
            for (const auto& [k, v] : p.heat_.channel_J) { auto& w = next.load_channel_J[k]; w = w + v; }
            next.wall_heat_J = next.wall_heat_J + Words(p.heat_.wall_heat_J); next.external_J = next.external_J + Words(p.heat_.external_J);
        }
        p.ledgers_ = std::move(next); p.ledgers_sha256_ = p.ledgers_.sha256();
        p.owner_ = lifetime_; p.generation_ = generation_; p.before_heat_ = heat_sha256_;
        p.after_heat_ = heat_sha256_of(heat_identity_, p.heat_);
        p.after_free_ = std::move(after); p.current_ = std::move(current);
        p.check_energy();
        return p;
    }
    bool quick_preflight(const ResidentFieldHeatProposal& p) const noexcept {
        return p.live_ && p.owner_.lock() == lifetime_ && p.generation_ == generation_ && p.token_ == token_ && thermal_->preflight(p.heat_);
    }
    void publish(ResidentFieldHeatProposal&& p) noexcept {
        if (!p.final_approved_ || !quick_preflight(p) || generation_ == UINT64_MAX) std::terminate();
        if (!thermal_->commit(std::move(p.heat_))) std::terminate();
        dev_->publish();
        using std::swap;
        swap(meta_, p.meta_); swap(free_, p.after_free_);
        free_stale_ = p.device_free_;   // v3.4 device path: the candidate free charge is the device slot (no host copy)
        swap(field_sha256_, p.after_); swap(heat_sha256_, p.after_heat_);
        swap(ledgers_, p.ledgers_); swap(ledgers_sha256_, p.ledgers_sha256_);
        ++generation_; ++token_; p.live_ = false; p.final_approved_ = false;
    }
    struct RestoreTag {};
    static ChargeLedger bootstrap_ledger(ChargeLedger q) {
        need(zero_words(q.gauss_admitted_C), "bootstrap admitted Gauss ledger must be zero");
        for (double x : q.divb_gross_T_m) need(x == 0, "bootstrap divB gross ledger must be zero");
        return q;
    }
    // Bootstrap (restored = false: exact census charge) or restore (restored = true: same net charge as Owner::restore).
    // v3.5 `sealed`: called after the upload and before the census check and the host/device admission (restore_as_genesis
    // checks the seals, the configuration transfer and the vacuum state there, in that order).
    ResidentOwner(std::shared_ptr<const Model> m, std::shared_ptr<gpu::Device> device, axial::State field, ChargeLedger charge,
                  const native::GapPicOptions& o, std::span<const native::ArenaRateEntry> census, gpu::ResidentOptions options, bool restored,
                  const std::function<void(ResidentOwner&, const axial::State&, const ChargeLedger&)>& sealed = {})
        : model_(checked_model(std::move(m))), device_(std::move(device)), geometry_(o),
          motion_(deferred_v3(with_forms(with_dt_max(make_motion_context(model_->mesh(), geometry_, model_->mesh_json(), model_->limits().maximum_carriers),
                                                       model_->dt_max()), model_->motion_forms()))),
          vacuum_(census.empty()) {
        heat_identity_ = model_->material_sha256();
        dev_ = std::make_unique<gpu::ResidentField>(device_, *model_, options);
        need(ResidentBackend(dev_->descriptor_sha256()).descriptor_sha256() == model_->solver().descriptor_sha256(),
             "resident kernels differ from the backend descriptor bound into the Model configuration");
        if (sealed) {
            meta_.time_s = field.time_s; meta_.revision = field.revision; meta_.conductor_declaration = charge.conductor_declaration;
            meta_.exterior_declaration = charge.exterior_declaration; meta_.residue_present = !charge.residue_C.empty();
            dev_->upload(0, field, charge);
            sealed(*this, field, charge);
        }
        const auto actual = charge_of_census(mesh(), geometry_, census, model_->limits().maximum_carriers);
        auto same_words = [](const WordVec& a, const WordVec& b) {
            if (a.size() != b.size()) return false;
            for (std::size_t i = 0; i < a.size(); ++i) if (a[i].hi != b[i].hi || a[i].lo != b[i].lo) return false;
            return true; };
        need(restored ? same_net_charge(charge.particles, actual) : (same_words(actual.net_C, charge.particles.net_C) && actual.gross_C == charge.particles.gross_C),
             "initial/restored net charge differs from actual canonical carriers");
        model_->admit(field, charge);   // host admission of the bootstrap state (sealed arithmetic)
        meta_.time_s = field.time_s; meta_.revision = field.revision; meta_.conductor_declaration = charge.conductor_declaration;
        meta_.exterior_declaration = charge.exterior_declaration; meta_.residue_present = !charge.residue_C.empty();
        free_ = charge.particles;
        if (!sealed) dev_->upload(0, field, charge);
        (void)dev_->admit(0, meta_);
        admit_census(census);
        if (!thermal_) thermal_ = std::make_unique<ThermalV3>(mesh(), model_->drive().start(), model_->thermal_declaration(), model_->dt_max());
        field_sha256_ = dev_->state_digest(0, api, configuration_sha256(), meta_);
        heat_sha256_ = heat_sha256_of(heat_identity_, *thermal_);
        ledgers_sha256_ = ledgers_.sha256();
        if (!restored) clock_guard();
    }
public:
    // Bootstrap from an accepted host state (same checks as Owner): the census charge must
    // equal the declared free charge and the admitted Gauss ledger must be zero; the state is
    // uploaded and admitted on the device.
    ResidentOwner(std::shared_ptr<const Model> m, std::shared_ptr<gpu::Device> device, axial::State field, ChargeLedger charge,
                  const native::GapPicOptions& o, std::span<const native::ArenaRateEntry> census, gpu::ResidentOptions options = {})
        : ResidentOwner(std::move(m), std::move(device), std::move(field), bootstrap_ledger(std::move(charge)), o, census, options, false) {}
    ResidentOwner(const ResidentOwner&) = delete;
    ResidentOwner& operator=(const ResidentOwner&) = delete;
    const Mesh& mesh() const noexcept { return model_->mesh(); }
    const Model& model() const noexcept { return *model_; }
    double time_s() const noexcept { return meta_.time_s; }
    std::uint64_t revision() const noexcept { return meta_.revision; }
    const std::string& field_sha256() const noexcept { return field_sha256_; }
    const std::string& heat_sha256() const noexcept { return heat_sha256_; }
    const std::string& configuration_sha256() const noexcept { return model_->configuration_sha256(); }
    const std::string& motion_mesh_sha256() const noexcept { return motion_.mesh_sha256; }
    const MotionContext& motion_context() const noexcept { return motion_; }
    const ThermalV3& thermal() const noexcept { return *thermal_; }
    const KineticLedgers& kinetic_ledgers() const noexcept { return ledgers_; }
    const std::string& kinetic_ledgers_sha256() const noexcept { return ledgers_sha256_; }
    const FreeCharge& free_charge() const { refresh_free(); return free_; }
    // v3.4: Owner-compatible accessors of the accepted state (host copies downloaded once per accepted generation;
    // O(mesh), for consumers that read the field on the host: source frames, controller, frames, probes).
    const axial::State& field() const {
        if (field_cache_gen_ != generation_ + 1) { dev_->download(0, field_cache_, charge_cache_, meta_); field_cache_gen_ = generation_ + 1; }
        return field_cache_;
    }
    const ChargeLedger& charge() const { (void)field(); return charge_cache_; }
    // v3.4 host-work mode and the device path engine.
    gpu::HostWork host_work() const noexcept { return dev_->options().host_work; }
    bool device_mode() const noexcept { return dev_->options().host_work == gpu::HostWork::device && model_->motion_forms().native_v7; }
    gpu::DeviceStep& device_step() {
        if (!step_) step_ = std::make_unique<gpu::DeviceStep>(device_, *model_, motion_, std::max<std::size_t>(model_->limits().maximum_carriers, 1));
        return *step_;
    }
    bool has_device_step() const noexcept { return bool(step_); }
    // Trial-field digest of host E/B arrays in the record's form (ResidentOptions::field_digest_form).
    std::string host_field_digest(const Vec& E, const Vec& B) const {
        if (dev_->options().field_digest_form == gpu::FieldDigestForm::blocks_v71) {
#if HALL_FORMS_V71
            return pf_v7::trial_field_blocks_sha256(E, B);
#else
            throw std::runtime_error("block-form field digests need the forms library of native v7.1");
#endif
        }
        return hall_driver::total_field_sha256(E, B);
    }
    const gpu::SlotMeta& accepted_meta() const noexcept { return meta_; }
    gpu::ResidentField& device_field() noexcept { return *dev_; }
    // Accepted state (download; checkpoints, comparisons).
    void accepted_state(axial::State& s, ChargeLedger& q) const { dev_->download(0, s, q, meta_); }
    void candidate_state(axial::State& s, ChargeLedger& q, const ResidentFieldHeatProposal& p) const {
        need(p.token_ == token_ && p.live_, "candidate slot no longer holds this proposal");
        dev_->download(1, s, q, p.meta_);
    }

    ResidentFieldHeatProposal prepare_vacuum(double dt, std::span<const SurfaceLoad> loads = {}) {
        need(vacuum_, "field-only publication requires an empty actual census");
        clock_guard();
        interval_guard(dt);
        const auto sigma = thermal_->conductivity();
        const auto sigma_sha = Model::conductivity_sha256(sigma);
        auto in = inputs(dt, sigma);
        ++token_;
        auto step = dev_->prepare(in, meta_);
        refresh_free();
        return finish(std::move(step), dt, sigma_sha, loads, nullptr, FreeCharge(free_), Vec());
    }
    bool preflight_field_heat(const ResidentFieldHeatProposal& p) const {
        try {
            if (!quick_preflight(p) || p.before_ != field_sha256_ || p.before_heat_ != heat_sha256_ || p.configuration_ != configuration_sha256() ||
                p.conductivity_ != conductivity_sha256())
                return false;
            clock_guard();
            // device recomputation: both admissions, every leaf of both slots (equal to the digested
            // leaves), a fresh host sample of leaves, and the midpoint from both slots
            if (!dev_->preflight_recheck(meta_, p.meta_)) return false;
            if (dev_->state_digest(0, api, configuration_sha256(), meta_) != field_sha256_ ||
                dev_->state_digest(1, api, configuration_sha256(), p.meta_) != p.after_ || dev_->midpoint_seal() != p.midpoint_seal_)
                return false;
            if (heat_sha256_of(heat_identity_, *thermal_) != heat_sha256_ || p.after_heat_ != heat_sha256_of(heat_identity_, p.heat_) ||
                p.end_s() != time_s() + p.dt_ || p.endpoint_revision() != revision() + 1 || p.ledgers_sha256_ != p.ledgers_.sha256())
                return false;
            for (const auto& body : p.heat_.bodies) if (body.state().time != p.end_s() || body.state().revision != p.endpoint_revision()) return false;
            return true;
        } catch (const std::exception&) { return false; }
    }
    template <class Boundary> bool preflight_vacuum_final(const ResidentFieldHeatProposal& p, const Boundary& b) const {
        p.final_approved_ = false;
        const bool ok = vacuum_ && preflight_field_heat(p) && b.field_sha256 == p.candidate_field_sha256() && b.heat_sha256 == p.candidate_heat_sha256() &&
                        b.time_s == p.end_s() && b.accepted_step == p.endpoint_revision();
        p.final_approved_ = ok;
        return ok;
    }
    void publish_vacuum_prevalidated(ResidentFieldHeatProposal&& p) noexcept { if (!vacuum_) std::terminate(); publish(std::move(p)); }

    bool preflight(const ResidentJointProposal& p, std::span<const native::ArenaRateEntry> census) const {
        if (p.device_) return preflight_device(p, census);
        try {
            if (!preflight_field_heat(p.field_heat_)) return false;
            const auto& m = p.motion(); const auto& b = m.binding();
            { std::vector<Newborn> nb; for (const auto& x : m.births()) nb.push_back({x.native_id, x.birth_time_s}); admit_census(residents(census, nb)); }
            // v3.2: native v7 forms: census v2 and the motion record v7 recomputed; else the v3 census/seal/record checks
            if (m.forms_v7() != motion_.forms_v7) return false;
            if (m.forms_v7()) {
                if (m.census_sha256() != census_v2(census, p.census_states_, motion_.maximum_carriers).sha256 || !digest(m.record_v7_sha256()) ||
                    MotionAssembly::record_v7(m) != m.record_v7_sha256())
                    return false;
            } else if (m.census_sha256() != canonical_census_sha256(census, motion_.maximum_carriers) || !digest(m.seal_sha256()) ||
                       !digest(m.record_v3_sha256()) || MotionAssembly::record_v3(m) != m.record_v3_sha256())
                return false;
            if (b.phase != hall_driver::ParticlePhase::full_time_position_and_proper_velocity ||
                b.accepted_step != revision() || b.begin_s != time_s() || b.end_s != p.end_s() || b.dt_s != p.field_heat_.dt_ ||
                b.mesh_sha256 != motion_.mesh_sha256 || b.material_sha256 != model_->material_sha256() || b.kinetic_policy_sha256 != motion_.policy_sha256 ||
                b.trial_field_sha256 != hall_driver::total_field_sha256(p.midE_, p.midB_) || m.current_A_m() != p.field_heat_.current_)
                return false;
            auto actual = charge_of_endpoints(mesh(), geometry_, m.endpoints(), motion_.maximum_carriers);
            if (!m.v2_identical() && ((!m.forms_v7() && MotionAssembly::sink_birth_digest(m).empty()) ||
                                      p.field_heat_.candidate_ledgers().sinks.count() != kinetic_ledgers().sinks.count() + m.sinks().size() ||
                                      p.field_heat_.candidate_ledgers().births.count() != kinetic_ledgers().births.count() + m.births().size()))
                return false;
            const auto& after = p.field_heat_.after_free_;
            if (actual.net_C.size() != after.net_C.size() || actual.gross_C != after.gross_C) return false;
            for (std::size_t i = 0; i < actual.net_C.size(); ++i) if (actual.net_C[i].hi != after.net_C[i].hi || actual.net_C[i].lo != after.net_C[i].lo) return false;
            return digest(p.source_) && std::isfinite(p.combined_balance_J());
        } catch (const std::exception&) { return false; }
    }
    // v3.4 device path: field/heat as v3.3 (device recheck, cached digests), census v2, census charge, record v7 section roots
    // and the endpoint free charge recomputed on the device (DeviceStep::preflight); O(1) binding checks; the trial
    // digest compared with the one computed for the binding and the stored midpoint compared bitwise with the trial buffers.
    bool preflight_device(const ResidentJointProposal& p, std::span<const native::ArenaRateEntry> census) const {
        try {
            if (!preflight_field_heat(p.field_heat_)) return false;
            const auto& m = p.motion(); const auto& b = m.binding();
            if (!m.forms_v7() || m.forms_v7() != motion_.forms_v7 || !step_) return false;
            if (!step_->preflight(m, census, *dev_, meta_)) return false;
            if (b.phase != hall_driver::ParticlePhase::full_time_position_and_proper_velocity ||
                b.accepted_step != revision() || b.begin_s != time_s() || b.end_s != p.end_s() || b.dt_s != p.field_heat_.dt_ ||
                b.mesh_sha256 != motion_.mesh_sha256 || b.material_sha256 != model_->material_sha256() || b.kinetic_policy_sha256 != motion_.policy_sha256 ||
                b.trial_field_sha256 != p.trial_digest_ || !dev_->midpoint_is_trial())
                return false;
            if (!p.midE_.empty() && b.trial_field_sha256 != host_field_digest(p.midE_, p.midB_)) return false;
            if (!m.v2_identical() && (p.field_heat_.candidate_ledgers().sinks.count() != kinetic_ledgers().sinks.count() + m.sinks().size() ||
                                      p.field_heat_.candidate_ledgers().births.count() != kinetic_ledgers().births.count() + m.births().size()))
                return false;
            return digest(p.source_) && std::isfinite(p.combined_balance_J());
        } catch (const std::exception&) { return false; }
    }
    template <class Boundary> bool preflight_final(const ResidentJointProposal& p, std::span<const native::ArenaRateEntry> census, const Boundary& b) const {
        p.field_heat_.final_approved_ = false;
        const auto& m = p.motion().binding();
        const bool ok = preflight(p, census) && b.field_sha256 == p.candidate_field_sha256() && b.heat_sha256 == p.candidate_heat_sha256() &&
                        b.source_sha256 == p.source_ && b.time_s == p.end_s() && b.accepted_step == p.endpoint_revision() && b.mesh_sha256 == m.mesh_sha256 &&
                        b.material_sha256 == m.material_sha256 && b.common_clock_sha256 == m.common_clock_sha256;
        p.field_heat_.final_approved_ = ok;
        return ok;
    }
    void publish_prevalidated(ResidentJointProposal&& p) noexcept { publish(std::move(p.field_heat_)); }
    // v3.2: checkpoint of the accepted state (field words, complete charge ledger incl. the admitted Gauss
    // ledger, heat slices, kinetic ledgers) with the resident identities; see resident_checkpoint below.
    void save(const std::filesystem::path& dir) const;
    static std::unique_ptr<ResidentOwner> restore(std::shared_ptr<const Model>, std::shared_ptr<gpu::Device>, const native::GapPicOptions&,
                                                  std::span<const native::ArenaRateEntry> census, const std::filesystem::path& dir,
                                                  gpu::ResidentOptions options = {});
    // v3.5: genesis at a vacuum pre-run checkpoint (v3.3, v3.4 or v3.5 resident checkpoint): seals checked under the
    // recorded configuration, configuration transfer (GenesisTransfer), vacuum state, admission in this Model; then the
    // revision is rebased to 0 with the clock kept and the provenance bound into the digests (README section 7).
    static std::unique_ptr<ResidentOwner> restore_as_genesis(std::shared_ptr<const Model> model, std::shared_ptr<gpu::Device> device,
                                                             const native::GapPicOptions& geometry, std::span<const native::ArenaRateEntry> census,
                                                             const std::filesystem::path& dir, gpu::ResidentOptions options = {},
                                                             GenesisTransfer transfer = {});
    const GenesisProvenance& genesis_provenance() const noexcept { return provenance_; }
};

// ------------------------------------------------------------------ resident checkpoint (v3.2; v3.5 schema)
// checkpoint.json: schema (api), checkpoint "hall-v3.5-resident-checkpoint-v1", configuration_sha256 and the
// operands that may differ on restore (dt_max_s, solver_name, solver_descriptor_sha256), motion_mesh_sha256,
// field (words, clock, revision), charge (v3 record), bank_residue_C, admitted_gauss_C, divb_gross_T_m (v3.4),
// field_sha256 (the resident block-Merkle state digest under configuration_sha256, with the genesis provenance of a
// rebased owner), heat_sha256 (+ heat/ slices), kinetic ledgers and their digest; v3.5: genesis_provenance (null or
// the record of GenesisProvenance), thermal_declaration (null or {sha256, patches}) and carrier_limits [cap, budget]
// (the operands a later restore_as_genesis recomputes the configuration from). Restore (v3.4 and v3.5 schemas): the
// recorded configuration must equal the Model's configuration, or Model::configuration_with(recorded dt_max, solver
// name, solver descriptor, motion forms) (same mesh, drive, material, rules, limits); the state is uploaded and
// admitted (host and device), its digest under the recorded configuration must equal field_sha256, and the heat and
// kinetic-ledger seals must match.
inline constexpr char resident_checkpoint_schema[] = "hall-v3.5-resident-checkpoint-v1";      // v3.5: + provenance, operands
inline constexpr char resident_checkpoint_schema_v34[] = "hall-v3.4-resident-checkpoint-v1";  // restore and genesis
inline constexpr char resident_checkpoint_schema_v33[] = "hall-v3.2-resident-checkpoint-v1";  // Hall v3.3 (genesis only)

inline void ResidentOwner::save(const std::filesystem::path& dir) const {
    clock_guard();
    axial::State s; ChargeLedger q; accepted_state(s, q);
    model_->admit(s, q);
    need(!std::filesystem::exists(dir), "checkpoint destination already exists");
    std::filesystem::create_directory(dir); thermal_->save(dir / "heat");
    Json j = Json::object();
    j["schema"] = api; j["checkpoint"] = resident_checkpoint_schema; j["configuration_sha256"] = configuration_sha256();
    j["dt_max_s"] = model_->dt_max(); j["solver_name"] = model_->solver().name(); j["solver_descriptor_sha256"] = model_->solver().descriptor_sha256();
    j["motion_forms_native_v7"] = model_->motion_forms().native_v7; j["neutral_stream_dt_max_s"] = model_->motion_forms().neutral_stream_dt_max_s;
    j["motion_mesh_sha256"] = motion_.mesh_sha256;
    j["field"] = state_record(s); j["charge"] = record(q); j["bank_residue_C"] = q.residue_C;
    j["admitted_gauss_C"] = checkpoint_detail::ledger_json(q.gauss_admitted_C);
    j["divb_gross_T_m"] = q.divb_gross_T_m;   // v3.4
    j["residue_present"] = meta_.residue_present;
    j["field_sha256"] = field_sha256_; j["heat_sha256"] = heat_sha256_;
    j["kinetic_ledgers"] = checkpoint_detail::ledgers_json(ledgers_); j["kinetic_ledgers_sha256"] = ledgers_sha256_;
    j["genesis_provenance"] = genesis_detail::provenance_json(provenance_);   // v3.5
    {   const auto& d = model_->thermal_declaration();
        if (d.empty()) j["thermal_declaration"] = nullptr;
        else { std::size_t n = 0; for (const auto& [k, sl] : d.slices) n += sl.patches.size();
               Json t = Json::object(); t["sha256"] = d.sha256(); t["patches"] = n; j["thermal_declaration"] = std::move(t); } }
    j["carrier_limits"] = Json::array({model_->limits().maximum_carriers, model_->limits().memory_budget_bytes});
    std::string text = j.dump(); drain(j);
    std::ofstream out(dir / "checkpoint.json", std::ios::binary); out << text << '\n'; out.close();
    need(bool(out), "checkpoint metadata write failed");
}
namespace genesis_detail {
// Parsed resident checkpoint (v3.3, v3.4 or v3.5 schema): integrity of form (schema, keys, byte cap, canonical words).
struct ParsedCheckpoint {
    std::string schema, recorded, solver_name, solver_descriptor, field_sha256, heat_sha256, ledgers_sha256, motion_mesh_sha256;
    double dt_max = 0; MotionForms forms; axial::State state; ChargeLedger ledger; KineticLedgers ledgers;
    GenesisProvenance provenance; std::optional<CarrierLimits> limits; bool v33 = false;
};
inline ParsedCheckpoint parse_checkpoint(const Mesh& g, const std::filesystem::path& dir, bool allow_v33) {
    using namespace checkpoint_detail;
    const auto file = dir / "checkpoint.json";
    const auto limit = std::uintmax_t(128) * (2 * std::uintmax_t(g.ne) + 2 * std::uintmax_t(g.nb) + 9 * std::uintmax_t(g.nn) + g.cells.size()) + 1048576;
    need(std::filesystem::is_regular_file(file) && std::filesystem::file_size(file) <= limit, "checkpoint byte cap");
    std::ifstream in(file, std::ios::binary); Json j = Json::parse(in);
    struct Guard { Json& j; ~Guard() { drain(j); } } guard{j};
    need(j.is_object() && j.contains("checkpoint") && j.at("checkpoint").is_string(), "checkpoint schema");
    ParsedCheckpoint c; c.schema = j.at("checkpoint").get<std::string>();
    c.v33 = c.schema == resident_checkpoint_schema_v33;
    if (c.schema == resident_checkpoint_schema)
        keys(j, {"schema", "checkpoint", "configuration_sha256", "dt_max_s", "solver_name", "solver_descriptor_sha256", "motion_forms_native_v7",
                 "neutral_stream_dt_max_s", "motion_mesh_sha256", "field", "charge", "bank_residue_C", "admitted_gauss_C", "divb_gross_T_m", "residue_present",
                 "field_sha256", "heat_sha256", "kinetic_ledgers", "kinetic_ledgers_sha256", "genesis_provenance", "thermal_declaration", "carrier_limits"});
    else if (c.schema == resident_checkpoint_schema_v34)
        keys(j, {"schema", "checkpoint", "configuration_sha256", "dt_max_s", "solver_name", "solver_descriptor_sha256", "motion_forms_native_v7",
                 "neutral_stream_dt_max_s", "motion_mesh_sha256", "field", "charge", "bank_residue_C", "admitted_gauss_C", "divb_gross_T_m", "residue_present",
                 "field_sha256", "heat_sha256", "kinetic_ledgers", "kinetic_ledgers_sha256"});
    else if (c.v33 && allow_v33)
        keys(j, {"schema", "checkpoint", "configuration_sha256", "dt_max_s", "solver_name", "solver_descriptor_sha256", "motion_forms_native_v7",
                 "neutral_stream_dt_max_s", "motion_mesh_sha256", "field", "charge", "bank_residue_C", "admitted_gauss_C", "residue_present",
                 "field_sha256", "heat_sha256", "kinetic_ledgers", "kinetic_ledgers_sha256"});
    else need(false, "checkpoint schema");
    need(j.at("schema") == api, "checkpoint schema");
    c.recorded = hash(j.at("configuration_sha256"));
    need(j.at("solver_name").is_string() && j.at("residue_present").is_boolean() && j.at("motion_forms_native_v7").is_boolean(), "checkpoint schema types");
    c.solver_name = j.at("solver_name").get<std::string>(); c.solver_descriptor = hash(j.at("solver_descriptor_sha256"));
    c.forms = MotionForms{j.at("motion_forms_native_v7").get<bool>(), number(j.at("neutral_stream_dt_max_s"))};
    c.dt_max = number(j.at("dt_max_s"));
    c.motion_mesh_sha256 = hash(j.at("motion_mesh_sha256"));
    c.state = field(j.at("field"), g); c.ledger = charge(j.at("charge"), g);
    {   const auto& rj = j.at("bank_residue_C"); need(rj.is_array() && (rj.empty() || rj.size() == std::size_t(g.nn)), "checkpoint bank residue dimensions");
        for (const auto& x : rj) { const double v = number(x); need(v >= 0, "checkpoint bank residue sign"); c.ledger.residue_C.push_back(v); } }
    need(j.at("residue_present").get<bool>() == !c.ledger.residue_C.empty(), "checkpoint bank residue presence");
    c.ledger.gauss_admitted_C = ledger_words(j.at("admitted_gauss_C"), g.nn);
    if (!c.v33) c.ledger.divb_gross_T_m = checkpoint_detail::divb_gross(j.at("divb_gross_T_m"), g.cells.size());   // v3.3: none (zero)
    c.field_sha256 = hash(j.at("field_sha256")); c.heat_sha256 = hash(j.at("heat_sha256"));
    c.ledgers = checkpoint_detail::ledgers_of(j.at("kinetic_ledgers")); c.ledgers_sha256 = hash(j.at("kinetic_ledgers_sha256"));
    if (c.schema == resident_checkpoint_schema) {
        c.provenance = provenance_of(j.at("genesis_provenance"));
        const auto& cl = j.at("carrier_limits");
        need(cl.is_array() && cl.size() == 2 && cl[0].is_number_unsigned() && cl[1].is_number_unsigned(), "checkpoint carrier limits");
        c.limits = CarrierLimits{cl[0].get<std::size_t>(), cl[1].get<std::uint64_t>()};
        const auto& td = j.at("thermal_declaration");
        need(td.is_null() || (td.is_object() && td.size() == 2 && td.contains("sha256") && td.contains("patches") && td.at("patches").is_number_unsigned()),
             "checkpoint thermal declaration record");
        if (!td.is_null()) (void)hash(td.at("sha256"));
    }
    return c;
}
} // namespace genesis_detail
inline std::unique_ptr<ResidentOwner> ResidentOwner::restore(std::shared_ptr<const Model> model, std::shared_ptr<gpu::Device> device,
                                                             const native::GapPicOptions& o, std::span<const native::ArenaRateEntry> census,
                                                             const std::filesystem::path& dir, gpu::ResidentOptions options) {
    model = checked_model(std::move(model));
    auto c = genesis_detail::parse_checkpoint(model->mesh(), dir, false);
    need(c.recorded == model->configuration_sha256() || c.recorded == model->configuration_with(c.dt_max, c.solver_name, c.solver_descriptor, c.forms),
         "checkpoint field/forcing/mesh configuration changed");
    need(c.motion_mesh_sha256 == model->motion_mesh_sha256(o), "checkpoint field/forcing/mesh configuration changed");
    auto candidate = std::unique_ptr<ResidentOwner>(new ResidentOwner(model, std::move(device), std::move(c.state), std::move(c.ledger), o, census, options, true));
    if (!c.provenance.empty()) {   // v3.5: a rebased owner keeps its provenance
        candidate->provenance_ = c.provenance; candidate->meta_.provenance = c.provenance.sha256();
        candidate->heat_identity_ = genesis_heat_identity(model->material_sha256(), candidate->meta_.provenance);
    }
    // seal under the recorded configuration (the uploaded words), then the identity under the restoring configuration
    need(candidate->dev_->state_digest(0, api, c.recorded, candidate->meta_) == c.field_sha256, "checkpoint field/charge seal mismatch");
    candidate->field_sha256_ = candidate->dev_->state_digest(0, api, candidate->configuration_sha256(), candidate->meta_);
    candidate->thermal_->load(dir / "heat"); candidate->clock_guard();
    candidate->heat_sha256_ = heat_sha256_of(candidate->heat_identity_, *candidate->thermal_);
    need(c.heat_sha256 == candidate->heat_sha256_, "checkpoint thermal seal mismatch");
    candidate->ledgers_ = std::move(c.ledgers); candidate->ledgers_sha256_ = candidate->ledgers_.sha256();
    need(c.ledgers_sha256 == candidate->ledgers_sha256_, "checkpoint kinetic ledger seal mismatch");
    return candidate;
}
// v3.5 genesis at a vacuum pre-run checkpoint. Order of the checks (any failure refuses; no owner is kept):
//  a. integrity: schema, keys, byte cap, canonical words (parse); the recorded configuration R identified as this
//     Model's composition with the checkpoint's operands (dt_max, solver record, motion forms, rule form, thermal
//     declaration, carrier limits; GenesisTransfer); the field seal (resident state digest of the uploaded words under
//     R, the v3.3 form for a v3.3 checkpoint), the thermal seal (slices loaded under R's declaration, digest with R's
//     material binding) and the kinetic-ledger seal;
//  b. transfer: declaration only adds heat-only patches, other carrier limits or native geometry only as allowed;
//  c. vacuum: zero free charge (net and gross) at every node, no bank residue, empty kinetic ledgers;
//  d. admission in this Model (census charge, host and device admission: absolute Gauss with the restored ledger, divB);
//  e. rebase: revision 0, clock kept; slice revisions 0 (temperatures, cumulative energies, watermarks, time kept;
//     zero inventory for added heat-only patches); kinetic ledgers empty; every field/charge word kept;
//  f. identity: state and heat digests under this Model's configuration, bound to GenesisProvenance::sha256().
inline std::unique_ptr<ResidentOwner> ResidentOwner::restore_as_genesis(std::shared_ptr<const Model> model, std::shared_ptr<gpu::Device> device,
                                                                        const native::GapPicOptions& o, std::span<const native::ArenaRateEntry> census,
                                                                        const std::filesystem::path& dir, gpu::ResidentOptions options,
                                                                        GenesisTransfer transfer) {
    using namespace genesis_detail;
    model = checked_model(std::move(model)); const auto& g = model->mesh();
    auto c = parse_checkpoint(g, dir, true);
    const ConfigurationForm form = c.v33 ? ConfigurationForm::v33 : ConfigurationForm::v34;
    // (a) the recorded configuration as a transfer of this Model
    std::vector<ThermalDeclaration> decls;
    if (transfer.recorded_declaration) decls.push_back(*transfer.recorded_declaration);
    else { decls.push_back(model->thermal_declaration()); if (!model->thermal_declaration().empty()) decls.push_back(ThermalDeclaration{}); }
    std::vector<CarrierLimits> lims{model->limits()};
    for (const auto& l : {transfer.recorded_limits, c.limits, std::optional<CarrierLimits>(CarrierLimits{})})
        if (l && std::none_of(lims.begin(), lims.end(), [&](const CarrierLimits& x) { return same_limits(x, *l); })) lims.push_back(*l);
    const ThermalDeclaration* dr = nullptr; const CarrierLimits* lr = nullptr;
    for (const auto& d : decls) {
        for (const auto& l : lims)
            if (model->configuration_transfer(c.dt_max, c.solver_name, c.solver_descriptor, c.forms, d, l, form) == c.recorded) { dr = &d; lr = &l; break; }
        if (dr) break;
    }
    need(dr != nullptr, "genesis checkpoint configuration is not a transfer of the restoring Model (mesh, drive, material law, seed, rules, guards or operands differ)");
    const std::string current_mesh = model->motion_mesh_sha256(o);
    std::optional<ThermalV3> recorded_heat;
    const std::string prior = c.provenance.empty() ? std::string() : c.provenance.sha256();
    auto sealed = [&](ResidentOwner& w, const axial::State&, const ChargeLedger& q) {
        // (a) seals under the recorded configuration
        w.meta_.provenance = prior;
        const std::string seal = c.v33 ? w.dev_->state_digest_v33(0, api, c.recorded, w.meta_) : w.dev_->state_digest(0, api, c.recorded, w.meta_);
        need(seal == c.field_sha256, "checkpoint field/charge seal mismatch");
        recorded_heat.emplace(g, model->drive().start(), *dr, c.dt_max);
        recorded_heat->load(dir / "heat");
        need(recorded_heat->time_s() == c.state.time_s && recorded_heat->revision() == c.state.revision, "accepted field/heat clock mismatch");
        need(heat_sha256_of(genesis_heat_identity(Model::material_sha256_of(*dr), prior), *recorded_heat) == c.heat_sha256, "checkpoint thermal seal mismatch");
        need(c.ledgers.sha256() == c.ledgers_sha256, "checkpoint kinetic ledger seal mismatch");
        // (b) transfer
        const bool same_declaration = dr->sha256() == model->thermal_declaration().sha256() && dr->empty() == model->thermal_declaration().empty();
        need(same_declaration || (transfer.thermal_declaration && declaration_extends(*dr, model->thermal_declaration())),
             "genesis thermal declaration may only add heat-only source patches");
        const bool vacuum = vacuum_charge(q) && !w.meta_.residue_present && c.ledgers.sha256() == KineticLedgers{}.sha256();
        need(same_limits(*lr, model->limits()) || (transfer.carrier_limits && vacuum), "genesis carrier limits differ");
        need(c.motion_mesh_sha256 == current_mesh || (transfer.native_geometry && vacuum), "genesis native geometry differs");
        // (c) vacuum state
        need(vacuum, "genesis checkpoint is not a vacuum state (free charge, bank residue or kinetic ledgers)");
    };
    // (d) admission in the restoring Model (census charge, host and device admission) by the restoring constructor
    auto owner = std::unique_ptr<ResidentOwner>(new ResidentOwner(model, std::move(device), c.state, c.ledger, o, census, options, true, sealed));
    // (e) rebase
    GenesisProvenance p;
    p.checkpoint_schema = c.schema; p.configuration_sha256 = c.recorded; p.field_sha256 = c.field_sha256; p.heat_sha256 = c.heat_sha256;
    p.kinetic_ledgers_sha256 = c.ledgers_sha256; p.motion_mesh_sha256 = c.motion_mesh_sha256; p.revision = c.state.revision; p.time_s = c.state.time_s;
    p.transfer_text = std::string("hall-v3.5-genesis-transfer-v1") +
        "; rule_form=" + (c.v33 ? "v3.3" : "v3.4") +
        "; dt_max_s=" + g17(c.dt_max) + "->" + g17(model->dt_max()) +
        "; solver=" + c.solver_name + "/" + c.solver_descriptor + "->" + model->solver().name() + "/" + model->solver().descriptor_sha256() +
        "; motion_forms=" + (c.forms.native_v7 ? "native-v7" : "native-v6") + "/" + g17(c.forms.neutral_stream_dt_max_s) + "->" +
            (model->motion_forms().native_v7 ? "native-v7" : "native-v6") + "/" + g17(model->motion_forms().neutral_stream_dt_max_s) +
        "; thermal_declaration=" + declaration_text(*dr) + "->" + declaration_text(model->thermal_declaration()) +
        "; carrier_limits=" + limits_text(*lr) + "->" + limits_text(model->limits()) +
        "; motion_mesh=" + c.motion_mesh_sha256 + "->" + current_mesh +
        "; divb_gross_ledger=" + (c.v33 ? "zero (v3.3 checkpoint)" : "kept") +
        "; prior_provenance=" + (prior.empty() ? std::string("none") : prior);
    owner->thermal_->install_genesis(*recorded_heat);
    owner->meta_.revision = 0;
    owner->ledgers_ = KineticLedgers{}; owner->ledgers_sha256_ = owner->ledgers_.sha256();
    // (f) identity
    owner->provenance_ = std::move(p); owner->meta_.provenance = owner->provenance_.sha256();
    owner->heat_identity_ = genesis_heat_identity(model->material_sha256(), owner->meta_.provenance);
    owner->field_sha256_ = owner->dev_->state_digest(0, api, owner->configuration_sha256(), owner->meta_);
    owner->heat_sha256_ = heat_sha256_of(owner->heat_identity_, *owner->thermal_);
    owner->clock_guard();
    return owner;
}

// Closed, source-free coupled proposal on the resident owner (same fixed point, checks and
// order as the v3 prepare_closed_proposal; the field prepare runs on the device).
inline ResidentJointProposal prepare_closed_proposal(ResidentOwner& owner, const native::GapPicOptions& o, std::span<const native::ArenaRateEntry> census,
                                                     hall_driver::MotionBinding common, const std::string& source, MotionProvider& provider,
                                                     const StepInputs& inputs) {
    // v3.4: the device path; a census above the configured carrier cap takes the v3.3 path (its refusal text and point)
    const std::size_t rows = inputs.device_rows ? inputs.device_count : census.size();
    if (owner.dev_->options().field_digest_form == gpu::FieldDigestForm::blocks_v71 && !owner.device_mode())
        throw std::runtime_error("field digest form blocks-v7.1 requires the device coupled step (host work device, native v7 motion forms)");
    if (owner.device_mode() && rows <= owner.model_->limits().maximum_carriers)
        return prepare_closed_proposal_device(owner, o, census, std::move(common), source, inputs);
    using clock = std::chrono::steady_clock;
    const auto since = [](clock::time_point a) { return std::chrono::duration<double>(clock::now() - a).count(); };
    const auto start = clock::now();
    ProposalTimings T; T.motion_provider = provider.name();
    owner.clock_guard(); owner.admit_census(residents(census, inputs.newborns));
    need(digest(source), "closed proposal requires unchanged source identity");
    T.admit_census_s = since(start);
    need(common.phase == hall_driver::ParticlePhase::full_time_position_and_proper_velocity && common.accepted_step == owner.revision() &&
         common.begin_s == owner.time_s() && common.end_s == common.begin_s + common.dt_s && std::isfinite(common.dt_s) && common.dt_s > 0 &&
         common.dt_s <= owner.model_->dt_max(), "closed proposal common clock mismatch");
    need(common.mesh_sha256 == owner.motion_.mesh_sha256 && owner.model_->motion_mesh_sha256(o) == owner.motion_.mesh_sha256 &&
         common.material_sha256 == owner.model_->material_sha256(), "closed proposal mesh/material identity mismatch");
    owner.interval_guard(common.dt_s);
    const auto sigma = owner.thermal_->conductivity(); const auto sigma_sha = Model::conductivity_sha256(sigma);
    Vec E, B;
    owner.dev_->accepted_rounded(E, B);
    owner.dev_->begin_fixed_point();
    provider.set_newborns(inputs.newborns, inputs.birth_scope); provider.set_census_states(inputs.census_states); provider.begin(owner.motion_, census);
    struct End { MotionProvider& p; ~End() { p.end(); } } end{provider};
    for (unsigned iteration = 0; iteration < 16; ++iteration) {
        auto binding = common; binding.trial_field_sha256 = hall_driver::total_field_sha256(E, B);
        provider.set_verified_field(&E, &B, binding.trial_field_sha256);
        auto t = clock::now(); FreeCharge after;
        auto motion = provider.motion(binding, E, B, after);
        need(motion.births().size() == inputs.newborns.size(), "motion provider did not move every newborn");
        T.motion_s += since(t); t = clock::now();
        const auto bank = bank_transfer(owner.mesh(), owner.geometry_, motion, owner.model_->exterior_support());
        const auto support = DepositSupportBuilder::build(owner.mesh(), owner.geometry_, census, motion);
        T.support_s += since(t); t = clock::now();
        auto in = owner.inputs(common.dt_s, sigma);
        in.plasma = &motion.current_A_m(); in.after = &after; in.support = &support.absolute_C(); in.bank = &bank;
        ++owner.token_;
        auto step = owner.dev_->prepare(in, owner.meta_);
        T.field_prepare_s += since(t);
        T.solve_s += step.solve.solve_s; T.residual_check_s += step.solve.check_s; T.last_solve = step.solve;
        T.continuity_ratio = std::max(T.continuity_ratio, step.continuity_ratio);
        if (owner.dev_->midpoint_unchanged()) {
            t = clock::now();
            if (motion.forms_v7()) MotionAssembly::seal_v7(motion); else { MotionAssembly::seal_v3(motion); MotionAssembly::seal_v2(motion); }
            T.motion_v2_seal_s = since(t);
            t = clock::now();
            const auto loads = step_loads(owner, motion, inputs.loads, common);
            Vec current = motion.current_A_m();
            auto fh = owner.finish(std::move(step), common.dt_s, sigma_sha, loads, &motion, std::move(after), std::move(current));
            T.heat_s = since(t); T.iterations = iteration + 1; T.total_s = since(start);
            ResidentJointProposal jp(std::move(fh), std::move(motion), source, iteration + 1, std::move(T), std::move(E), std::move(B));
            jp.census_states_ = provider.census_v2().states;
            return jp;
        }
        owner.dev_->midpoint_download(E, B);
    }
    throw std::runtime_error("closed proposal did not reach the identical authoritative midpoint within 16 iterations");
}
// v3.4: the coupled proposal on the device (host-work mode device, native v7 motion forms). Same fixed point,
// checks, refusal texts and order as the v3.3 path above; the census admission, the motion fold, the deposit
// support, the record v7 and the independent residual run on the device (DeviceStep, ResidentField); the host
// keeps O(1) checks, the bank transfer and the sink/birth records (O(sinks + newborns)), the trial-field digest of
// each iterate (streaming form: O(mesh) download and hash) and, with host_record_arrays, the record array downloads.
inline ResidentJointProposal prepare_closed_proposal_device(ResidentOwner& owner, const native::GapPicOptions& o, std::span<const native::ArenaRateEntry> census,
                                                            hall_driver::MotionBinding common, const std::string& source, const StepInputs& inputs) {
    using clock = std::chrono::steady_clock;
    const auto since = [](clock::time_point a) { return std::chrono::duration<double>(clock::now() - a).count(); };
    const auto start = clock::now();
    ProposalTimings T; T.motion_provider = "device-step-v34";
    auto& step = owner.device_step();
    step.reset_times();
    owner.clock_guard();
    {
        gpu::CensusSource src;
        src.rows = census; src.packet_keys = inputs.census_packet_keys; src.state_indices = inputs.census_state_indices; src.states = inputs.census_states;
        src.device_rows = inputs.device_rows; src.device_keys = inputs.device_keys; src.device_states = inputs.device_states; src.device_count = inputs.device_count;
        src.packet_of = inputs.packet_of;
        step.begin(src, inputs.newborns, common.begin_s, common.begin_s + common.dt_s, common.dt_s, *owner.dev_, owner.meta_);
    }
    need(digest(source), "closed proposal requires unchanged source identity");
    T.admit_census_s = since(start);
    need(common.phase == hall_driver::ParticlePhase::full_time_position_and_proper_velocity && common.accepted_step == owner.revision() &&
         common.begin_s == owner.time_s() && common.end_s == common.begin_s + common.dt_s && std::isfinite(common.dt_s) && common.dt_s > 0 &&
         common.dt_s <= owner.model_->dt_max(), "closed proposal common clock mismatch");
    need(common.mesh_sha256 == owner.motion_.mesh_sha256 && owner.model_->motion_mesh_sha256(o) == owner.motion_.mesh_sha256 &&
         common.material_sha256 == owner.model_->material_sha256(), "closed proposal mesh/material identity mismatch");
    owner.interval_guard(common.dt_s);
    const auto sigma = owner.thermal_->conductivity(); const auto sigma_sha = Model::conductivity_sha256(sigma);
    owner.dev_->begin_fixed_point();
    const auto trial = owner.dev_->trial_addresses();
    const bool scope = inputs.birth_scope || !inputs.newborns.empty();
    const bool host_arrays = owner.dev_->options().host_record_arrays;
    const auto cpu_s = [] { timespec ts{}; clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts); return double(ts.tv_sec) + ts.tv_nsec * 1e-9; };
    // The trial-field digest of the binding (O(mesh) download + streaming SHA-256) is computed once, for the accepted
    // iterate: the iterates before it are discarded, the device motion reads the trial buffers, not the digest, and the
    // fixed point leaves the trial buffers unchanged (ResidentField::midpoint_unchanged copies only when they differ).
    // The iterates carry a placeholder that is never sealed (v3.3 computed the digest per iterate; same value, same
    // decisions: the digest enters no check but its own format and the record header).
    static const std::string deferred(64, '0');
    const Vec no_field;
    for (unsigned iteration = 0; iteration < 16; ++iteration) {
        auto t = clock::now(); double c0 = cpu_s();
        auto binding = common;
        binding.trial_field_sha256 = deferred;
        const VerifiedField verified{&no_field, &no_field, deferred};
        auto it = step.motion(binding, trial[0], trial[1], verified, scope);
        auto& motion = it.receipt;
        need(motion.births().size() == inputs.newborns.size(), "motion provider did not move every newborn");
        T.motion_s += since(t); t = clock::now();
        c0 = cpu_s();
        const auto bank = gpu::bank_transfer_sparse(owner.mesh(), owner.geometry_, motion, owner.model_->exterior_support());
        step.support_check();
        T.support_s += since(t); T.support_cpu_s += cpu_s() - c0; t = clock::now(); c0 = cpu_s();
        auto in = owner.inputs(common.dt_s, sigma);
        in.plasma_dev = it.current; in.after_hi_dev = it.after_hi; in.after_lo_dev = it.after_lo; in.after_gross_dev = it.after_gross;
        in.support_dev = it.support; in.bank_sparse = &bank;
        ++owner.token_;
        auto fstep = owner.dev_->prepare(in, owner.meta_);
        const bool fixed = owner.dev_->midpoint_unchanged();
        T.field_prepare_s += since(t); T.field_prepare_cpu_s += cpu_s() - c0;
        T.solve_s += fstep.solve.solve_s; T.residual_check_s += fstep.solve.check_s; T.last_solve = fstep.solve;
        T.continuity_ratio = std::max(T.continuity_ratio, fstep.continuity_ratio);
        if (fixed) {
            t = clock::now(); c0 = cpu_s();
            const auto form = owner.dev_->options().field_digest_form;
            Vec E, B;
            if (form == gpu::FieldDigestForm::streaming_v7) {
                owner.dev_->trial_download(E, B);
                MotionAssembly::set_trial_field_sha256(motion, hall_driver::total_field_sha256(E, B));
            } else if (host_arrays) owner.dev_->trial_download(E, B);   // host copy for native's host admission only
            T.motion_trial_digest_s += since(t); T.trial_digest_cpu_s += cpu_s() - c0;
            t = clock::now();
            // block form: the binding's trial-field digest and the current's digest from the device block-list roots
            step.seal(motion, host_arrays, trial[0], trial[1], form);
            binding.trial_field_sha256 = motion.binding().trial_field_sha256;
            T.motion_v2_seal_s = since(t);
            t = clock::now(); c0 = cpu_s();
            const auto loads = step_loads(owner, motion, inputs.loads, common);
            auto fh = owner.finish(std::move(fstep), common.dt_s, sigma_sha, loads, &motion, FreeCharge{}, Vec());
            fh.device_free_ = true;
            T.heat_s = since(t); T.finish_cpu_s = cpu_s() - c0; T.iterations = iteration + 1; T.total_s = since(start);
            ResidentJointProposal jp(std::move(fh), std::move(motion), source, iteration + 1, std::move(T), host_arrays ? std::move(E) : Vec(),
                                     host_arrays ? std::move(B) : Vec());
            jp.census_states_ = step.census_v2().states;
            jp.device_ = true; jp.trial_digest_ = binding.trial_field_sha256;
            jp.seal_generation_ = step.record_view().seal_generation;
            return jp;
        }
    }
    throw std::runtime_error("closed proposal did not reach the identical authoritative midpoint within 16 iterations");
}
}  // namespace hall_corrective_v3
