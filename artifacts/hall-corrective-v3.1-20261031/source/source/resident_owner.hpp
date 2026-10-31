#pragma once
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

namespace hall_corrective_v3 {
class ResidentOwner;
class ResidentJointProposal;
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
    gpu::ResidentStep step_; gpu::SlotMeta meta_; ThermalTrialV3 heat_;
    std::weak_ptr<const OwnerLifetime> owner_; std::uint64_t generation_ = 0, token_ = 0;
    std::string before_, after_, configuration_, conductivity_, midpoint_seal_, before_heat_, after_heat_;
    FreeCharge after_free_; Vec current_; double dt_ = 0;
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
    ResidentFieldHeatProposal field_heat_; hall_driver::MotionReceipt motion_; std::string source_;
    Vec midE_, midB_;   // host copy of the accepted fixed-point midpoint (equal to the device midpoint)
    double combined_balance_ = 0; unsigned iterations_ = 0; ProposalTimings timings_; StepEnergy energy_;
    ResidentJointProposal(ResidentFieldHeatProposal&& p, hall_driver::MotionReceipt&& motion, std::string source, unsigned iterations, ProposalTimings timings,
                          Vec midE, Vec midB)
        : field_heat_(std::move(p)), motion_(std::move(motion)), source_(std::move(source)), midE_(std::move(midE)), midB_(std::move(midB)),
          iterations_(iterations), timings_(std::move(timings)) {
        // Same checks as JointProposal (root_api.hpp).
        const auto& r = field_heat_.step_.receipt; const auto& h = field_heat_.heat_;
        const bool births = !motion_.births().empty();
        const double work = motion_.budget().kinetic_change_J - r.plasma_work_J;
        need(std::isfinite(work) && std::abs(work) <= 1e-20 + 1e-10 * (std::abs(motion_.budget().kinetic_change_J) + std::abs(r.plasma_work_J)),
             "unchanged particle/Maxwell work bound");
        combined_balance_ = field_heat_.field_heat_balance_J() + work;
        double transfers = std::abs(r.field_delta_J) + std::abs(motion_.budget().kinetic_change_J) + std::abs(h.thermal_delta_J) + std::abs(r.outer_loss_J) +
                           std::abs(r.outer_supply_J) + std::abs(r.circuit_work_J);
        if (births || h.loads) transfers += std::abs(h.wall_heat_J) + h.external_heating_J + h.external_cooling_J;
        need(std::isfinite(combined_balance_) && std::abs(combined_balance_) <= 2e-20 + 1e-10 * transfers + h.roundoff_bound_J,
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
};

class ResidentOwner {
    friend ResidentJointProposal prepare_closed_proposal(ResidentOwner&, const native::GapPicOptions&, std::span<const native::ArenaRateEntry>,
                                                         hall_driver::MotionBinding, const std::string&, MotionProvider&, const StepInputs&);
    std::shared_ptr<const Model> model_; std::shared_ptr<gpu::Device> device_; const native::GapPicOptions geometry_; const MotionContext motion_;
    std::unique_ptr<gpu::ResidentField> dev_;
    gpu::SlotMeta meta_; FreeCharge free_;
    std::unique_ptr<ThermalV3> thermal_; KineticLedgers ledgers_; std::string ledgers_sha256_;
    const std::shared_ptr<const OwnerLifetime> lifetime_ = std::make_shared<const OwnerLifetime>();
    std::uint64_t generation_ = 0, token_ = 0; bool vacuum_ = false;
    std::string field_sha256_, heat_sha256_;
    static std::shared_ptr<const Model> checked_model(std::shared_ptr<const Model> m) {
        need(bool(m), "missing immutable model");
        need(m->solver().name() == "gpu1-resident-v1", "ResidentOwner requires a Model configured with the resident field backend");
        return m;
    }
    void clock_guard() const { need(meta_.time_s == thermal_->time_s() && meta_.revision == thermal_->revision(), "accepted field/heat clock mismatch"); }
    std::string conductivity_sha256() const { return Model::conductivity_sha256(thermal_->conductivity()); }
    void admit_census(std::span<const native::ArenaRateEntry> census) const {
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
        p.after_heat_ = heat_sha256_of(model_->material_sha256(), p.heat_);
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
        swap(field_sha256_, p.after_); swap(heat_sha256_, p.after_heat_);
        swap(ledgers_, p.ledgers_); swap(ledgers_sha256_, p.ledgers_sha256_);
        ++generation_; ++token_; p.live_ = false; p.final_approved_ = false;
    }
public:
    // Bootstrap from an accepted host state (same checks as Owner): the census charge must
    // equal the declared free charge; the state is uploaded and admitted on the device.
    ResidentOwner(std::shared_ptr<const Model> m, std::shared_ptr<gpu::Device> device, axial::State field, ChargeLedger charge,
                  const native::GapPicOptions& o, std::span<const native::ArenaRateEntry> census, gpu::ResidentOptions options = {})
        : model_(checked_model(std::move(m))), device_(std::move(device)), geometry_(o),
          motion_(deferred_v3(with_dt_max(make_motion_context(model_->mesh(), geometry_, model_->mesh_json(), model_->limits().maximum_carriers), model_->dt_max()))),
          vacuum_(census.empty()) {
        dev_ = std::make_unique<gpu::ResidentField>(device_, *model_, options);
        need(ResidentBackend(dev_->descriptor_sha256()).descriptor_sha256() == model_->solver().descriptor_sha256(),
             "resident kernels differ from the backend descriptor bound into the Model configuration");
        const auto actual = charge_of_census(mesh(), geometry_, census, model_->limits().maximum_carriers);
        auto same_words = [](const WordVec& a, const WordVec& b) {
            if (a.size() != b.size()) return false;
            for (std::size_t i = 0; i < a.size(); ++i) if (a[i].hi != b[i].hi || a[i].lo != b[i].lo) return false;
            return true; };
        need(same_words(actual.net_C, charge.particles.net_C) && actual.gross_C == charge.particles.gross_C,
             "initial/restored net charge differs from actual canonical carriers");
        model_->admit(field, charge);   // host admission of the bootstrap state (sealed arithmetic)
        meta_.time_s = field.time_s; meta_.revision = field.revision; meta_.conductor_declaration = charge.conductor_declaration;
        meta_.exterior_declaration = charge.exterior_declaration; meta_.residue_present = !charge.residue_C.empty();
        free_ = charge.particles;
        dev_->upload(0, field, charge);
        (void)dev_->admit(0, meta_);
        admit_census(census);
        thermal_ = std::make_unique<ThermalV3>(mesh(), model_->drive().start(), model_->thermal_declaration(), model_->dt_max());
        field_sha256_ = dev_->state_digest(0, api, configuration_sha256(), meta_);
        heat_sha256_ = heat_sha256_of(model_->material_sha256(), *thermal_);
        ledgers_sha256_ = ledgers_.sha256();
        clock_guard();
    }
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
    const FreeCharge& free_charge() const noexcept { return free_; }
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
            if (heat_sha256_of(model_->material_sha256(), *thermal_) != heat_sha256_ || p.after_heat_ != heat_sha256_of(model_->material_sha256(), p.heat_) ||
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
        try {
            if (!preflight_field_heat(p.field_heat_)) return false;
            const auto& m = p.motion(); const auto& b = m.binding();
            { std::vector<Newborn> nb; for (const auto& x : m.births()) nb.push_back({x.native_id, x.birth_time_s}); admit_census(residents(census, nb)); }
            if (m.census_sha256() != canonical_census_sha256(census, motion_.maximum_carriers) || !digest(m.seal_sha256()) || !digest(m.record_v3_sha256()) ||
                MotionAssembly::record_v3(m) != m.record_v3_sha256() || b.phase != hall_driver::ParticlePhase::full_time_position_and_proper_velocity ||
                b.accepted_step != revision() || b.begin_s != time_s() || b.end_s != p.end_s() || b.dt_s != p.field_heat_.dt_ ||
                b.mesh_sha256 != motion_.mesh_sha256 || b.material_sha256 != model_->material_sha256() || b.kinetic_policy_sha256 != motion_.policy_sha256 ||
                b.trial_field_sha256 != hall_driver::total_field_sha256(p.midE_, p.midB_) || m.current_A_m() != p.field_heat_.current_)
                return false;
            auto actual = charge_of_endpoints(mesh(), geometry_, m.endpoints(), motion_.maximum_carriers);
            if (!m.v2_identical() && (MotionAssembly::sink_birth_digest(m).empty() ||
                                      p.field_heat_.candidate_ledgers().sinks.count() != kinetic_ledgers().sinks.count() + m.sinks().size() ||
                                      p.field_heat_.candidate_ledgers().births.count() != kinetic_ledgers().births.count() + m.births().size()))
                return false;
            const auto& after = p.field_heat_.after_free_;
            if (actual.net_C.size() != after.net_C.size() || actual.gross_C != after.gross_C) return false;
            for (std::size_t i = 0; i < actual.net_C.size(); ++i) if (actual.net_C[i].hi != after.net_C[i].hi || actual.net_C[i].lo != after.net_C[i].lo) return false;
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
};

// Closed, source-free coupled proposal on the resident owner (same fixed point, checks and
// order as the v3 prepare_closed_proposal; the field prepare runs on the device).
inline ResidentJointProposal prepare_closed_proposal(ResidentOwner& owner, const native::GapPicOptions& o, std::span<const native::ArenaRateEntry> census,
                                                     hall_driver::MotionBinding common, const std::string& source, MotionProvider& provider,
                                                     const StepInputs& inputs) {
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
    provider.set_newborns(inputs.newborns, inputs.birth_scope); provider.begin(owner.motion_, census);
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
            t = clock::now(); MotionAssembly::seal_v3(motion); MotionAssembly::seal_v2(motion); T.motion_v2_seal_s = since(t);
            t = clock::now();
            const auto loads = step_loads(owner, motion, inputs.loads, common);
            Vec current = motion.current_A_m();
            auto fh = owner.finish(std::move(step), common.dt_s, sigma_sha, loads, &motion, std::move(after), std::move(current));
            T.heat_s = since(t); T.iterations = iteration + 1; T.total_s = since(start);
            return ResidentJointProposal(std::move(fh), std::move(motion), source, iteration + 1, std::move(T), std::move(E), std::move(B));
        }
        owner.dev_->midpoint_download(E, B);
    }
    throw std::runtime_error("closed proposal did not reach the identical authoritative midpoint within 16 iterations");
}
}  // namespace hall_corrective_v3
