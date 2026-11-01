#pragma once
// Hall v3.2 device-resident field step (namespace hall_corrective_v3::gpu); v3.1 with the admitted
// Gauss ledger words in each slot (see resident_field.cpp, field_contract.hpp).
//
// ResidentField keeps two field-state slots on the device (0 = accepted, 1 = candidate:
// E, B words and the complete charge ledger as structure-of-arrays) and runs the whole
// arithmetic of Model::prepare there with the deterministic kernels of resident_kernels.cu
// (physical state bitwise equal to the host arithmetic for the same increment) and the
// gpu1 CG kernels (cg_kernels.hpp) on the device right-hand side. The host keeps:
//   * the independent FP64 residual ||b - A x|| with the sealed apply arithmetic and the
//     sealed long-double norm (rhs and x downloaded once per solve);
//   * the thermal owner (cell Joule power of the metal cells downloaded once per step);
//   * the identities: block-Merkle digests of the device arrays (leaves hashed on the
//     device, root on the host, a random sample of leaves re-hashed on the host every
//     time, every leaf in qualification mode).
// The header exposes no CUDA types.
#include "gpu_backends.hpp"
#include "../control_receipt.hpp"
#include <array>
#include <map>
#include <random>

namespace hall_corrective_v3::gpu {
struct SparseBank;   // device_step_v34.hpp

// v3.4 host-work mode of the resident owner and its field step:
//   device: the coupled step's moved items on the device (device_step_v34.hpp; independent residual on the device);
//   v33:    the v3.3 host path of every moved item (verification mode; same decisions, texts, ledgers and state).
// Default: device, or HALL_HOST_WORK=v33 in the environment (qualification switch of the unchanged v3.3 suites).
enum class HostWork { v33, device };
HostWork default_host_work();
const char* host_work_name(HostWork);
// Form of the trial-field and current digests in the motion record (binding trial_field_sha256, header
// current_A_m_sha256): streaming_v7 = total-midpoint-E-B-binary64-v1 / motion-array-binary64-v1 (native v7, O(mesh)
// host hashing); blocks_v71 = total-midpoint-E-B-blocks-v1 / motion-array-blocks-v1 (native v7.1 API-071.txt 1.1:
// block-list roots of the device arrays, O(1) host; device coupled step only, forms library of native v7.1).
enum class FieldDigestForm { streaming_v7, blocks_v71 };
const char* field_digest_form_name(FieldDigestForm);
struct ResidentOptions {
    enum class Solve { device_cg, host_cpu_pcg };
    HostWork host_work = default_host_work();
    // v3.4 device path: download the accepted iterate's endpoint/path records and trial midpoint into the proposal
    // (native v7 host admission); false with a device admission (native v7.1 device lane).
    bool host_record_arrays = true;
    FieldDigestForm field_digest_form = FieldDigestForm::streaming_v7;
    Solve solve = Solve::device_cg;     // host_cpu_pcg: the sealed foil::solve on the downloaded rhs (parity mode)
    int verify_samples = 16;            // leaves re-hashed on the host per digest
    bool verify_all = false;            // qualification: every leaf re-hashed on the host
    std::uint64_t sample_seed = 0;      // 0: seeded from std::random_device
    // Diagnostic: per-node charge bound of each solve's residual, sum_{e at n} sqrt(M_e) |rho_e| / len_e
    // (rho = b - A x of the independent check), accumulated over accepted steps.
    bool trace_residual_charge = false;
};

// Per-step quantities for a step-size controller (control_receipt.hpp).
using ::hall_corrective_v3::ControlReceipt;

struct ResidentStep {
    axial::Receipt receipt;             // scalars and cell_joule_W (full cell vector); mid_E/solid/short arrays empty
    Constraints constraints;            // candidate absolute Gauss/divB maxima (guard values)
    double continuity_ratio = 0;
    SolveRecord solve;
    ControlReceipt control;
    bool residue_present = false;      // candidate ledger carries the bank residue array
    double residual_charge_step_max_C = 0, residual_charge_cumulative_max_C = 0;   // trace_residual_charge only
    std::map<std::string, double> host_ms;
};

// Metadata of a slot that is not on the device.
struct SlotMeta {
    double time_s = 0; std::uint64_t revision = 0;
    std::string conductor_declaration, exterior_declaration;
    bool residue_present = false;
};

struct ResidentInputs {
    double dt = 0, current = 0;
    const Vec* sigma = nullptr;                  // cell conductivity (required)
    const Vec* plasma = nullptr;                 // edge current (nullptr: zero)
    const FreeCharge* after = nullptr;           // endpoint free charge (nullptr: the accepted free charge)
    const Vec* support = nullptr;                // deposit support absolute_C (nullptr: zero)
    const BankTransfer* bank = nullptr;          // bank transfer (nullptr or empty: none)
    double drive_rel_change = 0, drive_knot_interval_s = 0;
    // v3.4 device inputs (device addresses on the session; read in place, no upload): plasma current [ne],
    // endpoint free charge words and gross [nn], deposit support [nn]; sparse bank transfer (nodes ascending).
    std::uint64_t plasma_dev = 0, after_hi_dev = 0, after_lo_dev = 0, after_gross_dev = 0, support_dev = 0;
    const SparseBank* bank_sparse = nullptr;
};

struct DigestStats { long digests = 0, leaves_hashed_on_host = 0, samples = 0; double host_ms = 0; };

class ResidentField {
    struct Impl;
    std::unique_ptr<Impl> p_;
public:
    ResidentField(std::shared_ptr<Device>, const Model&, ResidentOptions = {});
    ~ResidentField();
    ResidentField(const ResidentField&) = delete;
    ResidentField& operator=(const ResidentField&) = delete;

    // Descriptor of the resident step (kernel sources, NVRTC options and target, block size,
    // independent-residual rule, digest scheme); bound into the Model configuration.
    static std::string descriptor_text();
    const std::string& descriptor_sha256() const noexcept;
    std::string compile_json() const;
    const ResidentOptions& options() const noexcept;

    // Slot transfer (bootstrap, checkpoints, tests).
    void upload(int slot, const axial::State&, const ChargeLedger&);
    void download(int slot, axial::State&, ChargeLedger&, const SlotMeta&);

    // Absolute Gauss/divB admission of a slot (throws the sealed refusal texts in the sealed order).
    Constraints admit(int slot, const SlotMeta&, const FreeCharge* census_charge = nullptr);

    // Candidate (slot 1) from the accepted state (slot 0). Throws the sealed refusal texts.
    ResidentStep prepare(const ResidentInputs&, const SlotMeta& accepted);

    // Merkle identities. state_digest hashes the slot's arrays on the device and verifies
    // a sample (or all) of the leaves on the host.
    std::string state_digest(int slot, const std::string& api, const std::string& configuration, const SlotMeta&);
    std::string midpoint_seal();                 // candidate midpoint (mid E, mid B) role seal
    void midpoint_download(Vec& E, Vec& B);      // candidate midpoint (host copy)
    void accepted_rounded(Vec& E, Vec& B);       // Model::at(accepted)
    // Fixed point: true when the candidate midpoint equals the previous iterate's midpoint
    // bitwise; the candidate midpoint then becomes the previous one. begin_fixed_point sets
    // the previous midpoint to the rounded accepted words.
    void begin_fixed_point();
    bool midpoint_unchanged();
    // Preflight recomputation: midpoint from both slots equals the stored midpoint bitwise,
    // the slots' leaves equal the digested leaves, both slots admitted.
    bool preflight_recheck(const SlotMeta& accepted, const SlotMeta& candidate);

    // Publication: candidate becomes accepted (pointer swap, no allocation, no device work).
    void publish() noexcept;

    // ---- v3.4 device path
    // Free-charge words and gross of a slot (device addresses: net hi, net lo, gross).
    std::array<std::uint64_t, 3> free_addresses(int slot);
    // Admission of a slot with a free charge given by device addresses (admit_census on the device).
    Constraints admit_device(int slot, const SlotMeta&, std::uint64_t hi, std::uint64_t lo, std::uint64_t gross);
    // Trial midpoint of the fixed point (previous iterate's midpoint; the rounded accepted words at the start)
    // and the candidate midpoint: device addresses of E [ne] and B [nb].
    std::array<std::uint64_t, 2> trial_addresses();
    void trial_download(Vec& E, Vec& B);     // the trial midpoint (host copy)
    bool midpoint_is_trial();                // candidate midpoint bitwise equal to the trial buffers (comparison only)
    std::array<std::uint64_t, 2> midpoint_addresses();
    // Static node/edge tables: node->edges (ptr, edges, signs, ascending edge), edge lengths, edge end nodes.
    struct Tables { std::uint64_t ntp = 0, nte = 0, nts = 0, elen = 0, ea = 0, eb = 0; };
    Tables tables() const noexcept;
    // Independent residual of the last prepare (both modes): exact norms (identity rule hall-v3.4).
    static const char* residual_rule();

    DigestStats digest_stats() const noexcept;
    std::map<std::string, double> device_profile();   // last step device phase times (ms)
    std::string memory_json();
    // Increment history used by the field-accuracy indicator (accepted steps only).
    int history() const noexcept;
    // Diagnostic: per-node |D - L| (the checked absolute Gauss value) and bound of a slot (device admission
    // values), optionally |D| (uncorrected) and the signed admitted Gauss ledger (hi + lo); returns the
    // admission constraints. Also the accumulated residual-charge bound of the accepted steps (trace_residual_charge).
    Constraints gauss_profile(int slot, const SlotMeta&, Vec& defect, Vec& bound, Vec* uncorrected = nullptr, Vec* ledger = nullptr,
                              Vec* divb = nullptr, Vec* divb_bound = nullptr,    // v3.3: per-cell |divB| and bound (cell index order)
                              Vec* divb_floor = nullptr);                        // v3.4: per-cell 8 eps^2 G_c
    const Vec& residual_charge_cumulative() const noexcept;
    // Accumulated per-node |Gauss increment| (|Q1 - Q0 - (div plasma + div metal + div outer) dt| of each accepted step).
    void gauss_increment_cumulative(Vec& out);
};

// Host form of the sealed matrix-free apply with preallocated workspace (bitwise equal to
// the lambda of Model::prepare), and the independent residual check of record.
struct HostApply {
    const Mesh& g; Vec sqrtM, E, B, K;
    explicit HostApply(const Mesh&);
    void apply(const Vec& S, double dt, const Vec& v, Vec& out);
};
// Conduction S from cell conductivity using only the metal cells (bitwise equal to Mesh::conduction).
void conduction_metal(const Mesh&, const std::vector<int>& metal_cells, const Vec& sigma, Vec& S);

}  // namespace hall_corrective_v3::gpu
