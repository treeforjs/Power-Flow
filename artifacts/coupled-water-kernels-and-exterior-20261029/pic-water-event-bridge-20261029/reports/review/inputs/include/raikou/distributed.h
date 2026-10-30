#pragma once

#include "raikou/backend.h"
#include <cstdint>
#include <memory>
#include <span>
#include <string>
#include <vector>

namespace raikou::backend {
using GlobalIndex = std::uint64_t;

// Ranks own consecutive, disjoint intervals in communicator rank order. Empty
// intervals are allowed. Global IDs remain stable when the rank count changes.
struct RowPartition {
    GlobalIndex global_rows{};
    GlobalIndex row_begin{};
    std::size_t local_rows{};
};
RowPartition balanced_partition(GlobalIndex global_rows, int rank, int ranks);
RowPartition distributed_partition(GlobalIndex global_rows, const ExecutionOptions& = {});

struct DistributedCsrMatrix {
    GlobalIndex global_rows{};
    GlobalIndex row_begin{};
    std::vector<std::size_t> offsets;
    std::vector<GlobalIndex> columns;
    std::vector<double> values;
    [[nodiscard]] std::size_t local_rows() const noexcept { return offsets.empty() ? 0 : offsets.size() - 1; }
    [[nodiscard]] RowPartition partition() const noexcept { return {global_rows, row_begin, local_rows()}; }
    void validate_local() const;
};

// Ghosts are unique and ordered by global ID. Packed sends are ordered by peer
// rank, then the requesting peer's ghost order. Neither list contains a global
// vector. The only collective replication during setup is O(number of ranks)
// ownership metadata; applications exchange only requested halo entries.
struct HaloLayout {
    std::vector<GlobalIndex> ghost_global_ids;
    std::vector<std::size_t> send_local_indices;
};

class DistributedHalo {
public:
    DistributedHalo(RowPartition, std::span<const GlobalIndex> requested_global_ids,
                    const ExecutionOptions& = {});
    [[nodiscard]] RowPartition partition() const noexcept;
    [[nodiscard]] const HaloLayout& layout() const noexcept;
    [[nodiscard]] std::size_t compact_index(GlobalIndex) const;
    [[nodiscard]] std::vector<double> exchange(std::span<const double> local_values) const;
    [[nodiscard]] std::vector<double> exchange_packed_halo(std::span<const double> packed_send) const;
    [[nodiscard]] double global_sum(double) const;
    [[nodiscard]] double global_max(double) const;
    // Sum of values on lower communicator ranks (zero on rank zero).
    [[nodiscard]] double exclusive_sum(double) const;
    // Every rank must call this, even when its error is empty. A local failure
    // is reported on all ranks before the next communication phase starts.
    void collective_require(const std::string& local_error) const;
    [[nodiscard]] int rank() const noexcept;
    [[nodiscard]] int ranks() const noexcept;
    [[nodiscard]] int communicator_f() const noexcept;
private:
    struct Impl;
    std::shared_ptr<Impl> impl_;
};

struct CompactCsrView {
    std::size_t rows{};
    std::size_t input_columns{};
    std::span<const std::size_t> offsets;
    std::span<const int> columns;
    std::span<const double> values;
};

struct DistributedLinearResult {
    std::vector<double> local_solution;
    GlobalIndex row_begin{};
    GlobalIndex global_rows{};
    bool converged{false};
    int iterations{};
    double residual{};
    // Upper bound on local vector entries reserved by the Krylov solver.
    // Hessenberg/rotation scalars and matrix/preconditioner entries are separate.
    std::size_t local_working_vector_entries{};
};

// The operator owns its local rows and duplicates the supplied communicator.
// Destroy it before MPI_Finalize to release MPI resources. Collective calls
// must occur in the same order on every rank, on MPI-permitted host threads.
// Local invocation serialization prevents simultaneous uses of its halo tags.
class DistributedLinearOperator {
public:
    explicit DistributedLinearOperator(DistributedCsrMatrix, const ExecutionOptions& = {});
    [[nodiscard]] const DistributedCsrMatrix& matrix() const noexcept;
    [[nodiscard]] CompactCsrView compact_matrix() const noexcept;
    [[nodiscard]] const DistributedHalo& halo() const noexcept;
    [[nodiscard]] const ExecutionOptions& options() const noexcept;
    [[nodiscard]] std::vector<double> apply(std::span<const double> local_values) const;
    [[nodiscard]] DistributedLinearResult solve(std::span<const double> local_rhs) const;
private:
    struct Impl;
    std::shared_ptr<Impl> impl_;
};

DistributedLinearResult solve_distributed(DistributedCsrMatrix, std::span<const double> local_rhs,
                                        const ExecutionOptions& = {});
// CPU matrix-free FGMRES for nonlocal scalar-moment operators. The callback
// returns owned rows and must obey the halo's collective ordering/preflight
// rules. Its sparse preconditioning matrix stores only owned rows.
DistributedLinearResult solve_distributed_matrix_free(const DistributedHalo&,
    const DistributedCsrMatrix& preconditioning_matrix, const CsrOperator& local_apply,
    std::span<const double> local_rhs, const ExecutionOptions& = {});
// Explicit compatibility boundary: gathers a completed local vector once. No
// distributed operator or Krylov iteration uses this routine internally.
std::vector<double> gather_distributed(const DistributedHalo&, std::span<const double> local_values);
DistributedCsrMatrix partition_replicated(const CsrMatrix&, RowPartition);
}
