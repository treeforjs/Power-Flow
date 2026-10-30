#pragma once
#include "raikou/json.h"
#include <cstddef>
#include <functional>
#include <span>
#include <string>
#include <vector>

namespace raikou::backend {
struct CsrMatrix {
    std::size_t rows{};
    std::vector<std::size_t> offsets;
    std::vector<int> columns;
    std::vector<double> values;
    void validate() const;
};
struct ExecutionOptions {
    std::string device{"cpu"};
    std::size_t threads{};
    bool mpi{false};
    // MPI_Comm_c2f supplied by an embedding host; -1 selects the launcher team.
    int communicator_f{-1};
    double tolerance{1e-10};
    int max_iterations{1000};
    int restart{40};
    // Distributed CPU solves also support block_ilu0 (rank-local ILU(0)).
    std::string preconditioner{"jacobi"};
    int cuda_device{-1}; // -1 preserves the caller's current device; explicit per-rank ordinal otherwise.
};
struct LinearResult {
    std::vector<double> solution;
    bool converged{false};
    int iterations{};
    double residual{};
};
// The host owns MPI initialization. MPI-enabled calls require an initialized runtime.
void parallel_for(std::size_t n, std::size_t threads, const std::function<void(std::size_t)>& action);
using CsrOperator = std::function<std::vector<double>(std::span<const double>)>;
// Owns its matrix/device storage and duplicates an MPI communicator. Each MPI
// application remains collective. This compatibility API receives/returns full
// vectors; DistributedLinearOperator exposes owned rows and halo-only matvecs.
CsrOperator make_csr_operator(const CsrMatrix&, const ExecutionOptions& = {});
std::vector<double> csr_multiply(const CsrMatrix&, std::span<const double>, const ExecutionOptions& = {});
LinearResult solve_csr(const CsrMatrix&, std::span<const double>, const ExecutionOptions& = {});
// Independent jobs are partitioned across ranks and local worker threads. The
// callback must not call MPI. Every rank supplies the same ordered job list.
Json map_jobs(const Json& jobs, const std::function<Json(const Json&)>&,
              const ExecutionOptions& = {});
// Collective request distribution; only the root's value is used. In serial
// mode this returns the supplied value. Intended for launcher input streams.
Json broadcast_json(const Json& value, int root = 0, const ExecutionOptions& = {});
Json capabilities_json();
Json linear_solve_json(const Json&);
bool cuda_available();
int mpi_rank();
int mpi_size();
class Runtime {
public:
    explicit Runtime(bool initialize_mpi);
    ~Runtime();
    Runtime(const Runtime&) = delete;
    Runtime& operator=(const Runtime&) = delete;
private:
    bool owns_mpi_{false};
};
} // namespace raikou::backend
