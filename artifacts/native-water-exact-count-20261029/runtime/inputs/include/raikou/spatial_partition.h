#pragma once
#include "raikou/spatial_mesh.h"
#include "raikou/distributed.h"
#include <limits>
#include <map>
#include <optional>

namespace raikou::spatial::partition {
using GlobalIndex = backend::GlobalIndex;
inline constexpr GlobalIndex no_dof = std::numeric_limits<GlobalIndex>::max();
struct CellInterval { std::size_t begin = 0, count = 0; };

// Geometry and offset metadata may be replicated; no state or matrix is stored.
// The base block is species, phi, photons, temperatures. Each cell then owns
// sigma/face-potential pairs for its dielectric walls, in global wall order.
// The optional circuit unknown follows the last cell's final wall pair.
class DofLayout {
public:
    DofLayout(const mesh::Geometry&, std::size_t species, std::size_t photons,
        std::size_t thermal, std::span<const unsigned char> dielectric_walls,
        bool circuit, const backend::ExecutionOptions& = {},
        std::optional<CellInterval> owned_cells = {});
    [[nodiscard]] std::size_t cells() const noexcept { return offsets_.size() - 1; }
    [[nodiscard]] std::size_t components() const noexcept { return components_; }
    [[nodiscard]] std::size_t phi_component() const noexcept { return species_; }
    [[nodiscard]] std::size_t photon_component(std::size_t group) const;
    [[nodiscard]] std::size_t thermal_component(std::size_t field) const;
    [[nodiscard]] GlobalIndex cell_dof(std::size_t cell, std::size_t component) const;
    [[nodiscard]] GlobalIndex cell_begin(std::size_t cell) const;
    [[nodiscard]] GlobalIndex cell_end(std::size_t cell) const;
    [[nodiscard]] GlobalIndex sigma_dof(std::size_t wall) const;
    [[nodiscard]] GlobalIndex dielectric_dof(std::size_t wall) const;
    [[nodiscard]] bool is_dielectric(std::size_t wall) const;
    [[nodiscard]] GlobalIndex circuit_dof() const noexcept { return circuit_; }
    [[nodiscard]] CellInterval owned_cells() const noexcept { return owned_; }
    [[nodiscard]] backend::RowPartition row_partition() const noexcept { return rows_; }
    [[nodiscard]] bool owns_cell(std::size_t cell) const noexcept;
    [[nodiscard]] bool owns_dof(GlobalIndex dof) const noexcept;
    [[nodiscard]] const backend::DistributedHalo& control_halo() const { return *ownership_; }
    [[nodiscard]] backend::DistributedHalo make_halo(std::span<const GlobalIndex> requested) const;
private:
    std::size_t species_ = 0, photons_ = 0, thermal_ = 0, components_ = 0;
    std::vector<GlobalIndex> offsets_, sigma_;
    GlobalIndex circuit_ = no_dof;
    CellInterval owned_;
    backend::RowPartition rows_;
    backend::ExecutionOptions options_;
    std::optional<backend::DistributedHalo> ownership_;
};

struct StencilDependencies {
    std::span<const mesh::ScalarStencil> cells;
    std::span<const mesh::ScalarStencil> faces;
    std::span<const mesh::ScalarStencil> walls;
    std::span<const std::array<mesh::ScalarStencil, 3>> gradients;
    // Empty selects every base component. Boundary slots request both wall
    // unknowns where present, and the circuit unknown where enabled.
    std::span<const std::size_t> components;
};
// Includes direct face neighbors plus the supplied expanded cell/face/wall
// patches touching owned cells. Second-ring and symbolic wall dependencies are
// retained. Explicit extras support nonlocal circuit/other physics dependencies.
// Returns sorted unique nonowned DOFs; validation is collective.
[[nodiscard]] std::vector<GlobalIndex> requested_ghost_ids(const DofLayout&,
    const mesh::Geometry&, const StencilDependencies& = {},
    std::span<const GlobalIndex> extras = {});

class OwnedValues {
public:
    explicit OwnedValues(backend::DistributedHalo, double initial_value = 0);
    OwnedValues(backend::DistributedHalo, std::vector<double> owned_values);
    void exchange(); // collective; validates local values on every rank
    [[nodiscard]] double value(GlobalIndex) const;
    [[nodiscard]] std::span<double> local_values(); // marks halo stale
    [[nodiscard]] std::span<const double> local_values() const noexcept { return local_; }
    [[nodiscard]] const backend::DistributedHalo& halo() const noexcept { return halo_; }
    [[nodiscard]] std::size_t storage_entries() const noexcept { return local_.size() + ghosts_.size(); }
    // Includes old/new receive arrays plus packed sends during exchange.
    [[nodiscard]] std::size_t exchange_working_entries_bound() const noexcept;
private:
    backend::DistributedHalo halo_;
    std::vector<double> local_, ghosts_;
    bool current_ = false;
};

struct LocalSystem {
    backend::DistributedCsrMatrix matrix;
    std::vector<double> residual;
};
// Assembly never accumulates remote rows. A globally owned circuit row needs
// an explicit reduction or assembly by its owner with all dependencies present.
// add_* contain no collectives: local failures are retained and rethrown on all
// ranks by finish(), which every rank must call. finish consumes this assembly.
class LocalAssembly {
public:
    explicit LocalAssembly(backend::DistributedHalo);
    void add_residual(GlobalIndex row, double value);
    void add_matrix(GlobalIndex row, GlobalIndex column, double value);
    [[nodiscard]] LocalSystem finish();
    [[nodiscard]] std::size_t residual_storage_entries() const noexcept { return residual_.size(); }
    [[nodiscard]] std::size_t matrix_storage_entries() const noexcept { return entries_; }
private:
    backend::DistributedHalo halo_;
    std::vector<double> residual_;
    std::vector<std::map<GlobalIndex, double>> rows_;
    std::string error_;
    std::size_t entries_ = 0;
    bool finished_ = false;
    [[nodiscard]] std::optional<std::size_t> owned_row(GlobalIndex);
};
} // namespace raikou::spatial::partition
