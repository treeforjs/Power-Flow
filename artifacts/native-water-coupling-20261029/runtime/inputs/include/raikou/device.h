#pragma once
#include "raikou/backend.h"
#include <limits>
#include <memory>
#include <span>

namespace raikou::backend {

// Compact CSR: local owned columns come first, followed by halo columns.
// Row i's diagonal is column i. Arrays are copied at executor construction.
struct DeviceCsrView {
    std::size_t rows{}, columns{};
    std::span<const std::size_t> offsets;
    std::span<const int> indices;
    std::span<const double> values;
};
struct DeviceTransferMetrics {
    std::size_t host_to_device_bytes{}, device_to_host_bytes{}, allocated_bytes{};
};
class DeviceVector {
public:
    virtual ~DeviceVector() = default;
    [[nodiscard]] virtual std::size_t size() const = 0;
};
class DeviceIndexPlan {
public:
    virtual ~DeviceIndexPlan() = default;
    [[nodiscard]] virtual std::size_t size() const = 0;
};
// Vectors and index plans belong to one executor. Bulk data stays on its device;
// only explicit upload/download/gather/scatter and scalar reductions transfer.
class DeviceExecutor {
public:
    virtual ~DeviceExecutor() = default;
    virtual std::unique_ptr<DeviceVector> allocate(std::size_t size) = 0;
    virtual std::unique_ptr<DeviceIndexPlan> index_plan(std::span<const std::size_t>) = 0;
    virtual void upload(std::span<const double>, DeviceVector&) = 0;
    virtual std::vector<double> download(const DeviceVector&) = 0;
    virtual std::vector<double> gather(const DeviceVector&, const DeviceIndexPlan&) = 0;
    virtual void scatter(std::span<const double>, const DeviceIndexPlan&, DeviceVector&) = 0;
    virtual void fill(DeviceVector&, double) = 0;
    virtual void copy(const DeviceVector&, DeviceVector&, std::size_t source_begin = 0,
                      std::size_t destination_begin = 0,
                      std::size_t count = std::numeric_limits<std::size_t>::max()) = 0;
    virtual void scale(double, DeviceVector&) = 0;
    virtual void axpy(double, const DeviceVector&, DeviceVector&) = 0;
    virtual double dot(const DeviceVector&, const DeviceVector&) = 0;
    virtual double norm(const DeviceVector&) = 0;
    virtual void apply(const DeviceVector&, DeviceVector&) = 0;
    virtual void precondition(const DeviceVector&, DeviceVector&) = 0;
    [[nodiscard]] virtual DeviceTransferMetrics metrics() const = 0;
};

// -1 uses the calling thread's current device; explicit ordinals are process-local.
// block_ilu0 factors the owned block during host setup and applies triangular
// factors on the device; jacobi uses only its diagonal. Neither transfers
// iterative vectors to the host for preconditioning.
std::unique_ptr<DeviceExecutor> make_cuda_executor(const DeviceCsrView&,int device_id=-1,
                                                 const std::string& preconditioner="jacobi");
struct ResidentLinearResult {
    LinearResult linear;
    DeviceTransferMetrics transfers;
};
// The matrix, Krylov basis, preconditioner and iterates remain device-resident.
// Small Hessenberg/Givens coefficients are host control data, not field vectors.
ResidentLinearResult solve_cuda_resident(const CsrMatrix&, std::span<const double>,
                                        const ExecutionOptions& = {});
} // namespace raikou::backend
