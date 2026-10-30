#pragma once
#include <stdexcept>
namespace raikou_pic_parallel {
class CudssPoisson {
public:
    explicit CudssPoisson(const mfem::SparseMatrix&) {
        throw std::runtime_error("GPU execution is disabled in this validation build");
    }
    void Solve(const mfem::Vector&,mfem::Vector&) {
        throw std::runtime_error("GPU execution is disabled in this validation build");
    }
};
}
