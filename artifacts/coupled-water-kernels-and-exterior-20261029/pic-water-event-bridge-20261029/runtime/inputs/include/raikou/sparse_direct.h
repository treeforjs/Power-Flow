#pragma once
#include "raikou/backend.h"
#include <span>
#include <string>
#include <vector>

namespace raikou::backend {

// Whether this build carries the sparse direct solver (Eigen's SparseLU, found at configure time).
[[nodiscard]] bool sparse_direct_available();

struct DirectResult {
    std::vector<double> solution;
    bool solved{false};
    std::string message;
};
// One serial LU factorization and solve with a column approximate-minimum-degree ordering.
// A singular or non-finite system returns solved=false with the factorization's message.
[[nodiscard]] DirectResult solve_sparse_direct(const CsrMatrix& matrix, std::span<const double> rhs);

} // namespace raikou::backend
