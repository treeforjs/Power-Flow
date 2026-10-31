#pragma once
// Source owner v6.1: exact fast evaluation of the sealed thermal-field emission kernel.
//
//   electrode_heat_v6_1::emission::evaluate(T, F, config)
//
// returns a surface_emission::Distribution that is bit-identical, field by field and node
// by node, to the sealed surface_emission::evaluate(T, F, config) (source kernel of the Hall
// heat-v2 binding, compiled from its sealed source). It never approximates a value:
//
//  1. Retarding or zero field (F <= 0, T > 0): the sealed evaluate uses F only to store
//     extraction_field_V_m (no barrier, no Schottky lowering, quadrature on [phi, phi +
//     tail_widths kT]). The first call for a (T bits, config) runs the sealed evaluate; later
//     calls copy that Distribution and store their own F.
//  2. Extracting field (F >= 1 V/m, T > 0): the sealed quadrature is re-evaluated with the
//     same expressions in the same order. Two changes do not alter any bit:
//     a. The Gauss-Legendre rules and the action-node tables (cos, sin^2, weights and their
//        products) are computed once per order by the same code.
//     b. Energy nodes whose quadrature weight exp(l - largest) is exactly 0 in binary64 are
//        not evaluated. A node is dropped only when a rigorous upper bound U on its log
//        weight l satisfies U < M - 800, M being the exact log weight of an evaluated node
//        (so l - largest < -800 and exp() underflows to +0 below -745.14). The bound uses
//        the kernel's own monotonicity: every action term w s^2 / sqrt(mid + half cos) is at
//        least w s^2 / sqrt(max(mid + half, mid - half)) (correctly rounded operations are
//        monotone), the long-double sum of nonnegative terms is monotone, softplus(G) >= G,
//        and the supply and panel-width logs are bounded per quadrature cell. Dropped nodes
//        add +0 to the long-double normalization sum and have p == 0 in the moment loop
//        (they are not stored in the sealed node list either). The sealed guard on a dropped
//        node ("Nonfinite Schottky-Nordheim action") is shown to pass by a finite upper bound
//        of its action; when it cannot be shown, the node is evaluated exactly.
//     c. The remaining action sums run four energy nodes at a time: the binary64 terms are
//        formed with vector sqrt/div/mul/add (IEEE correctly rounded, no contraction) and
//        each node's long-double sum is accumulated in the sealed term order.
//  3. T == 0 runs the sealed evaluate.
// A self-check per configuration compares the fast result with the sealed evaluate on
// reference inputs before the fast path is used; on any difference the configuration is
// served by the sealed evaluate only (statistics().sealed_fallback_configs).
#include "surface_emission.hpp"
#include <cstdint>

namespace electrode_heat_v6_1::emission {
struct Statistics {
  std::uint64_t calls = 0, sealed_calls = 0, retarding_reused = 0, retarding_sealed = 0, extracting_fast = 0;
  std::uint64_t nodes = 0, nodes_exact = 0, nodes_bounded_zero = 0, self_checks = 0, sealed_fallback_configs = 0;
};
surface_emission::Distribution evaluate(double T_K, double signed_extraction_field_V_m, const surface_emission::Config& config);
// The sealed evaluate on the same inputs, for comparisons.
surface_emission::Distribution evaluate_sealed(double T_K, double signed_extraction_field_V_m, const surface_emission::Config& config);
bool same_distribution(const surface_emission::Distribution& a, const surface_emission::Distribution& b) noexcept;
Statistics statistics() noexcept;
void reset_statistics() noexcept;
// Test switch: false makes evaluate() call the sealed evaluate for every input.
void set_fast_path(bool enabled) noexcept;
}  // namespace electrode_heat_v6_1::emission
