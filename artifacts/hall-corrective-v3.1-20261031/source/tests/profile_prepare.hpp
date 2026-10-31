#pragma once
// Phase-instrumented copy of Model::prepare (Hall v3 field_model.hpp, unchanged
// arithmetic in the same order) used only for profiling. The result of the
// instrumented copy is compared word for word with Model::prepare by the
// profiler, so the phase split refers to exactly the arithmetic of record.
#include <chrono>
#include <map>
#include <string>

namespace hall_corrective_v3 {
struct PhaseClock {
    std::map<std::string, double> ms;
    std::chrono::steady_clock::time_point t = std::chrono::steady_clock::now();
    void lap(const char* name) {
        const auto now = std::chrono::steady_clock::now();
        ms[name] += std::chrono::duration<double, std::milli>(now - t).count();
        t = now;
    }
};
// Same members as FieldCandidate (whose members are private to Model/Owner).
struct ProfiledCandidate {
    axial::State field_; ChargeLedger charge_; Midpoint midpoint_{{}, {}}; axial::Receipt receipt_; Constraints constraints_;
    std::string before_, configuration_, conductivity_; Vec current_; double dt_ = 0; SolveRecord solve_; double continuity_ratio_ = 0; bool bank_transfer_ = false;
    const axial::State& field() const noexcept { return field_; }
    const ChargeLedger& charge() const noexcept { return charge_; }
    const axial::Receipt& receipt() const noexcept { return receipt_; }
    const Midpoint& midpoint() const noexcept { return midpoint_; }
};
struct ModelTestAccess {
    // Same body as Model::prepare (v3) with lap() calls between phases.
    static ProfiledCandidate prepare(const Model& m, PhaseClock& T, const axial::State& old, const ChargeLedger& before, const std::string& before_sha256,
                                  const Vec& plasma, const FreeCharge& after, const DepositSupport& support, const Vec& sigma,
                                  const std::string& sigma_sha256, double dt, const BankTransfer* bank = nullptr) {
        const auto& g = m.mesh_;
        T.t = std::chrono::steady_clock::now();
        m.admit(old, before);
        T.lap("prepare.admit_old(Gauss+divB)");
        free_charge_guard(g, after);
        need(std::isfinite(dt) && dt > 0 && dt <= 1e-14 && old.time_s + dt > old.time_s && old.time_s + dt <= m.drive_.end() && old.revision < UINT64_MAX - 1, "authoritative field interval");
        need(plasma.size() == std::size_t(g.ne) && axial::finite(plasma), "field current shape/finiteness");
        need(support.absolute_C().size() == std::size_t(g.nn), "deposit support shape");
        ProfiledCandidate out;
        const bool transfer = bank && !bank->empty;
        out.bank_transfer_ = transfer;
        auto divergence = incidence(g, plasma);
        for (int n = 0; n < g.nn; ++n) {
            const auto defect = transfer ? after.net_C[n] + bank->conductor_C[n] + bank->exterior_C[n] - before.particles.net_C[n] - divergence[n] * dt
                                         : after.net_C[n] - before.particles.net_C[n] - divergence[n] * dt;
            const long double scale = before.particles.gross_C[n] + after.gross_C[n] + magnitude(divergence[n] * dt) + (transfer ? (long double)bank->gross_C[n] : 0.L);
            const long double tolerance = 1e-12L * scale + (long double)support.absolute_C()[n];
            const long double size = std::abs(defect.wide());
            need(size <= tolerance, "actual particle endpoint charge/current mismatch");
            if (size > 0) out.continuity_ratio_ = std::max(out.continuity_ratio_, double(size / tolerance));
        }
        T.lap("prepare.free_charge_guard+continuity");
        const double current = m.drive_.current(old.time_s + .5 * dt);
        auto S = g.conduction(sigma);
        Vec diagonal(g.ne, 1), sqrtM(g.ne);
        for (int e = 0; e < g.ne; ++e) { sqrtM[e] = std::sqrt(g.M[e]); if (!g.fixed_E[e]) diagonal[e] += .5 * dt * (S[e] + g.Y[e]) / g.M[e]; }
        for (int b = 0; b < g.nb; ++b) for (int k = 0; k < 4; ++k) { int e = g.C[b].index[k]; if (!g.fixed_E[e]) diagonal[e] += .25 * dt * dt * g.H[b] * std::pow(g.C[b].coefficient[k], 2) / g.M[e]; }
        T.lap("prepare.conduction+diagonal");
        auto apply = [&](const Vec& v) {
            Vec E(g.ne);
            for (int e = 0; e < g.ne; ++e) if (!g.fixed_E[e]) E[e] = v[e] / sqrtM[e];
            auto K = g.curl_adjoint(g.curl(E));
            Vec result(g.ne);
            for (int e = 0; e < g.ne; ++e) result[e] = g.fixed_E[e] ? v[e] : v[e] + (.5 * dt * (S[e] + g.Y[e]) * E[e] + .25 * dt * dt * K[e]) / sqrtM[e];
            return result;
        };
        const auto oldE = words(old.E, old.E_low), oldB = words(old.B, old.B_low);
        auto curl = m.curl_words(oldE);
        auto estimatedB = oldB;
        for (int b = 0; b < g.nb; ++b) estimatedB[b] = estimatedB[b] - Words(curl[b]) * (.5 * dt);
        auto ampere = m.ampere_words(estimatedB);
        Vec rhs(g.ne);
        for (int e = 0; e < g.ne; ++e)
            if (!g.fixed_E[e]) rhs[e] = ((Words(ampere[e]) - Words(plasma[e]) - Words(current) * m.initial_.feed[e] - oldE[e] * (S[e] + g.Y[e]) +
                                         Words(m.initial_.background_force[e])) * dt / sqrtM[e]).rounded();
        T.lap("prepare.rhs(curl_words+ampere_words)");
        auto& sr = out.solve_;
        sr.backend = m.solver_->name();
        const Apply closure = std::cref(apply);
        const auto increment = m.solver_->solve(closure, rhs, 200, SolveContext{g, dt, S, sqrtM, diagonal}, sr);
        T.lap("prepare.solve");
        need(increment.size() == rhs.size() && axial::finite(increment), "field solve backend returned a misshaped or nonfinite increment");
        {
            const auto Ax = apply(increment);
            Vec r(rhs.size());
            for (std::size_t i = 0; i < r.size(); ++i) r[i] = rhs[i] - Ax[i];
            sr.rhs_norm = foil::norm(rhs); sr.target = std::max(1e-30, 2e-13 * sr.rhs_norm); sr.independent_residual = foil::norm(r);
            need(sr.independent_residual <= sr.target, "independent FP64 field residual check failed");
        }
        T.lap("prepare.independent_residual");
        out.field_ = old; out.charge_ = before; out.charge_.particles = after; out.dt_ = dt;
        out.before_ = before_sha256; out.configuration_ = m.configuration_; out.conductivity_ = sigma_sha256; out.current_ = plasma;
        T.lap("prepare.copy_state");
        auto& r = out.receipt_;
        r.current_A = current; r.mid_E.resize(g.ne); r.solid_current.resize(g.ne); r.short_current.resize(g.ne);
        WordVec nextE(g.ne), nextB(g.nb);
        Vec midB(g.nb);
        for (int e = 0; e < g.ne; ++e) {
            nextE[e] = oldE[e] + Words(increment[e] / sqrtM[e]); out.field_.E[e] = nextE[e].hi; out.field_.E_low[e] = nextE[e].lo;
            r.mid_E[e] = ((oldE[e] + nextE[e]) * .5).rounded();
        }
        const auto cmid = g.curl(r.mid_E);
        for (int b = 0; b < g.nb; ++b) {
            nextB[b] = oldB[b] - Words(cmid[b]) * dt; out.field_.B[b] = nextB[b].hi; out.field_.B_low[b] = nextB[b].lo;
            midB[b] = ((oldB[b] + nextB[b]) * .5).rounded();
        }
        T.lap("prepare.update(E,B,midpoint)");
        out.midpoint_ = Midpoint(r.mid_E, midB);
        T.lap("prepare.midpoint_role_seal(hash)");
        ampere = g.curl_adjoint(midB);
        Vec metal_current(g.ne), exterior_current(g.ne), total(g.ne), displacement(g.ne);
        for (int e = 0; e < g.ne; ++e) {
            const auto delta = nextE[e] - oldE[e]; const auto true_mid = (nextE[e] + oldE[e]) * .5;
            displacement[e] = (delta * (g.M[e] / dt)).rounded(); r.solid_current[e] = S[e] * r.mid_E[e];
            exterior_current[e] = g.Y[e] * r.mid_E[e] - m.initial_.background_force[e];
            metal_current[e] = r.solid_current[e] + current * m.initial_.feed[e];
            total[e] = metal_current[e] + plasma[e] + exterior_current[e] + displacement[e]; const double residual = ampere[e] - total[e];
            if (g.fixed_E[e]) { r.short_current[e] = residual; metal_current[e] += residual; total[e] += residual; }
            else r.ampere_residual_A_m = std::max(r.ampere_residual_A_m, std::abs(residual));
            r.field_delta_J += double(g.M[e] * true_mid.wide() * delta.wide());
            r.plasma_work_J += dt * r.mid_E[e] * plasma[e]; r.joule_J += dt * r.mid_E[e] * r.solid_current[e];
            r.circuit_work_J -= dt * current * m.initial_.feed[e] * r.mid_E[e]; r.outer_loss_J += dt * g.Y[e] * r.mid_E[e] * r.mid_E[e];
            r.outer_supply_J += dt * m.initial_.background_force[e] * r.mid_E[e];
        }
        for (int b = 0; b < g.nb; ++b) r.field_delta_J += double(g.H[b] * ((nextB[b] + oldB[b]) * .5).wide() * (nextB[b] - oldB[b]).wide());
        T.lap("prepare.ampere+currents+energy_sums");
        const auto metal_delta = incidence(g, metal_current), outer_delta = incidence(g, exterior_current);
        if (transfer) {
            for (int n = 0; n < g.nn; ++n) {
                if (bank->conductor_C[n].hi != 0 || bank->conductor_C[n].lo != 0) out.charge_.conductor_C[n] = out.charge_.conductor_C[n] + bank->conductor_C[n];
                if (bank->exterior_C[n].hi != 0 || bank->exterior_C[n].lo != 0) out.charge_.exterior_C[n] = out.charge_.exterior_C[n] + bank->exterior_C[n];
            }
            bool any = false; for (double x : bank->residue_C) any = any || x != 0;
            if (any) {
                if (out.charge_.residue_C.empty()) out.charge_.residue_C.assign(g.nn, 0.);
                for (int n = 0; n < g.nn; ++n) if (bank->residue_C[n] != 0) {
                    const long double x = (long double)out.charge_.residue_C[n] + bank->residue_C[n];
                    double v = double(x); if ((long double)v < x) v = std::nextafter(v, INFINITY); out.charge_.residue_C[n] = v;
                }
            }
        }
        for (int n = 0; n < g.nn; ++n) {
            out.charge_.conductor_C[n] = out.charge_.conductor_C[n] + metal_delta[n] * dt;
            out.charge_.exterior_C[n] = out.charge_.exterior_C[n] + outer_delta[n] * dt;
        }
        T.lap("prepare.ledger(incidence words)");
        r.terminal_V = -axial::dot(m.initial_.feed, r.mid_E); r.cell_joule_W = g.cell_joule_W(sigma, r.mid_E);
        T.lap("prepare.cell_joule+terminal");
        r.balance_J = r.field_delta_J + r.plasma_work_J + r.joule_J + r.outer_loss_J - r.outer_supply_J - r.circuit_work_J;
        r.kcl_residual_A = axial::maximum(g.gradient_transpose(total));
        T.lap("prepare.kcl");
        const auto Q0 = field_charge(g, old), Q1 = field_charge(g, out.field_);
        double activityQ = 0;
        for (int n = 0; n < g.nn; ++n) {
            const auto change = Q1[n] - Q0[n]; const auto expected = (divergence[n] + metal_delta[n] + outer_delta[n]) * dt;
            r.gauss_increment_C = std::max(r.gauss_increment_C, double(std::abs((change - expected).wide()))); activityQ = std::max(activityQ, double(magnitude(change)));
        }
        need(r.gauss_increment_C <= 1e-22 + 2e-10 * activityQ, "existing field Gauss increment bound");
        T.lap("prepare.gauss_increment");
        const int cut = g.kshort + g.nz_foil / 2 - 1; r.region_conduction_A = g.upper_region_inflow(r.solid_current, cut);
        r.region_plasma_A = g.upper_region_inflow(plasma, cut); r.region_displacement_A = g.upper_region_inflow(displacement, cut);
        r.region_boundary_A = g.upper_region_inflow(exterior_current, cut);
        T.lap("prepare.region_inflows");
        const double activity = std::abs(r.field_delta_J) + std::abs(r.plasma_work_J) + std::abs(r.joule_J) + std::abs(r.outer_loss_J) + std::abs(r.outer_supply_J) + std::abs(r.circuit_work_J);
        need(std::isfinite(r.balance_J) && std::abs(r.balance_J) <= 1e-20 + 1e-10 * activity, "unchanged field/current/Joule energy bound");
        out.field_.time_s = old.time_s + dt; ++out.field_.revision;
        out.constraints_ = m.admit(out.field_, out.charge_);
        T.lap("prepare.admit_candidate(Gauss+divB)");
        return out;
    }
};
}  // namespace hall_corrective_v3
