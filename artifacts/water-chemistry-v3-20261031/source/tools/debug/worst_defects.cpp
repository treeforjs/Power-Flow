// Prints the events with the largest relative closure defect (development aid).
#include "chem_host.hpp"
#include <iostream>
#include <random>
using namespace pfc;
int main(int, char** argv) {
    const Blob b = load_blob(argv[1]);
    StepConfig cfg{0x5EED, 1, 5e-11, 7.0, PFC_PURPOSE_MCC, 4096, 64, 0};
    std::mt19937_64 g(1);
    std::uniform_real_distribution<double> U(0, 1);
    std::vector<CellTargets> cells(1);
    CellTargets& ct = cells[0];
    ct = CellTargets{};
    pfc_u32 targets[3] = {species_index(b, "H2O"), species_index(b, "H3O+"), 0xFFFFFFFFu};
    ct.n_targets = 2;
    for (int k = 0; k < 2; ++k) { ct.species[k] = targets[k]; ct.n[k] = k ? 1e19 : 1e23; ct.s[k] = 400; }
    for (int k = 2; k < PFC_MAX_TARGETS; ++k) ct.species[k] = 0xFFFFFFFFu;
    std::vector<CarrierIn> cs;
    for (int i = 0; i < 200000; ++i) {
        CarrierIn c{}; c.id = i + 1; c.cell = 0; c.weight = 1e6;
        c.species = (i % 2) ? 0 : species_index(b, "H2O+");
        const double m = b.view.species[c.species].mass_eV;
        const double E = 1e-3 * std::exp(std::log((c.species == 0 ? 1e6 : 1e3) / 1e-3) * U(g));
        const V3 u = u_from_p(m, scale(iso_dir(U(g), U(g)), std::sqrt(E * (E + 2 * m))));
        c.u[0] = u.x; c.u[1] = u.y; c.u[2] = u.z; cs.push_back(c);
    }
    StepResult r; cpu_step(b, cfg, cs, cells, r);
    std::vector<std::pair<double, const EventRec*>> w;
    for (std::size_t i = 0; i < cs.size(); ++i) for (pfc_u32 k = 0; k < r.outs[i].n_events; ++k) {
        const EventRec& e = r.events[i * 64 + k];
        const double d = event_defect(e);
        w.push_back({std::fabs(d) / event_tolerance(e) * 0x1p-40, &e});
    }
    std::sort(w.begin(), w.end(), [](auto& a, auto& b) { return a.first > b.first; });
    std::map<std::string, std::pair<int,double>> byp;
    for (auto& [rr, e] : w) { auto& x = byp[b.process_names[e->process]]; x.first++; x.second = std::max(x.second, rr); }
    for (auto& [n, x] : byp) std::cout << n << " events " << x.first << " max rel " << x.second << "\n";
    for (int i = 0; i < 6 && i < (int)w.size(); ++i) {
        const EventRec& e = *w[i].second;
        std::cout.precision(17);
        std::cout << "rel " << w[i].first << " " << b.process_names[e.process] << " E_tf " << e.E_tf << " ke_in " << e.ke_in_p << " + " << e.ke_in_t << " q " << e.q_eV << " int " << e.e_internal << " defect " << event_defect(e) << " outs";
        for (pfc_u32 j = 0; j < e.n_out; ++j) std::cout << " " << b.species_names[e.out_species[j]] << ":" << e.ke_out[j];
        std::cout << "\n";
    }
}
