// harness.hpp -- chemistry-driven qualification harness for gas target owner v1.
//
// A population of carriers (one carrier per packet) on the native-v8
// collision grid runs complete chemistry v3 stages with the C3 kernel of
// chemistry v3.1.1 (headers and blob read in place, pinned) and native v8's
// resolution semantics (API-008 3.3-3.5): CellTargets from exact moments of
// the target carriers, prefix-rule debits of target-consuming events,
// survivor/consumption, births, reservoir deposits and the conservative
// parcel release of v3.1.1, debit allocation in ascending native ID with
// exact remainders. It produces exactly the hook inputs native v8.x produces
// (INTEGRATION.txt 1): recoil records, debit and release flows, and calls the
// owner's feedback before the kernel. Between the stages carriers drift
// ballistically (zero field); leaving the box is a sink.
#pragma once
#include "chem_host.hpp"
#include "pf_gas_v1/gas_owner.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <map>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace gh {
using pf_gas_v1::Exact;
using pfg::u32;
using pfg::u64;

static_assert(sizeof(pfg::ChemEvent) == sizeof(pfc::EventRec), "ChemEvent layout");
static_assert(offsetof(pfg::ChemEvent, out_u) == offsetof(pfc::EventRec, out_u), "ChemEvent layout");
static_assert(offsetof(pfg::ChemEvent, e_internal) == offsetof(pfc::EventRec, e_internal), "ChemEvent layout");
static_assert(offsetof(pfg::ChemEvent, out_dest) == offsetof(pfc::EventRec, out_dest), "ChemEvent layout");
static_assert(sizeof(pfg::CellTargetsAbi) == sizeof(pfc::CellTargets), "CellTargets layout");
static_assert(offsetof(pfg::CellTargetsAbi, s) == offsetof(pfc::CellTargets, s), "CellTargets layout");
static_assert(offsetof(pfg::CellTargetsAbi, n_targets) == offsetof(pfc::CellTargets, n_targets), "CellTargets layout");

inline void need(bool ok, const std::string &m) {
  if (!ok)
    throw std::runtime_error("harness: " + m);
}
inline bool negative(const Exact &x) { return x.normalised()[pfg::kDigits - 1] < 0; }

// Internal-energy carrier product of a process kind (chem_event.h
// event_kinematics, v3.1 kinds 1-6 and v3.2 kinds 7-8).
inline u32 internal_carrier(u32 kind, u32 n_out) {
  switch (kind) {
  case 2: return 1;                                  // EXCITE: target
  case 3: return n_out - 1;                          // DISSOC: last fragment
  case 4: return 2;                                  // IONIZE: ion
  case 5: return 0;                                  // CAPTURE: first product
  case 6: return n_out > 2 ? n_out - 1 : pfg::kNone; // CHARGE_X (v3.2 fragments)
  case 7: return 2;                                  // ION_IONIZE (v3.2)
  case 8: return n_out >= 2 ? n_out - 2 : pfg::kNone;// BREAKUP (v3.2)
  default: return pfg::kNone;                        // ELASTIC
  }
}

// Owner configuration from the chemistry blob (native registry conventions:
// state "pfc-v3/<name>", mass_kg = (mass_eV * e) / (c * c)).
inline pf_gas_v1::Configuration owner_configuration(const pfc::Blob &b, const pf_gas_v1::Grid &g,
                                                    const std::vector<pf_gas_v1::Species> &extra = {}) {
  pf_gas_v1::Configuration c;
  const double e = 1.602176634e-19, cc = 299792458.0, c2 = cc * cc;
  std::vector<u32> target(b.n_species, 0);
  for (u32 k = 0; k < b.n_groups; ++k)
    target[b.view.groups[k].target] = 1;
  for (u32 s = 0; s < b.n_species; ++s) {
    pf_gas_v1::Species x;
    x.state = "pfc-v3/" + b.species_names[s];
    x.mass_kg = (b.view.species[s].mass_eV * e) / c2;
    x.charge = b.view.species[s].charge;
    x.chemistry_target = target[s] != 0;
    c.species.push_back(x);
  }
  for (const auto &x : extra)
    c.species.push_back(x);
  for (u32 p = 0; p < b.n_proc; ++p) {
    const auto &pr = b.view.procs[p];
    c.processes.push_back({b.view.groups[pr.group].projectile, pr.n_out, internal_carrier(pr.kind, pr.n_out), pr.kind});
  }
  c.chemistry_blob_sha256 = b.sha256;
  c.chemistry_registry_sha256 = "harness-registry/" + b.sha256;
  c.grid = g;
  return c;
}

inline pf_gas_v1::Grid test_grid(u32 nx = 4, u32 ny = 4, u32 nz = 2) {
  pf_gas_v1::Grid g;
  g.low = {0.0, 0.0, 0.0};
  g.high = {1e-3, 1e-3, 2e-3};
  g.n = {nx, ny, nz};
  g.volume = ((g.high[0] - g.low[0]) / nx) * ((g.high[1] - g.low[1]) / ny) * ((g.high[2] - g.low[2]) / nz);
  return g;
}

struct Carrier {
  u64 id;
  u32 species;
  double w;
  double x[3], u[3];
};

struct StageResult {
  std::vector<pfg::RecoilRecord> records;
  std::vector<pfg::FlowRecord> flows;
  std::vector<pfg::CellTargetsAbi> cells_after_feedback;
  std::vector<u32> atlas;
  std::string records_sha256, cells_sha256;
  u64 events = 0, refused_debits = 0, births = 0, releases = 0, consumed = 0;
  pf_gas_v1::FeedbackReceipt feedback;
};

class World {
public:
  const pfc::Blob *blob = nullptr;
  pf_gas_v1::Configuration cfg;
  std::vector<Carrier> pop;   // ascending id
  u64 next_id = 1;
  std::map<pfc::InventoryKey, pfc::Reservoir> res;
  u64 seed = 0x5EED5EED20261101ull;
  double w0 = 1e9;
  double trunc = 7.0;
  u32 max_events = 64, max_candidates = 4096;
  std::vector<u32> projectile, target;

  World(const pfc::Blob &b, const pf_gas_v1::Configuration &c) : blob(&b), cfg(c) {
    projectile.assign(b.n_species, 0);
    target.assign(b.n_species, 0);
    for (u32 k = 0; k < b.n_groups; ++k) {
      projectile[b.view.groups[k].projectile] = 1;
      target[b.view.groups[k].target] = 1;
    }
  }
  double mass(u32 s) const { return cfg.species[s].mass_kg; }
  void add(u32 species, double w, const double *x, const double *u) {
    Carrier c;
    c.id = next_id++;
    c.species = species;
    c.w = w;
    std::memcpy(c.x, x, 24);
    std::memcpy(c.u, u, 24);
    pop.push_back(c);
  }
  u32 species_index(const std::string &state) const {
    for (u32 s = 0; s < cfg.species.size(); ++s)
      if (cfg.species[s].state == state)
        return s;
    throw std::runtime_error("harness: no species " + state);
  }
  // exact population sums (all species): K (J), P (kg m/s)
  void totals(Exact &K, Exact P[3]) const {
    for (const auto &c : pop) {
      K.add_product(c.w, pfg::kinetic(mass(c.species), c.u));
      for (int a = 0; a < 3; ++a)
        P[a].add_triple(c.w, mass(c.species), c.u[a]);
    }
  }
  std::vector<pf_gas_v1::CensusRow> rows() const {
    std::vector<pf_gas_v1::CensusRow> r;
    r.reserve(pop.size());
    for (const auto &c : pop) {
      pf_gas_v1::CensusRow x;
      x.native_id = c.id;
      x.species = c.species;
      std::memcpy(x.position_m, c.x, 24);
      std::memcpy(x.proper_velocity_m_s, c.u, 24);
      x.weight = c.w;
      r.push_back(x);
    }
    return r;
  }
  std::string population_sha256() const {
    pfc::Sha256 h;
    h.text("gas-harness/population");
    for (const auto &c : pop) {
      h.pod(c.id);
      h.pod(c.species);
      h.pod(c.w);
      h.update(c.x, 24);
      h.update(c.u, 24);
    }
    return h.hex();
  }

  // One chemistry stage (API-008 3.3-3.5). owner may be null (field-target
  // reference without an owner); feedback is called when an owner is given.
  StageResult stage(pf_gas_v1::Owner *owner, u64 step, u32 stage_no, double dt_stage, bool take = true) {
    const auto &b = *blob;
    StageResult out;
    const u32 none = 0xFFFFFFFFu;
    // 1. cells
    std::vector<u32> row_cell(pop.size(), none);
    std::vector<u32> occupied;
    for (std::size_t i = 0; i < pop.size(); ++i) {
      const u32 c = cfg.grid.locate(pop[i].x);
      need(c != none, "carrier outside the grid");
      row_cell[i] = c;
      if (pop[i].species < b.n_species && projectile[pop[i].species])
        occupied.push_back(c);
    }
    std::sort(occupied.begin(), occupied.end());
    occupied.erase(std::unique(occupied.begin(), occupied.end()), occupied.end());
    auto compact = [&](u32 atlas) -> u32 {
      auto it = std::lower_bound(occupied.begin(), occupied.end(), atlas);
      return (it != occupied.end() && *it == atlas) ? u32(it - occupied.begin()) : none;
    };
    struct T {
      u32 cell, species, row;
    };
    std::vector<T> targets;
    for (std::size_t i = 0; i < pop.size(); ++i) {
      const u32 sp = pop[i].species;
      if (sp >= b.n_species || !target[sp])
        continue;
      const u32 c = compact(row_cell[i]);
      if (c != none)
        targets.push_back({c, sp, u32(i)});
    }
    std::sort(targets.begin(), targets.end(), [](const T &a, const T &c) {
      return a.cell != c.cell ? a.cell < c.cell : (a.species != c.species ? a.species < c.species : a.row < c.row);
    });
    std::vector<pfc::CellTargets> cells(occupied.size());
    for (auto &c : cells)
      std::memset(&c, 0, sizeof c);
    std::map<std::pair<u32, u32>, Exact> inventory;
    for (std::size_t i = 0; i < targets.size();) {
      std::size_t j = i;
      while (j < targets.size() && targets[j].cell == targets[i].cell && targets[j].species == targets[i].species)
        ++j;
      Exact sw, su[3], s2;
      for (std::size_t k = i; k < j; ++k) {
        const auto &r = pop[targets[k].row];
        sw.add_product(r.w, 1.0);
        for (int a = 0; a < 3; ++a)
          su[a].add_product(r.w, r.u[a]);
      }
      const double W = sw.value();
      double drift[3];
      for (int a = 0; a < 3; ++a)
        drift[a] = su[a].value() / W;
      for (std::size_t k = i; k < j; ++k) {
        const auto &r = pop[targets[k].row];
        const double t0 = r.u[0] - drift[0], t1 = r.u[1] - drift[1], t2 = r.u[2] - drift[2];
        s2.add_product(r.w, t0 * t0 + t1 * t1 + t2 * t2);
      }
      auto &cell = cells[targets[i].cell];
      need(cell.n_targets < PFC_MAX_TARGETS, "more than 12 targets");
      const u32 slot = cell.n_targets++;
      cell.species[slot] = targets[i].species;
      cell.n[slot] = W / cfg.grid.volume;
      for (int a = 0; a < 3; ++a)
        cell.drift[slot][a] = drift[a];
      cell.s[slot] = std::sqrt((s2.value() / W) / 3.0);
      inventory[{targets[i].cell, targets[i].species}] = sw;
      i = j;
    }
    out.atlas = occupied;
    if (owner)
      out.feedback = owner->feedback(step, stage_no, occupied,
                                     std::span<pfg::CellTargetsAbi>(reinterpret_cast<pfg::CellTargetsAbi *>(cells.data()),
                                                                    cells.size()));
    out.cells_after_feedback.assign(reinterpret_cast<const pfg::CellTargetsAbi *>(cells.data()),
                                    reinterpret_cast<const pfg::CellTargetsAbi *>(cells.data()) + cells.size());
    {
      pfc::Sha256 h;
      h.text("gas-harness/cells");
      h.update(cells.data(), cells.size() * sizeof(pfc::CellTargets));
      h.update(occupied.data(), occupied.size() * 4);
      out.cells_sha256 = h.hex();
    }
    // 2. kernel (C3 CPU reference), carriers ascending ID (= population order)
    pfc::StepConfig sc;
    std::memset(&sc, 0, sizeof sc);
    sc.seed = seed;
    sc.step = step;
    sc.dt_s = dt_stage;
    sc.trunc_sigmas = trunc;
    sc.purpose = stage_no == 1 ? PFC_PURPOSE_MCC : PFC_PURPOSE_MCC_SECOND;
    sc.max_candidates = max_candidates;
    sc.max_events = max_events;
    struct Active {
      std::size_t row;
      u32 cell;
      pfc::CarrierOut o;
      std::vector<pfc::EventRec> ev;
    };
    std::vector<Active> act;
    std::vector<pfc::EventRec> scratch(max_events);
    pfc::Sha256 rh;
    rh.text("gas-harness/records");
    for (std::size_t i = 0; i < pop.size(); ++i) {
      const auto &c = pop[i];
      if (c.species >= b.n_species || !projectile[c.species])
        continue;
      pfc::CarrierIn in;
      std::memset(&in, 0, sizeof in);
      in.id = c.id;
      for (int a = 0; a < 3; ++a)
        in.u[a] = c.u[a];
      in.weight = c.w;
      in.species = c.species;
      in.cell = compact(row_cell[i]);
      pfc::CarrierOut o;
      std::memset(&o, 0, sizeof o);
      std::memset(scratch.data(), 0, scratch.size() * sizeof(pfc::EventRec));
      pfc::mcc_carrier(b.view, b.cover, sc, in, cells[in.cell], scratch.data(), &o);
      need(o.status == PFC_ST_OK || o.status == PFC_ST_CONSUMED, "carrier status " + std::to_string(o.status));
      rh.pod(o);
      for (u32 k = 0; k < o.n_events; ++k)
        rh.pod(scratch[k]);
      if (o.n_events)
        act.push_back({i, in.cell, o, std::vector<pfc::EventRec>(scratch.begin(), scratch.begin() + o.n_events)});
    }
    out.records_sha256 = rh.hex();
    // 3. debits (order: cell, carrier ID; prefix rule; target-consuming events)
    auto consumes = [](const pfc::EventRec &e) {
      for (u32 p = 0; p < e.n_out; ++p)
        if (e.out_dest[p] == PFC_DEST_GAS && e.out_species[p] == e.target_species)
          return false;
      return true;
    };
    std::vector<std::size_t> order(act.size());
    for (std::size_t i = 0; i < order.size(); ++i)
      order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t c) {
      return act[a].cell != act[c].cell ? act[a].cell < act[c].cell : pop[act[a].row].id < pop[act[c].row].id;
    });
    std::map<std::pair<u32, u32>, Exact> left, debited;
    std::vector<u32> admitted(act.size(), 0);
    for (auto i : order) {
      u32 k = 0;
      for (; k < act[i].ev.size(); ++k) {
        const auto &e = act[i].ev[k];
        if (!target[e.target_species] || !consumes(e))
          continue;
        const std::pair<u32, u32> key{e.cell, e.target_species};
        auto it = left.find(key);
        if (it == left.end())
          it = left.emplace(key, inventory.count(key) ? inventory[key] : Exact()).first;
        Exact rem = it->second;
        rem.add_product(e.weight, 1.0, 1);
        if (negative(rem)) {
          ++out.refused_debits;
          break;
        }
        it->second = rem;
        debited[key].add_product(e.weight, 1.0);
      }
      admitted[i] = k;
    }
    // 4. outcomes (ascending carrier ID = act order)
    std::set<std::size_t> consumed_rows;
    std::vector<Carrier> births;
    for (std::size_t i = 0; i < act.size(); ++i) {
      const u32 k = admitted[i];
      if (!k)
        continue;
      Carrier &c = pop[act[i].row];
      const u32 atlas = occupied[act[i].cell];
      const auto &last = act[i].ev[k - 1];
      int surv = -1;
      for (u32 p = 0; p < last.n_out; ++p)
        if (last.out_dest[p] == PFC_DEST_SURVIVOR)
          surv = int(p);
      for (u32 j = 0; j < k; ++j) {
        const auto &e = act[i].ev[j];
        pfg::RecoilRecord r;
        std::memset(&r, 0, sizeof r);
        r.accepted_step = step;
        r.stage = stage_no;
        r.grid_cell = atlas;
        std::memcpy(&r.ev, &e, sizeof e);
        out.records.push_back(r);
        ++out.events;
        for (u32 p = 0; p < e.n_out; ++p)
          if (e.out_dest[p] == PFC_DEST_BIRTH) {
            Carrier nb;
            nb.id = 0;
            nb.species = e.out_species[p];
            nb.w = e.weight;
            std::memcpy(nb.x, c.x, 24);
            for (int a = 0; a < 3; ++a)
              nb.u[a] = e.out_u[p][a];
            births.push_back(nb);
          }
        pfc::EventRec em = e;
        em.cell = atlas;
        pfc::deposit_reservoir(b, res, em, c.x);
      }
      if (surv >= 0)
        for (int a = 0; a < 3; ++a)
          c.u[a] = last.out_u[surv][a];
      else {
        consumed_rows.insert(act[i].row);
        ++out.consumed;
      }
    }
    // 5. debit allocation: ascending native ID among the cell's target carriers
    std::vector<Carrier> remainders;
    std::set<std::size_t> removed = consumed_rows;
    for (auto &[key, d] : debited) {
      Exact left_d = d;
      for (const auto &t : targets) {
        if (t.cell != key.first || t.species != key.second)
          continue;
        if (left_d.zero())
          break;
        if (consumed_rows.count(t.row))
          continue;
        Carrier &c = pop[t.row];
        const u32 atlas = occupied[key.first];
        Exact diff = left_d;
        diff.add_product(c.w, 1.0, 1);
        pfg::FlowRecord f{};
        f.kind = pfg::kFlowDebit;
        f.negative = 1;
        f.grid_cell = atlas;
        f.species = c.species;
        f.native_id = c.id;
        f.weight = c.w;
        std::memcpy(f.u, c.u, 24);
        out.flows.push_back(f);
        removed.insert(t.row);
        if (!negative(diff)) {
          left_d = diff;
        } else {
          // remainder R = c - left_d as binary64 pieces (first keeps the ID)
          Exact R;
          R.add_product(c.w, 1.0);
          R.merge(left_d, 1);
          left_d = Exact();
          bool first = true;
          while (!R.zero()) {
            double r = R.value();
            Exact piece;
            piece.add_product(r, 1.0);
            Exact rest = R;
            rest.merge(piece, 1);
            if (negative(rest)) {
              r = std::nextafter(r, 0.0);
              piece = Exact();
              piece.add_product(r, 1.0);
              rest = R;
              rest.merge(piece, 1);
            }
            need(r > 0 && !negative(rest), "remainder split");
            Carrier rc = c;
            rc.w = r;
            rc.id = first ? c.id : 0;
            first = false;
            remainders.push_back(rc);
            pfg::FlowRecord g = f;
            g.negative = 0;
            g.weight = r;
            g.native_id = rc.id;
            out.flows.push_back(g);
            R = rest;
          }
        }
      }
      need(left_d.zero(), "debit conflict");
    }
    // 6. apply: remove consumed/debited rows, append births, remainders
    std::vector<Carrier> next;
    next.reserve(pop.size() + births.size() + remainders.size());
    for (std::size_t i = 0; i < pop.size(); ++i)
      if (!removed.count(i))
        next.push_back(pop[i]);
    for (auto &nb : births) {
      nb.id = next_id++;
      next.push_back(nb);
    }
    out.births = births.size();
    for (auto &rc : remainders) {
      if (rc.id == 0)
        rc.id = next_id++;
      next.push_back(rc);
    }
    // 7. stage 2: conservative parcel release (v3.1.1)
    if (stage_no == 2) {
      const auto rel = pfc::release_reservoirs(b, res, w0, seed, step);
      for (const auto &p : rel) {
        Carrier rc;
        rc.id = next_id++;
        rc.species = p.species;
        rc.w = p.weight;
        std::memcpy(rc.x, p.x, 24);
        std::memcpy(rc.u, p.u, 24);
        next.push_back(rc);
        pfg::FlowRecord f{};
        f.kind = pfg::kFlowRelease;
        f.negative = 0;
        f.grid_cell = p.cell;
        f.species = p.species;
        f.native_id = rc.id;
        f.weight = p.weight;
        std::memcpy(f.u, p.u, 24);
        out.flows.push_back(f);
      }
      out.releases = rel.size();
    }
    std::sort(next.begin(), next.end(), [](const Carrier &a, const Carrier &c) { return a.id < c.id; });
    pop = std::move(next);
    if (owner && take)
      owner->take(step, stage_no, out.records, out.flows);
    return out;
  }

  // Ballistic drift for dt (zero field); carriers leaving the box are sunk.
  // Returns the sink flows (all species) and the sunk IDs.
  std::vector<pfg::FlowRecord> drift(double dt, std::vector<u64> &sunk) {
    std::vector<pfg::FlowRecord> flows;
    std::vector<Carrier> next;
    next.reserve(pop.size());
    const double c2 = pf7::c_squared;
    for (auto c : pop) {
      const double n = c.u[0] * c.u[0] + c.u[1] * c.u[1] + c.u[2] * c.u[2];
      const double g = std::sqrt(1 + n / c2);
      const u32 from = cfg.grid.locate(c.x);
      double x[3];
      for (int a = 0; a < 3; ++a)
        x[a] = c.x[a] + (c.u[a] / g) * dt;
      if (cfg.grid.locate(x) == 0xFFFFFFFFu) {
        pfg::FlowRecord f{};
        f.kind = pfg::kFlowSink;
        f.negative = 1;
        f.grid_cell = from;
        f.species = c.species;
        f.native_id = c.id;
        f.weight = c.w;
        std::memcpy(f.u, c.u, 24);
        flows.push_back(f);
        sunk.push_back(c.id);
        continue;
      }
      std::memcpy(c.x, x, 24);
      next.push_back(c);
    }
    pop = std::move(next);
    return flows;
  }
};

// Thermal (Maxwellian) proper velocity components at temperature T for mass m.
inline void maxwell(std::mt19937_64 &rng, double m, double T, double *u) {
  std::normal_distribution<double> N(0.0, std::sqrt(1.380649e-23 * T / m));
  for (int a = 0; a < 3; ++a)
    u[a] = N(rng);
}
// Electron/ion with kinetic energy E (eV), isotropic direction.
inline void energetic(std::mt19937_64 &rng, double m_kg, double E_eV, double *u) {
  std::uniform_real_distribution<double> U(0, 1);
  const double c = 299792458.0, mc2 = m_kg * c * c, K = E_eV * 1.602176634e-19;
  const double gamma = 1 + K / mc2;
  const double un = c * std::sqrt(gamma * gamma - 1);
  const double ct = 2 * U(rng) - 1, ph = 6.283185307179586 * U(rng), st = std::sqrt(1 - ct * ct);
  u[0] = un * st * std::cos(ph);
  u[1] = un * st * std::sin(ph);
  u[2] = un * ct;
}

}  // namespace gh
