// gas_core.h -- gas target owner v1, portable core.
//
// Compiled unchanged by the host C++ compiler (C++20, -ffp-contract=off,
// -fno-fast-math) and by NVRTC (--fmad=false --prec-div=true
// --prec-sqrt=true --ftz=false). Exact quantities are sums of exact products
// of binary64 values held in the native v7 forms superaccumulator
// (forms_core.h: 134 radix-2^32 digits, offset 2^-2148), so every ledger is
// independent of the order and partition of its terms. The only binary64
// roundings are the per-term operands (kinetic energy of one carrier, the eV
// to J conversion of one Q or internal energy, the internal-energy momentum
// of one product, a transport share) and the sampling representation of the
// feedback; each is a fixed sequence of IEEE + - * / sqrt.
//
// Definitions: README.txt section 2 (state), 3 (records), 4 (feedback),
// 5 (flows and transport).
#pragma once
#include "forms_core.h"

namespace pfg {

typedef pf7_u64 u64;
typedef pf7_i64 i64;
typedef pf7_u32 u32;
#define PFG_HD PF7_HD

constexpr int kDigits = pf7::acc_digits;     // 134
constexpr double kEv = 1.602176634e-19;      // J per eV (exact SI value, binary64)
constexpr double kC = 299792458.0;
constexpr u32 kNone = 0xFFFFFFFFu;

// chemistry v3 product destinations (chem_tables.h)
constexpr u32 kDestSurvivor = 1, kDestBirth = 2, kDestReservoir = 3, kDestGas = 4;
constexpr u32 kMaxOut = 6;
constexpr u32 kMaxTargets = 12;
constexpr u32 kSinks = 6;

// ------------------------------------------------------------------ records
// Chemistry v3 event record: byte layout of pfc::EventRec (chem_event.h,
// 384 bytes; asserted against the chemistry header in the adapters).
struct ChemEvent {
  u64 carrier;
  u32 candidate, process, cell, target_species;
  u32 n_out, cover_class, internal_sink, flags;
  double weight;
  double time_s;
  double u_in[3];
  double ut[3];
  double E_tf;
  double ke_in_p, ke_in_t;
  double q_eV, e_internal;
  u32 out_species[kMaxOut];
  u32 out_dest[kMaxOut];
  double out_u[kMaxOut][3];
  double ke_out[kMaxOut];
};
// Target recoil record (API-008 3.8 amendment, INTEGRATION.txt 1.1): one per
// admitted chemistry event, canonical order (projectile native ID, candidate).
struct RecoilRecord {
  u64 accepted_step;
  u32 stage;       // 1 or 2
  u32 grid_cell;   // native-v8/collision-grid-v1 cell ID of the event (atlas ID)
  ChemEvent ev;
};
// Population flow of an owner species (INTEGRATION.txt 1.2). negative = 1:
// the carrier leaves the population (debited packet before the debit, sink);
// 0: it enters (debit remainder, birth, release, transaction output).
// internal_J: internal energy carried by the macro (J, total, not per
// particle; native v8.1 FullTimeBirth::internal_energy_J convention).
constexpr u32 kFlowDebit = 0, kFlowRelease = 1, kFlowBirth = 2, kFlowSink = 3, kFlowTransaction = 4;
constexpr u32 kFlowKinds = 5;
struct FlowRecord {
  u32 kind, negative, grid_cell, species;
  u64 native_id;
  double weight;
  double u[3];
  double internal_J;
};
// Transport of the deposit shares (carrier_shares mode): a carrier of a gas
// species that was in grid cell from_cell at the start of the motion and is
// in to_cell after it (kNone: sunk at a wall). cell_count = RN(exact count of
// the species' carriers in from_cell at the start of the motion).
struct MoverRecord {
  u64 native_id;
  u32 species, from_cell, to_cell, pad;
  double weight, cell_count;
};
// Chemistry v3 CellTargets byte layout (536 bytes).
struct CellTargetsAbi {
  double n[kMaxTargets];
  double drift[kMaxTargets][3];
  double s[kMaxTargets];
  u32 species[kMaxTargets];
  u32 n_targets, pad;
};

// ------------------------------------------------------------------ tables
constexpr u32 kSpGas = 1u;        // neutral: deposits held per (cell, species)
constexpr u32 kSpTarget = 2u;     // chemistry target species
constexpr u32 kSpFeedback = 4u;   // neutral target: eligible for sampling feedback
struct SpeciesInfo {              // 16 bytes
  double mass_kg;
  int charge;
  u32 flags;
};
struct ProcessInfo {              // 16 bytes
  u32 projectile, n_out, internal_carrier, kind;   // internal_carrier: product index or kNone
};

// Entry quantities (per (cell, gas species)).
constexpr u32 kQE = 0, kQPX = 1, kQPY = 2, kQPZ = 3, kQU = 4, kEntryQ = 5;
// Fixed global accumulators.
enum : u32 {
  G_W = 0,           // sum of event weights (physical events)
  G_W_GAS = 1,       // events returning their target to the gas
  G_W_CONSUME = 2,   // target-consuming events
  G_EIN_P = 3,       // projectile kinetic energy in (J)
  G_EIN_T = 4,       // sampled target kinetic energy in (J)
  G_EOUT = 5,        // 4: product kinetic energy out by destination (survivor, birth, reservoir, gas)
  G_Q = 9,           // reference-energy change (J)
  G_I_GAS = 10,      // internal energy held by the gas (J)
  G_I_UN = 11,       // 6: internal energy not held by the gas, by sink category (J)
  G_DEF_E = 17,      // event energy defect (J): in - out - Q - internal
  G_PIN_P = 18,      // 3: projectile momentum in (kg m/s)
  G_PIN_T = 21,      // 3: sampled target momentum in
  G_POUT = 24,       // 12: product momentum out by destination
  G_PINT_GAS = 36,   // 3: momentum of gas-held internal energy
  G_PINT_UN = 39,    // 3: momentum of internal energy not held by the gas
  G_DEF_P = 42,      // 3: event momentum defect
  G_RES_E = 45,      // reservoir mirror: kinetic energy (J), RESERVOIR products - releases
  G_RES_P = 46,      // 3: reservoir mirror momentum
  G_MOVED = 49,      // 5: deposit moved by transport (E, P3, U)
  G_FIXED = 54
};
struct Layout {
  u32 S, NP;   // species, processes
  PFG_HD u32 chg(u32 s, u32 q) const { return G_FIXED + s * 4u + q; }               // charged-target substitution (E, P3)
  PFG_HD u32 wall(u32 s, u32 q) const { return G_FIXED + 4u * S + s * 5u + q; }     // transport shares to walls (E, P3, U)
  PFG_HD u32 flow(u32 kind, u32 s, u32 q) const {                                   // flows (count, E, P3, internal)
    return G_FIXED + 9u * S + (kind * S + s) * 6u + q;
  }
  PFG_HD u32 proc(u32 p, u32 q) const { return G_FIXED + 39u * S + p * 5u + q; }    // per process (W, E_gas, P_gas3)
  PFG_HD u32 total() const { return G_FIXED + 39u * S + 5u * NP; }
};

// Counters (u64; part of the state digest).
enum : u32 {
  C_STAGES = 0, C_RECORDS, C_FLOWS, C_MOVERS, C_FEEDBACK_SLOTS, C_FEEDBACK_CLAMPS, C_CLOSURE_E, C_CLOSURE_P,
  C_STATUS, C_MAXREL_E, C_MAXREL_P, C_TRANSPORTS, C_COUNTERS
};
// status bits (C_STATUS)
constexpr u64 kStBadInput = 1, kStBelowGrid = 2, kStOverflow = 4, kStNoEntry = 8, kStInternal = 16, kStSpecies = 32;

struct Tables {
  const SpeciesInfo *species;
  u32 n_species;
  const ProcessInfo *procs;
  u32 n_proc;
  const u64 *keys;   // ascending (cell << 32 | species)
  u64 n_keys;
  Layout L;
};

// ------------------------------------------------------------------ exact terms
// An exact value placed into the digit grid: digits d .. d+n-1 receive the
// 32-bit chunks c[0..n-1] with sign neg.
struct Chunks {
  u32 c[6];
  int d, n, neg, zero, bad;   // bad: 1 nonfinite factor, 2 below the grid, 4 above the grid
};
PFG_HD Chunks chunks_product(double a, double b) {
  const pf7::Term t = pf7::product_term(a, b);
  Chunks r;
  for (int k = 0; k < 6; ++k)
    r.c[k] = 0;
  for (int k = 0; k < 5; ++k)
    r.c[k] = t.c[k];
  r.d = t.d;
  r.n = 5;
  r.neg = t.neg;
  r.zero = t.zero;
  r.bad = t.bad ? 1 : 0;
  return r;
}
PFG_HD Chunks chunks_value(double a) { return chunks_product(a, 1.0); }
// Exact a*b*c (three binary64 factors; W m u). The 159-bit integer product of
// the significands is placed at grid position ea + eb + ec + 2148; a value
// with nonzero bits below 2^-2148 or at/above bit 32*133 is refused (bad).
PFG_HD Chunks chunks_triple(double a, double b, double c) {
  Chunks r;
  for (int k = 0; k < 6; ++k)
    r.c[k] = 0;
  r.d = 0;
  r.n = 6;
  r.neg = 0;
  r.zero = 1;
  r.bad = 0;
  pf7_u64 ma = 0, mb = 0, mc = 0;
  int ea = 0, eb = 0, ec = 0, sa = 0, sb = 0, sc = 0;
  const int fa = pf7::decompose(a, ma, ea, sa), fb = pf7::decompose(b, mb, eb, sb), fc = pf7::decompose(c, mc, ec, sc);
  if (!fa || !fb || !fc) {
    r.bad = 1;
    return r;
  }
  if (ma == 0 || mb == 0 || mc == 0)
    return r;
  // P = ma mb (< 2^106) as (h1:l1); Q = P mc (< 2^159) as (w2:w1:w0).
  const pf7_u64 l1 = ma * mb, h1 = pf7::mul_hi(ma, mb);
  const pf7_u64 a0 = l1 * mc, a1 = pf7::mul_hi(l1, mc);
  const pf7_u64 b0 = h1 * mc, b1 = pf7::mul_hi(h1, mc);
  pf7_u64 w0 = a0, w1 = a1 + b0, w2 = b1 + (w1 < a1 ? 1ull : 0ull);
  int pos = ea + eb + ec + pf7::acc_offset;
  if (pos < 0) {
    const int s = -pos;
    if (s >= 192) {
      r.bad = 2;
      return r;
    }
    // bits shifted out must be zero
    pf7_u64 lost = 0;
    if (s >= 128) {
      lost = w0 | w1 | (s > 128 ? (w2 & ((1ull << (s - 128)) - 1)) : 0ull);
      w0 = s - 128 < 64 ? (w2 >> (s - 128)) : 0;
      w1 = 0;
      w2 = 0;
    } else if (s >= 64) {
      lost = w0 | (s > 64 ? (w1 & ((1ull << (s - 64)) - 1)) : 0ull);
      const int t = s - 64;
      w0 = t ? ((w1 >> t) | (w2 << (64 - t))) : w1;
      w1 = t ? (w2 >> t) : w2;
      w2 = 0;
    } else if (s > 0) {
      lost = w0 & ((1ull << s) - 1);
      w0 = (w0 >> s) | (w1 << (64 - s));
      w1 = (w1 >> s) | (w2 << (64 - s));
      w2 = w2 >> s;
    }
    if (lost) {
      r.bad = 2;
      return r;
    }
    pos = 0;
  }
  const int sh = pos & 31;
  r.d = pos >> 5;
  if (r.d + 6 > kDigits - 1) {
    r.bad = 4;
    return r;
  }
  const pf7_u64 v0 = w0 << sh;
  const pf7_u64 v1 = sh ? ((w1 << sh) | (w0 >> (64 - sh))) : w1;
  const pf7_u64 v2 = sh ? ((w2 << sh) | (w1 >> (64 - sh))) : w2;   // w2 < 2^31: no overflow
  r.c[0] = (u32)(v0 & 0xffffffffull);
  r.c[1] = (u32)(v0 >> 32);
  r.c[2] = (u32)(v1 & 0xffffffffull);
  r.c[3] = (u32)(v1 >> 32);
  r.c[4] = (u32)(v2 & 0xffffffffull);
  r.c[5] = (u32)(v2 >> 32);
  r.neg = sa ^ sb ^ sc;
  r.zero = 0;
  return r;
}
// Host/sequential accumulation (the device uses atomic adds of the same chunks).
PFG_HD void add_chunks(i64 *digit, const Chunks &t, int negate) {
  if (t.zero || t.bad)
    return;
  const int neg = t.neg ^ negate;
  for (int k = 0; k < t.n; ++k)
    digit[t.d + k] += neg ? -(i64)t.c[k] : (i64)t.c[k];
}
PFG_HD u64 bad_status(const Chunks &t) {
  if (!t.bad)
    return 0;
  return t.bad == 1 ? kStBadInput : (t.bad == 2 ? kStBelowGrid : kStOverflow);
}

// ------------------------------------------------------------------ physics helpers
// Kinetic energy (J): native pic_kinetic_energy_J (forms kinetic_energy).
PFG_HD double kinetic(double m, const double *u) { return pf7::kinetic_energy(m, u[0], u[1], u[2]); }
// Momentum of the internal energy I_lab (J) carried by a product with proper
// velocity u: the internal rest energy I = I_lab / gamma moves with the
// product, p = (I / c^2) u. Fixed operation order.
PFG_HD void internal_momentum(double I_lab, const double *u, double *p) {
  const double n = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
  const double g = pf7::fsqrt(1 + n / pf7::c_squared);
  const double f = (I_lab / g) / pf7::c_squared;
  p[0] = f * u[0];
  p[1] = f * u[1];
  p[2] = f * u[2];
}
PFG_HD double fabs_(double x) { return x < 0 ? -x : x; }
PFG_HD int finite3(const double *u) { return pf7::finite(u[0]) && pf7::finite(u[1]) && pf7::finite(u[2]); }

PFG_HD u64 entry_key(u32 cell, u32 species) { return ((u64)cell << 32) | (u64)species; }
// Index of key in the ascending table, or -1.
PFG_HD i64 find_key(const u64 *keys, u64 n, u64 key) {
  u64 lo = 0, hi = n;
  while (lo < hi) {
    const u64 mid = lo + ((hi - lo) >> 1);
    if (keys[mid] < key)
      lo = mid + 1;
    else
      hi = mid;
  }
  return (lo < n && keys[lo] == key) ? (i64)lo : -1;
}

// ------------------------------------------------------------------ record application
// Sink interface (host: plain adds; device: atomics):
//   void entry(i64 e, u32 q, const Chunks &, int negate);
//   void global(u32 g, const Chunks &, int negate);
//   void status(u64 bits);
//   void closure(double relE, double relP);   // per-event closure (counters)
// Gas side of a species: an owner entry (gas species) or the charged-target
// substitution ledger.
template <class Sink>
PFG_HD void put_side(Sink &S, const Tables &T, u32 cell, u32 species, u32 q, const Chunks &c, int negate) {
  if (T.species[species].flags & kSpGas) {
    const i64 e = find_key(T.keys, T.n_keys, entry_key(cell, species));
    if (e < 0) {
      S.status(kStNoEntry);
      return;
    }
    S.entry(e, q, c, negate);
  } else if (q < 4)
    S.global(T.L.chg(species, q), c, negate);
  else
    S.status(kStSpecies);   // internal energy is held only by gas species
}
template <class Sink> PFG_HD void put_global(Sink &S, u32 g, const Chunks &c, int negate) {
  if (c.bad)
    S.status(bad_status(c));
  S.global(g, c, negate);
}

// One admitted chemistry event (README 3). Per-event closure: the energy and
// momentum defects are evaluated in binary64 in a fixed order, relative to
// the event's energy/momentum scale (the exact defects are global ledgers).
template <class Sink> PFG_HD void apply_record(Sink &S, const Tables &T, const RecoilRecord &r) {
  const ChemEvent &ev = r.ev;
  const double W = ev.weight;
  if (ev.process >= T.n_proc || ev.target_species >= T.n_species || ev.n_out < 1 || ev.n_out > kMaxOut ||
      !(W > 0) || !pf7::finite(W) || !finite3(ev.u_in) || !finite3(ev.ut) || !pf7::finite(ev.q_eV) ||
      !pf7::finite(ev.e_internal) || ev.internal_sink >= kSinks) {
    S.status(kStBadInput);
    return;
  }
  const ProcessInfo &pr = T.procs[ev.process];
  const u32 sp = pr.projectile, st = ev.target_species, cell = r.grid_cell;
  if (sp >= T.n_species || pr.n_out != ev.n_out) {
    S.status(kStBadInput);
    return;
  }
  for (u32 i = 0; i < ev.n_out; ++i)
    if (ev.out_species[i] >= T.n_species || ev.out_dest[i] < kDestSurvivor || ev.out_dest[i] > kDestGas ||
        !finite3(ev.out_u[i])) {
      S.status(kStBadInput);
      return;
    }
  const double mp = T.species[sp].mass_kg, mt = T.species[st].mass_kg;
  const double Kp = kinetic(mp, ev.u_in), Kt = kinetic(mt, ev.ut);
  int consuming = 1;
  for (u32 i = 0; i < ev.n_out; ++i)
    if (ev.out_dest[i] == kDestGas)
      consuming = 0;
  put_global(S, G_W, chunks_value(W), 0);
  put_global(S, consuming ? G_W_CONSUME : G_W_GAS, chunks_value(W), 0);
  put_global(S, T.L.proc(ev.process, 0), chunks_value(W), 0);
  // in
  const Chunks cKp = chunks_product(W, Kp), cKt = chunks_product(W, Kt);
  put_global(S, G_EIN_P, cKp, 0);
  put_global(S, G_EIN_T, cKt, 0);
  put_global(S, G_DEF_E, cKp, 0);
  put_global(S, G_DEF_E, cKt, 0);
  double sumE = Kp + Kt, scaleE = fabs_(Kp) + fabs_(Kt);
  double sumP[3], scaleP = 0;
  for (u32 a = 0; a < 3; ++a) {
    const Chunks pp = chunks_triple(W, mp, ev.u_in[a]), pt = chunks_triple(W, mt, ev.ut[a]);
    put_global(S, G_PIN_P + a, pp, 0);
    put_global(S, G_PIN_T + a, pt, 0);
    put_global(S, G_DEF_P + a, pp, 0);
    put_global(S, G_DEF_P + a, pt, 0);
    const double xp = mp * ev.u_in[a], xt = mt * ev.ut[a];
    sumP[a] = xp + xt;
    scaleP = scaleP + fabs_(xp) + fabs_(xt);
  }
  // the sampled target leaves the gas (consumed, or replaced by its GAS product)
  put_side(S, T, cell, st, kQE, cKt, 1);
  for (u32 a = 0; a < 3; ++a)
    put_side(S, T, cell, st, kQPX + a, chunks_triple(W, mt, ev.ut[a]), 1);
  if (T.species[st].flags & kSpGas) {
    put_global(S, T.L.proc(ev.process, 1), cKt, 1);
    for (u32 a = 0; a < 3; ++a)
      put_global(S, T.L.proc(ev.process, 2 + a), chunks_triple(W, mt, ev.ut[a]), 1);
  }
  // out
  for (u32 i = 0; i < ev.n_out; ++i) {
    const u32 so = ev.out_species[i], dest = ev.out_dest[i];
    const double mo = T.species[so].mass_kg;
    const double Ko = kinetic(mo, ev.out_u[i]);
    const Chunks cKo = chunks_product(W, Ko);
    put_global(S, G_EOUT + (dest - 1), cKo, 0);
    put_global(S, G_DEF_E, cKo, 1);
    sumE = sumE - Ko;
    scaleE = scaleE + fabs_(Ko);
    for (u32 a = 0; a < 3; ++a) {
      const Chunks po = chunks_triple(W, mo, ev.out_u[i][a]);
      put_global(S, G_POUT + 3 * (dest - 1) + a, po, 0);
      put_global(S, G_DEF_P + a, po, 1);
      const double xo = mo * ev.out_u[i][a];
      sumP[a] = sumP[a] - xo;
      scaleP = scaleP + fabs_(xo);
      if (dest == kDestGas) {
        put_side(S, T, cell, so, kQPX + a, po, 0);
        if (T.species[so].flags & kSpGas)
          put_global(S, T.L.proc(ev.process, 2 + a), po, 0);
      }
      if (dest == kDestReservoir)
        put_global(S, G_RES_P + a, po, 0);
    }
    if (dest == kDestGas) {
      put_side(S, T, cell, so, kQE, cKo, 0);
      if (T.species[so].flags & kSpGas)
        put_global(S, T.L.proc(ev.process, 1), cKo, 0);
    }
    if (dest == kDestReservoir)
      put_global(S, G_RES_E, cKo, 0);
  }
  // reference energy
  const double qJ = ev.q_eV * kEv;
  const Chunks cQ = chunks_product(W, qJ);
  put_global(S, G_Q, cQ, 0);
  put_global(S, G_DEF_E, cQ, 1);
  sumE = sumE - qJ;
  scaleE = scaleE + fabs_(qJ);
  // internal energy and its momentum
  if (ev.e_internal != 0.0) {
    const u32 k = pr.internal_carrier;
    if (k >= ev.n_out) {
      S.status(kStInternal);
    } else {
      const double IJ = ev.e_internal * kEv;
      double pint[3];
      internal_momentum(IJ, ev.out_u[k], pint);
      const Chunks cI = chunks_product(W, IJ);
      const u32 sk = ev.out_species[k], dk = ev.out_dest[k];
      const int gas_held = (dk == kDestGas || dk == kDestReservoir) && (T.species[sk].flags & kSpGas);
      put_global(S, gas_held ? G_I_GAS : G_I_UN + ev.internal_sink, cI, 0);
      put_global(S, G_DEF_E, cI, 1);
      if (gas_held) {
        put_side(S, T, cell, sk, kQU, cI, 0);
        put_global(S, T.L.proc(ev.process, 1), cI, 0);
      }
      for (u32 a = 0; a < 3; ++a) {
        const Chunks cp = chunks_product(W, pint[a]);
        put_global(S, (gas_held ? G_PINT_GAS : G_PINT_UN) + a, cp, 0);
        put_global(S, G_DEF_P + a, cp, 1);
        if (gas_held) {
          put_side(S, T, cell, sk, kQPX + a, cp, 0);
          put_global(S, T.L.proc(ev.process, 2 + a), cp, 0);
        }
        sumP[a] = sumP[a] - pint[a];
        scaleP = scaleP + fabs_(pint[a]);
      }
      sumE = sumE - IJ;
      scaleE = scaleE + fabs_(IJ);
    }
  }
  // per-event closure (relative to the event scale; per physical event)
  const double relE = fabs_(sumE) / (scaleE + 1e-300);
  double relP = 0;
  for (u32 a = 0; a < 3; ++a) {
    const double x = fabs_(sumP[a]) / (scaleP + 1e-300);
    relP = x > relP ? x : relP;
  }
  S.closure(relE, relP);
}

// Flow (README 5.1): population ledger of (kind, species); a debit of a gas
// species returns its actual carriers to the gas deposit (the sampled target
// was removed by the record), a debit of a charged target to its substitution
// ledger; debited internal energy (native v8.1 release at seal) goes to the
// gas internal energy.
template <class Sink> PFG_HD void apply_flow(Sink &S, const Tables &T, const FlowRecord &f) {
  if (f.kind >= kFlowKinds || f.species >= T.n_species || !(f.weight > 0) || !pf7::finite(f.weight) ||
      !finite3(f.u) || !pf7::finite(f.internal_J) || f.negative > 1) {
    S.status(kStBadInput);
    return;
  }
  const double m = T.species[f.species].mass_kg;
  const double K = kinetic(m, f.u);
  const int neg = (int)f.negative;
  const Layout &L = T.L;
  put_global(S, L.flow(f.kind, f.species, 0), chunks_value(f.weight), neg);
  const Chunks cK = chunks_product(f.weight, K);
  put_global(S, L.flow(f.kind, f.species, 1), cK, neg);
  for (u32 a = 0; a < 3; ++a)
    put_global(S, L.flow(f.kind, f.species, 2 + a), chunks_triple(f.weight, m, f.u[a]), neg);
  put_global(S, L.flow(f.kind, f.species, 5), chunks_value(f.internal_J), neg);
  if (f.kind == kFlowDebit) {
    // deposit effect = -(population effect)
    put_side(S, T, f.grid_cell, f.species, kQE, cK, 1 - neg);
    for (u32 a = 0; a < 3; ++a)
      put_side(S, T, f.grid_cell, f.species, kQPX + a, chunks_triple(f.weight, m, f.u[a]), 1 - neg);
    if (neg && f.internal_J != 0.0 && (T.species[f.species].flags & kSpGas))
      put_side(S, T, f.grid_cell, f.species, kQU, chunks_value(f.internal_J), 0);
  }
}

// Transport share (README 5.2). value: rounded entry quantities at the start
// of the transport (5 per entry).
template <class Sink> PFG_HD void apply_mover(Sink &S, const Tables &T, const MoverRecord &m, const double *value) {
  if (m.species >= T.n_species || !(m.weight > 0) || !(m.cell_count >= m.weight) || !pf7::finite(m.cell_count)) {
    S.status(kStBadInput);
    return;
  }
  if (!(T.species[m.species].flags & kSpGas))
    return;
  const i64 e = find_key(T.keys, T.n_keys, entry_key(m.from_cell, m.species));
  if (e < 0)
    return;   // nothing held in the source cell
  i64 t = -1;
  if (m.to_cell != kNone) {
    t = find_key(T.keys, T.n_keys, entry_key(m.to_cell, m.species));
    if (t < 0) {
      S.status(kStNoEntry);
      return;
    }
  }
  const double f = m.weight / m.cell_count;
  for (u32 q = 0; q < kEntryQ; ++q) {
    const double x = value[(u64)e * kEntryQ + q] * f;
    if (x == 0.0)
      continue;
    const Chunks c = chunks_value(x);
    S.entry(e, q, c, 1);
    if (t >= 0)
      S.entry(t, q, c, 0);
    else
      S.global(T.L.wall(m.species, q), c, 0);
    S.global(G_MOVED + q, c, 0);
  }
}

// ------------------------------------------------------------------ feedback
// Sampling representation of the local gas state (README 4). Slot k of a cell:
// N = n V; shift d = P_d / (N m); thermal energy per unit mass changes by
// E_d / (N m) - drift.d - |d|^2 / 2 (the bulk-flow energy of the shift is
// removed); s'^2 = s^2 + (2/3) of that, clamped at 0. A slot whose deposit is
// exactly zero is left untouched (bit-identical to feedback off).
// value: E, Px, Py, Pz, U of the entry. Returns 1 if the slot changed, and sets
// *clamped when s'^2 < 0.
PFG_HD int feedback_slot(double n, double V, double m, const double *value, double *drift, double *s, int *clamped) {
  *clamped = 0;
  const double E = value[kQE], P0 = value[kQPX], P1 = value[kQPY], P2 = value[kQPZ];
  if (E == 0.0 && P0 == 0.0 && P1 == 0.0 && P2 == 0.0)
    return 0;
  const double N = n * V;
  if (!(N > 0) || !(m > 0))
    return 0;
  const double Nm = N * m;
  const double d0 = P0 / Nm, d1 = P1 / Nm, d2 = P2 / Nm;
  const double cross = (drift[0] * d0 + drift[1] * d1) + drift[2] * d2;
  const double half = ((d0 * d0 + d1 * d1) + d2 * d2) / 2;
  const double de = (E / Nm - cross) - half;
  const double s2 = (*s) * (*s) + (2.0 / 3.0) * de;
  drift[0] = drift[0] + d0;
  drift[1] = drift[1] + d1;
  drift[2] = drift[2] + d2;
  if (s2 > 0)
    *s = pf7::fsqrt(s2);
  else {
    *s = 0.0;
    *clamped = 1;
  }
  return 1;
}

// ------------------------------------------------------------------ digests
// Entry leaf: SHA256(text "gas-target-owner-v1/entry", u64 key, 5 x 134
// normalised digits as u64). Global leaf: SHA256(text
// "gas-target-owner-v1/global", u64 index, 134 digits). digit must be
// normalised (pf7::normalise).
PFG_HD void entry_leaf(u64 key, const i64 *digit5, pf7_u32 *out) {
  pf7::Sha s;
  pf7::sha_init(s);
  const char tag[] = "gas-target-owner-v1/entry";
  pf7::sha_text(s, tag, (pf7_u32)(sizeof(tag) - 1));
  pf7::sha_u64(s, key);
  pf7::sha_u64_array(s, (const pf7_u64 *)digit5, (pf7_u64)kEntryQ * kDigits);
  pf7::sha_finish(s, out);
}
PFG_HD void global_leaf(u64 index, const i64 *digit, pf7_u32 *out) {
  pf7::Sha s;
  pf7::sha_init(s);
  const char tag[] = "gas-target-owner-v1/global";
  pf7::sha_text(s, tag, (pf7_u32)(sizeof(tag) - 1));
  pf7::sha_u64(s, index);
  pf7::sha_u64_array(s, (const pf7_u64 *)digit, (pf7_u64)kDigits);
  pf7::sha_finish(s, out);
}
PFG_HD int digits_zero(const i64 *digit, int n) {
  for (int k = 0; k < n; ++k)
    if (digit[k])
      return 0;
  return 1;
}

}  // namespace pfg
