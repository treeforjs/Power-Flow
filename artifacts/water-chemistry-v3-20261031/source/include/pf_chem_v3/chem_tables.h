// chem_tables.h -- immutable chemistry v3 tables (binary blob "PFCHEMV3") and
// their evaluation, shared by host and NVRTC.
//
// The blob is produced by tools/build_tables.py from data/tables/
// chemistry-v3-spec.json and the extracted source tables; its SHA-256 is the
// table identity. All records below are plain little-endian structures with
// the exact sizes asserted on the host.
//
// Cross sections. Processes are grouped by (projectile species, target
// species). A group owns one ascending energy grid E_0 < ... < E_{n-1} (eV,
// projectile kinetic energy in the target rest frame) that contains every
// breakpoint of every source table of its processes, every threshold and
// every coverage edge. sigma_k(E) is stored at every grid point; between grid
// points it is linear in E and sigma (the declared interpolation of the
// tabulated sources, so for tabulated sources the evaluation equals the source
// law up to rounding; analytic models are sampled densely and their
// interpolation error is recorded by the builder). Outside [E_0, E_{n-1}]
// every sampled process is zero (coverage classes say whether that zero is
// known or unresolved).
//
// Majorant. rmq holds a sparse table (levels l = 0..L-1) of upper bounds of
// F(E) = sum_k sigma_k(E) v(E) over grid intervals [E_i, E_{i+1}]: level l,
// entry i bounds intervals i .. i + 2^l - 1. v is the projectile speed in the
// target frame (m/s). Bounds are computed by the builder from the exact
// piecewise-linear sigma and the exact speed law, then multiplied by
// (1 + 2^-20).
#pragma once
#include "det_math.h"

namespace pfc {

// Species kinds.
#define PFC_KIND_ELECTRON 0u
#define PFC_KIND_CATION 1u
#define PFC_KIND_ANION 2u
#define PFC_KIND_NEUTRAL 3u

// Process kinds.
#define PFC_P_ELASTIC 1u        // A + B -> A + B (angular law)
#define PFC_P_EXCITE 2u         // A + B -> A + B(internal loss)
#define PFC_P_DISSOC 3u         // A + B -> A + B*(loss) -> A + chain(fragments)
#define PFC_P_IONIZE 4u         // e + B -> e + e + Xv -> chain(ion, fragments)
#define PFC_P_CAPTURE 5u        // A + B -> C*(K_cm - Q) -> chain(products)
#define PFC_P_CHARGE_X 6u       // A+ + B -> A(fast, forward) + B+(slow)

// Angular laws (scatter kinds).
#define PFC_ANG_ISO 0u
#define PFC_ANG_SCREENED 1u     // xi(E) row of the group
#define PFC_ANG_FORWARD 2u      // CM angle 0 for the first product

// Product destinations.
#define PFC_DEST_SURVIVOR 1u    // the projectile carrier continues (same id)
#define PFC_DEST_BIRTH 2u       // new charged carrier (electron or ion)
#define PFC_DEST_RESERVOIR 3u   // neutral product to the quantized-release reservoir
#define PFC_DEST_GAS 4u         // target returned to its gas inventory (recoil)

// Internal-energy sink categories.
#define PFC_SINK_NONE 0u
#define PFC_SINK_ROT 1u
#define PFC_SINK_VIB 2u
#define PFC_SINK_ELEC 3u          // non-dissociative electronic excitation
#define PFC_SINK_ION_INTERNAL 4u  // ion left in an excited state (orbital hole)
#define PFC_SINK_FRAGMENT 5u      // unpartitioned fragment excitation
#define PFC_SINK_COUNT 6u

// Coverage classes.
#define PFC_COVER_ZERO 0u         // known zero (below threshold)
#define PFC_COVER_DATA 1u         // measured or recommended data
#define PFC_COVER_MODEL 2u        // model-derived values (tagged)
#define PFC_COVER_UNRESOLVED 3u   // no data: sampled as zero, exposure counted

#define PFC_MAX_OUT 6
#define PFC_MAX_ORB 4
#define PFC_MAX_COVER 6
#define PFC_MAX_REGIONS 32
#define PFC_MAX_TARGETS 12

struct SpeciesRec {          // 64 bytes
    double mass_eV;          // rest energy
    double ref_eV;           // 0 K formation enthalpy (chemical reference energy)
    int charge;              // elementary charges
    pfc_u32 kind;
    pfc_u32 nuclei[6];       // H, O, Fe, Cr, Ni, spare
    pfc_u32 pad[4];
};
struct GroupRec {            // 64 bytes
    pfc_u32 projectile, target;
    pfc_u32 grid_off, grid_n;
    pfc_u32 proc_off, proc_n;      // sampled processes (contiguous in the process table)
    pfc_u32 sigma_off;             // proc j (local), point i: sigma[sigma_off + j*grid_n + i]
    pfc_u32 bound_off, bound_n;    // unresolved-bound rows: bound[bound_off + r*grid_n + i]
    pfc_u32 rmq_off, rmq_levels;   // rmq[rmq_off + l*(grid_n-1) + i]
    pfc_u32 xi_off;                // xi[xi_off + i], or 0xFFFFFFFF
    pfc_u32 pad[4];
};
struct ProcessRec {          // 512 bytes
    pfc_u32 kind, group, angular, tag;      // tag: 1 measured, 2 recommended, 3 calculated, 4 model_derived
    pfc_u32 n_out, internal_sink, n_chain, n_orb;
    pfc_u32 out_species[PFC_MAX_OUT];       // final particles, in record order
    pfc_u32 out_dest[PFC_MAX_OUT];
    pfc_u32 n_cover, flags, pad0, pad1;
    double loss_eV;          // EXCITE: internal energy; DISSOC: state energy
    double q_eV;             // reference-energy change of the full event (products - reactants)
    double kinetic_fraction; // fraction of the chain's available energy given to kinetic energy
    double orb_B[PFC_MAX_ORB], orb_U[PFC_MAX_ORB], orb_N[PFC_MAX_ORB], orb_w[PFC_MAX_ORB];
    double cover_lo[PFC_MAX_COVER], cover_hi[PFC_MAX_COVER];
    pfc_u32 cover_class[PFC_MAX_COVER], cover_region[PFC_MAX_COVER];
    double chain_q_eV;       // reference energy of the chain parent above its products (eV, >= 0 enters K)
    double bohr_m, rydberg_eV, relativistic_from_eV;
    double pad2[11];
};

// Pointers into one blob (host memory or device memory).
struct TableView {
    pfc_u32 n_species, n_groups, n_proc, n_regions;
    const SpeciesRec* species;
    const GroupRec* groups;
    const ProcessRec* procs;
    const double* grid;
    const double* sigma;
    const double* bound;
    const double* rmq;
    const double* xi;
    const pfc_u32* group_of;    // [projectile * n_species + target] -> group index, 0xFFFFFFFF if none
};
PFC_FN pfc_u32 group_index(const TableView& t, pfc_u32 projectile, pfc_u32 target) {
    if (projectile >= t.n_species || target >= t.n_species) return 0xFFFFFFFFu;
    return t.group_of[projectile * t.n_species + target];
}

// Index i with grid[i] <= E < grid[i+1] for E inside [grid[0], grid[n-1]];
// returns -1 below, n-1 at or above the last point.
PFC_FN int grid_locate(const double* g, pfc_u32 n, double E) {
    if (!(E >= g[0])) return -1;
    if (E >= g[n - 1]) return (int)n - 1;
    pfc_u32 lo = 0, hi = n - 1;   // g[lo] <= E < g[hi]
    while (hi - lo > 1u) {
        const pfc_u32 mid = lo + ((hi - lo) >> 1);
        if (E < g[mid]) hi = mid;
        else lo = mid;
    }
    return (int)lo;
}
// Linear interpolation of row r (grid_n values) at E with located index i.
PFC_FN double row_eval(const double* row, const double* g, pfc_u32 n, int i, double E) {
    if (i < 0) return 0.0;
    if (i >= (int)n - 1) return (E == g[n - 1]) ? row[n - 1] : 0.0;
    const double x0 = g[i], x1 = g[i + 1], y0 = row[i], y1 = row[i + 1];
    return y0 + (y1 - y0) * ((E - x0) / (x1 - x0));
}
PFC_FN double sigma_eval(const TableView& t, const GroupRec& gr, pfc_u32 local, int i, double E) {
    return row_eval(t.sigma + gr.sigma_off + (pfc_u64)local * gr.grid_n, t.grid + gr.grid_off, gr.grid_n, i, E);
}
PFC_FN double bound_eval(const TableView& t, const GroupRec& gr, pfc_u32 r, int i, double E) {
    return row_eval(t.bound + gr.bound_off + (pfc_u64)r * gr.grid_n, t.grid + gr.grid_off, gr.grid_n, i, E);
}
PFC_FN double xi_eval(const TableView& t, const GroupRec& gr, int i, double E) {
    if (gr.xi_off == 0xFFFFFFFFu) return 0.0;
    const double* g = t.grid + gr.grid_off;
    const double* row = t.xi + gr.xi_off;
    if (i < 0) return row[0];
    if (i >= (int)gr.grid_n - 1) return row[gr.grid_n - 1];
    return row_eval(row, g, gr.grid_n, i, E);
}
PFC_FN int floor_log2(pfc_u32 x) {
    int r = -1;
    while (x) {
        x >>= 1;
        ++r;
    }
    return r;
}
// Upper bound of F(E) = sum_k sigma_k(E) v(E) over E in [Elo, Ehi] (target frame).
PFC_FN double majorant(const TableView& t, const GroupRec& gr, double Elo, double Ehi) {
    const double* g = t.grid + gr.grid_off;
    const pfc_u32 n = gr.grid_n;
    if (n < 2u) return 0.0;
    if (!(Ehi >= g[0]) || Elo > g[n - 1]) return 0.0;
    int a = grid_locate(g, n, Elo < g[0] ? g[0] : Elo);
    int b = grid_locate(g, n, Ehi);
    if (a < 0) a = 0;
    if (a > (int)n - 2) a = (int)n - 2;
    if (b > (int)n - 2) b = (int)n - 2;
    if (b < a) b = a;
    const pfc_u32 len = (pfc_u32)(b - a + 1);
    const int l = floor_log2(len);
    const double* r = t.rmq + gr.rmq_off + (pfc_u64)l * (n - 1);
    const double x = r[a], y = r[b - (1 << l) + 1];
    return x > y ? x : y;
}
// Coverage class of process p at energy E and its region index.
PFC_FN pfc_u32 coverage_class(const ProcessRec& p, double E, pfc_u32* region) {
    for (pfc_u32 c = 0; c < p.n_cover; ++c) {
        if (E >= p.cover_lo[c] && E < p.cover_hi[c]) {
            *region = p.cover_region[c];
            return p.cover_class[c];
        }
    }
    *region = 0xFFFFFFFFu;
    return PFC_COVER_UNRESOLVED;
}

}  // namespace pfc
