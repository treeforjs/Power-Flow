// merge_split.cpp -- CPU reference of the merge/split transaction (see merge_split.hpp).
#include "merge_split.hpp"

#include <openssl/evp.h>

#include <algorithm>
#include <bit>
#include <cstring>
#include <limits>
#include <numeric>
#include <queue>
#include <random>
#include <string_view>

namespace merge_split {
namespace {

void need(bool ok, const char* message) {
    if (!ok) throw Error(message);
}

std::uint64_t splitmix64(std::uint64_t x) {
    x += 0x9E3779B97F4A7C15ull;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ull;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBull;
    return x ^ (x >> 31);
}

// Canonical little-endian encoder of the digest input.
struct Encoder {
    std::string s;
    void u64(std::uint64_t v) {
        for (int i = 0; i < 8; ++i) s.push_back(char((v >> (8 * i)) & 0xff));
    }
    void u32(std::uint32_t v) {
        for (int i = 0; i < 4; ++i) s.push_back(char((v >> (8 * i)) & 0xff));
    }
    void f64(double x) { u64(std::bit_cast<std::uint64_t>(x)); }
    void text(std::string_view t) {
        u64(t.size());
        s.append(t);
    }
    void exact(const Exact& e) {
        for (auto x : e.d) u64(std::uint64_t(x));
    }
};

// Nonnegative magnitude of an exact value: digits in [0, 2^32).
struct Magnitude {
    std::array<std::uint64_t, MS_ACC_DIGITS> d{};
    int top() const {
        for (int i = MS_ACC_DIGITS - 1; i >= 0; --i)
            if (d[i]) return 32 * i + 31 - std::countl_zero(std::uint32_t(d[i]));
        return -1;
    }
    int bit(int pos) const { return pos < 0 ? 0 : int((d[pos >> 5] >> (pos & 31)) & 1u); }
    bool any_below(int pos) const {   // bits [0, pos)
        if (pos <= 0) return false;
        const int i = pos >> 5;
        for (int j = 0; j < i; ++j)
            if (d[j]) return true;
        const std::uint64_t mask = (std::uint64_t(1) << (pos & 31)) - 1;
        return (d[i] & mask) != 0;
    }
    std::uint64_t bits(int lo, int count) const {
        std::uint64_t r = 0;
        for (int k = 0; k < count; ++k) r |= std::uint64_t(bit(lo + k)) << k;
        return r;
    }
};

Magnitude magnitude(const Exact& e, bool& negative) {
    negative = e.negative();
    MsAcc a;
    for (int i = 0; i < MS_ACC_DIGITS; ++i) a.d[i] = negative ? -e.d[i] : e.d[i];
    ms_acc_normalize(&a);
    Magnitude m;
    for (int i = 0; i < MS_ACC_DIGITS; ++i) {
        need(a.d[i] >= 0 && a.d[i] < (std::int64_t(1) << 32), "exact value outside the superaccumulator range");
        m.d[i] = std::uint64_t(a.d[i]);
    }
    return m;
}

double round_magnitude(const Magnitude& m) {
    const int H = m.top();
    if (H < 0) return 0;
    if (H <= 52) return std::bit_cast<double>(m.bits(0, H + 1));
    std::uint64_t M = m.bits(H - 52, 53);
    const bool half = m.bit(H - 53), sticky = m.any_below(H - 53);
    if (half && (sticky || (M & 1))) ++M;
    const std::uint64_t b = (std::uint64_t(H - 52) << 52) + M;
    need((b >> 52) < 2047, "exact value overflows binary64");
    return std::bit_cast<double>(b);
}

// Quantity layout per species (shared with the GPU ledger kernel).
enum Q { W_IN = 0, W_OUT = 1, P_IN = 2, P_OUT = 5, E_IN = 8, E_OUT = 9, NQ = 10 };

void add_terms(MsAcc* q, bool out, double w, const double u[3], double ke, double m) {
    ms_acc_add(&q[out ? W_OUT : W_IN], w);
    for (int d = 0; d < 3; ++d) ms_acc_add_product3(&q[(out ? P_OUT : P_IN) + d], w, u[d], m);
    ms_acc_add_product2(&q[out ? E_OUT : E_IN], w, ke);
}

struct Rec {
    std::uint32_t cell, species, bin, index;
    std::uint64_t id;
    double ke;
};
bool key_less(const Rec& a, const Rec& b) {
    if (a.cell != b.cell) return a.cell < b.cell;
    if (a.species != b.species) return a.species < b.species;
    if (a.bin != b.bin) return a.bin < b.bin;
    return a.id < b.id;
}

void encode_policy(Encoder& e, const MsParams& P) {
    e.text("merge-split-policy-v1");
    e.u32(P.n_species);
    e.u32(P.nsub);
    for (auto v : {P.nx, P.ny, P.nz, P.kz_lo, P.kz_hi, P.cell_bits, P.sub[0], P.sub[1], P.sub[2]}) e.u32(v);
    for (int d = 0; d < 3; ++d) e.f64(P.shift[d]);
    for (int d = 0; d < 3; ++d) e.f64(P.origin[d]);
    for (int d = 0; d < 3; ++d) e.f64(P.h[d]);
    for (unsigned s = 0; s < P.n_species; ++s) {
        const auto& sp = P.species[s];
        for (double x : {sp.mass_kg, sp.charge_C, sp.mc2_J, sp.split_min_weight}) e.f64(x);
        e.u32(sp.n_thresh);
        for (unsigned i = 0; i < sp.n_thresh; ++i) e.f64(sp.e_thresh[i]);
        for (auto v : {sp.merge_above, sp.merge_target, sp.split_below, sp.split_target, sp.k_min, sp.k_max, sp.split_levels})
            e.u32(v);
    }
    e.text("philox4x32-10/merge-split-stream-v1");
}

}  // namespace

// ------------------------------------------------------------------ Carriers / Exact
void Carriers::reserve(std::size_t n) {
    id.reserve(n);
    species.reserve(n);
    flags.reserve(n);
    for (int d = 0; d < 3; ++d) { x[d].reserve(n); u[d].reserve(n); }
    w.reserve(n);
}
void Carriers::push(std::uint64_t i, std::uint32_t s, std::uint32_t f, const double xx[3], const double uu[3], double ww) {
    id.push_back(i);
    species.push_back(s);
    flags.push_back(f);
    for (int d = 0; d < 3; ++d) { x[d].push_back(xx[d]); u[d].push_back(uu[d]); }
    w.push_back(ww);
}

Exact Exact::of(MsAcc a) {
    ms_acc_normalize(&a);
    Exact e;
    for (int i = 0; i < MS_ACC_DIGITS; ++i) e.d[i] = a.d[i];
    need(e.d[MS_ACC_DIGITS - 1] > -(std::int64_t(1) << 32) && e.d[MS_ACC_DIGITS - 1] < (std::int64_t(1) << 32),
         "exact value outside the superaccumulator range");
    return e;
}
Exact Exact::operator-(const Exact& o) const {
    MsAcc a;
    for (int i = 0; i < MS_ACC_DIGITS; ++i) a.d[i] = d[i] - o.d[i];
    return of(a);
}
Exact Exact::operator+(const Exact& o) const {
    MsAcc a;
    for (int i = 0; i < MS_ACC_DIGITS; ++i) a.d[i] = d[i] + o.d[i];
    return of(a);
}
bool Exact::zero() const noexcept {
    for (auto x : d)
        if (x) return false;
    return true;
}
double Exact::rounded() const {
    bool neg = false;
    const double r = round_magnitude(magnitude(*this, neg));
    return neg ? -r : r;
}
std::vector<double> Exact::expansion() const {
    std::vector<double> out;
    Exact v = *this;
    while (!v.zero()) {
        need(out.size() < 48, "exact expansion length");
        const double x = v.rounded();
        out.push_back(x);
        MsAcc a;
        for (int i = 0; i < MS_ACC_DIGITS; ++i) a.d[i] = v.d[i];
        ms_acc_add(&a, -x);
        v = of(a);
    }
    return out;
}
std::string Exact::hex() const {
    bool neg = false;
    const Magnitude m = magnitude(*this, neg);
    static const char* digits = "0123456789abcdef";
    std::string s;
    bool started = false;
    for (int i = MS_ACC_DIGITS - 1; i >= 0; --i)
        for (int k = 7; k >= 0; --k) {
            const int nib = int((m.d[i] >> (4 * k)) & 15u);
            if (nib || started) { s.push_back(digits[nib]); started = true; }
        }
    if (!started) return "0";
    return (neg ? "-0x" : "0x") + s;
}

// ------------------------------------------------------------------ parameters
void validate(const MsParams& P) {
    need(P.n_species >= 1 && P.n_species <= MS_MAX_SPECIES, "species count");
    need(P.nsub <= 16, "direction subdivision (0: no direction bins)");
    for (int d = 0; d < 3; ++d) need(P.sub[d] >= 1 && P.sub[d] <= 64, "spatial sub-bins per cell axis");
    need(P.nx >= 1 && P.ny >= 1 && P.nz >= 1 && P.kz_lo < P.kz_hi && P.kz_hi <= P.nz, "mesh counts");
    const std::uint64_t cells = std::uint64_t(P.nx) * P.ny * P.nz;
    need(cells * P.n_species < (std::uint64_t(1) << 32), "cell x species key range");
    need(P.cell_bits >= 1 && P.cell_bits <= 32 && (cells - 1) < (std::uint64_t(1) << P.cell_bits), "cell key bits");
    for (int d = 0; d < 3; ++d)
        need(std::isfinite(P.shift[d]) && std::isfinite(P.origin[d]) && std::isfinite(P.h[d]) && P.h[d] > 0, "mesh geometry");
    need(P.step < (std::uint64_t(1) << 40), "step outside the Philox layout domain");
    need(P.base_id > 0, "base native ID must be positive");
    for (unsigned s = 0; s < P.n_species; ++s) {
        const auto& sp = P.species[s];
        need(std::isfinite(sp.mass_kg) && sp.mass_kg > 0 && std::isfinite(sp.charge_C), "species mass/charge");
        need(sp.mc2_J == sp.mass_kg * MS_C_LIGHT * MS_C_LIGHT, "species m c^2 must be fl(fl(m c) c)");
        need(sp.n_thresh <= MS_MAX_THRESH, "energy bin count");
        for (unsigned i = 0; i < sp.n_thresh; ++i)
            need(std::isfinite(sp.e_thresh[i]) && sp.e_thresh[i] > 0 && (i == 0 || sp.e_thresh[i] > sp.e_thresh[i - 1]),
                 "energy bin edges ascending, positive, finite");
        need(std::uint64_t(P.sub[0]) * P.sub[1] * P.sub[2] * (sp.n_thresh + 1) * (P.nsub ? 6 * P.nsub * P.nsub : 1) <
                 (std::uint64_t(1) << 24),
             "bin index range");
        need(sp.k_min >= 3 && sp.k_max <= MS_MAX_GROUP && sp.k_min <= sp.k_max, "group size bounds");
        need(sp.merge_target >= 2 && sp.merge_target <= sp.merge_above, "merge target");
        need(std::uint64_t(sp.split_below) <= std::uint64_t(sp.merge_above) + 1, "split and merge ranges overlap");
        need(sp.split_levels >= 1 && sp.split_levels <= MS_MAX_SPLIT_LEVELS, "split levels");
        need(std::isfinite(sp.split_min_weight) && sp.split_min_weight > 0, "split minimum weight");
    }
}

std::string sha256_hex(const std::string& s) {
    unsigned char md[EVP_MAX_MD_SIZE];
    unsigned len = 0;
    need(EVP_Digest(s.data(), s.size(), md, &len, EVP_sha256(), nullptr) == 1 && len == 32, "sha256");
    static const char* hexd = "0123456789abcdef";
    std::string out;
    for (unsigned i = 0; i < len; ++i) { out.push_back(hexd[md[i] >> 4]); out.push_back(hexd[md[i] & 15]); }
    return out;
}

std::string params_sha256(const MsParams& P) {
    Encoder e;
    e.text("merge-split-transaction-params-v1");
    e.u64(P.seed);
    e.u64(P.step);
    e.u64(P.base_id);
    encode_policy(e, P);
    return sha256_hex(e.s);
}
std::string policy_sha256(const MsParams& P) {
    Encoder e;
    encode_policy(e, P);
    return sha256_hex(e.s);
}

MsParams hall_mesh_params(int n, int axial, double side, int air, double depth) {
    need(n >= 9 && n <= 500 && axial >= 2 && air >= 1 && side > .0054, "Hall mesh specification");
    MsParams P{};
    const double dx = side / n, dy = side / n, dz = depth / axial;
    const double zmin = -depth / 2 - air * dz;
    P.nx = unsigned(n);
    P.ny = unsigned(n);
    P.nz = unsigned(axial + 2 * air);
    P.kz_lo = unsigned(air);
    P.kz_hi = unsigned(air + axial);
    const std::uint64_t cells = std::uint64_t(P.nx) * P.ny * P.nz;
    P.cell_bits = 1;
    while ((cells - 1) >> P.cell_bits) ++P.cell_bits;
    P.shift[0] = -side / 2;
    P.shift[1] = -side / 2;
    P.shift[2] = 0;
    P.origin[0] = -side / 2;
    P.origin[1] = -side / 2;
    P.origin[2] = zmin;
    P.h[0] = dx;
    P.h[1] = dy;
    P.h[2] = dz;
    P.nsub = 1;
    P.sub[0] = P.sub[1] = P.sub[2] = 1;
    return P;
}
void set_species(MsParams& P, unsigned s, double mass_kg, double charge_C) {
    need(s < MS_MAX_SPECIES, "species index");
    auto& sp = P.species[s];
    sp = MsSpecies{};
    sp.mass_kg = mass_kg;
    sp.charge_C = charge_C;
    sp.mc2_J = mass_kg * MS_C_LIGHT * MS_C_LIGHT;
    if (P.n_species < s + 1) P.n_species = s + 1;
}

// ------------------------------------------------------------------ transaction
Result transact(const MsParams& P, const Carriers& C, const PartitionPlan& plan) {
    validate(P);
    const std::size_t n = C.size();
    need(C.species.size() == n && C.flags.size() == n && C.w.size() == n, "carrier array sizes");
    for (int d = 0; d < 3; ++d) need(C.x[d].size() == n && C.u[d].size() == n, "carrier array sizes");
    need(n < 0xFFFFFFF0u, "carrier count exceeds the 32-bit index range");
    need(plan.parts >= 1 && plan.parts <= 64, "partition count");
    const unsigned parts = plan.parts;
    Result R;
    R.params_sha256 = params_sha256(P);
    R.counters.carriers = n;
    {
        std::vector<std::uint64_t> ids(C.id);
        std::sort(ids.begin(), ids.end());
        for (std::size_t i = 0; i < n; ++i) need(ids[i] != 0 && (i == 0 || ids[i] != ids[i - 1]), "native IDs must be distinct and nonzero");
    }
    // 1. partitions: contiguous ID blocks, or hashed owners with shuffled local order
    std::vector<std::vector<std::uint32_t>> part(parts);
    std::vector<std::uint8_t> owner(n);
    if (!plan.shuffled) {
        std::vector<std::uint32_t> by_id(n);
        std::iota(by_id.begin(), by_id.end(), 0u);
        std::sort(by_id.begin(), by_id.end(), [&](std::uint32_t a, std::uint32_t b) { return C.id[a] < C.id[b]; });
        for (std::size_t r = 0; r < n; ++r) {
            const unsigned p = unsigned((std::uint64_t(r) * parts) / std::max<std::size_t>(n, 1));
            part[p].push_back(by_id[r]);
            owner[by_id[r]] = std::uint8_t(p);
        }
    } else {
        for (std::size_t i = 0; i < n; ++i) {
            const unsigned p = unsigned(splitmix64(C.id[i] ^ plan.shuffle_seed) % parts);
            part[p].push_back(std::uint32_t(i));
            owner[i] = std::uint8_t(p);
        }
        for (unsigned p = 0; p < parts; ++p) {
            std::mt19937_64 rng(plan.shuffle_seed + 977 * p + 1);
            std::shuffle(part[p].begin(), part[p].end(), rng);
        }
    }
    // 2. local classification and local canonical sort
    std::vector<std::vector<Rec>> local(parts);
    for (unsigned p = 0; p < parts; ++p) {
        local[p].reserve(part[p].size());
        for (auto i : part[p]) {
            Rec r{};
            r.index = i;
            r.id = C.id[i];
            r.species = C.species[i];
            need(ms_classify(&P, C.species[i], C.x[0][i], C.x[1][i], C.x[2][i], C.u[0][i], C.u[1][i], C.u[2][i], C.w[i],
                             &r.cell, &r.bin, &r.ke) == 1,
                 "invalid census carrier (species, finiteness, weight or position outside the kinetic cells)");
            local[p].push_back(r);
        }
        std::sort(local[p].begin(), local[p].end(), key_less);
    }
    // 3. k-way merge of the sorted runs (the coordinator's serial order)
    std::vector<Rec> canon;
    canon.reserve(n);
    {
        auto greater = [&](const std::pair<unsigned, std::size_t>& a, const std::pair<unsigned, std::size_t>& b) {
            return key_less(local[b.first][b.second], local[a.first][a.second]);
        };
        std::priority_queue<std::pair<unsigned, std::size_t>, std::vector<std::pair<unsigned, std::size_t>>, decltype(greater)> heap(greater);
        for (unsigned p = 0; p < parts; ++p)
            if (!local[p].empty()) heap.push({p, 0});
        while (!heap.empty()) {
            auto [p, k] = heap.top();
            heap.pop();
            canon.push_back(local[p][k]);
            if (k + 1 < local[p].size()) heap.push({p, k + 1});
        }
        for (std::size_t i = 1; i < canon.size(); ++i) need(key_less(canon[i - 1], canon[i]), "canonical order");
    }
    local.clear();
    local.shrink_to_fit();
    // 4. sorted structure of arrays
    std::vector<std::uint64_t> sid(n);
    std::vector<std::uint32_t> ssp(n), sfl(n), sbin(n);
    std::array<std::vector<double>, 3> sx, su;
    std::vector<double> sw(n), ske(n);
    for (int d = 0; d < 3; ++d) { sx[d].resize(n); su[d].resize(n); }
    for (std::size_t s = 0; s < n; ++s) {
        const auto i = canon[s].index;
        sid[s] = C.id[i];
        ssp[s] = C.species[i];
        sfl[s] = C.flags[i];
        sbin[s] = canon[s].bin;
        for (int d = 0; d < 3; ++d) { sx[d][s] = C.x[d][i]; su[d][s] = C.u[d][i]; }
        sw[s] = C.w[i];
        ske[s] = canon[s].ke;
    }
    MsView v{sid.data(), ssp.data(), sfl.data(), sbin.data(), {sx[0].data(), sx[1].data(), sx[2].data()},
             {su[0].data(), su[1].data(), su[2].data()}, sw.data(), ske.data()};
    // 5. (cell, species) segments and group formation
    std::vector<MsGroup> tmp;
    for (std::size_t b = 0; b < n;) {
        std::size_t e = b + 1;
        while (e < n && canon[e].cell == canon[b].cell && canon[e].species == canon[b].species) ++e;
        ++R.counters.segments;
        tmp.resize(e - b);
        std::uint32_t refused = 0;
        const auto count = ms_form_segment(&P, v, std::uint32_t(b), std::uint32_t(e), canon[b].cell, canon[b].species, tmp.data(), &refused);
        R.groups.insert(R.groups.end(), tmp.begin(), tmp.begin() + count);
        R.counters.formation_refused += refused;
        b = e;
    }
    const std::size_t G = R.groups.size();
    R.counters.groups = G;
    R.merge.assign(G, MsMergeResult{});
    R.group_first_id.resize(G);
    // 6. kinematics; group g is computed by partition g mod parts (pure function of its inputs)
    for (unsigned p = 0; p < parts; ++p)
        for (std::size_t g = p; g < G; g += parts) {
            auto& gr = R.groups[g];
            std::uint32_t gf = 0;
            if (gr.kind == MS_KIND_MERGE) {
                gr.status = ms_merge_group(&P, v, &gr, &R.merge[g], &gf);
                gr.n_out = gr.status == MS_OK ? 2u : 0u;
            } else {
                gr.status = ms_split_group(&P, v, &gr);
                gr.n_out = gr.status == MS_OK ? (1u << P.species[gr.species].split_levels) : 0u;
            }
            gr.gflags = gf;
            if (gr.status != MS_OK) R.merge[g] = MsMergeResult{};
        }
    // 7. exclusive scans: output ranks and relocation record ranks
    std::uint64_t outs = 0, relocs = 0;
    for (std::size_t g = 0; g < G; ++g) {
        auto& gr = R.groups[g];
        need(outs <= 0xFFFFFFFFull && relocs <= 0xFFFFFFFFull, "output count exceeds the 32-bit index range");
        gr.out_offset = std::uint32_t(outs);
        gr.reloc_offset = std::uint32_t(relocs);
        outs += gr.n_out;
        if (gr.kind == MS_KIND_MERGE && gr.status == MS_OK) relocs += gr.k;
        R.group_first_id[g] = sid[gr.start];
    }
    detail::fill_group_counters(R, P);
    // 8. emission
    R.outputs.reserve(outs);
    for (int d = 0; d < 3; ++d) { R.reloc.start[d].reserve(relocs); R.reloc.D[d].reserve(relocs); }
    std::vector<std::uint8_t> consumed(n, 0);
    for (std::size_t g = 0; g < G; ++g) {
        const auto& gr = R.groups[g];
        if (gr.status != MS_OK) continue;
        const auto& sp = P.species[gr.species];
        const std::uint64_t base = P.base_id + gr.out_offset;
        if (gr.kind == MS_KIND_MERGE) {
            const auto& m = R.merge[g];
            const double Wh = gr.W * 0.5;
            R.outputs.push(base, gr.species, MS_FLAG_MERGE | MS_FLAG_SPLIT, m.c, m.ua, Wh);
            R.outputs.push(base + 1, gr.species, MS_FLAG_MERGE | MS_FLAG_SPLIT, m.c, m.ub, Wh);
            double cf[3];
            for (int d = 0; d < 3; ++d) cf[d] = m.c[d] + P.shift[d];
            for (std::uint32_t s = gr.start; s < gr.end; ++s) {
                if (!(sfl[s] & MS_FLAG_MERGE)) continue;
                consumed[s] = 1;
                R.reloc.id.push_back(sid[s]);
                R.reloc.group.push_back(std::uint32_t(g));
                for (int d = 0; d < 3; ++d) {
                    const double st = sx[d][s] + P.shift[d];
                    R.reloc.start[d].push_back(st);
                    R.reloc.D[d].push_back(cf[d] - st);
                }
                R.reloc.charge.push_back(sw[s] * sp.charge_C);
            }
        } else {
            const std::uint32_t s = gr.start;
            consumed[s] = 1;
            const std::uint32_t clones = 1u << sp.split_levels;
            const double wm = sw[s] * (1.0 / double(clones));
            const double xx[3] = {sx[0][s], sx[1][s], sx[2][s]}, uu[3] = {su[0][s], su[1][s], su[2][s]};
            for (std::uint32_t j = 0; j < clones; ++j) R.outputs.push(base + j, gr.species, MS_FLAG_MERGE | MS_FLAG_SPLIT, xx, uu, wm);
        }
    }
    need(R.outputs.size() == outs && R.reloc.size() == relocs, "emission count");
    for (std::size_t s = 0; s < n; ++s)
        if (consumed[s]) { R.removed.push_back(sid[s]); ++R.counters.inputs_consumed; }
    std::sort(R.removed.begin(), R.removed.end());
    // 9. exact ledger: partial accumulators per partition, combined exactly
    const unsigned S = P.n_species;
    std::vector<MsAcc> partial(std::size_t(parts) * S * NQ);
    for (auto& a : partial) ms_acc_zero(&a);
    auto acc = [&](unsigned p, unsigned s) { return &partial[(std::size_t(p) * S + s) * NQ]; };
    for (std::size_t s = 0; s < n; ++s) {
        if (!consumed[s]) continue;
        const double uu[3] = {su[0][s], su[1][s], su[2][s]};
        add_terms(acc(owner[canon[s].index], ssp[s]), false, sw[s], uu, ske[s], P.species[ssp[s]].mass_kg);
    }
    for (std::size_t g = 0; g < G; ++g) {
        const auto& gr = R.groups[g];
        if (gr.status != MS_OK) continue;
        const double m = P.species[gr.species].mass_kg;
        for (std::uint32_t j = 0; j < gr.n_out; ++j) {
            const std::size_t o = gr.out_offset + j;
            const double uu[3] = {R.outputs.u[0][o], R.outputs.u[1][o], R.outputs.u[2][o]};
            add_terms(acc(unsigned(g % parts), gr.species), true, R.outputs.w[o], uu, ms_kinetic_energy(m, uu[0], uu[1], uu[2]), m);
        }
    }
    std::vector<Exact> totals(std::size_t(S) * NQ);
    for (unsigned s = 0; s < S; ++s)
        for (int q = 0; q < NQ; ++q) {
            MsAcc a;
            ms_acc_zero(&a);
            for (unsigned p = 0; p < parts; ++p) {
                const MsAcc& b = acc(p, s)[q];
                for (int i = 0; i < MS_ACC_DIGITS; ++i) a.d[i] += b.d[i];
            }
            totals[std::size_t(s) * NQ + q] = Exact::of(a);
        }
    detail::fill_ledger(R, P, totals);
    R.digest_sha256 = sha256_hex(encode(R));
    return R;
}

namespace detail {
void fill_group_counters(Result& R, const MsParams& P) {
    auto& c = R.counters;
    c.groups = R.groups.size();
    c.merge_groups = c.split_groups = c.accepted_merge = c.accepted_split = c.centroid_fallback = c.isotropic = 0;
    c.status = {};
    std::uint64_t outs = 0, relocs = 0;
    for (const auto& gr : R.groups) {
        if (gr.kind == MS_KIND_MERGE) ++c.merge_groups;
        else ++c.split_groups;
        need(gr.status < MS_STATUS_KINDS, "group status");
        ++c.status[gr.status];
        if (gr.status == MS_OK) (gr.kind == MS_KIND_MERGE ? c.accepted_merge : c.accepted_split)++;
        if (gr.gflags & MS_GFLAG_CENTROID_FALLBACK) ++c.centroid_fallback;
        if (gr.gflags & MS_GFLAG_ISOTROPIC) ++c.isotropic;
        outs += gr.n_out;
        if (gr.kind == MS_KIND_MERGE && gr.status == MS_OK) relocs += gr.k;
    }
    need(outs <= ~std::uint64_t(0) - P.base_id, "native ID exhaustion");
    R.next_id_after = P.base_id + outs;
    c.outputs = outs;
    c.relocations = relocs;
}

void fill_ledger(Result& R, const MsParams& P, const std::vector<Exact>& totals) {
    const unsigned S = P.n_species;
    need(totals.size() == std::size_t(S) * NQ, "ledger totals shape");
    R.ledger.assign(S, SpeciesLedger{});
    for (unsigned s = 0; s < S; ++s) {
        const Exact* tot = &totals[std::size_t(s) * NQ];
        auto& L = R.ledger[s];
        L.w_in = tot[W_IN];
        L.w_out = tot[W_OUT];
        L.e_in = tot[E_IN];
        L.e_out = tot[E_OUT];
        L.e_res = L.e_out - L.e_in;
        L.e_res_rn = L.e_res.rounded();
        L.e_res_expansion = L.e_res.expansion();
        for (int d = 0; d < 3; ++d) {
            L.p_in[d] = tot[P_IN + d];
            L.p_out[d] = tot[P_OUT + d];
            L.p_res[d] = L.p_out[d] - L.p_in[d];
            L.p_res_rn[d] = L.p_res[d].rounded();
            L.p_res_expansion[d] = L.p_res[d].expansion();
        }
        L.count_equal = L.w_in == L.w_out;
        need(L.count_equal, "exact particle count changed by the transaction");
    }
}
}  // namespace detail

std::string encode(const Result& R) {
    Encoder e;
    e.text("merge-split-transaction-v1");
    e.text(R.params_sha256);
    e.u64(R.next_id_after);
    const auto& c = R.counters;
    for (auto x : {c.carriers, c.segments, c.groups, c.merge_groups, c.split_groups, c.accepted_merge, c.accepted_split,
                   c.formation_refused, c.centroid_fallback, c.isotropic, c.inputs_consumed, c.outputs, c.relocations})
        e.u64(x);
    for (auto x : c.status) e.u64(x);
    e.text("groups");
    e.u64(R.groups.size());
    for (std::size_t g = 0; g < R.groups.size(); ++g) {
        const auto& gr = R.groups[g];
        e.f64(gr.W);
        for (auto x : {gr.start, gr.end, gr.k, gr.kind, gr.cell, gr.species, gr.bin, gr.index, gr.status, gr.n_out,
                       gr.out_offset, gr.reloc_offset, gr.gflags})
            e.u32(x);
        e.u64(R.group_first_id[g]);
        const auto& m = R.merge[g];
        for (int d = 0; d < 3; ++d) e.f64(m.c[d]);
        for (int d = 0; d < 3; ++d) e.f64(m.ua[d]);
        for (int d = 0; d < 3; ++d) e.f64(m.ub[d]);
    }
    e.text("outputs");
    const auto& o = R.outputs;
    e.u64(o.size());
    for (std::size_t i = 0; i < o.size(); ++i) {
        e.u64(o.id[i]);
        e.u32(o.species[i]);
        e.u32(o.flags[i]);
        for (int d = 0; d < 3; ++d) e.f64(o.x[d][i]);
        for (int d = 0; d < 3; ++d) e.f64(o.u[d][i]);
        e.f64(o.w[i]);
    }
    e.text("removed");
    e.u64(R.removed.size());
    for (auto x : R.removed) e.u64(x);
    e.text("relocation");
    const auto& r = R.reloc;
    e.u64(r.size());
    for (std::size_t i = 0; i < r.size(); ++i) {
        e.u64(r.id[i]);
        e.u32(r.group[i]);
        for (int d = 0; d < 3; ++d) e.f64(r.start[d][i]);
        for (int d = 0; d < 3; ++d) e.f64(r.D[d][i]);
        e.f64(r.charge[i]);
    }
    e.text("ledger");
    e.u64(R.ledger.size());
    for (const auto& L : R.ledger) {
        e.exact(L.w_in);
        e.exact(L.w_out);
        for (int d = 0; d < 3; ++d) e.exact(L.p_in[d]);
        for (int d = 0; d < 3; ++d) e.exact(L.p_out[d]);
        e.exact(L.e_in);
        e.exact(L.e_out);
        for (int d = 0; d < 3; ++d) e.exact(L.p_res[d]);
        e.exact(L.e_res);
        for (int d = 0; d < 3; ++d) e.f64(L.p_res_rn[d]);
        e.f64(L.e_res_rn);
        e.u32(L.count_equal ? 1u : 0u);
    }
    return e.s;
}

Carriers apply(const Carriers& C, const Result& R) {
    Carriers out;
    out.reserve(C.size() - R.removed.size() + R.outputs.size());
    for (std::size_t i = 0; i < C.size(); ++i) {
        if (std::binary_search(R.removed.begin(), R.removed.end(), C.id[i])) continue;
        const double xx[3] = {C.x[0][i], C.x[1][i], C.x[2][i]}, uu[3] = {C.u[0][i], C.u[1][i], C.u[2][i]};
        out.push(C.id[i], C.species[i], C.flags[i], xx, uu, C.w[i]);
    }
    const auto& o = R.outputs;
    for (std::size_t i = 0; i < o.size(); ++i) {
        const double xx[3] = {o.x[0][i], o.x[1][i], o.x[2][i]}, uu[3] = {o.u[0][i], o.u[1][i], o.u[2][i]};
        out.push(o.id[i], o.species[i], o.flags[i], xx, uu, o.w[i]);
    }
    return out;
}

std::vector<Exact> exact_weights(const MsParams& P, const Carriers& C) {
    std::vector<MsAcc> a(P.n_species);
    for (auto& x : a) ms_acc_zero(&x);
    for (std::size_t i = 0; i < C.size(); ++i) {
        need(C.species[i] < P.n_species, "species index");
        ms_acc_add(&a[C.species[i]], C.w[i]);
    }
    std::vector<Exact> out;
    for (auto& x : a) out.push_back(Exact::of(x));
    return out;
}

}  // namespace merge_split
