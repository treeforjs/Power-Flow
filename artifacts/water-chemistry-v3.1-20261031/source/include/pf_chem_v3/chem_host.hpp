// chem_host.hpp -- host side of chemistry v3: table blob loading, exact
// ledgers, the CPU reference step, canonical record order and digests,
// tracked-target debit resolution and the quantized neutral release.
//
// Exact arithmetic. ExactSum holds an integer multiple of 2^-1074 in 36
// 64-bit two's-complement limbs (range beyond +-2^1100), so sums of binary64
// values and exact products of two binary64 values (Dekker two-product, valid
// while no partial product under- or overflows; checked) are exact and
// independent of summation order. Ledger totals are therefore identical for
// any partition or processing order of the same records.
#pragma once
#include "chem_event.h"

#include <openssl/evp.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace pfc {

static_assert(sizeof(SpeciesRec) == 64, "SpeciesRec layout");
static_assert(sizeof(GroupRec) == 64, "GroupRec layout");
static_assert(sizeof(ProcessRec) == 512, "ProcessRec layout");
static_assert(sizeof(EventRec) == 384, "EventRec layout");
static_assert(sizeof(CarrierIn) == 48, "CarrierIn layout");
static_assert(sizeof(CellTargets) == 536, "CellTargets layout");
static_assert(sizeof(StepConfig) == 48, "StepConfig layout");
static_assert(sizeof(CarrierOut) == 384, "CarrierOut layout");

struct Error : std::runtime_error {
    using std::runtime_error::runtime_error;
};
inline void require(bool ok, const std::string& what) {
    if (!ok) throw Error(what);
}

// ------------------------------------------------------------------ SHA-256
inline std::string sha256_hex(const void* data, std::size_t n) {
    unsigned char md[32];
    unsigned int len = 0;
    EVP_MD_CTX* c = EVP_MD_CTX_new();
    require(c && EVP_DigestInit_ex(c, EVP_sha256(), nullptr) == 1 && EVP_DigestUpdate(c, data, n) == 1 &&
                EVP_DigestFinal_ex(c, md, &len) == 1,
            "sha256");
    EVP_MD_CTX_free(c);
    static const char* hx = "0123456789abcdef";
    std::string s;
    for (unsigned i = 0; i < len; ++i) {
        s += hx[md[i] >> 4];
        s += hx[md[i] & 15];
    }
    return s;
}
class Sha256 {
public:
    Sha256() : c_(EVP_MD_CTX_new()) { require(c_ && EVP_DigestInit_ex(c_, EVP_sha256(), nullptr) == 1, "sha256 init"); }
    ~Sha256() { EVP_MD_CTX_free(c_); }
    Sha256(const Sha256&) = delete;
    Sha256& operator=(const Sha256&) = delete;
    void update(const void* d, std::size_t n) { require(EVP_DigestUpdate(c_, d, n) == 1, "sha256 update"); }
    template <class T> void pod(const T& v) { update(&v, sizeof(T)); }
    void text(const std::string& s) {
        const std::uint64_t n = s.size();
        pod(n);
        update(s.data(), s.size());
    }
    std::string hex() {
        unsigned char md[32];
        unsigned int len = 0;
        require(EVP_DigestFinal_ex(c_, md, &len) == 1, "sha256 final");
        static const char* hx = "0123456789abcdef";
        std::string s;
        for (unsigned i = 0; i < len; ++i) {
            s += hx[md[i] >> 4];
            s += hx[md[i] & 15];
        }
        return s;
    }

private:
    EVP_MD_CTX* c_;
};

// ------------------------------------------------------------------ ExactSum
class ExactSum {
public:
    static constexpr int kLimbs = 36;   // bit 0 = 2^-1074
    ExactSum() { limb_.fill(0); }
    void add(double x) {
        require(std::isfinite(x), "ExactSum: nonfinite term");
        if (x == 0.0) return;
        const std::uint64_t b = std::bit_cast<std::uint64_t>(x);
        const bool neg = b >> 63;
        const int ex = int((b >> 52) & 0x7FF);
        std::uint64_t m = b & 0xFFFFFFFFFFFFFull;
        int pos;
        if (ex) {
            m |= 1ull << 52;
            pos = ex - 1;   // x = m * 2^(ex - 1075) = m * 2^(pos - 1074)
        } else {
            pos = 0;
        }
        add_shifted(m, pos, neg);
    }
    // Exact a*b (Dekker); refuses if a partial product could have underflowed.
    void add_product(double a, double b) {
        if (a == 0.0 || b == 0.0) return;
        const double p = a * b;
        require(std::isfinite(p), "ExactSum: product overflow");
        require(std::fabs(p) > 0x1p-900, "ExactSum: product too small for an exact two-product");
        double ah, al, bh, bl;
        split(a, ah, al);
        split(b, bh, bl);
        const double err = ((ah * bh - p) + ah * bl + al * bh) + al * bl;
        add(p);
        add(err);
    }
    // a + b + carry with carry out (portable add-with-carry)
    static std::uint64_t adc(std::uint64_t a, std::uint64_t b, std::uint64_t& carry) {
        const std::uint64_t s1 = a + b;
        const std::uint64_t c1 = s1 < a ? 1u : 0u;
        const std::uint64_t s2 = s1 + carry;
        const std::uint64_t c2 = s2 < s1 ? 1u : 0u;
        carry = c1 | c2;
        return s2;
    }
    void add(const ExactSum& o) {
        std::uint64_t carry = 0;
        for (int i = 0; i < kLimbs; ++i) limb_[i] = adc(limb_[i], o.limb_[i], carry);
    }
    void negate() {
        for (auto& l : limb_) l = ~l;
        std::uint64_t carry = 1;
        for (int i = 0; i < kLimbs; ++i) limb_[i] = adc(limb_[i], 0, carry);
    }
    void sub(const ExactSum& o) {
        ExactSum n = o;
        n.negate();
        add(n);
    }
    bool is_zero() const {
        for (auto l : limb_)
            if (l) return false;
        return true;
    }
    bool negative() const { return limb_[kLimbs - 1] >> 63; }
    bool operator==(const ExactSum& o) const { return limb_ == o.limb_; }
    // Nearest binary64 (within one ulp; exact when representable).
    double approx() const {
        ExactSum a = *this;
        const bool neg = a.negative();
        if (neg) a.negate();
        double r = 0.0;
        for (int i = kLimbs - 1; i >= 0; --i) {
            if (!a.limb_[i]) continue;
            const double hi = std::ldexp(double(a.limb_[i] >> 32), 64 * i + 32 - 1074);
            const double lo = std::ldexp(double(a.limb_[i] & 0xFFFFFFFFull), 64 * i - 1074);
            r += hi;
            r += lo;
        }
        return neg ? -r : r;
    }
    // Canonical hexadecimal encoding (two's complement, most significant limb first).
    std::string hex() const {
        static const char* hx = "0123456789abcdef";
        std::string s;
        for (int i = kLimbs - 1; i >= 0; --i)
            for (int k = 60; k >= 0; k -= 4) s += hx[(limb_[i] >> k) & 15];
        return s;
    }
    const std::array<std::uint64_t, kLimbs>& limbs() const { return limb_; }

private:
    static void split(double a, double& hi, double& lo) {
        const double c = 134217729.0 * a;   // 2^27 + 1
        hi = c - (c - a);
        lo = a - hi;
    }
    void add_shifted(std::uint64_t m, int pos, bool neg) {
        const int li = pos >> 6, sh = pos & 63;
        std::array<std::uint64_t, 2> v{m << sh, sh ? (m >> (64 - sh)) : 0};
        if (!neg) {
            std::uint64_t carry = 0;
            for (int i = li; i < kLimbs; ++i) {
                const std::uint64_t add = (i - li) < 2 ? v[i - li] : 0;
                if (add == 0 && carry == 0 && i - li >= 2) break;
                limb_[i] = adc(limb_[i], add, carry);
            }
        } else {
            std::uint64_t borrow = 0;
            for (int i = li; i < kLimbs; ++i) {
                const std::uint64_t sv = (i - li) < 2 ? v[i - li] : 0;
                if (sv == 0 && borrow == 0 && i - li >= 2) break;
                const std::uint64_t a = limb_[i];
                const std::uint64_t d = a - sv - borrow;
                borrow = (a < sv || (a - sv) < borrow) ? 1 : 0;
                limb_[i] = d;
            }
        }
    }
    std::array<std::uint64_t, kLimbs> limb_;
};

// ------------------------------------------------------------------- Blob
struct Blob {
    std::vector<std::uint64_t> storage;   // 8-byte aligned bytes
    std::size_t bytes = 0;
    std::string sha256;
    struct Section {
        std::uint64_t offset = 0, size = 0;
    };
    std::map<std::string, Section> sections;
    std::vector<std::string> species_names, process_names, process_tags, region_names;
    TableView view{};
    CoverView cover{};
    pfc_u32 n_species = 0, n_groups = 0, n_proc = 0, n_regions = 0;
    const unsigned char* base() const { return reinterpret_cast<const unsigned char*>(storage.data()); }
    template <class T> const T* at(const std::string& name, std::size_t* count = nullptr) const {
        auto it = sections.find(name);
        require(it != sections.end(), "blob section missing: " + name);
        require(it->second.offset % 8 == 0 && it->second.offset + it->second.size <= bytes, "blob section bounds: " + name);
        require(it->second.size % sizeof(T) == 0, "blob section size: " + name);
        if (count) *count = it->second.size / sizeof(T);
        return reinterpret_cast<const T*>(base() + it->second.offset);
    }
    // Device views: pointers rebased onto a device copy of the same bytes.
    TableView view_at(std::uint64_t device_base) const {
        TableView v = view;
        auto rebase = [&](const void* p) { return device_base + (std::uint64_t)((const unsigned char*)p - base()); };
        v.species = (const SpeciesRec*)rebase(view.species);
        v.groups = (const GroupRec*)rebase(view.groups);
        v.procs = (const ProcessRec*)rebase(view.procs);
        v.grid = (const double*)rebase(view.grid);
        v.sigma = (const double*)rebase(view.sigma);
        v.bound = (const double*)rebase(view.bound);
        v.rmq = (const double*)rebase(view.rmq);
        v.xi = (const double*)rebase(view.xi);
        v.group_of = (const pfc_u32*)rebase(view.group_of);
        return v;
    }
    CoverView cover_at(std::uint64_t device_base) const {
        CoverView c = cover;
        auto rebase = [&](const void* p) { return device_base + (std::uint64_t)((const unsigned char*)p - base()); };
        c.bin_region = (const pfc_u32*)rebase(cover.bin_region);
        c.bin_row = (const pfc_u32*)rebase(cover.bin_row);
        c.bin_off = (const pfc_u32*)rebase(cover.bin_off);
        return c;
    }
};

inline std::vector<std::string> split_strings(const char* p, std::size_t n) {
    std::vector<std::string> out;
    std::string cur;
    for (std::size_t i = 0; i < n; ++i) {
        if (p[i] == '\n') {
            out.push_back(cur);
            cur.clear();
        } else cur += p[i];
    }
    return out;
}

// Blob format v1: "PFCHEMV3" | u32 version=1 | u32 n_sections | n x {char name[24];
// u64 offset; u64 size} | data. Offsets are multiples of 64.
inline Blob load_blob(const std::string& path, const std::string& expected_sha256 = "") {
    std::ifstream f(path, std::ios::binary);
    require(bool(f), "cannot open table blob " + path);
    f.seekg(0, std::ios::end);
    const std::size_t n = std::size_t(f.tellg());
    f.seekg(0);
    Blob b;
    b.bytes = n;
    b.storage.assign((n + 7) / 8, 0);
    f.read(reinterpret_cast<char*>(b.storage.data()), std::streamsize(n));
    require(bool(f), "short read of table blob");
    b.sha256 = sha256_hex(b.base(), n);
    if (!expected_sha256.empty()) require(b.sha256 == expected_sha256, "table blob sha256 mismatch: " + b.sha256);
    require(n >= 16 && std::memcmp(b.base(), "PFCHEMV3", 8) == 0, "table blob magic");
    std::uint32_t version, ns;
    std::memcpy(&version, b.base() + 8, 4);
    std::memcpy(&ns, b.base() + 12, 4);
    require(version == 1 && ns < 64 && 16 + ns * 40ull <= n, "table blob header");
    for (std::uint32_t i = 0; i < ns; ++i) {
        const unsigned char* e = b.base() + 16 + 40ull * i;
        char name[25] = {0};
        std::memcpy(name, e, 24);
        Blob::Section s;
        std::memcpy(&s.offset, e + 24, 8);
        std::memcpy(&s.size, e + 32, 8);
        require(s.offset % 64 == 0 && s.offset <= n && s.size <= n - s.offset, std::string("section bounds ") + name);
        b.sections[name] = s;
    }
    std::size_t c;
    const pfc_u32* meta = b.at<pfc_u32>("meta", &c);
    require(c >= 4, "meta");
    b.n_species = meta[0];
    b.n_groups = meta[1];
    b.n_proc = meta[2];
    b.n_regions = meta[3];
    TableView& v = b.view;
    v.n_species = b.n_species;
    v.n_groups = b.n_groups;
    v.n_proc = b.n_proc;
    v.n_regions = b.n_regions;
    std::size_t ns_, ng, np, ngrid, nsig, nb, nr, nxi, ngo;
    v.species = b.at<SpeciesRec>("species", &ns_);
    v.groups = b.at<GroupRec>("groups", &ng);
    v.procs = b.at<ProcessRec>("procs", &np);
    v.grid = b.at<double>("grid", &ngrid);
    v.sigma = b.at<double>("sigma", &nsig);
    v.bound = b.at<double>("bound", &nb);
    v.rmq = b.at<double>("rmq", &nr);
    v.xi = b.at<double>("xi", &nxi);
    v.group_of = b.at<pfc_u32>("group_of", &ngo);
    require(ns_ == b.n_species && ng == b.n_groups && np == b.n_proc && ngo == std::size_t(b.n_species) * b.n_species,
            "table blob counts");
    std::size_t nbr, nbrow, nboff;
    b.cover.bin_region = b.at<pfc_u32>("bin_region", &nbr);
    b.cover.bin_row = b.at<pfc_u32>("bin_row", &nbrow);
    b.cover.bin_off = b.at<pfc_u32>("bin_off", &nboff);
    require(nbr == nbrow && nboff == b.n_groups, "coverage bins");
    // structural validation (every offset in bounds)
    for (pfc_u32 g = 0; g < b.n_groups; ++g) {
        const GroupRec& gr = v.groups[g];
        require(gr.projectile < b.n_species && gr.target < b.n_species, "group species");
        require(gr.grid_n >= 2 && gr.grid_off + std::size_t(gr.grid_n) <= ngrid, "group grid");
        for (pfc_u32 i = 0; i + 1 < gr.grid_n; ++i) require(v.grid[gr.grid_off + i] < v.grid[gr.grid_off + i + 1], "grid ascending");
        require(gr.proc_off + std::size_t(gr.proc_n) <= np, "group processes");
        require(gr.sigma_off + std::size_t(gr.proc_n) * gr.grid_n <= nsig, "group sigma");
        require(gr.bound_off + std::size_t(gr.bound_n) * gr.grid_n <= nb, "group bound");
        require(gr.rmq_levels >= 1 && gr.rmq_off + std::size_t(gr.rmq_levels) * (gr.grid_n - 1) <= nr &&
                    (1u << (gr.rmq_levels - 1)) <= gr.grid_n - 1,
                "group rmq");
        require(gr.xi_off == 0xFFFFFFFFu || gr.xi_off + std::size_t(gr.grid_n) <= nxi, "group xi");
        require(b.cover.bin_off[g] + std::size_t(gr.grid_n) + 1 <= nbr, "group bins");
        require(group_index(v, gr.projectile, gr.target) == g, "group map");
        for (pfc_u32 i = 0; i <= gr.grid_n; ++i) {
            const pfc_u32 row = b.cover.bin_row[b.cover.bin_off[g] + i];
            require(row == 0xFFFFFFFFu || row < gr.bound_n, "bin bound row");
            const pfc_u32 reg = b.cover.bin_region[b.cover.bin_off[g] + i];
            require(reg == 0xFFFFFFFFu || reg < b.n_regions, "bin region");
        }
    }
    for (pfc_u32 p = 0; p < b.n_proc; ++p) {
        const ProcessRec& pr = v.procs[p];
        require(pr.group < b.n_groups && pr.n_out >= 1 && pr.n_out <= PFC_MAX_OUT && pr.n_orb <= PFC_MAX_ORB && pr.n_cover <= PFC_MAX_COVER,
                "process record");
        const GroupRec& gr = v.groups[pr.group];
        require(p >= gr.proc_off && p < gr.proc_off + gr.proc_n, "process group range");
        for (pfc_u32 i = 0; i < pr.n_out; ++i) require(pr.out_species[i] < b.n_species, "process product species");
    }
    std::size_t nstr;
    const char* str = b.at<char>("strings", &nstr);
    const auto lines = split_strings(str, nstr);
    // strings: species names, then process names, process tags, region names
    require(lines.size() == std::size_t(b.n_species) + 2 * b.n_proc + b.n_regions, "strings section");
    b.species_names.assign(lines.begin(), lines.begin() + b.n_species);
    b.process_names.assign(lines.begin() + b.n_species, lines.begin() + b.n_species + b.n_proc);
    b.process_tags.assign(lines.begin() + b.n_species + b.n_proc, lines.begin() + b.n_species + 2 * b.n_proc);
    b.region_names.assign(lines.begin() + b.n_species + 2 * b.n_proc, lines.end());
    return b;
}

inline pfc_u32 species_index(const Blob& b, const std::string& name) {
    for (pfc_u32 i = 0; i < b.n_species; ++i)
        if (b.species_names[i] == name) return i;
    throw Error("unknown species " + name);
}

// ------------------------------------------------------------ CPU reference
struct StepResult {
    std::vector<CarrierOut> outs;
    std::vector<EventRec> events;   // carrier-major: carrier i owns events[i*max_events ..)
};
inline void cpu_step(const Blob& b, const StepConfig& cfg, const std::vector<CarrierIn>& carriers,
                     const std::vector<CellTargets>& cells, StepResult& r, std::size_t begin = 0,
                     std::size_t end = std::size_t(-1)) {
    end = std::min(end, carriers.size());
    r.outs.assign(carriers.size(), CarrierOut{});
    r.events.assign(carriers.size() * cfg.max_events, EventRec{});
    for (std::size_t i = begin; i < end; ++i) {
        const CarrierIn& c = carriers[i];
        if (c.cell >= cells.size()) {
            r.outs[i].id = c.id;
            r.outs[i].status = PFC_ST_BAD_INPUT;
            continue;
        }
        mcc_carrier(b.view, b.cover, cfg, c, cells[c.cell], r.events.data() + i * cfg.max_events, &r.outs[i]);
    }
}

// Canonical record stream: carriers in ascending id, each carrier's events in
// candidate order. The digest covers every byte of every record, so it is the
// bit-identity witness for CPU/GPU and partition comparisons.
inline std::string step_digest(const std::vector<CarrierOut>& outs, const std::vector<EventRec>& events, pfc_u32 max_events,
                               const std::vector<CarrierIn>& carriers) {
    std::vector<std::size_t> order(carriers.size());
    for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return carriers[a].id < carriers[b].id; });
    Sha256 h;
    h.text("pfc-chem-v3-step-records-v1");
    for (auto i : order) {
        h.pod(outs[i]);
        for (pfc_u32 k = 0; k < outs[i].n_events; ++k) h.pod(events[i * max_events + k]);
    }
    return h.hex();
}

// --------------------------------------------------------------- ledgers
struct ProcessTally {
    std::uint64_t events = 0;
    ExactSum weight;
};
struct Ledger {
    // energy (eV, physical: per-particle energy times weight, exact)
    ExactSum ke_in_projectile, ke_in_target;
    ExactSum ke_out_survivor, ke_out_birth, ke_out_reservoir, ke_out_gas;
    ExactSum q_reference;
    std::array<ExactSum, PFC_SINK_COUNT> internal;
    ExactSum defect;   // = in - out - q - internal (exact residual of rounding)
    // counts (physical particles, exact) per species
    std::vector<ExactSum> created, destroyed;
    ExactSum charge_created, charge_destroyed;   // elementary charges
    std::array<ExactSum, 6> nuclei_created, nuclei_destroyed;
    std::vector<ProcessTally> process;
    std::array<std::uint64_t, 4> events_by_cover{};
    // coverage exposure per region (exact)
    std::vector<ExactSum> region_exposure, region_bound;
    // per-event closure bound check
    std::uint64_t events = 0, closure_violations = 0;
    double max_rel_defect = 0.0;
    std::uint64_t carriers = 0, candidates = 0;
    std::map<pfc_u32, std::uint64_t> status_counts;
    double nu_dt_max = 0.0;
    void init(const Blob& b) {
        created.assign(b.n_species, ExactSum{});
        destroyed.assign(b.n_species, ExactSum{});
        process.assign(b.n_proc, ProcessTally{});
        region_exposure.assign(b.n_regions, ExactSum{});
        region_bound.assign(b.n_regions, ExactSum{});
    }
};

// Per-event exact closure: in = out + q + internal + defect_e, with
// |defect_e| <= tol(in, q). Returns the per-particle defect approximation.
inline double event_defect(const EventRec& e, ExactSum* exact = nullptr) {
    ExactSum s;
    s.add(e.ke_in_p);
    s.add(e.ke_in_t);
    for (pfc_u32 i = 0; i < e.n_out; ++i) s.add(-e.ke_out[i]);
    s.add(-e.q_eV);
    s.add(-e.e_internal);
    if (exact) *exact = s;
    return s.approx();
}
// Tolerance of the per-event closure: 2^-40 of the largest energy scale of the
// event (kinetic in, |q|, internal, and the kinetic energy of every product).
inline double event_tolerance(const EventRec& e) {
    double scale = std::fabs(e.ke_in_p) + std::fabs(e.ke_in_t) + std::fabs(e.q_eV) + std::fabs(e.e_internal);
    for (pfc_u32 i = 0; i < e.n_out; ++i) scale += std::fabs(e.ke_out[i]);
    return 0x1p-40 * scale + 1e-300;
}

inline void ledger_add_event(const Blob& b, Ledger& L, const EventRec& e) {
    const double W = e.weight;
    const ProcessRec& pr = b.view.procs[e.process];
    L.events++;
    L.process[e.process].events++;
    L.process[e.process].weight.add(W);
    if (e.cover_class < 4) L.events_by_cover[e.cover_class]++;
    L.ke_in_projectile.add_product(W, e.ke_in_p);
    L.ke_in_target.add_product(W, e.ke_in_t);
    for (pfc_u32 i = 0; i < e.n_out; ++i) {
        switch (e.out_dest[i]) {
            case PFC_DEST_SURVIVOR: L.ke_out_survivor.add_product(W, e.ke_out[i]); break;
            case PFC_DEST_BIRTH: L.ke_out_birth.add_product(W, e.ke_out[i]); break;
            case PFC_DEST_RESERVOIR: L.ke_out_reservoir.add_product(W, e.ke_out[i]); break;
            case PFC_DEST_GAS: L.ke_out_gas.add_product(W, e.ke_out[i]); break;
            default: throw Error("event product without destination");
        }
    }
    L.q_reference.add_product(W, e.q_eV);
    require(e.internal_sink < PFC_SINK_COUNT, "internal sink category");
    L.internal[e.internal_sink].add_product(W, e.e_internal);
    // per-event defect, exact and bounded
    ExactSum d;
    const double da = event_defect(e, &d);
    const double tol = event_tolerance(e);
    if (!(std::fabs(da) <= tol)) L.closure_violations++;
    double scale = std::fabs(e.ke_in_p) + std::fabs(e.ke_in_t) + std::fabs(e.q_eV) + 1e-300;
    L.max_rel_defect = std::max(L.max_rel_defect, std::fabs(da) / scale);
    // counts: projectile and target destroyed; outputs created (survivor: the
    // projectile continues -- counted as destroyed and created so that species
    // balances are uniform)
    const pfc_u32 sp = [&] {
        // projectile species = group projectile
        return b.view.groups[pr.group].projectile;
    }();
    const pfc_u32 st = e.target_species;
    auto destroy = [&](pfc_u32 s) {
        L.destroyed[s].add(W);
        const SpeciesRec& r = b.view.species[s];
        L.charge_destroyed.add_product(W, double(r.charge));
        for (int k = 0; k < 6; ++k)
            if (r.nuclei[k]) L.nuclei_destroyed[k].add_product(W, double(r.nuclei[k]));
    };
    auto create = [&](pfc_u32 s) {
        L.created[s].add(W);
        const SpeciesRec& r = b.view.species[s];
        L.charge_created.add_product(W, double(r.charge));
        for (int k = 0; k < 6; ++k)
            if (r.nuclei[k]) L.nuclei_created[k].add_product(W, double(r.nuclei[k]));
    };
    destroy(sp);
    destroy(st);
    for (pfc_u32 i = 0; i < e.n_out; ++i) create(e.out_species[i]);
}

inline void ledger_add_carrier(Ledger& L, const CarrierOut& o) {
    L.carriers++;
    L.candidates += o.n_candidates;
    L.status_counts[o.status]++;
    L.nu_dt_max = std::max(L.nu_dt_max, o.nu_dt_max);
    for (pfc_u32 i = 0; i < o.n_cover; ++i) {
        require(o.cover_region[i] < L.region_exposure.size(), "coverage region index");
        L.region_exposure[o.cover_region[i]].add(o.cover_exposure[i]);
        L.region_bound[o.cover_region[i]].add(o.cover_bound[i]);
    }
}

inline void ledger_finish(Ledger& L) {
    ExactSum d = L.ke_in_projectile;
    d.add(L.ke_in_target);
    for (auto* x : {&L.ke_out_survivor, &L.ke_out_birth, &L.ke_out_reservoir, &L.ke_out_gas, &L.q_reference}) d.sub(*x);
    for (auto& x : L.internal) d.sub(x);
    L.defect = d;
}

// Event-level charge and nuclei conservation (integer stoichiometry).
inline bool event_conserves(const Blob& b, const EventRec& e) {
    const ProcessRec& pr = b.view.procs[e.process];
    const SpeciesRec& P = b.view.species[b.view.groups[pr.group].projectile];
    const SpeciesRec& T = b.view.species[e.target_species];
    long q = P.charge + T.charge;
    long nuc[6];
    for (int k = 0; k < 6; ++k) nuc[k] = long(P.nuclei[k]) + long(T.nuclei[k]);
    for (pfc_u32 i = 0; i < e.n_out; ++i) {
        const SpeciesRec& r = b.view.species[e.out_species[i]];
        q -= r.charge;
        for (int k = 0; k < 6; ++k) nuc[k] -= long(r.nuclei[k]);
    }
    for (int k = 0; k < 6; ++k)
        if (nuc[k]) return false;
    return q == 0;
}

// -------------------------------------------------- tracked-target debits
// Events against tracked target inventories are admitted in canonical order
// (cell, carrier id, candidate). An event whose target inventory (exact count at
// step begin minus earlier admitted debits) is below its weight is refused; the
// carrier's later events in the step are then dropped as well and the carrier
// keeps its state from before the refused event. Background targets never
// deplete.
struct InventoryKey {
    pfc_u32 cell, species;
    bool operator<(const InventoryKey& o) const { return cell != o.cell ? cell < o.cell : species < o.species; }
};
struct DebitResult {
    std::vector<std::uint8_t> admitted;     // per event slot (carrier-major), 1 admitted
    std::vector<pfc_u32> admitted_events;   // per carrier: number of admitted leading events
    std::uint64_t refused = 0;
    std::map<InventoryKey, ExactSum> debited;
};
inline DebitResult resolve_debits(const Blob& b, const std::vector<CarrierIn>& carriers, const StepResult& r, pfc_u32 max_events,
                                  const std::map<InventoryKey, ExactSum>& inventory) {
    (void)b;
    DebitResult d;
    d.admitted.assign(r.events.size(), 0);
    d.admitted_events.assign(carriers.size(), 0);
    std::vector<std::size_t> order(carriers.size());
    for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t c) {
        return carriers[a].cell != carriers[c].cell ? carriers[a].cell < carriers[c].cell : carriers[a].id < carriers[c].id;
    });
    std::map<InventoryKey, ExactSum> left = inventory;   // only tracked targets appear here
    for (auto i : order) {
        pfc_u32 k = 0;
        for (; k < r.outs[i].n_events; ++k) {
            const EventRec& e = r.events[i * max_events + k];
            const InventoryKey key{e.cell, e.target_species};
            auto it = left.find(key);
            if (it != left.end()) {
                ExactSum rem = it->second;
                ExactSum w;
                w.add(e.weight);
                rem.sub(w);
                if (rem.negative()) {
                    d.refused++;
                    break;
                }
                it->second = rem;
                d.debited[key].add(e.weight);
            }
            d.admitted[i * max_events + k] = 1;
        }
        d.admitted_events[i] = k;
    }
    return d;
}

// ------------------------------------------- quantized neutral release
// Reservoir per (cell, neutral species): exact count C, exact momentum P (eV,
// physical: W p), exact kinetic energy K (eV, physical). Deposits add the
// products of admitted events. Release: x = C/w0 (binary64 of the exact count
// divided by w0), n = floor(x) + Bernoulli(x - floor(x)) on the chemistry
// release stream, capped at floor(C/w0) exactly; each released packet has
// weight w0 and the proper velocity of a product of this step chosen with
// probability proportional to its weight (keyed pick), or the reservoir mean
// velocity when the step deposited none. C, P and K are debited exactly by the
// released packets; their signed remainders stay in the reservoir.
struct ReservoirEntry {
    double weight;
    double u[3];
    double x[3];   // position of the event (the projectile carrier's position), m
};
struct Reservoir {
    ExactSum count;
    std::array<ExactSum, 3> momentum;   // eV (W * m * u / c)
    ExactSum kinetic;                   // eV (W * KE)
    std::vector<ReservoirEntry> step_entries;
    double last_x[3] = {0.0, 0.0, 0.0}; // position of the last deposit (used when a step deposits none)
};
struct ReleasedPacket {
    pfc_u32 cell, species, index;
    double weight;
    double u[3];
    double x[3];   // position of the picked entry (reservoir centroid when the step deposited none)
};
inline double word_bernoulli(const ChemKey& k, pfc_u32 event, double p) {
    struct Words {
        ChemKey k;
        pfc_u32 ev;
        pfa_u64 operator()(int j) const { return chem_word(k, ev, pfc_u32(j)); }
    } w{k, event};
    return pfa::bernoulli_digits(w, p);
}
inline std::vector<ReleasedPacket> release_reservoirs(const Blob& b, std::map<InventoryKey, Reservoir>& res, double w0,
                                                      pfc_u64 seed, pfc_u64 step) {
    require(w0 > 0.0 && std::isfinite(w0), "release weight");
    std::vector<ReleasedPacket> out;
    for (auto& [key, R] : res) {
        const double C = R.count.approx();
        if (!(C > 0.0)) {
            R.step_entries.clear();
            continue;
        }
        const double x = C / w0;
        require(x < 0x1p52, "release count too large");
        const double fl = std::floor(x);
        ChemKey ck{seed, step, (pfa_u64(key.species) << 32) | key.cell, PFC_PURPOSE_RELEASE_COUNT};
        pfa_u64 n = pfa_u64(fl) + pfa_u64(word_bernoulli(ck, 0, x - fl));
        // cap at the exact stock: n*w0 <= C
        while (n > 0) {
            ExactSum need;
            need.add_product(double(n), w0);
            ExactSum rem = R.count;
            rem.sub(need);
            if (!rem.negative()) break;
            --n;
        }
        const double m = b.view.species[key.species].mass_eV;
        double total_w = 0.0;
        for (auto& e : R.step_entries) total_w += e.weight;
        for (pfa_u64 k = 0; k < n; ++k) {
            V3 u;
            double xp[3] = {R.last_x[0], R.last_x[1], R.last_x[2]};
            if (!R.step_entries.empty()) {
                ChemKey pk{seed, step, (pfa_u64(key.species) << 32) | key.cell, PFC_PURPOSE_RELEASE_PICK};
                const double y = chem_u01(pk, pfc_u32(k), 0) * total_w;
                double cum = 0.0;
                std::size_t j = R.step_entries.size() - 1;
                for (std::size_t q = 0; q < R.step_entries.size(); ++q) {
                    cum += R.step_entries[q].weight;
                    if (y < cum) {
                        j = q;
                        break;
                    }
                }
                u = v3(R.step_entries[j].u[0], R.step_entries[j].u[1], R.step_entries[j].u[2]);
                for (int a = 0; a < 3; ++a) xp[a] = R.step_entries[j].x[a];
            } else {
                const double c = PFC_C_M_S;
                u = v3(R.momentum[0].approx() / C / m * c, R.momentum[1].approx() / C / m * c, R.momentum[2].approx() / C / m * c);
            }
            const V3 p = p_from_u(m, u);
            R.count.add(-w0);
            ExactSum t;
            for (int a = 0; a < 3; ++a) {
                const double comp = a == 0 ? p.x : (a == 1 ? p.y : p.z);
                ExactSum s;
                s.add_product(w0, comp);
                R.momentum[a].sub(s);
            }
            ExactSum ke;
            ke.add_product(w0, ke_from_u(m, u));
            R.kinetic.sub(ke);
            out.push_back(ReleasedPacket{key.cell, key.species, pfc_u32(k), w0, {u.x, u.y, u.z}, {xp[0], xp[1], xp[2]}});
        }
        if (!R.step_entries.empty())
            for (int a = 0; a < 3; ++a) R.last_x[a] = R.step_entries.back().x[a];
        R.step_entries.clear();
    }
    return out;
}
// pos: the projectile carrier's position at the event (m); nullptr records the origin.
inline void deposit_reservoir(const Blob& b, std::map<InventoryKey, Reservoir>& res, const EventRec& e, const double* pos = nullptr) {
    for (pfc_u32 i = 0; i < e.n_out; ++i) {
        if (e.out_dest[i] != PFC_DEST_RESERVOIR) continue;
        const pfc_u32 s = e.out_species[i];
        Reservoir& R = res[InventoryKey{e.cell, s}];
        const double m = b.view.species[s].mass_eV;
        const V3 u = v3(e.out_u[i][0], e.out_u[i][1], e.out_u[i][2]);
        const V3 p = p_from_u(m, u);
        R.count.add(e.weight);
        R.momentum[0].add_product(e.weight, p.x);
        R.momentum[1].add_product(e.weight, p.y);
        R.momentum[2].add_product(e.weight, p.z);
        R.kinetic.add_product(e.weight, e.ke_out[i]);
        R.step_entries.push_back(ReservoirEntry{e.weight, {u.x, u.y, u.z}, {pos ? pos[0] : 0.0, pos ? pos[1] : 0.0, pos ? pos[2] : 0.0}});
    }
}

}  // namespace pfc
