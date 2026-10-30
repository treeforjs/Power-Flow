// encoding.hpp -- canonical binary encoding of named arrays ("HGPENC1").
// Layout: 8-byte magic "HGPENC1\0", u64 entry count, then per entry:
//   u64 name length, name bytes, u8 type ('d' binary64, 'q' u64, 'i' i32, 's' UTF-8 bytes),
//   u64 element count, elements little-endian.
// Used for fixture inputs and for every output compared by tests/compare.py.
#pragma once
#include <bit>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <map>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace hgp_enc {

struct Entry {
    char type = 'd';
    std::vector<unsigned char> bytes;   // raw little-endian payload
    std::uint64_t count = 0;
};

class Writer {
public:
    void f64(const std::string& name, std::span<const double> v) { put(name, 'd', v.data(), v.size(), 8); }
    void u64(const std::string& name, std::span<const std::uint64_t> v) { put(name, 'q', v.data(), v.size(), 8); }
    void i32(const std::string& name, std::span<const std::int32_t> v) { put(name, 'i', v.data(), v.size(), 4); }
    void text(const std::string& name, const std::string& s) { put(name, 's', s.data(), s.size(), 1); }
    std::string bytes() const {
        std::string out("HGPENC1\0", 8);
        word(out, order_.size());
        for (const auto& name : order_) {
            const auto& e = entries_.at(name);
            word(out, name.size());
            out.append(name);
            out.push_back(e.type);
            word(out, e.count);
            out.append(reinterpret_cast<const char*>(e.bytes.data()), e.bytes.size());
        }
        return out;
    }
    void save(const std::string& path) const {
        const auto b = bytes();
        std::ofstream out(path, std::ios::binary | std::ios::trunc);
        out.write(b.data(), std::streamsize(b.size()));
        if (!out) throw std::runtime_error("write failed: " + path);
    }

private:
    static_assert(std::endian::native == std::endian::little, "little-endian host required");
    void put(const std::string& name, char type, const void* data, std::size_t count, std::size_t size) {
        if (entries_.count(name)) throw std::runtime_error("duplicate encoding entry " + name);
        Entry e;
        e.type = type;
        e.count = count;
        e.bytes.resize(count * size);
        if (count) std::memcpy(e.bytes.data(), data, count * size);
        entries_[name] = std::move(e);
        order_.push_back(name);
    }
    static void word(std::string& out, std::uint64_t v) { out.append(reinterpret_cast<const char*>(&v), 8); }
    std::map<std::string, Entry> entries_;
    std::vector<std::string> order_;
};

class Reader {
public:
    explicit Reader(const std::string& path) {
        std::ifstream in(path, std::ios::binary);
        if (!in) throw std::runtime_error("cannot read " + path);
        char magic[8];
        in.read(magic, 8);
        if (!in || std::memcmp(magic, "HGPENC1", 8) != 0) throw std::runtime_error("bad encoding magic in " + path);
        const auto entries = word(in);
        for (std::uint64_t k = 0; k < entries; ++k) {
            std::string name(word(in), '\0');
            in.read(name.data(), std::streamsize(name.size()));
            Entry e;
            e.type = char(in.get());
            e.count = word(in);
            const std::size_t size = e.type == 'i' ? 4 : e.type == 's' ? 1 : 8;
            e.bytes.resize(e.count * size);
            in.read(reinterpret_cast<char*>(e.bytes.data()), std::streamsize(e.bytes.size()));
            if (!in) throw std::runtime_error("truncated encoding " + path);
            map_[name] = std::move(e);
        }
    }
    bool has(const std::string& n) const { return map_.count(n) != 0; }
    std::vector<double> f64(const std::string& n) const { return take<double>(n, 'd'); }
    std::vector<std::uint64_t> u64(const std::string& n) const { return take<std::uint64_t>(n, 'q'); }
    std::vector<std::int32_t> i32(const std::string& n) const { return take<std::int32_t>(n, 'i'); }
    std::string text(const std::string& n) const {
        const auto& e = at(n, 's');
        return std::string(e.bytes.begin(), e.bytes.end());
    }

private:
    const Entry& at(const std::string& n, char type) const {
        auto it = map_.find(n);
        if (it == map_.end()) throw std::runtime_error("missing encoding entry " + n);
        if (it->second.type != type) throw std::runtime_error("encoding entry type mismatch " + n);
        return it->second;
    }
    template<class T> std::vector<T> take(const std::string& n, char type) const {
        const auto& e = at(n, type);
        std::vector<T> v(e.count);
        if (e.count) std::memcpy(v.data(), e.bytes.data(), e.bytes.size());
        return v;
    }
    static std::uint64_t word(std::ifstream& in) { std::uint64_t v = 0; in.read(reinterpret_cast<char*>(&v), 8); return v; }
    std::map<std::string, Entry> map_;
};

// Field names of the per-carrier diagnostics (same order as hgp::DiagReal / hgp::DiagInt).
inline const char* const diag_real_names[] = {"moved_dx", "moved_dy", "moved_dz", "moved_ux", "moved_uy", "moved_uz",
                                              "hit_fraction", "flight_s", "gyro", "displacement_cells", "kinetic_fraction", "kinetic_plane"};
inline const char* const diag_int_names[] = {"moved_valid", "hit_surface", "kinetic_valid", "kinetic_sink", "kinetic_axis"};

// Common record-level writer: SoA endpoint/path arrays are written interleaved per
// carrier (x, y, z), exactly as the reference GapParticle / CarrierPath vectors.
struct RecordArrays {
    std::vector<std::uint64_t> id;
    std::vector<double> end_pos, end_u, charge, mass, weight;   // interleaved 3N for pos/u
    std::vector<double> start, disp, flight;                    // interleaved 3N for start/disp
    std::vector<double> current, initial_charge, endpoint_charge;
    std::vector<double> budget;                                 // 7 values
    std::vector<double> net_hi, net_lo, gross;                  // free charge of the endpoints
    std::string refusal, charge_refusal;
    std::int32_t refusal_rank = -1, charge_refusal_rank = -1;
    void write(Writer& w, bool include_rank) const {
        w.text("refusal.motion", refusal);
        w.text("refusal.charge", charge_refusal);
        if (include_rank) {
            const std::int32_t r[2] = {refusal_rank, charge_refusal_rank};
            w.i32("refusal.rank", r);
        }
        if (!refusal.empty()) return;
        w.u64("endpoint.native_id", id);
        w.f64("endpoint.position_m", end_pos);
        w.f64("endpoint.proper_velocity_m_s", end_u);
        w.f64("endpoint.charge_C", charge);
        w.f64("endpoint.mass_kg", mass);
        w.f64("endpoint.weight", weight);
        w.u64("path.native_id", id);
        w.f64("path.field_start_m", start);
        w.f64("path.field_displacement_m", disp);
        w.f64("path.flight_s", flight);
        w.f64("current_A_m", current);
        w.f64("initial_charge_C", initial_charge);
        w.f64("endpoint_charge_C", endpoint_charge);
        w.f64("budget", budget);
        if (!charge_refusal.empty()) return;
        w.f64("free_charge.net_hi", net_hi);
        w.f64("free_charge.net_lo", net_lo);
        w.f64("free_charge.gross_C", gross);
    }
};

}  // namespace hgp_enc
