// canonical.hpp -- canonical binary encoding and SHA-256 for the arrival owner v1.
//
// Primitive rules (the same rules as native v5 include/raikou/mfem_pic_canonical.h):
//  * u64: 8 bytes little-endian; i64: two's complement as u64.
//  * f64: the IEEE-754 binary64 bit pattern as u64 (exact; -0.0 != +0.0).
//  * bool: one byte 0x00 or 0x01 (any other byte refuses on decode).
//  * u8: one byte.
//  * text: u64 byte length, then the bytes (no terminator).
//  * every record starts with its version tag encoded as text.
//  * sequences: u64 element count, then the elements in a stated canonical order.
//  * an exact count: its canonical positive-carrier decomposition (u64 n + n f64),
//    as native v5 encodes ExactParticleCount.
// Digests are SHA-256 over these bytes, lower-case hex. No JSON is hashed.
#pragma once
#include <openssl/evp.h>
#include <bit>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace pf_arrival {

struct EncodingError : std::invalid_argument {
    using std::invalid_argument::invalid_argument;
};

inline std::string hex_bytes(const unsigned char* p, std::size_t n) {
    static const char* d = "0123456789abcdef";
    std::string s(2 * n, '0');
    for (std::size_t i = 0; i < n; ++i) {
        s[2 * i] = d[p[i] >> 4];
        s[2 * i + 1] = d[p[i] & 15];
    }
    return s;
}

inline std::string sha256_hex(std::string_view bytes) {
    unsigned char out[EVP_MAX_MD_SIZE];
    unsigned int n = 0;
    if (EVP_Digest(bytes.data(), bytes.size(), out, &n, EVP_sha256(), nullptr) != 1 || n != 32)
        throw std::runtime_error("SHA-256 evaluation failed");
    return hex_bytes(out, n);
}

class Writer {
public:
    std::string bytes;
    void u8(std::uint8_t v) { bytes.push_back(char(v)); }
    void u64(std::uint64_t v) {
        char b[8];
        for (int i = 0; i < 8; ++i) b[i] = char((v >> (8 * i)) & 255u);
        bytes.append(b, 8);
    }
    void i64(std::int64_t v) { u64(static_cast<std::uint64_t>(v)); }
    void f64(double x) { u64(std::bit_cast<std::uint64_t>(x)); }
    void boolean(bool b) { u8(b ? 1 : 0); }
    void text(std::string_view s) {
        u64(s.size());
        bytes.append(s.data(), s.size());
    }
    void tag(std::string_view s) { text(s); }
    void f64s(const std::vector<double>& v) {
        u64(v.size());
        for (double x : v) f64(x);
    }
    std::string sha256() const { return sha256_hex(bytes); }
};

// Strict reader over caller-owned bytes; every malformed input refuses.
class Reader {
public:
    explicit Reader(std::string_view b) : b_(b) {}
    std::uint8_t u8() { return static_cast<std::uint8_t>(take(1)[0]); }
    std::uint64_t u64() {
        const auto s = take(8);
        std::uint64_t v = 0;
        for (int i = 0; i < 8; ++i) v |= std::uint64_t(static_cast<unsigned char>(s[i])) << (8 * i);
        return v;
    }
    std::int64_t i64() { return static_cast<std::int64_t>(u64()); }
    double f64() { return std::bit_cast<double>(u64()); }
    bool boolean() {
        const auto v = u8();
        if (v > 1) throw EncodingError("canonical bool byte other than 0/1");
        return v == 1;
    }
    std::string text(std::size_t cap = 1u << 20) {
        const auto n = u64();
        if (n > cap) throw EncodingError("canonical text above its cap");
        return std::string(take(n));
    }
    void expect_tag(std::string_view t) {
        if (text(256) != t) throw EncodingError("canonical record tag mismatch");
    }
    std::vector<double> f64s(std::size_t cap) {
        const auto n = u64();
        if (n > cap) throw EncodingError("canonical sequence above its cap");
        std::vector<double> v;
        v.reserve(n);
        for (std::uint64_t i = 0; i < n; ++i) v.push_back(f64());
        return v;
    }
    std::size_t position() const noexcept { return pos_; }
    bool done() const noexcept { return pos_ == b_.size(); }
    void end() const {
        if (!done()) throw EncodingError("trailing bytes after canonical record");
    }

private:
    std::string_view take(std::size_t n) {
        if (n > b_.size() - pos_) throw EncodingError("canonical record truncated");
        auto s = b_.substr(pos_, n);
        pos_ += n;
        return s;
    }
    std::string_view b_;
    std::size_t pos_ = 0;
};

}  // namespace pf_arrival
