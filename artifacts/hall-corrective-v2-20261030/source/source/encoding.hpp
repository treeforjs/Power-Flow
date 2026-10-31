#pragma once
// Canonical binary encoding and streaming SHA-256 for Hall corrective v2.
//
// Encoding rules (identical to the native v5 convention, include/raikou/
// mfem_pic_canonical.h): u64 = 8 bytes little-endian; f64 = the IEEE-754
// binary64 bit pattern as u64 (exact; -0 and +0 differ); bool = one byte 0/1;
// text = u64 byte length + bytes; a sequence = u64 element count + elements in
// the stated canonical order. Every record starts with its version tag (text).
// The digest is SHA-256 over these bytes, written as lower-case hex.
//
// The encoder streams into the hash through a fixed staging buffer: no JSON
// value, no growing byte string and no heap allocation per field. The only
// allocations are the OpenSSL digest context (released in the destructor) and
// the returned hex string.
#include <openssl/evp.h>
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstring>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace hall_corrective_v2 {
static_assert(std::endian::native == std::endian::little, "canonical encoding assumes a little-endian host");
static_assert(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559, "binary64 required");

class Digest {
    EVP_MD_CTX* ctx_ = nullptr;
    std::array<unsigned char, 16384> buffer_{};
    std::size_t used_ = 0;
    bool finished_ = false;
    void flush() {
        if (used_ && EVP_DigestUpdate(ctx_, buffer_.data(), used_) != 1) throw std::runtime_error("SHA-256 update failed");
        used_ = 0;
    }
    void raw(const void* data, std::size_t bytes) {
        if (finished_) throw std::logic_error("digest already finished");
        const auto* p = static_cast<const unsigned char*>(data);
        while (bytes) {
            if (used_ == buffer_.size()) flush();
            const std::size_t n = std::min(bytes, buffer_.size() - used_);
            std::memcpy(buffer_.data() + used_, p, n);
            used_ += n; p += n; bytes -= n;
        }
    }
public:
    explicit Digest(std::string_view tag) {
        ctx_ = EVP_MD_CTX_new();
        if (!ctx_) throw std::runtime_error("SHA-256 context allocation failed");
        if (EVP_DigestInit_ex(ctx_, EVP_sha256(), nullptr) != 1) { EVP_MD_CTX_free(ctx_); ctx_ = nullptr; throw std::runtime_error("SHA-256 init failed"); }
        text(tag);
    }
    ~Digest() { EVP_MD_CTX_free(ctx_); }
    Digest(const Digest&) = delete;
    Digest& operator=(const Digest&) = delete;
    void u64(std::uint64_t v) { raw(&v, 8); }
    void i64(std::int64_t v) { u64(static_cast<std::uint64_t>(v)); }
    void f64(double v) { u64(std::bit_cast<std::uint64_t>(v)); }
    void boolean(bool v) { const unsigned char b = v ? 1 : 0; raw(&b, 1); }
    void text(std::string_view s) { u64(s.size()); raw(s.data(), s.size()); }
    void f64s(std::span<const double> v) { u64(v.size()); raw(v.data(), v.size() * 8); }
    void i32s_as_i64(std::span<const int> v) { u64(v.size()); for (int x : v) i64(x); }
    void bools(const std::vector<bool>& v) { u64(v.size()); for (bool x : v) boolean(x); }
    // Lower-case hex SHA-256 of everything written so far. The digest is closed.
    std::string hex() {
        flush();
        unsigned char out[EVP_MAX_MD_SIZE];
        unsigned int n = 0;
        if (EVP_DigestFinal_ex(ctx_, out, &n) != 1 || n != 32) throw std::runtime_error("SHA-256 final failed");
        finished_ = true;
        static constexpr char hexdigits[] = "0123456789abcdef";
        std::string s(64, '0');
        for (unsigned i = 0; i < 32; ++i) { s[2 * i] = hexdigits[out[i] >> 4]; s[2 * i + 1] = hexdigits[out[i] & 15]; }
        return s;
    }
};

// SHA-256 of a byte string, as lower-case hex (no tag; used for file bytes
// and for digests that must equal the sealed text digests).
inline std::string bytes_sha256(std::string_view bytes) {
    EVP_MD_CTX* ctx = EVP_MD_CTX_new();
    if (!ctx) throw std::runtime_error("SHA-256 context allocation failed");
    unsigned char out[EVP_MAX_MD_SIZE];
    unsigned int n = 0;
    const bool ok = EVP_DigestInit_ex(ctx, EVP_sha256(), nullptr) == 1 && EVP_DigestUpdate(ctx, bytes.data(), bytes.size()) == 1 &&
                    EVP_DigestFinal_ex(ctx, out, &n) == 1 && n == 32;
    EVP_MD_CTX_free(ctx);
    if (!ok) throw std::runtime_error("SHA-256 failed");
    static constexpr char hexdigits[] = "0123456789abcdef";
    std::string s(64, '0');
    for (unsigned i = 0; i < 32; ++i) { s[2 * i] = hexdigits[out[i] >> 4]; s[2 * i + 1] = hexdigits[out[i] & 15]; }
    return s;
}
}  // namespace hall_corrective_v2
