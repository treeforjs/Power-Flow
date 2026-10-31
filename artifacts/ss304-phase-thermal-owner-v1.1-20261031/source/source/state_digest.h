// state_digest.h -- block-Merkle digest of the owner state (host and NVRTC), in the layout of the native v7 forms
// (API-007 section 1.2): records in index order; chunk j = 32 records, c_j = SHA256(text "<tag>/chunk", u64 j,
// u64 m_j, record bytes); leaf k = 1024 records, l_k = SHA256(text "<tag>/block", u64 k, u64 n_k, u64 q_k, chunk
// digests); tree: node(a, b) = SHA256(text "native-v5/merkle-node", a, b), pairs left to right, an odd last node
// carried up, no leaves: SHA256(text "native-v5/merkle-empty"); root = SHA256(text "<tag>/root", u64 1024, u64 n,
// u64 L, tree). Encoding: u64 little-endian 8 bytes; text = u64 byte length then the bytes; SHA-256 FIPS 180-4.
// Sections: "ss304-phase-v1.1/cells" (record: E lo, E hi, then N lo/hi for the six elements, then flags as u64:
// 15 words) and "ss304-phase-v1.1/face-columns" (record: face mass words lo/hi for the six elements: 12 words).
#pragma once
#include "det_math.h"

#if !PQ_DEVICE
#include <array>
#include <vector>
#endif

struct pq_sha { pq_u32 h[8]; unsigned char buf[64]; pq_u64 len; int nb; };

PQ_HD pq_u32 pq_rotr(pq_u32 x, int n) { return (x >> n) | (x << (32 - n)); }
PQ_HD void pq_sha_compress(pq_sha& s) {
    const pq_u32 K[64] = {
        0x428a2f98u,0x71374491u,0xb5c0fbcfu,0xe9b5dba5u,0x3956c25bu,0x59f111f1u,0x923f82a4u,0xab1c5ed5u,0xd807aa98u,0x12835b01u,0x243185beu,0x550c7dc3u,
        0x72be5d74u,0x80deb1feu,0x9bdc06a7u,0xc19bf174u,0xe49b69c1u,0xefbe4786u,0x0fc19dc6u,0x240ca1ccu,0x2de92c6fu,0x4a7484aau,0x5cb0a9dcu,0x76f988dau,
        0x983e5152u,0xa831c66du,0xb00327c8u,0xbf597fc7u,0xc6e00bf3u,0xd5a79147u,0x06ca6351u,0x14292967u,0x27b70a85u,0x2e1b2138u,0x4d2c6dfcu,0x53380d13u,
        0x650a7354u,0x766a0abbu,0x81c2c92eu,0x92722c85u,0xa2bfe8a1u,0xa81a664bu,0xc24b8b70u,0xc76c51a3u,0xd192e819u,0xd6990624u,0xf40e3585u,0x106aa070u,
        0x19a4c116u,0x1e376c08u,0x2748774cu,0x34b0bcb5u,0x391c0cb3u,0x4ed8aa4au,0x5b9cca4fu,0x682e6ff3u,0x748f82eeu,0x78a5636fu,0x84c87814u,0x8cc70208u,
        0x90befffau,0xa4506cebu,0xbef9a3f7u,0xc67178f2u};
    pq_u32 w[64];
    for (int i = 0; i < 16; ++i) w[i] = ((pq_u32)s.buf[4 * i] << 24) | ((pq_u32)s.buf[4 * i + 1] << 16) | ((pq_u32)s.buf[4 * i + 2] << 8) | (pq_u32)s.buf[4 * i + 3];
    for (int i = 16; i < 64; ++i) {
        const pq_u32 s0 = pq_rotr(w[i - 15], 7) ^ pq_rotr(w[i - 15], 18) ^ (w[i - 15] >> 3);
        const pq_u32 s1 = pq_rotr(w[i - 2], 17) ^ pq_rotr(w[i - 2], 19) ^ (w[i - 2] >> 10);
        w[i] = w[i - 16] + s0 + w[i - 7] + s1;
    }
    pq_u32 a = s.h[0], b = s.h[1], c = s.h[2], d = s.h[3], e = s.h[4], f = s.h[5], g = s.h[6], h = s.h[7];
    for (int i = 0; i < 64; ++i) {
        const pq_u32 S1 = pq_rotr(e, 6) ^ pq_rotr(e, 11) ^ pq_rotr(e, 25);
        const pq_u32 ch = (e & f) ^ (~e & g);
        const pq_u32 t1 = h + S1 + ch + K[i] + w[i];
        const pq_u32 S0 = pq_rotr(a, 2) ^ pq_rotr(a, 13) ^ pq_rotr(a, 22);
        const pq_u32 mj = (a & b) ^ (a & c) ^ (b & c);
        const pq_u32 t2 = S0 + mj;
        h = g; g = f; f = e; e = d + t1; d = c; c = b; b = a; a = t1 + t2;
    }
    s.h[0] += a; s.h[1] += b; s.h[2] += c; s.h[3] += d; s.h[4] += e; s.h[5] += f; s.h[6] += g; s.h[7] += h;
}
PQ_HD void pq_sha_init(pq_sha& s) {
    s.h[0] = 0x6a09e667u; s.h[1] = 0xbb67ae85u; s.h[2] = 0x3c6ef372u; s.h[3] = 0xa54ff53au;
    s.h[4] = 0x510e527fu; s.h[5] = 0x9b05688cu; s.h[6] = 0x1f83d9abu; s.h[7] = 0x5be0cd19u;
    s.len = 0; s.nb = 0;
}
PQ_HD void pq_sha_byte(pq_sha& s, pq_u32 v) {
    s.buf[s.nb++] = (unsigned char)(v & 0xffu); s.len += 1;
    if (s.nb == 64) { pq_sha_compress(s); s.nb = 0; }
}
PQ_HD void pq_sha_u64(pq_sha& s, pq_u64 v) { for (int i = 0; i < 8; ++i) pq_sha_byte(s, (pq_u32)(v >> (8 * i))); }
PQ_HD void pq_sha_text(pq_sha& s, const char* t) {
    pq_u64 n = 0; while (t[n]) ++n;
    pq_sha_u64(s, n);
    for (pq_u64 i = 0; i < n; ++i) pq_sha_byte(s, (pq_u32)(unsigned char)t[i]);
}
PQ_HD void pq_sha_digest(pq_sha& s, const pq_u32* d) { for (int i = 0; i < 8; ++i) for (int k = 3; k >= 0; --k) pq_sha_byte(s, d[i] >> (8 * k)); }
PQ_HD void pq_sha_finish(pq_sha& s, pq_u32* out) {
    const pq_u64 bits = s.len * 8;
    pq_sha_byte(s, 0x80u);
    while (s.nb != 56) pq_sha_byte(s, 0u);
    for (int k = 7; k >= 0; --k) pq_sha_byte(s, (pq_u32)(bits >> (8 * k)));
    for (int i = 0; i < 8; ++i) out[i] = s.h[i];
}

#define PQ_SEC_CELLS 0
#define PQ_SEC_FC 1
PQ_HD const char* pq_tag(int section, int level) {
    // level 0 chunk, 1 block, 2 root
    if (section == PQ_SEC_CELLS) return level == 0 ? "ss304-phase-v1.1/cells/chunk" : level == 1 ? "ss304-phase-v1.1/cells/block" : "ss304-phase-v1.1/cells/root";
    return level == 0 ? "ss304-phase-v1.1/face-columns/chunk" : level == 1 ? "ss304-phase-v1.1/face-columns/block" : "ss304-phase-v1.1/face-columns/root";
}
PQ_HD void pq_record_words(int section, const pq_i128* E, const pq_i128* N, const int* flags, const pq_i128* fcm, pq_u64 r, pq_u64* w) {
    if (section == PQ_SEC_CELLS) {
        w[0] = E[r].lo; w[1] = (pq_u64)E[r].hi;
        for (int el = 0; el < 6; ++el) { w[2 + 2 * el] = N[6 * r + el].lo; w[3 + 2 * el] = (pq_u64)N[6 * r + el].hi; }
        w[14] = (pq_u64)(pq_i64)flags[r];
    } else {
        for (int el = 0; el < 6; ++el) { w[2 * el] = fcm[6 * r + el].lo; w[2 * el + 1] = (pq_u64)fcm[6 * r + el].hi; }
    }
}
PQ_HD int pq_record_count(int section) { return section == PQ_SEC_CELLS ? 15 : 12; }
// chunk j of n records
PQ_HD void pq_chunk_digest(int section, pq_u64 j, pq_u64 n, const pq_i128* E, const pq_i128* N, const int* flags, const pq_i128* fcm, pq_u32* out) {
    pq_sha s; pq_sha_init(s);
    const pq_u64 r0 = 32 * j, r1 = r0 + 32 < n ? r0 + 32 : n;
    pq_sha_text(s, pq_tag(section, 0)); pq_sha_u64(s, j); pq_sha_u64(s, r1 - r0);
    pq_u64 w[15];
    for (pq_u64 r = r0; r < r1; ++r) { pq_record_words(section, E, N, flags, fcm, r, w); for (int q = 0; q < pq_record_count(section); ++q) pq_sha_u64(s, w[q]); }
    pq_sha_finish(s, out);
}
// leaf k from its chunk digests (chunks [32k, 32k + q_k))
PQ_HD void pq_leaf_digest(int section, pq_u64 k, pq_u64 n, const pq_u32* chunks, pq_u32* out) {
    const pq_u64 r0 = 1024 * k, r1 = r0 + 1024 < n ? r0 + 1024 : n, nk = r1 - r0, qk = (nk + 31) / 32;
    pq_sha s; pq_sha_init(s);
    pq_sha_text(s, pq_tag(section, 1)); pq_sha_u64(s, k); pq_sha_u64(s, nk); pq_sha_u64(s, qk);
    for (pq_u64 q = 0; q < qk; ++q) pq_sha_digest(s, chunks + 8 * (32 * k + q));
    pq_sha_finish(s, out);
}
PQ_HD void pq_merkle_node(const pq_u32* a, const pq_u32* b, pq_u32* out) {
    pq_sha s; pq_sha_init(s); pq_sha_text(s, "native-v5/merkle-node"); pq_sha_digest(s, a); pq_sha_digest(s, b); pq_sha_finish(s, out);
}
PQ_HD void pq_merkle_empty(pq_u32* out) { pq_sha s; pq_sha_init(s); pq_sha_text(s, "native-v5/merkle-empty"); pq_sha_finish(s, out); }

#if !PQ_DEVICE
struct StateRoots { std::array<pq_u32, 8> cells{}, fcm{}; std::size_t cell_leaves = 0, fc_leaves = 0; };
// Root of a section from its leaf digests (host; identical for CPU and GPU engines).
inline std::array<pq_u32, 8> pq_section_root(int section, pq_u64 n, std::vector<std::array<pq_u32, 8>> level) {
    const pq_u64 L = level.size();
    std::array<pq_u32, 8> tree{};
    if (level.empty()) pq_merkle_empty(tree.data());
    else {
        while (level.size() > 1) {
            std::vector<std::array<pq_u32, 8>> up;
            for (std::size_t i = 0; i + 1 < level.size(); i += 2) { std::array<pq_u32, 8> o; pq_merkle_node(level[i].data(), level[i + 1].data(), o.data()); up.push_back(o); }
            if (level.size() % 2) up.push_back(level.back());
            level.swap(up);
        }
        tree = level[0];
    }
    pq_sha s; pq_sha_init(s);
    pq_sha_text(s, pq_tag(section, 2)); pq_sha_u64(s, 1024); pq_sha_u64(s, n); pq_sha_u64(s, L); pq_sha_digest(s, tree.data());
    std::array<pq_u32, 8> root; pq_sha_finish(s, root.data());
    return root;
}
inline std::vector<std::array<pq_u32, 8>> pq_host_leaves(int section, pq_u64 n, const pq_i128* E, const pq_i128* N, const int* flags, const pq_i128* fcm) {
    const pq_u64 nc = (n + 31) / 32, nl = (n + 1023) / 1024;
    std::vector<pq_u32> chunks(8 * nc);
    for (pq_u64 j = 0; j < nc; ++j) pq_chunk_digest(section, j, n, E, N, flags, fcm, chunks.data() + 8 * j);
    std::vector<std::array<pq_u32, 8>> leaves(nl);
    for (pq_u64 k = 0; k < nl; ++k) pq_leaf_digest(section, k, n, chunks.data(), leaves[k].data());
    return leaves;
}
inline StateRoots host_state_roots(const pq_i128* E, const pq_i128* N, const int* flags, int n_cells, const pq_i128* fcm, int n_fc,
                                   std::vector<std::array<pq_u32, 8>>* cell_leaves) {
    StateRoots r;
    auto lc = pq_host_leaves(PQ_SEC_CELLS, pq_u64(n_cells), E, N, flags, fcm);
    auto lf = pq_host_leaves(PQ_SEC_FC, pq_u64(n_fc), E, N, flags, fcm);
    r.cell_leaves = lc.size(); r.fc_leaves = lf.size();
    if (cell_leaves) *cell_leaves = lc;
    r.cells = pq_section_root(PQ_SEC_CELLS, pq_u64(n_cells), lc);
    r.fcm = pq_section_root(PQ_SEC_FC, pq_u64(n_fc), lf);
    return r;
}
#endif
