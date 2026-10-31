// chem_rng.h -- chemistry v3 random-number stream layout on the arrival owner
// v1 Philox4x32-10 core (pic-arrival-owner-v1-20261030/include/pf_arrival/
// philox_core.h, read in place, sha256 pinned by the build).
//
// Chemistry stream layout v3 ("pfc-chem-stream-v3"):
//   key      k0 = run seed bits 0..31, k1 = run seed bits 32..63 (the same run
//            seed as the arrival owner)
//   counter  c0 = (event << 8) | block      event < 2^24, block < 2^8
//            c1 = carrier id bits 0..31
//            c2 = carrier id bits 32..47 (bits 0..15) | purpose (bits 16..23)
//                 | step bits 32..39 (bits 24..31)
//            c3 = step bits 0..31
//   domain   carrier id < 2^48, step < 2^40, event < 2^24, block < 2^8,
//            purpose in the chemistry range 0x30..0x3F (production) or
//            0xF0..0xFF (qualification fixtures).
// The arrival owner layout v1 puts its purpose tags 1..3 in the same bit
// field, so the two layouts never share a counter: the joint map is
// injective. "Carrier id" is the native carrier id of the projectile (purpose
// MCC) or the reservoir key (purpose RELEASE: cell index < 2^32 in bits 0..31,
// species index in bits 32..47). "Event" is the candidate index of the
// projectile within the step (null and real candidates both count), so every
// draw is a function of (seed, step, carrier, candidate, word) only: any
// partition of carriers over ranks, threads or devices draws the same numbers.
// Word j of an event (j < 512) is R(j) of the block sequence, as in v1.
#pragma once
#include "det_math.h"
#include "philox_core.h"

namespace pfc {

#define PFC_PURPOSE_MCC 0x30u            // collision stage 1 of a step (first half-step exposure)
#define PFC_PURPOSE_RELEASE_COUNT 0x31u
#define PFC_PURPOSE_RELEASE_PICK 0x32u
#define PFC_PURPOSE_MCC_SECOND 0x34u     // collision stage 2 of the same step (second half-step exposure)
#define PFC_PURPOSE_FIXTURE_SWARM 0xF1u
#define PFC_PURPOSE_FIXTURE_TEST 0xF2u
#define PFC_CARRIER_LIMIT 0x1000000000000ull   // 2^48
#define PFC_EVENT_LIMIT 0x1000000u              // 2^24
#define PFC_EVENT_WORDS 512u

// Word slots of one MCC candidate event (fixed; unused slots are skipped).
#define PFC_W_TIME 0u          // waiting time to this candidate
#define PFC_W_TARGET 1u        // target species selection
#define PFC_W_CHANNEL 2u       // real/null decision and channel selection
#define PFC_W_TVEL 4u          // 4 words: target velocity (Box-Muller pairs)
#define PFC_W_KIN 8u           // kinematics words 8..63

struct ChemKey {
    pfa_u64 seed;
    pfa_u64 step;      // accepted step index, < 2^40
    pfa_u64 carrier;   // < 2^48
    pfa_u32 purpose;
};

PFC_FN pfa::Block chem_block(const ChemKey& k, pfa_u32 event, pfa_u32 block) {
    pfa::Key2 key;
    key.v[0] = (pfa_u32)k.seed;
    key.v[1] = (pfa_u32)(k.seed >> 32);
    pfa::Block c;
    c.v[0] = (event << 8) | (block & 0xFFu);
    c.v[1] = (pfa_u32)k.carrier;
    c.v[2] = (pfa_u32)((k.carrier >> 32) & 0xFFFFull) | ((k.purpose & 0xFFu) << 16) | ((pfa_u32)((k.step >> 32) & 0xFFull) << 24);
    c.v[3] = (pfa_u32)k.step;
    return pfa::philox4x32_10(c, key);
}
PFC_FN pfa_u64 chem_word(const ChemKey& k, pfa_u32 event, pfa_u32 j) {
    const pfa::Block b = chem_block(k, event, j >> 1);
    return (j & 1u) ? pfa::high_word(b) : pfa::low_word(b);
}
// Open uniform in (0,1) from word j of event (arrival u01_open conversion).
PFC_FN double chem_u01(const ChemKey& k, pfa_u32 event, pfa_u32 j) {
    return pfa::u01_open(chem_word(k, event, j));
}
// Domain check (host and device): returns 1 iff the key and event are inside
// the injective domain.
PFC_FN int chem_key_valid(const ChemKey& k, pfa_u32 event) {
    if (k.carrier >= PFC_CARRIER_LIMIT) return 0;
    if (k.step >= PFA_STEP_LIMIT) return 0;
    if (event >= PFC_EVENT_LIMIT) return 0;
    const pfa_u32 p = k.purpose;
    if (!((p >= 0x30u && p <= 0x3Fu) || (p >= 0xF0u && p <= 0xFFu))) return 0;
    return 1;
}

}  // namespace pfc
