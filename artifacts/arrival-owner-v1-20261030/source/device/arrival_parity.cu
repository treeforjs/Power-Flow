// arrival_parity.cu -- NVRTC program: the arrival owner v1 RNG layer on the device.
// Compiled by the shared GPU layer (hall-gpu-resident-v1 common/) with
// --fmad=false --prec-div=true --prec-sqrt=true --ftz=false, from the same
// philox_core.h bytes the host compiles. One thread per key record.
#include "philox_core.h"

struct PfaKeyRecord {
    pfa_u64 seed, step;
    pfa_u32 stream, candidate, purpose, block;
};

extern "C" __global__ void pfa_parity(const PfaKeyRecord* keys, const double* x, const double* p, unsigned n, pfa_u32* raw, double* uniforms,
                                      double* electron, double* neutral, pfa_u64* rounded, unsigned* bernoulli) {
    const unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    const PfaKeyRecord k = keys[i];
    pfa::StreamKey s;
    s.seed = k.seed;
    s.step = k.step;
    s.stream = k.stream;
    s.candidate = k.candidate;
    s.purpose = k.purpose;
    const pfa::Block b = pfa::stream_block(s, k.block);
    for (int j = 0; j < 4; ++j) raw[4ull * i + j] = b.v[j];
    uniforms[3ull * i] = pfa::u01_open(pfa::low_word(b));
    uniforms[3ull * i + 1] = pfa::u01_open(pfa::high_word(b));
    uniforms[3ull * i + 2] = pfa::u01_half_open(pfa::low_word(b));
    pfa::electron_marks(s, electron + 4ull * i);
    pfa::neutral_marks(s, neutral + 6ull * i);
    rounded[i] = pfa::stochastic_round(s, x[i]);
    bernoulli[i] = (unsigned)pfa::bernoulli_exact(s, p[i], k.block);
}
