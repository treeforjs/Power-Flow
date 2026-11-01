// gpu_program.h -- main translation unit of the owner's NVRTC program (shared by gpu_engine.cpp and
// tests/nvrtc_archs.cpp): the phase_kernels.h / state_digest.h kernels plus the Merkle chunk and leaf kernels.
#pragma once
static const char kPxMain[] =
    "#include \"phase_kernels.h\"\n"
    "#include \"state_digest.h\"\n"
    "extern \"C\" __global__ void k_chunks(int section, px_u64 n, px_digest_src g, px_u32* out) {\n"
    "    const px_u64 j = (px_u64)blockIdx.x * blockDim.x + threadIdx.x; if (j < (n + 31) / 32) px_chunk_digest(section, j, n, g, out + 8 * j);\n"
    "}\n"
    "extern \"C\" __global__ void k_leaves(int section, px_u64 n, const px_u32* chunks, px_u32* out) {\n"
    "    const px_u64 k = (px_u64)blockIdx.x * blockDim.x + threadIdx.x; if (k < (n + 1023) / 1024) px_leaf_digest(section, k, n, chunks, out + 8 * k);\n"
    "}\n";
