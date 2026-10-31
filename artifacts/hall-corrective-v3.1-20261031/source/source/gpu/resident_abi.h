// resident_abi.h -- layout shared by the host engine (resident_field.cpp) and the
// NVRTC device code (resident_kernels.cu) of the Hall v3.1 device-resident field step.
// Plain data only: device addresses as 64-bit integers, sizes, scalars.
#ifndef HALL_V31_RESIDENT_ABI_H
#define HALL_V31_RESIDENT_ABI_H

typedef unsigned long long rf_u64;
typedef unsigned int rf_u32;

#define RF_TPB 256
#define RF_RBLK 256          /* blocks of the first reduction stage (fixed tree) */
#define RF_ROWS 32           /* reduction rows per launch */
#define RF_BLOCK_BYTES 2048  /* Merkle leaf payload size (bytes) */
#define RF_MAX_ARRAYS 16
/* edge terms */
#define ET_FDE 0
#define ET_PW 1
#define ET_JOULE 2
#define ET_CIRCUIT 3
#define ET_OLOSS 4
#define ET_OSUPPLY 5
#define ET_MDE2 6
#define ET_MTHIRD2 7
#define ET_ME2 8
#define ET_AMPRES 9
#define ET_TERMINAL 10
#define ET_REG_SOLID 11
#define ET_REG_PLASMA 12
#define ET_REG_DISP 13
#define ET_REG_EXT 14
#define ET_COUNT 15
/* face terms */
#define FT_FDB 0
#define FT_COUNT 1
/* node values */
#define NV_CONT 0
#define NV_KCL 1
#define NV_GINC 2
#define NV_GACT 3
#define NV_OMEGAP 4
#define NV_DEF0 5
#define NV_BND0 6
#define NV_DEF1 7
#define NV_BND1 8
#define NV_ACT0 9
#define NV_ACT1 10
#define NV_COUNT 11
/* cell values */
#define CV_BMAG 0
#define CV_DIV0 1
#define CV_DBND0 2
#define CV_DIV1 3
#define CV_DBND1 4
#define CV_COUNT 5

/* One field-state slot (accepted or candidate): words of E and B and the charge
   ledger as structure-of-arrays (hi and lo words in separate arrays). */
struct RfSlot {
    rf_u64 E, El, B, Bl;
    rf_u64 fnet_hi, fnet_lo, fgross, cond_hi, cond_lo, ext_hi, ext_lo, residue;
};

/* Reduction rows: sums (double-word fixed tree), maxima (exact), first index (exact). */
struct RfRows {
    rf_u64 ptr[RF_ROWS];
    int n[RF_ROWS];
    int kind[RF_ROWS];       /* 0 sum (double terms), 1 max (non-negative doubles), 2 min (u32) */
    int count;
    int pad;
};

/* Merkle leaf table: arrays of the canonical record, each cut into RF_BLOCK_BYTES blocks. */
struct RfLeafTable {
    rf_u64 ptr[RF_MAX_ARRAYS];
    rf_u64 bytes[RF_MAX_ARRAYS];
    rf_u64 tag[RF_MAX_ARRAYS];        /* array tag number written into every leaf */
    int first_leaf[RF_MAX_ARRAYS + 1];
    int arrays;
    int pad;
    rf_u64 out;                       /* 32 bytes per leaf */
};

struct RfCtx {
    /* mesh (static) */
    int ne, nb, nn, ncell, nx, ny, nz, kshort, kfeed, pad0;
    double dx, dy, dz, initial_gauss_bound;
    rf_u64 fe, fc, H, ctp, ctf, ctc, pc2, sqrtM, M, Y, fixed, ea, eb, elen;
    rf_u64 ecp, ecc, qv, ntp, nte, nts, cedge, cvol, feed, bg, node_body_or_short, exterior_support, cell_body, region_grad;
    /* step scalars */
    double dt, current, sentinel0;
    int transfer, residue_any, pad1, pad2;
    /* slots */
    struct RfSlot old_, new_;
    /* per-step inputs (uploaded) */
    rf_u64 sigma, plasma, after_hi, after_lo, after_gross, support, bank_cond_hi, bank_cond_lo, bank_ext_hi, bank_ext_lo, bank_gross, bank_residue;
    /* operator and solve */
    rf_u64 S, a, diag, rhs, x;
    /* work */
    rf_u64 HBh, HBl, midE, midB, HmidB, cmid, solid, short_cur, metal, exterior, total, disp;
    rf_u64 div_hi, div_lo, Q0_hi, Q0_lo, Q1_hi, Q1_lo;
    rf_u64 joule_cell;
    /* per-element guard terms/values (strided: term k of edge e at eterm[k*ne + e], etc.) */
    rf_u64 eterm, fterm, nval, cval;
    /* u32 codes */
    rf_u64 code_cont, code_joule, code_scratch;
    rf_u64 metal_cells, joule_metal;   /* int index list of metal cells, gathered Joule power */
    int nmetal, pad4;
    /* previous increments (field-accuracy indicator) */
    rf_u64 inc0, inc1, inc2;
    int history, pad3;
};

#endif
