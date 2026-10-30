// push_abi.hpp -- kernel parameter blocks shared by the NVRTC program and the host
// launcher. Plain structs of ints, doubles and pointers: identical layout on both
// sides (LP64). Host code stores device addresses in the pointer members.
#pragma once
#include "push_core.hpp"

namespace hgp {

// Per-carrier diagnostic layout (structure of arrays, n entries per field).
enum DiagReal : int {
    diag_moved_dx = 0, diag_moved_dy, diag_moved_dz, diag_moved_ux, diag_moved_uy, diag_moved_uz,
    diag_hit_fraction, diag_flight_s, diag_gyro, diag_displacement_cells, diag_kinetic_fraction, diag_kinetic_plane,
    diag_real_fields
};
enum DiagInt : int {
    diag_moved_valid = 0, diag_hit_surface, diag_kinetic_valid, diag_kinetic_sink, diag_kinetic_axis, diag_int_fields
};

struct CarrierArgs {
    Mesh g;
    Policy p;
    double off[3];   // native -> field translation
    double dt;
    hgp_u32 n, node_sentinel;
    int diagnostics;
    int pad_;
    const double* x; const double* y; const double* z;
    const double* ux; const double* uy; const double* uz;
    const double* q; const double* m; const double* w;
    const double* E; const double* B;
    int* code1; int* code2; int* code3;
    hgp_u32* fail1; hgp_u32* fail2; hgp_u32* fail3;     // rank if failing, else 0xffffffff
    double* endpoint;    // [6][n]: x y z ux uy uz (native frame)
    double* start;       // [3][n]: field-frame start of the canonical path
    double* disp;        // [3][n]: canonical field displacement D'
    double* dke;         // [n]: KE(endpoint) - KE(start), binary64
    double* gyro; double* dcells;   // [n]
    hgp_u32* edge_count; // [n]: 24 * path segments (0 when the carrier fails)
    hgp_u32* start_key; hgp_u32* start_pay;   // [8n]
    double* start_omega; double* start_init; double* start_gross;   // [8n]
    hgp_u32* end_key; hgp_u32* end_pay;       // [8n]
    double* end_charge; double* end_hi; double* end_lo;              // [8n]
    double* diag_real;   // [diag_real_fields][n] or null
    int* diag_int;       // [diag_int_fields][n] or null
};

struct EmitArgs {
    Mesh g;
    double dt;
    hgp_u32 n, pad_;
    const double* q; const double* w;
    const double* start; const double* disp;   // [3][n]
    const hgp_u32* edge_count; const hgp_u32* edge_offset;
    hgp_u32* key; hgp_u32* pay; double* value;
};

struct NodeReduceArgs {
    hgp_u32 nn, pad_;
    const hgp_u32* pay; const hgp_u32* begin; const hgp_u32* end;
    const double* v0; const double* v1; const double* v2;   // per-record values
    double* out0; double* out1; double* out2;              // per-node results
};

struct EndNodeArgs {
    hgp_u32 nn, pad_;
    const hgp_u32* pay; const hgp_u32* begin; const hgp_u32* end;
    const double* charge; const double* hi; const double* lo;
    double* out_charge; double* out_hi; double* out_lo; double* out_gross;
};

struct FinalizeArgs {
    const hgp_u32* min_fail;   // [3]
    const int* code1; const int* code2; const int* code3;
    const hgp_u32* edge_total;
    hgp_u32* status;           // [8]: r1 c1 r2 c2 r3 c3 edge_total spare
};

}  // namespace hgp
