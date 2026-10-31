// push_abi_v3.hpp -- kernel parameter blocks of the Hall v3 motion kernels: the v1
// blocks of hall-gpu-resident-v1-20261030/push/push_abi.hpp (unchanged, included)
// plus the carrier block with the sink geometry and sink outputs, and the endpoint
// node reduction that excludes sink records from the free-charge census.
#pragma once
#include "push_abi.hpp"
#include "sink_core.hpp"

namespace hgp {

struct CarrierArgsV3 {
    CarrierArgs v1;            // every v1 operand and output, same meaning
    SinkGeometry sg;           // kinetic policy, translation, sink node support
    int* sink_kind;            // [n]: KineticSink of a truncated carrier, -1 otherwise
    int* sink_axis;            // [n]
    double* sink_real;         // [7][n]: flight, discarded flight, fraction, plane, segment x, y, z
    const double* flight;      // [n]: carrier flight F (dt; newborn end - tau), null = dt for all
    const int* face;           // [n]: foil of a newborn face start, -1 otherwise; null = none
    const int* newborn;        // [n]: 1 for a newborn of the step; null = none
};

struct EndNodeArgsV3 {
    hgp_u32 nn, pad_;
    const hgp_u32* pay; const hgp_u32* begin; const hgp_u32* end;
    const double* charge; const double* hi; const double* lo;
    const int* sink_kind;      // per carrier; records of sink carriers count in the endpoint charge only
    double* out_charge; double* out_hi; double* out_lo; double* out_gross;
};

}  // namespace hgp
