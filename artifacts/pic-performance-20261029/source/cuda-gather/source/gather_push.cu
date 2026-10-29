#include "gather_push_core.hpp"
// Tables contain immutable field snapshots for one already-selected midpoint.
// No std::function, MFEM object, source state, or particle ownership crosses here.
extern "C" __global__ void gather_push_kernel(gather_push::Tables tables,
    const gather_push::Inputs* inputs,gather_push::Outputs* outputs,unsigned int n) {
    const unsigned int i=blockIdx.x*blockDim.x+threadIdx.x;
    if(i<n) outputs[i]=gather_push::evaluate(tables,inputs[i]);
}
