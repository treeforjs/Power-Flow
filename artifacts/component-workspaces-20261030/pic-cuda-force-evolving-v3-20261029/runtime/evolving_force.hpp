#pragma once
#include "cuda_force_runtime.hpp"
#include <functional>
namespace pic_evolving {
using namespace raikou::mfem_backend;
enum class Mode {Off,Shadow,OptimizedCandidate};
struct NativeSelection {
    const mfem::GridFunction *psi=nullptr,*az0=nullptr,*az1=nullptr,*rate0=nullptr,*rate1=nullptr;
    const std::vector<int>* elements=nullptr;
    const std::vector<double> *xs=nullptr,*ys=nullptr;
    std::uint64_t layout_epoch=0,unit_epoch=0,lower_frame=0,upper_frame=0;
    double origin=0,absolute_midpoint=0,lower_time=0,upper_time=0,fraction=0,voltage=0;
};
using Provider=std::function<NativeSelection(double relative_midpoint)>;
// Must run before any application CUDA operation. Only parent GPU7 dispatch,
// explicit bounded replay token, one1ps window and this owned prefix are allowed.
void candidate_environment();
void validate_output(const std::filesystem::path&);
void bound_restart(double origin,double relative_start,double end,std::uint64_t step);
void configure(GapPic&,Mode,const std::filesystem::path& ptx,const std::filesystem::path& output,Provider);
Mode requested_mode();
bool active(const GapPic*);
bool shadow(const GapPic*);
void begin_force(const GapPic*,double relative_midpoint,double dt);
void reference_ready(const GapPic*);
void end_force(const GapPic*);
void force(GapPic*,std::span<GapParticle>,const mfem::GridFunction& phi,const GapPicOptions&,
           double relative_midpoint,double dt,std::span<const gather_push::Outputs> independent_reference);
void observe_step(const GapPic&,std::uint64_t step,double dt);
void observe_sources(const std::array<double,4>& spans,
                    const std::array<std::span<const double>,4>& integral,
                    const std::array<std::span<const double>,4>& temperature);
void shutdown(const char* original_error=nullptr);
}
