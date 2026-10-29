#pragma once
#include "raikou/mfem_gap_pic.h"
#include <filesystem>
#include <iosfwd>

namespace pic_capture {
// A diagnostic snapshot writer, NOT a field evaluator. It reads the native
// provider's already-selected state after the original force loop has run.
using NativeWriter=void(*)(const void*,const std::filesystem::path&,double,std::ostream&);
void configure(const std::filesystem::path&,const mfem::GridFunction& unit_psi,
               double global_origin,const void* provider,NativeWriter);
void begin(const raikou::mfem_backend::GapPic&,double dt);
void record(const raikou::mfem_backend::GapParticle& pushed,
            const raikou::mfem_backend::PicVector& E,const raikou::mfem_backend::PicVector& B);
void finish(const raikou::mfem_backend::GapPic&,double dt);
}
