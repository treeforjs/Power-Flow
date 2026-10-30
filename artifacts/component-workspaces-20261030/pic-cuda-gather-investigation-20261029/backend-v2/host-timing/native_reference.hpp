#pragma once
#include "batch_snapshot.hpp"
#include <functional>
#include <memory>

namespace native_reference {
using namespace raikou::mfem_backend;
struct Fields;
struct Voltage;
struct RectField {
    int nx,ny;double width,height;
    mfem::Mesh mesh;
    mfem::H1_FECollection collection;
    mfem::FiniteElementSpace space;
    mfem::GridFunction field;
    Q1ScalarGather gather;
    RectField(const pic_batch::HostSnapshot&,const std::vector<double>&);
    PicVector electric(double,double);
    double value(double,double);
};
// Uses separate self/unit meshes, the original Fields/Voltage state and three
// std::function callbacks, matching the native candidate's un-fused boundary.
class NativeReference {
public:
    NativeReference(const pic_batch::HostSnapshot&,const std::filesystem::path& fixture,const std::filesystem::path& voltage_csv);
    ~NativeReference();
    void force_update(std::span<GapParticle>);
    void generate(std::span<const GapParticle>,std::span<gather_push::Outputs>);
private:
    pic_batch::Parameters p_;
    RectField self_,unit_;
    std::unique_ptr<Fields> native_;
    std::unique_ptr<Voltage> voltage_;
    double origin_,relative_midpoint_,absolute_midpoint_;
    std::function<PicVector(double,double,double)> electric_,magnetic_;
    std::function<PicVector(double,double,double,double)> axial_;
    std::pair<PicVector,PicVector> fields(const GapParticle&);
};
}
