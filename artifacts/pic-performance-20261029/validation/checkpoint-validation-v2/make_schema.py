import pathlib,json
out=pathlib.Path(__file__).parent
options=[("nx","i32"),("ny","i32")]+[(n,"f64")for n in ["width_m","gap_m","depth_m","lower_potential_V","upper_potential_V"]]+[("magnetic_field_T["+str(k)+"]","f64")for k in range(3)]+[("distributed_fields","bool"),("finite_foils","bool")]+[(foil+"."+n,"f64")for foil in ["lower_foil","upper_foil"]for n in ["left_m","right_m","bottom_m","top_m"]]+[("finite_axial_extent","bool"),("short_z_m","f64"),("feed_z_m","f64")]
diag="time_s particle_charge_C deposited_charge_C lower_field_charge_C upper_field_charge_C gauss_defect_C particle_charge_defect_C field_energy_J particle_kinetic_energy_J lower_absorbed_charge_C upper_absorbed_charge_C lower_absorbed_energy_J upper_absorbed_energy_J lower_injected_charge_C upper_injected_charge_C injected_energy_J authored_charge_C authored_kinetic_energy_J maximum_nodal_continuity_defect_C escaped_charge_C escaped_kinetic_energy_J short_absorbed_charge_C short_absorbed_energy_J axial_escaped_charge_C poisson_residual_C_per_m".split()
particle=[("id","u64")]+[(name+"["+str(k)+"]","f64")for name in ["position_m","proper_velocity_m_s"]for k in range(3)]+[(n,"f64")for n in ["charge_C","mass_kg","weight"]]
lines=["#include \"raikou/mfem_gap_pic.h\"","#include <cstddef>","#include <bit>","#include <iostream>","#include <type_traits>","using namespace raikou::mfem_backend;","int main(){","static_assert(sizeof(double)==8 && sizeof(int)==4 && sizeof(bool)==1);","static_assert(std::is_standard_layout_v<GapPicOptions> && std::is_standard_layout_v<GapPicDiagnostics> && std::is_standard_layout_v<GapParticle>);",r"std::cout<<\"{\\\"endian\\\":\\\"\"<<(std::endian::native==std::endian::little?\"little\":\"big\")<<\"\\\",\\\"options_size\\\":\"<<sizeof(GapPicOptions)<<\",\\\"diagnostics_size\\\":\"<<sizeof(GapPicDiagnostics)<<\",\\\"particle_size\\\":\"<<sizeof(GapParticle)<<\",\";".replace("\\\"","\\\"")]
# Use ordinary quoted C++ text, keeping output as strict JSON.
lines=lines[:-1]+["std::cout<<\"{\\\"endian\\\":\\\"\"<<(std::endian::native==std::endian::little?\"little\":\"big\")<<\"\\\",\\\"options_size\\\":\"<<sizeof(GapPicOptions)<<\",\\\"diagnostics_size\\\":\"<<sizeof(GapPicDiagnostics)<<\",\\\"particle_size\\\":\"<<sizeof(GapParticle)<<\",\";"]
for kind,typ,fields in [("options","GapPicOptions",options),("diagnostics","GapPicDiagnostics",[(n,"f64")for n in diag]),("particle","GapParticle",particle)]:
    lines.append("std::cout<<\"\\\""+kind+"\\\":[\";")
    for i,(name,fmt) in enumerate(fields):
        if "." in name:
            outer,inner=name.split(".");offset="offsetof("+typ+","+outer+")+offsetof(GapFoil,"+inner+")"
        elif "[" in name:
            outer,index=name[:-1].split("[");offset="offsetof("+typ+","+outer+")+"+index+"*sizeof(double)"
        else:offset="offsetof("+typ+","+name+")"
        prefix=("," if i else "")+"{\"name\":\""+name+"\",\"kind\":\""+fmt+"\",\"offset\":"
        lines.append("std::cout<<"+json.dumps(prefix)+"<<("+offset+")<<\"}\";")
    lines.append("std::cout<<\"]"+("}" if kind=="particle" else ",")+"\";")
lines.append("}")
(out/"schema_probe.cpp").write_text("\n".join(lines)+"\n")
(out/"progress.json").write_text(json.dumps({"status":"preparing_abi_schema","cpu":6,"nice":15,"gpu_work":False},indent=2)+"\n")
