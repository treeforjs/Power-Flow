#include "reduced.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(_mhd_core, m) {
    m.doc() = "C++ kernels for the mhdlab research prototype";

    py::class_<mhd_core::HalfRangeMoments>(m, "HalfRangeMoments")
        .def_readonly("number_flux", &mhd_core::HalfRangeMoments::number_flux)
        .def_readonly("mean_normal_speed", &mhd_core::HalfRangeMoments::mean_normal_speed)
        .def_readonly("mean_speed", &mhd_core::HalfRangeMoments::mean_speed);

    m.def("half_range_maxwell_moments", &mhd_core::half_range_maxwell_moments,
          py::arg("temperature_k"), py::arg("mass_kg"));

    m.def("laplacian2d", &mhd_core::laplacian2d,
          py::arg("values"), py::arg("ny"), py::arg("nx"), py::arg("dx"), py::arg("dy"));
}
