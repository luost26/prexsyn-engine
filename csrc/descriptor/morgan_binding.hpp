#pragma once

#include <pybind11/pybind11.h>

#include "base.hpp"
#include "morgan.hpp"

namespace prexsyn::descriptor {

namespace py = pybind11;

inline void def_submodule_morgan(py::module &m) {
    py::class_<MorganFingerprint, MoleculeDescriptor, py::smart_holder>(m, "MorganFingerprint")
        .def_static("ecfp4", &MorganFingerprint::ecfp4)
        .def_static("fcfp4", &MorganFingerprint::fcfp4);
}

} // namespace prexsyn::descriptor
