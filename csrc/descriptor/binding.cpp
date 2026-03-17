#include "binding.hpp"

#include <cstdint>
#include <stdexcept>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <pybind11/detail/common.h>
#include <pybind11/detail/using_smart_holder.h>

#include "descriptor.hpp"

#include "morgan_binding.hpp"
#include "synthesis_binding.hpp"

namespace py = pybind11;
using namespace prexsyn;
using namespace prexsyn::descriptor;

static auto convert_dtype(const DataType::T &d) {
    switch (d) {
    case DataType::float32:
        return py::dtype::of<float>();
    case DataType::int64:
        return py::dtype::of<std::int64_t>();
    case DataType::bool8:
        return py::dtype::of<bool>();
    default:
        throw std::runtime_error("Unsupported data type");
    }
};

void def_module_descriptor(pybind11::module &m) {
    py::class_<MoleculeDescriptor, py::smart_holder>(m, "_MoleculeDescriptor")
        .def("dtype", [](const MoleculeDescriptor &desc) { return convert_dtype(desc.dtype()); })
        .def("size", &MoleculeDescriptor::size)
        .def("num_elements", &MoleculeDescriptor::num_elements)
        .def("size_in_bytes", &MoleculeDescriptor::size_in_bytes);

    py::class_<SynthesisDescriptor, py::smart_holder>(m, "_SynthesisDescriptor")
        .def("dtype", [](const SynthesisDescriptor &desc) { return convert_dtype(desc.dtype()); })
        .def("size", &SynthesisDescriptor::size)
        .def("num_elements", &SynthesisDescriptor::num_elements)
        .def("size_in_bytes", &SynthesisDescriptor::size_in_bytes);

    def_submodule_morgan(m);
    def_submodule_synthesis(m);
}
