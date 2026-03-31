#include "bind.hpp"

#include <cstddef>
#include <span>
#include <vector>

#include <omp.h>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <pybind11/detail/common.h>
#include <pybind11/detail/using_smart_holder.h>

#include "../chemistry/chemistry.hpp"
#include "../chemspace/chemspace.hpp"
#include "../utility/data_type_bind.hpp"
#include "descriptor.hpp"

#include "morgan_bind.hpp"
#include "synthesis_bind.hpp"

namespace py = pybind11;
using namespace prexsyn;
using namespace prexsyn::descriptor;

template <typename DescriptorType, typename InputType>
static auto get_numpy_array(const DescriptorType &desc, const InputType &input) {
    auto dtype = data_type_to_numpy_dtype(desc.dtype());
    auto shape = desc.size();
    auto arr = py::array(dtype, shape);

    auto span =
        std::span<std::byte>(reinterpret_cast<std::byte *>(arr.mutable_data()), arr.nbytes());
    desc(input, span);
    return arr;
}

template <typename DescriptorType, typename InputType>
static auto get_batched_numpy_array(const DescriptorType &desc,
                                    const std::vector<InputType> &inputs) {
    auto dtype = data_type_to_numpy_dtype(desc.dtype());
    auto batched_shape = desc.size();
    batched_shape.insert(batched_shape.begin(), inputs.size());
    auto arr = py::array(dtype, batched_shape);

    auto *base_ptr = reinterpret_cast<std::byte *>(arr.mutable_data());

#pragma omp parallel for
    for (size_t i = 0; i < inputs.size(); ++i) {
        auto offset = i * desc.size_in_bytes();
        auto span_i = std::span<std::byte>(base_ptr + offset, desc.size_in_bytes());
        desc(inputs[i], span_i);
    }
    return arr;
}

void def_module_descriptor(pybind11::module &m) {
    py::class_<MoleculeDescriptor, py::smart_holder>(m, "_MoleculeDescriptor")
        .def("dtype",
             [](const MoleculeDescriptor &desc) { return data_type_to_numpy_dtype(desc.dtype()); })
        .def("size", &MoleculeDescriptor::size)
        .def("num_elements", &MoleculeDescriptor::num_elements)
        .def("size_in_bytes", &MoleculeDescriptor::size_in_bytes)
        .def("__call__", &get_numpy_array<MoleculeDescriptor, Molecule>)
        .def("__call__", &get_batched_numpy_array<MoleculeDescriptor, Molecule>);

    py::class_<SynthesisDescriptor, py::smart_holder>(m, "_SynthesisDescriptor")
        .def("dtype",
             [](const SynthesisDescriptor &desc) { return data_type_to_numpy_dtype(desc.dtype()); })
        .def("size", &SynthesisDescriptor::size)
        .def("num_elements", &SynthesisDescriptor::num_elements)
        .def("size_in_bytes", &SynthesisDescriptor::size_in_bytes)
        .def("__call__", &get_numpy_array<SynthesisDescriptor, chemspace::ChemicalSpaceSynthesis>)
        .def("__call__",
             &get_batched_numpy_array<SynthesisDescriptor, chemspace::ChemicalSpaceSynthesis>);

    def_submodule_morgan(m);
    def_submodule_synthesis(m);
}
