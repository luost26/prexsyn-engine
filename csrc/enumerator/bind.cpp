#include "bind.hpp"

#include <cstddef>
#include <map>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include <pybind11/cast.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

#include "../chemspace/chemspace.hpp"
#include "../descriptor/descriptor.hpp"
#include "../utility/data_type_bind.hpp"
#include "enumerator.hpp"

namespace py = pybind11;
using namespace prexsyn;
using namespace prexsyn::enumerator;

void def_module_enumerator(pybind11::module &m) {
    py::class_<EnumeratorConfig, py::smart_holder>(m, "EnumeratorConfig")
        .def(py::init<>())
        .def_readwrite("max_building_blocks", &EnumeratorConfig::max_building_blocks)
        .def_readwrite("heavy_atom_limit", &EnumeratorConfig::heavy_atom_limit);

    py::class_<RandomEnumerator, py::smart_holder>(m, "RandomEnumerator")
        .def(py::init<std::shared_ptr<chemspace::ChemicalSpace>, const RandomEnumerator::Config &,
                      std::optional<size_t>>(),
             py::arg("chemical_space"), py::arg("config") = default_config,
             py::arg("random_seed") = std::nullopt)
        .def("next", &RandomEnumerator::next)
        .def("next_with_product", &RandomEnumerator::next_with_product);
}
