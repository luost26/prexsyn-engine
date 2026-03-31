#include "bind.hpp"

#include <cstddef>
#include <memory>
#include <optional>

#include <pybind11/cast.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

#include "../chemspace/chemspace.hpp"
#include "enumerator.hpp"

namespace py = pybind11;
using namespace prexsyn;
using namespace prexsyn::enumerator;

void def_module_enumerator(pybind11::module &m) {
    py::class_<EnumeratorConfig, py::smart_holder>(m, "EnumeratorConfig")
        .def(py::init<>())
        .def_readwrite("max_building_blocks", &EnumeratorConfig::max_building_blocks)
        .def_readwrite("heavy_atom_limit", &EnumeratorConfig::heavy_atom_limit)
        .def_readwrite("selectivity_cutoff", &EnumeratorConfig::selectivity_cutoff)
        .def_readwrite("max_outcomes_per_reaction", &EnumeratorConfig::max_outcomes_per_reaction);

    py::class_<RandomEnumerator, py::smart_holder>(m, "RandomEnumerator")
        .def(py::init<std::shared_ptr<chemspace::ChemicalSpace>, const RandomEnumerator::Config &,
                      std::optional<size_t>>(),
             py::arg("chemical_space"), py::arg("config") = kDefaultEnumeratorConfig,
             py::arg("random_seed") = std::nullopt)
        .def("next", &RandomEnumerator::next)
        .def("next_with_product", &RandomEnumerator::next_with_product);
}
