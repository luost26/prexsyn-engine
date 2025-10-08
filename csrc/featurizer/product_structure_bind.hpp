#pragma once

#include "product_structure.hpp"

#include <boost/python.hpp>

namespace py = boost::python;
using namespace prexsyn_engine;

inline void product_structure_bind() {
    py::class_<ProductStructureFeaturizerOption>(
        "ProductStructureFeaturizerOption")
        .def_readwrite("fp_types", &ProductStructureFeaturizerOption::fp_types)
        .def_readwrite(
            "embedder_name_template",
            &ProductStructureFeaturizerOption::embedder_name_template)
        .def_readwrite("scaffold", &ProductStructureFeaturizerOption::scaffold)
        .def_readwrite("scaffold_fp_type",
                       &ProductStructureFeaturizerOption::scaffold_fp_type)
        .def_readwrite(
            "scaffold_embedder_name",
            &ProductStructureFeaturizerOption::scaffold_embedder_name)
        .def_readwrite("num_fragments",
                       &ProductStructureFeaturizerOption::num_fragments)
        .def_readwrite("fragment_fp_type",
                       &ProductStructureFeaturizerOption::fragment_fp_type)
        .def_readwrite(
            "fragment_embedder_name",
            &ProductStructureFeaturizerOption::fragment_embedder_name);

    py::class_<ProductStructureFeaturizer>(
        "ProductStructureFeaturizer",
        py::init<ProductStructureFeaturizerOption>(
            (py::arg("option") = ProductStructureFeaturizerOption())))
        .def("__call__", &ProductStructureFeaturizer::operator(),
             (py::arg("synthesis"), py::arg("builder")));

    py::register_ptr_to_python<std::shared_ptr<ProductStructureFeaturizer>>();
}
