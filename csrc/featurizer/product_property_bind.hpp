#pragma once

#include "product_property.hpp"

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

namespace py = boost::python;
using namespace synthesis_backend;

inline void product_property_bind() {
    py::class_<std::vector<std::string>>("RDKitPropertyList")
        .def(py::vector_indexing_suite<std::vector<std::string>>());
    py::class_<ProductRDKitPropertyFeaturizerOption>(
        "ProductRDKitPropertyFeaturizerOption")
        .def_readwrite("name", &ProductRDKitPropertyFeaturizerOption::name)
        .def_readwrite(
            "num_evaluated_properties",
            &ProductRDKitPropertyFeaturizerOption::num_evaluated_properties)
        .def_readwrite(
            "rdkit_property_index_offset",
            &ProductRDKitPropertyFeaturizerOption::rdkit_property_index_offset)
        .def_readwrite("rdkit_properties",
                       &ProductRDKitPropertyFeaturizerOption::rdkit_properties);

    py::class_<ProductRDKitPropertyFeaturizer>(
        "ProductRDKitPropertyFeaturizer",
        py::init<ProductRDKitPropertyFeaturizerOption>(
            (py::arg("option") = ProductRDKitPropertyFeaturizerOption())))
        .def("__call__", &ProductRDKitPropertyFeaturizer::operator(),
             (py::arg("synthesis"), py::arg("builder")))
        .def("max_property_index",
             &ProductRDKitPropertyFeaturizer::max_property_index);

    py::register_ptr_to_python<
        std::shared_ptr<ProductRDKitPropertyFeaturizer>>();
}
