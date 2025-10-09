#pragma once

#include "rdkit_descriptors.hpp"

#include <boost/python.hpp>

namespace py = boost::python;
using namespace prexsyn_engine;

BOOST_PYTHON_MODULE(featurizer__rdkit_descriptors) {
    py::class_<RDKitDescriptorsFeaturizer,
               std::shared_ptr<RDKitDescriptorsFeaturizer>,
               py::bases<Featurizer>>(
        "RDKitDescriptorsFeaturizer",
        py::init<std::string, unsigned int, std::vector<std::string>>(
            (py::arg("name") = "rdkit_descriptors",
             py::arg("num_evaluated_descriptors") = 4,
             py::arg("descriptor_names") = SUPPORTED_RDKIT_DESCRIPTORS)))
        .def_readonly("name", &RDKitDescriptorsFeaturizer::name)
        .def_readonly("num_evaluated_descriptors",
                      &RDKitDescriptorsFeaturizer::num_evaluated_descriptors)
        .def_readonly("descriptor_names",
                      &RDKitDescriptorsFeaturizer::descriptor_names)
        .def("max_property_index",
             &RDKitDescriptorsFeaturizer::max_property_index)
        .def("__call__", &RDKitDescriptorsFeaturizer::operator(),
             (py::arg("synthesis"), py::arg("builder")));
}
