#include <boost/python.hpp>

#include "substructures.hpp"

namespace py = boost::python;

using namespace prexsyn_engine;
BOOST_PYTHON_MODULE(substructures) {
    py::class_<BRICSFragmentsFeaturizer,
               std::shared_ptr<BRICSFragmentsFeaturizer>,
               py::bases<Featurizer>>(
        "BRICSFragmentsFeaturizer",
        py::init<std::string, std::string, unsigned int>(
            (py::arg("name"), py::arg("fp_type"),
             py::arg("max_num_fragments") = 8)))
        .def("__call__", &BRICSFragmentsFeaturizer::operator(),
             (py::arg("synthesis"), py::arg("builder")));

    py::implicitly_convertible<std::shared_ptr<BRICSFragmentsFeaturizer>,
                               std::shared_ptr<Featurizer>>();
}
