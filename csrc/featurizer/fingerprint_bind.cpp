
#include <boost/python.hpp>

#include "fingerprint.hpp"

namespace py = boost::python;

using namespace prexsyn_engine;
BOOST_PYTHON_MODULE(featurizer__fingerprint) {
    py::class_<FingerprintFeaturizer, std::shared_ptr<FingerprintFeaturizer>,
               py::bases<Featurizer>>(
        "FingerprintFeaturizer", py::init<std::string, std::string>(
                                     (py::arg("name"), py::arg("fp_type"))))
        .def("__call__", &FingerprintFeaturizer::operator(),
             (py::arg("synthesis"), py::arg("builder")))
        .def_readonly("name", &FingerprintFeaturizer::name);
}
