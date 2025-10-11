#include <boost/python.hpp>

#include "scaffold.hpp"

namespace py = boost::python;

using namespace prexsyn_engine;
BOOST_PYTHON_MODULE(scaffold) {
    py::class_<MurckoScaffoldFeaturizer,
               std::shared_ptr<MurckoScaffoldFeaturizer>,
               py::bases<Featurizer>>(
        "MurckoScaffoldFeaturizer", py::init<std::string, std::string>(
                                        (py::arg("name"), py::arg("fp_type"))))
        .def("__call__", &MurckoScaffoldFeaturizer::operator(),
             (py::arg("synthesis"), py::arg("builder")))
        .def_readonly("name", &MurckoScaffoldFeaturizer::name);

    py::implicitly_convertible<std::shared_ptr<MurckoScaffoldFeaturizer>,
                               std::shared_ptr<Featurizer>>();
}
