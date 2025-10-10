
#include <boost/python/class.hpp>
#include <boost/python/import.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <memory>

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

#include "base.hpp"

namespace py = boost::python;
namespace np = boost::python::numpy;

using namespace prexsyn_engine;
BOOST_PYTHON_MODULE(base) {
    np::initialize();

    std::ignore =
        py::class_<Featurizer, std::shared_ptr<Featurizer>, boost::noncopyable>(
            "Featurizer", py::no_init);

    auto featurizer_set_class =
        py::class_<FeaturizerSet, std::shared_ptr<FeaturizerSet>>(
            "FeaturizerSet", py::init<>())
            .def("__call__", &FeaturizerSet::operator(),
                 (py::arg("synthesis"), py::arg("builder")))
            .def("add", &FeaturizerSet::add, py::return_internal_reference<>(),
                 (py::arg("featurizer")));
}
