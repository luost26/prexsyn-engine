#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

#include "pydict_builder.hpp"

namespace py = boost::python;
namespace np = boost::python::numpy;

using namespace prexsyn_engine;
BOOST_PYTHON_MODULE(feature) {
    np::initialize();

    py::class_<PyDictBuilder, std::shared_ptr<PyDictBuilder>,
               py::bases<FeatureBuilder>>("PyDictBuilder", py::init<>())
        .def("get", &PyDictBuilder::get)
        .def("erase_type", &PyDictBuilder::erase_type,
             py::return_internal_reference<>());
}
