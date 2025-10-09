#include <boost/python.hpp>

#include "pydict_builder.hpp"

namespace py = boost::python;

using namespace prexsyn_engine;
BOOST_PYTHON_MODULE(feature) {
    py::class_<PyDictBuilder, std::shared_ptr<PyDictBuilder>>("PyDictBuilder",
                                                              py::init<>())
        .def("get", &PyDictBuilder::get)
        .def("erase_type", &PyDictBuilder::erase_type,
             py::return_internal_reference<>());
}
