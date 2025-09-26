#include "reaction_list.hpp"

#include <boost/python.hpp>

namespace py = boost::python;
using namespace synthesis_backend;

BOOST_PYTHON_MODULE(reaction_list) {
    py::class_<ReactionList>("ReactionList", py::no_init)
        .def("__getitem__", &ReactionList::get, (py::arg("index")))
        .def("__len__", &ReactionList::size)
        .def("save", &ReactionList::save, (py::arg("path")))
        .def("load", &ReactionList::load, (py::arg("path")),
             py::return_value_policy<py::manage_new_object>())
        .staticmethod("load")
        .def("from_txt", &ReactionList::from_txt, (py::arg("path")),
             py::return_value_policy<py::manage_new_object>())
        .staticmethod("from_txt");
    py::register_ptr_to_python<std::shared_ptr<ReactionList>>();
}
