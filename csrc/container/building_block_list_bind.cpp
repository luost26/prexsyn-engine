#include "building_block_list.hpp"

#include <boost/python.hpp>

namespace py = boost::python;
using namespace synthesis_backend;

BOOST_PYTHON_MODULE(building_block_list) {
    py::class_<BuildingBlockPreprocessingOption>(
        "BuildingBlockPreprocessingOption", py::init<>())
        .def_readwrite("remove_Hs",
                       &BuildingBlockPreprocessingOption::remove_Hs)
        .def_readwrite(
            "largest_fragment_only",
            &BuildingBlockPreprocessingOption::largest_fragment_only);

    py::class_<BuildingBlockList>("BuildingBlockList", py::no_init)
        .def("__getitem__", &BuildingBlockList::get, (py::arg("index")))
        .def("__len__", &BuildingBlockList::size)
        .def("save", &BuildingBlockList::save, (py::arg("path")))
        .def("load", &BuildingBlockList::load, (py::arg("path")),
             py::return_value_policy<py::manage_new_object>())
        .staticmethod("load")
        .def("from_sdf", &BuildingBlockList::from_sdf,
             (py::arg("path"),
              py::arg("option") = BuildingBlockPreprocessingOption()),
             py::return_value_policy<py::manage_new_object>())
        .staticmethod("from_sdf");
    py::register_ptr_to_python<std::shared_ptr<BuildingBlockList>>();
}
