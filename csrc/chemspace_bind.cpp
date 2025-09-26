#include "chemspace.hpp"

#include <memory>

#include <boost/python.hpp>
#include <random>

namespace py = boost::python;
using namespace synthesis_backend;

BOOST_PYTHON_MODULE(chemspace) {
    py::class_<ChemicalSpaceDefinition>("ChemicalSpaceDefinition", py::no_init)
        .def("get_primary_building_blocks",
             &ChemicalSpaceDefinition::get_primary_building_blocks,
             py::return_value_policy<py::return_by_value>())
        .def("get_secondary_building_blocks",
             &ChemicalSpaceDefinition::get_secondary_building_blocks,
             py::return_value_policy<py::return_by_value>())
        .def("get_reactions", &ChemicalSpaceDefinition::get_reactions,
             py::return_value_policy<py::return_by_value>())
        .def("save", &ChemicalSpaceDefinition::save, (py::arg("dirpath")));
    py::register_ptr_to_python<std::shared_ptr<ChemicalSpaceDefinition>>();

    py::class_<ChemicalSpaceDefinitionBuilder, boost::noncopyable>(
        "ChemicalSpaceDefinitionBuilder", py::init<>())
        .def("building_blocks_from_sdf",
             &ChemicalSpaceDefinitionBuilder::building_blocks_from_sdf,
             (py::arg("path"),
              py::arg("option") = BuildingBlockPreprocessingOption()),
             py::return_internal_reference())
        .def("building_blocks_from_cache",
             &ChemicalSpaceDefinitionBuilder::building_blocks_from_cache,
             (py::arg("path")), py::return_internal_reference())
        .def("reactions_from_txt",
             &ChemicalSpaceDefinitionBuilder::reactions_from_txt,
             (py::arg("path")), py::return_internal_reference())
        .def("reactions_from_cache",
             &ChemicalSpaceDefinitionBuilder::reactions_from_cache,
             (py::arg("path")), py::return_internal_reference())
        .def("secondary_building_blocks_from_single_reaction",
             &ChemicalSpaceDefinitionBuilder::
                 secondary_building_blocks_from_single_reaction,
             py::return_internal_reference())
        .def("secondary_building_blocks_from_cache",
             &ChemicalSpaceDefinitionBuilder::
                 secondary_building_blocks_from_cache,
             (py::arg("path")), py::return_internal_reference())
        .def("no_secondary_building_blocks",
             &ChemicalSpaceDefinitionBuilder::no_secondary_building_blocks,
             py::return_internal_reference())
        .def("build_primary_index",
             &ChemicalSpaceDefinitionBuilder::build_primary_index,
             py::return_internal_reference())
        .def("build_secondary_index",
             &ChemicalSpaceDefinitionBuilder::build_secondary_index,
             py::return_internal_reference())
        .def("all_from_cache", &ChemicalSpaceDefinitionBuilder::all_from_cache,
             (py::arg("dirpath")), py::return_internal_reference())
        .def("build", &ChemicalSpaceDefinitionBuilder::build,
             py::return_value_policy<py::return_by_value>() /* shared_ptr */)
        .def("count_secondary_building_blocks",
             &ChemicalSpaceDefinitionBuilder::count_secondary_building_blocks)
        .def("get_secondary_building_block",
             &ChemicalSpaceDefinitionBuilder::get_secondary_building_block,
             (py::arg("index")));

    py::class_<SynthesisGeneratorOption>("SynthesisGeneratorOption",
                                         py::init<>())
        .def_readwrite("num_reactions_cutoff",
                       &SynthesisGeneratorOption::num_reactions_cutoff)
        .def_readwrite("num_product_atoms_cutoff",
                       &SynthesisGeneratorOption::num_product_atoms_cutoff);

    py::class_<SynthesisGenerator>(
        "SynthesisGenerator",
        py::init<std::shared_ptr<ChemicalSpaceDefinition>,
                 SynthesisGeneratorOption, std::mt19937::result_type>(
            (py::arg("csd"), py::arg("option") = SynthesisGeneratorOption(),
             py::arg("seed") = std::random_device{}())))
        .def_readwrite("option", &SynthesisGenerator::option)
        .def("next", &SynthesisGenerator::next);
    py::register_ptr_to_python<std::shared_ptr<SynthesisGenerator>>();
}
