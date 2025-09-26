#include "types.hpp"

#include <filesystem>

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

namespace py = boost::python;
using namespace synthesis_backend;

boost::python::list mol_set_to_list(const MolSet &set) {
    boost::python::list result;
    for (const auto &mol : set) {
        result.append(mol);
    }
    return result;
}

BOOST_PYTHON_MODULE(types) {
    py::class_<std::filesystem::path>("Path").def(py::init<std::string>());
    py::implicitly_convertible<std::string, std::filesystem::path>();

    auto rdchem_module = py::import("rdkit.Chem.rdchem");
    auto rdChemReactions_module = py::import("rdkit.Chem.rdChemReactions");
    py::scope().attr("Mol") = rdchem_module.attr("Mol");
    py::scope().attr("Reaction") =
        rdChemReactions_module.attr("ChemicalReaction");

    py::register_ptr_to_python<Reaction_sptr>();
    py::class_<MolVector>("MolVector")
        .def(py::vector_indexing_suite<MolVector, true>());
    py::class_<ReactionVector>("ReactionVector")
        .def(py::vector_indexing_suite<ReactionVector, true>());

    py::class_<MolSet>("MolSet")
        .def("to_list", &mol_set_to_list)
        .def("__len__", &MolSet::size);
}
