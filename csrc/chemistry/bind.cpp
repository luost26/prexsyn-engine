#include "bind.hpp"

#include <cstddef>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

#include <pybind11/detail/common.h>
#include <pybind11/detail/using_smart_holder.h>
#include <pyerrors.h>

#include "chemistry.hpp"

namespace py = pybind11;
using namespace prexsyn;

static void def_molecule(py::module &m) {
    py::class_<Molecule, py::smart_holder>(m, "Molecule")
        .def_static("from_smiles", &Molecule::from_smiles)
        .def_static("from_rdkit_mol",
                    [](const py::object &rdk_mol) -> std::shared_ptr<Molecule> {
                        py::tuple init_args = rdk_mol.attr("__getinitargs__")();
                        py::bytes rdk_pickle = init_args[0];
                        std::string pickle_str(rdk_pickle);
                        return Molecule::from_rdkit_pickle(pickle_str);
                    })
        .def("smiles", &Molecule::smiles, py::return_value_policy::copy)
        .def("to_rdkit_mol",
             [](const Molecule &mol) -> py::object {
                 py::bytes rdk_pickle(mol.rdkit_pickle());
                 py::module rdchem = py::module::import("rdkit.Chem.rdchem");
                 py::object rdk_mol = rdchem.attr("Mol").attr("__new__")(rdchem.attr("Mol"));
                 rdk_mol.attr("__init__")(rdk_pickle);
                 return rdk_mol;
             })
        .def("num_heavy_atoms", &Molecule::num_heavy_atoms)
        .def("largest_fragment", &Molecule::largest_fragment)
        .def(py::pickle([](const Molecule &mol) { return py::bytes(mol.serialize()); },
                        [](const py::bytes &pickle) {
                            std::string pickle_str(pickle);
                            return Molecule::deserialize(pickle_str);
                        }))
        .def("__repr__", [](const Molecule &mol) { return "<Molecule " + mol.smiles() + ">"; });

    py::register_exception<MoleculeError>(m, "MoleculeError", PyExc_RuntimeError);
}

static void def_reaction(py::module &m) {

    py::class_<ReactionOutcome, py::smart_holder>(m, "ReactionOutcome")
        .def_readonly("products", &ReactionOutcome::products)
        .def("empty", &ReactionOutcome::empty)
        .def("num_products", &ReactionOutcome::num_products)
        .def("main_product", &ReactionOutcome::main_product);

    py::class_<ReactionOutcomeWithReactantAssignment, ReactionOutcome, py::smart_holder>(
        m, "ReactionOutcomeWithReactantAssignment")
        .def_readonly("reactant_names", &ReactionOutcomeWithReactantAssignment::reactant_names);

    py::class_<Reaction, py::smart_holder>(m, "Reaction")
        .def_static("from_smarts",
                    py::overload_cast<const std::string &, const std::vector<std::string> &>(
                        &Reaction::from_smarts),
                    py::arg("smarts"), py::arg("reactant_names"))
        .def_static("from_smarts", py::overload_cast<const std::string &>(&Reaction::from_smarts))
        .def("num_reactants", &Reaction::num_reactants)
        .def("reactant_names", &Reaction::reactant_names,
             py::return_value_policy::reference_internal)
        .def("match_reactants",
             [](const Reaction &rxn, std::shared_ptr<Molecule> &mol) {
                 if (!mol) {
                     throw ReactionError("Molecule pointer is null");
                 }
                 auto matches = rxn.match_reactants(*mol);
                 py::list py_matches;
                 for (const auto &match : matches) {
                     py::dict py_match;
                     py_match["index"] = match.index;
                     py_match["name"] = std::string(match.name);
                     py_match["count"] = match.count;
                     py_matches.append(py_match);
                 }
                 return py_matches;
             })
        .def("apply",
             py::overload_cast<const std::map<std::string, std::shared_ptr<Molecule>> &, bool>(
                 &Reaction::apply, py::const_),
             py::arg("reactants"), py::arg("ignore_errors") = false)
        .def("apply",
             py::overload_cast<const std::vector<std::shared_ptr<Molecule>> &, bool>(
                 &Reaction::apply, py::const_),
             py::arg("reactants"), py::arg("ignore_errors") = false)
        .def(py::pickle([](const Reaction &rxn) { return py::bytes(rxn.serialize()); },
                        [](const py::bytes &pickle) {
                            std::string pickle_str(pickle);
                            return Reaction::deserialize(pickle_str);
                        }))
        .def("__repr__", [](const Reaction &rxn) {
            std::string repr = "<Reaction(";
            for (size_t i = 0; i < rxn.num_reactants(); ++i) {
                repr += rxn.reactant_names()[i];
                if (i < rxn.num_reactants() - 1) {
                    repr += ", ";
                }
            }
            repr += ")>";
            return repr;
        });

    py::register_exception<ReactionError>(m, "ReactionError", PyExc_RuntimeError);
}

static void def_synthesis(py::module &m) {
    py::class_<SynthesisNode::PrecursorMolecule>(m, "PrecursorMolecule")
        .def_readonly("precursor_index", &SynthesisNode::PrecursorMolecule::precursor_index)
        .def_readonly("reactant_name", &SynthesisNode::PrecursorMolecule::reactant_name)
        .def_readonly("precursor_node", &SynthesisNode::PrecursorMolecule::precursor_node)
        .def_readonly("item_index", &SynthesisNode::PrecursorMolecule::item_index)
        .def_readonly("molecule", &SynthesisNode::PrecursorMolecule::molecule);

    py::class_<SynthesisNode, py::smart_holder>(m, "SynthesisNode")
        .def("index", &SynthesisNode::index)
        .def("size", &SynthesisNode::size)
        .def("precursor_nodes", &SynthesisNode::precursor_nodes,
             py::return_value_policy::reference_internal)
        .def("at", &SynthesisNode::at, py::arg("i"), py::return_value_policy::reference_internal)
        .def("precursors", &SynthesisNode::precursors, py::arg("index"),
             py::return_value_policy::reference_internal);

    py::class_<Synthesis, py::smart_holder>(m, "Synthesis")
        .def(py::init<>())
        .def("nodes", &Synthesis::nodes, py::return_value_policy::reference_internal)
        .def("stack_size", &Synthesis::stack_size)
        .def("stack_top", &Synthesis::stack_top, py::arg("i") = 0,
             py::return_value_policy::reference_internal)
        .def(
            "push_molecule",
            [](Synthesis &s, const std::shared_ptr<Molecule> &mol) { s.push(mol); },
            py::arg("molecule"))
        .def(
            "push_reaction",
            [](Synthesis &s, const std::shared_ptr<Reaction> &rxn,
               std::optional<size_t> max_outcomes) { s.push(rxn, max_outcomes); },
            py::arg("reaction"), py::arg("max_outcomes"))
        .def("undo", &Synthesis::undo);

    py::register_exception<SynthesisError>(m, "SynthesisError", PyExc_RuntimeError);
}

void def_module_chemistry(pybind11::module &m) {
    def_molecule(m);
    def_reaction(m);
    def_synthesis(m);
}
