#include "synthesis.hpp"

#include <stdexcept>
#include <variant>

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "utils/pickling.hpp"

namespace py = boost::python;
using namespace synthesis_backend;

py::object postfix_notation_getitem(const PostfixNotation &pfn, size_t index) {
    if (index >= pfn.size()) {
        throw std::out_of_range(
            "Index out of range in PostfixNotation::getitem: " +
            std::to_string(index) + " >= " + std::to_string(pfn.size()));
    }
    auto item = pfn[index];
    if (std::holds_alternative<Mol_sptr>(item)) {
        return py::object(std::get<Mol_sptr>(item));
    } else {
        return py::object(std::get<Reaction_sptr>(item));
    }
}

py::list postfix_notation_to_list(const PostfixNotation &pfn) {
    py::list result;
    for (size_t index = 0; index < pfn.size(); ++index) {
        result.append(postfix_notation_getitem(pfn, index));
    }
    return result;
}

py::tuple synthesis_vector_get_all_products(const SynthesisVector &syntheses) {
    py::list all_products;
    py::list synthesis_indices;
    for (size_t i = 0; i < syntheses.size(); ++i) {
        const auto &synthesis = syntheses[i];
        if (synthesis->stack_size() == 0) {
            continue;
        }
        const auto &top_molset = synthesis->top();
        for (const auto &mol : top_molset) {
            all_products.append(mol);
            synthesis_indices.append(i);
        }
    }
    return py::make_tuple(all_products, synthesis_indices);
}

BOOST_PYTHON_MODULE(synthesis) {
    py::enum_<PostfixNotation::ItemType>("PostfixNotationItemType")
        .value("Molecule", PostfixNotation::ItemType::Molecule)
        .value("Reaction", PostfixNotation::ItemType::Reaction)
        .export_values();

    py::class_<PostfixNotation>("PostfixNotation", py::init<>())
        .def("type", &PostfixNotation::type, (py::arg("index")))
        .def("append_mol",
             (void (PostfixNotation::*)(
                 const Mol_sptr &))&PostfixNotation::append,
             (py::arg("mol")))
        .def("append_reaction",
             (void (PostfixNotation::*)(
                 const Reaction_sptr &))&PostfixNotation::append,
             (py::arg("reaction")))
        .def("extend", &PostfixNotation::extend, (py::arg("other")))
        .def("__getitem__", &postfix_notation_getitem, ((py::arg("index"))))
        .def("__len__", &PostfixNotation::size)
        .def("to_list", &postfix_notation_to_list)
        .def("count_reactions", &PostfixNotation::count_reactions)
        .def("count_building_blocks", &PostfixNotation::count_building_blocks)
        .def_pickle(pickle_suite<PostfixNotation>());

    py::class_<push_reaction_exception> push_reaction_exception_class(
        "PushReactionException", py::no_init);
    push_reaction_exception_class.def("what", &push_reaction_exception::what);
    py::register_exception_translator<push_reaction_exception>(
        [&](const push_reaction_exception &e) {
            PyErr_SetString(push_reaction_exception_class.ptr(), e.what());
        });

    py::class_<std::vector<MolSet>>("MolSetVector")
        .def(py::vector_indexing_suite<std::vector<MolSet>, true>());

    py::class_<Synthesis>("Synthesis", py::init<>())
        .def("get_postfix_notation", &Synthesis::get_postfix_notation,
             py::return_internal_reference())
        .def("get_stack", &Synthesis::get_stack,
             py::return_internal_reference())
        .def("push_mol",
             (void (Synthesis::*)(const Mol_sptr &))&Synthesis::push,
             (py::arg("mol")))
        .def("push_reaction",
             (void (Synthesis::*)(const Reaction_sptr &,
                                  size_t))&Synthesis::push,
             (py::arg("reaction"), py::arg("max_products") = 8))
        .def("push_synthesis",
             (void (Synthesis::*)(const Synthesis &))&Synthesis::push,
             (py::arg("synthesis")))
        .def("top",
             (MolSet (Synthesis::*)(std::vector<MolSet>::size_type) const) &
                 Synthesis::top,
             (py::arg("index") = 0))
        .def("stack_size", &Synthesis::stack_size)
        .def("count_reactions", &Synthesis::count_reactions)
        .def("count_building_blocks", &Synthesis::count_building_blocks)
        .def_pickle(pickle_suite<Synthesis>());
    py::register_ptr_to_python<Synthesis_sptr>();
    py::class_<SynthesisVector>("SynthesisVector")
        .def(py::vector_indexing_suite<SynthesisVector, true>())
        .def("get_all_products", &synthesis_vector_get_all_products);
}
