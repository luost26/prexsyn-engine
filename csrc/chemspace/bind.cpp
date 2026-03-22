#include "bind.hpp"

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>

#include <pybind11/native_enum.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

#include <pybind11/detail/common.h>
#include <pybind11/detail/using_smart_holder.h>
#include <pyerrors.h>

#include "chemspace.hpp"

namespace py = pybind11;
using namespace prexsyn::chemspace;

template <typename T>
static void serialize_to_file(const T &obj, const std::filesystem::path &path) {
    std::ofstream ofs(path, std::ios::binary);
    if (!ofs) {
        throw std::runtime_error("failed to open file for writing: " + path.string());
    }
    obj.serialize(ofs);
}

template <typename T>
static std::unique_ptr<T> deserialize_from_file(const std::filesystem::path &path) {
    std::ifstream ifs(path, std::ios::binary);
    if (!ifs) {
        throw std::runtime_error("failed to open file for reading: " + path.string());
    }
    return T::deserialize(ifs);
}

static void def_postfix_notation(py::module &m) {
    py::native_enum<PostfixNotation::Token::Type>(m, "PostfixNotationTokenType", "enum.Enum")
        .value("BuildingBlock", PostfixNotation::Token::Type::BuildingBlock)
        .value("Reaction", PostfixNotation::Token::Type::Reaction)
        .export_values()
        .finalize();

    py::class_<PostfixNotation::Token>(m, "PostfixNotationToken")
        .def(py::init<>())
        .def_readwrite("index", &PostfixNotation::Token::index)
        .def_readwrite("type", &PostfixNotation::Token::type);

    py::class_<PostfixNotation, py::smart_holder>(m, "PostfixNotation")
        .def(py::init<>())
        .def("tokens", &PostfixNotation::tokens, py::return_value_policy::reference_internal)
        .def("append", &PostfixNotation::append, py::arg("index"), py::arg("type"))
        .def("extend", &PostfixNotation::extend, py::arg("tokens"))
        .def("pop_back", &PostfixNotation::pop_back)
        .def("__repr__", [](const PostfixNotation &pn) {
            std::ostringstream oss;
            oss << pn;
            return oss.str();
        });
}

static void def_bb_lib(py::module &m) {
    py::class_<BuildingBlockEntry>(m, "BuildingBlockEntry")
        .def(py::init<>())
        .def_readwrite("molecule", &BuildingBlockEntry::molecule)
        .def_readwrite("identifier", &BuildingBlockEntry::identifier)
        .def_readwrite("labels", &BuildingBlockEntry::labels);

    py::class_<BuildingBlockItem, BuildingBlockEntry>(m, "BuildingBlockItem")
        .def(py::init<>())
        .def_readonly("index", &BuildingBlockItem::index);

    py::class_<BuildingBlockLibrary, py::smart_holder>(m, "BuildingBlockLibrary")
        .def(py::init<>())
        .def("size", &BuildingBlockLibrary::size)
        .def("get",
             py::overload_cast<BuildingBlockLibrary::Index>(&BuildingBlockLibrary::get, py::const_),
             py::arg("index"), py::return_value_policy::reference_internal)
        .def("get", py::overload_cast<const std::string &>(&BuildingBlockLibrary::get, py::const_),
             py::arg("identifier"), py::return_value_policy::reference_internal)
        .def("add", &BuildingBlockLibrary::add, py::arg("entry"))
        .def("serialize", &serialize_to_file<BuildingBlockLibrary>, py::arg("path"))
        .def_static("deserialize", &deserialize_from_file<BuildingBlockLibrary>, py::arg("path"))
        .def("__len__", &BuildingBlockLibrary::size)
        .def("__getitem__",
             py::overload_cast<BuildingBlockLibrary::Index>(&BuildingBlockLibrary::get, py::const_),
             py::return_value_policy::reference_internal);

    py::class_<BuildingBlockPreprocessor>(m, "BuildingBlockPreprocessor")
        .def(py::init<>())
        .def_readwrite("largest_fragment_only", &BuildingBlockPreprocessor::largest_fragment_only);

    py::class_<BuildingBlockCSVConfig>(m, "BuildingBlockCSVConfig")
        .def(py::init<>())
        .def_readwrite("identifier_column", &BuildingBlockCSVConfig::identifier_column)
        .def_readwrite("smiles_column", &BuildingBlockCSVConfig::smiles_column);

    py::register_exception<BuildingBlockLibraryError>(m, "BuildingBlockLibraryError",
                                                      PyExc_RuntimeError);

    m.def("bb_lib_from_sdf", &bb_lib_from_sdf, py::arg("path"),
          py::arg("preprocessor") = BuildingBlockPreprocessor{});
    m.def("bb_lib_from_csv", &bb_lib_from_csv, py::arg("path"),
          py::arg("config") = BuildingBlockCSVConfig{},
          py::arg("preprocessor") = BuildingBlockPreprocessor{});
}

static void def_rxn_lib(py::module &m) {
    py::class_<ReactionEntry>(m, "ReactionEntry")
        .def(py::init<>())
        .def_readwrite("reaction", &ReactionEntry::reaction)
        .def_readwrite("name", &ReactionEntry::name);

    py::class_<ReactionItem, ReactionEntry>(m, "ReactionItem")
        .def(py::init<>())
        .def_readonly("index", &ReactionItem::index);

    py::class_<ReactionLibrary::Match>(m, "ReactionMatch")
        .def_readonly("reaction_index", &ReactionLibrary::Match::reaction_index)
        .def_property_readonly(
            "reaction_name",
            [](const ReactionLibrary::Match &match) { return std::string(match.reaction_name); })
        .def_readonly("reactant_index", &ReactionLibrary::Match::reactant_index)
        .def_property_readonly(
            "reactant_name",
            [](const ReactionLibrary::Match &match) { return std::string(match.reactant_name); })
        .def_readonly("count", &ReactionLibrary::Match::count);

    py::class_<ReactionLibrary, py::smart_holder>(m, "ReactionLibrary")
        .def(py::init<>())
        .def("size", &ReactionLibrary::size)
        .def("get", py::overload_cast<ReactionLibrary::Index>(&ReactionLibrary::get, py::const_),
             py::arg("index"), py::return_value_policy::reference_internal)
        .def("get", py::overload_cast<const std::string &>(&ReactionLibrary::get, py::const_),
             py::arg("name"), py::return_value_policy::reference_internal)
        .def("add", &ReactionLibrary::add, py::arg("entry"))
        .def("match_reactants", &ReactionLibrary::match_reactants, py::arg("molecule"))
        .def("serialize", &serialize_to_file<ReactionLibrary>, py::arg("path"))
        .def_static("deserialize", &deserialize_from_file<ReactionLibrary>, py::arg("path"))
        .def("__len__", &ReactionLibrary::size)
        .def("__getitem__",
             py::overload_cast<ReactionLibrary::Index>(&ReactionLibrary::get, py::const_),
             py::return_value_policy::reference_internal);

    m.def("rxn_lib_from_plain_text", &rxn_lib_from_plain_text, py::arg("path"),
          py::arg("ignore_errors") = false);
}

static void def_int_lib(py::module &m) {
    py::class_<IntermediateEntry>(m, "IntermediateEntry")
        .def(py::init<>())
        .def_readwrite("postfix_notation", &IntermediateEntry::postfix_notation)
        .def_readwrite("molecule", &IntermediateEntry::molecule)
        .def_readwrite("identifier", &IntermediateEntry::identifier);

    py::class_<IntermediateItem, IntermediateEntry>(m, "IntermediateItem")
        .def(py::init<>())
        .def_readonly("index", &IntermediateItem::index);

    py::class_<IntermediateLibrary, py::smart_holder>(m, "IntermediateLibrary")
        .def(py::init<>())
        .def("size", &IntermediateLibrary::size)
        .def("get",
             py::overload_cast<IntermediateLibrary::Index>(&IntermediateLibrary::get, py::const_),
             py::arg("index"), py::return_value_policy::reference_internal)
        .def("get", py::overload_cast<const std::string &>(&IntermediateLibrary::get, py::const_),
             py::arg("identifier"), py::return_value_policy::reference_internal)
        .def("add", &IntermediateLibrary::add, py::arg("entry"))
        .def("clear", &IntermediateLibrary::clear)
        .def("serialize", &serialize_to_file<IntermediateLibrary>, py::arg("path"))
        .def_static("deserialize", &deserialize_from_file<IntermediateLibrary>, py::arg("path"))
        .def("__len__", &IntermediateLibrary::size)
        .def("__getitem__",
             py::overload_cast<IntermediateLibrary::Index>(&IntermediateLibrary::get, py::const_),
             py::return_value_policy::reference_internal);
}

static void def_chemical_space(py::module &m) {
    py::class_<ReactantMatchingConfig>(m, "ReactantMatchingConfig")
        .def(py::init<>())
        .def_readwrite("selectivity_cutoff", &ReactantMatchingConfig::selectivity_cutoff)
        .def("check", &ReactantMatchingConfig::check, py::arg("match"));

    py::class_<ReactantLists, py::smart_holder>(m, "ReactantLists")
        .def(py::init<>())
        .def("get", &ReactantLists::get, py::arg("reaction_index"), py::arg("reactant_index"),
             py::return_value_policy::reference_internal)
        .def("num_matches", &ReactantLists::num_matches);

    py::class_<ChemicalSpace::PeekStats>(m, "ChemicalSpacePeekStats")
        .def(py::init<>())
        .def_readonly("num_reactions", &ChemicalSpace::PeekStats::num_reactions)
        .def_readonly("num_building_blocks", &ChemicalSpace::PeekStats::num_building_blocks)
        .def_readonly("num_intermediates", &ChemicalSpace::PeekStats::num_intermediates);

    py::class_<ChemicalSpace, py::smart_holder>(m, "ChemicalSpace")
        .def(py::init<std::unique_ptr<BuildingBlockLibrary>, std::unique_ptr<ReactionLibrary>,
                      std::unique_ptr<IntermediateLibrary>, const ReactantMatchingConfig &>(),
             py::arg("bb_lib"), py::arg("rxn_lib"), py::arg("int_lib"),
             py::arg("matching_config") = ReactantMatchingConfig{})
        .def("serialize", &serialize_to_file<ChemicalSpace>, py::arg("path"))
        .def_static("deserialize", &deserialize_from_file<ChemicalSpace>, py::arg("path"))
        .def_static("peek",
                    [](const py::bytes &data) {
                        std::string raw(data);
                        std::stringstream ss(raw);
                        return ChemicalSpace::peek(ss);
                    })
        .def_static("peek",
                    [](const std::filesystem::path &path) {
                        std::ifstream ifs(path, std::ios::binary);
                        if (!ifs) {
                            throw std::runtime_error("failed to open file for reading: " +
                                                     path.string());
                        }
                        return ChemicalSpace::peek(ifs);
                    })
        .def("bb_lib", &ChemicalSpace::bb_lib, py::return_value_policy::reference_internal)
        .def("rxn_lib", &ChemicalSpace::rxn_lib, py::return_value_policy::reference_internal)
        .def("int_lib", &ChemicalSpace::int_lib, py::return_value_policy::reference_internal)
        .def("reactant_matching_config", &ChemicalSpace::reactant_matching_config,
             py::return_value_policy::reference_internal)
        .def("building_block_reactant_lists", &ChemicalSpace::building_block_reactant_lists,
             py::return_value_policy::reference_internal)
        .def("intermediate_reactant_lists", &ChemicalSpace::intermediate_reactant_lists,
             py::return_value_policy::reference_internal)
        .def("generate_intermediates", &ChemicalSpace::generate_intermediates)
        .def("build_reactant_lists_for_building_blocks",
             &ChemicalSpace::build_reactant_lists_for_building_blocks)
        .def("build_reactant_lists_for_intermediates",
             &ChemicalSpace::build_reactant_lists_for_intermediates)
        .def("print_reactant_lists",
             [](const ChemicalSpace &chemspace) {
                 std::ostringstream oss;
                 chemspace.print_reactant_lists(oss);
                 return oss.str();
             })
        .def("new_synthesis", &ChemicalSpace::new_synthesis);
}

static void def_chemical_space_synthesis(py::module &m) {
    py::class_<ChemicalSpaceSynthesis::Result>(m, "SynthesisResult")
        .def(py::init<>())
        .def_readonly("is_ok", &ChemicalSpaceSynthesis::Result::is_ok)
        .def_readonly("message", &ChemicalSpaceSynthesis::Result::message)
        .def_static("ok", &ChemicalSpaceSynthesis::Result::ok)
        .def_static("error", &ChemicalSpaceSynthesis::Result::error, py::arg("message"))
        .def("__bool__", [](const ChemicalSpaceSynthesis::Result &result) {
            return static_cast<bool>(result);
        });

    py::class_<ChemicalSpaceSynthesis, py::smart_holder>(m, "Synthesis")
        .def("postfix_notation", &ChemicalSpaceSynthesis::postfix_notation,
             py::return_value_policy::reference_internal)
        .def("synthesis", &ChemicalSpaceSynthesis::synthesis,
             py::return_value_policy::reference_internal)
        .def("count_building_blocks", &ChemicalSpaceSynthesis::count_building_blocks)
        .def("count_reactions", &ChemicalSpaceSynthesis::count_reactions)
        .def("products", &ChemicalSpaceSynthesis::products)
        .def("add_building_block",
             py::overload_cast<BuildingBlockLibrary::Index>(
                 &ChemicalSpaceSynthesis::add_building_block),
             py::arg("index"))
        .def("add_building_block",
             py::overload_cast<const std::string &>(&ChemicalSpaceSynthesis::add_building_block),
             py::arg("identifier"))
        .def("add_reaction",
             py::overload_cast<ReactionLibrary::Index, std::optional<size_t>>(
                 &ChemicalSpaceSynthesis::add_reaction),
             py::arg("index"), py::arg("max_outcomes"))
        .def("add_reaction",
             py::overload_cast<const std::string &, std::optional<size_t>>(
                 &ChemicalSpaceSynthesis::add_reaction),
             py::arg("name"), py::arg("max_outcomes"))
        .def("add_postfix_notation", &ChemicalSpaceSynthesis::add_postfix_notation,
             py::arg("postfix_notation"), py::arg("max_outcomes"))
        .def("add_intermediate",
             py::overload_cast<IntermediateLibrary::Index, std::optional<size_t>>(
                 &ChemicalSpaceSynthesis::add_intermediate),
             py::arg("index"), py::arg("max_outcomes"))
        .def("add_intermediate",
             py::overload_cast<const std::string &, std::optional<size_t>>(
                 &ChemicalSpaceSynthesis::add_intermediate),
             py::arg("identifier"), py::arg("max_outcomes"))
        .def("undo", &ChemicalSpaceSynthesis::undo);
}

void def_module_chemspace(pybind11::module &m) {
    def_postfix_notation(m);
    def_bb_lib(m);
    def_rxn_lib(m);
    def_int_lib(m);
    def_chemical_space(m);
    def_chemical_space_synthesis(m);
}
