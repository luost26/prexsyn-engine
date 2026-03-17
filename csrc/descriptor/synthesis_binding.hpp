#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>

#include "base.hpp"
#include "synthesis.hpp"

namespace prexsyn::descriptor {

namespace py = pybind11;

inline void def_submodule_synthesis(py::module_ &m) {
    py::class_<TokenDef, py::smart_holder>(m, "TokenDef")
        .def(py::init<>())
        .def_static("from_dict",
                    [](const py::dict &d) {
                        TokenDef token_def;
                        if (d.contains("pad"))
                            token_def.pad = d["pad"].cast<std::int64_t>();
                        if (d.contains("end"))
                            token_def.end = d["end"].cast<std::int64_t>();
                        if (d.contains("start"))
                            token_def.start = d["start"].cast<std::int64_t>();
                        if (d.contains("bb"))
                            token_def.bb = d["bb"].cast<std::int64_t>();
                        if (d.contains("rxn"))
                            token_def.rxn = d["rxn"].cast<std::int64_t>();
                        return token_def;
                    })
        .def_readwrite("pad", &TokenDef::pad)
        .def_readwrite("end", &TokenDef::end)
        .def_readwrite("start", &TokenDef::start)
        .def_readwrite("bb", &TokenDef::bb)
        .def_readwrite("rxn", &TokenDef::rxn);

    py::class_<SynthesisPostfixNotation, SynthesisDescriptor, py::smart_holder>(
        m, "SynthesisPostfixNotation")
        .def_static("create",
                    py::overload_cast<const TokenDef &, size_t>(&SynthesisPostfixNotation::create),
                    py::arg("token_def"), py::arg("max_length") = 16)
        .def_static("create", py::overload_cast<size_t>(&SynthesisPostfixNotation::create),
                    py::arg("max_length") = 16);
}

} // namespace prexsyn::descriptor
