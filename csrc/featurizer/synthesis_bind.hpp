#pragma once

#include "synthesis.hpp"

#include <boost/python.hpp>

namespace py = boost::python;
using namespace synthesis_backend;

inline void synthesis_bind() {
    py::class_<PostfixNotationTokenDef>("PostfixNotationTokenDef")
        .def_readwrite("PAD", &PostfixNotationTokenDef::PAD)
        .def_readwrite("END", &PostfixNotationTokenDef::END)
        .def_readwrite("START", &PostfixNotationTokenDef::START)
        .def_readwrite("BB", &PostfixNotationTokenDef::BB)
        .def_readwrite("RXN", &PostfixNotationTokenDef::RXN);

    py::class_<PostfixNotationFeaturizerOption>(
        "PostfixNotationFeaturizerOption")
        .def_readwrite("length", &PostfixNotationFeaturizerOption::length)
        .def_readwrite("token_def",
                       &PostfixNotationFeaturizerOption::token_def);

    py::class_<PostfixNotationFeaturizer>(
        "PostfixNotationFeaturizer",
        py::init<PostfixNotationFeaturizerOption>(
            (py::arg("option") = PostfixNotationFeaturizerOption())))
        .def_readonly("option", &PostfixNotationFeaturizer::option)
        .def("__call__", &PostfixNotationFeaturizer::operator(),
             (py::arg("synthesis"), py::arg("builder")));

    py::register_ptr_to_python<std::shared_ptr<PostfixNotationFeaturizer>>();
}
