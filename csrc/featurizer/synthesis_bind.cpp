#pragma once

#include "synthesis.hpp"

#include <boost/python.hpp>

namespace py = boost::python;
using namespace prexsyn_engine;

BOOST_PYTHON_MODULE(featurizer__synthesis) {
    py::class_<PostfixNotationTokenDef>("PostfixNotationTokenDef")
        .def_readwrite("PAD", &PostfixNotationTokenDef::PAD)
        .def_readwrite("END", &PostfixNotationTokenDef::END)
        .def_readwrite("START", &PostfixNotationTokenDef::START)
        .def_readwrite("BB", &PostfixNotationTokenDef::BB)
        .def_readwrite("RXN", &PostfixNotationTokenDef::RXN);

    py::class_<PostfixNotationFeaturizer,
               std::shared_ptr<PostfixNotationFeaturizer>,
               py::bases<Featurizer>>(
        "PostfixNotationFeaturizer",
        py::init<unsigned int, PostfixNotationTokenDef>(
            (py::arg("max_length") = 16,
             py::arg("token_def") = PostfixNotationTokenDef())))
        .def_readonly("max_length", &PostfixNotationFeaturizer::max_length)
        .def_readonly("token_def", &PostfixNotationFeaturizer::token_def)
        .def("__call__", &PostfixNotationFeaturizer::operator(),
             (py::arg("synthesis"), py::arg("builder")));
}
