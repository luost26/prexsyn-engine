#include "synthesis.hpp"

#include <boost/python.hpp>

namespace py = boost::python;
using namespace prexsyn_engine;

BOOST_PYTHON_MODULE(synthesis) {
    py::class_<PostfixNotationTokenDef>(
        "PostfixNotationTokenDef",
        py::init<int, int, int, int, int>(
            (py::arg("pad") = DEFAULT_PAD, py::arg("end") = DEFAULT_END,
             py::arg("start") = DEFAULT_START, py::arg("bb") = DEFAULT_BB,
             py::arg("rxn") = DEFAULT_RXN)))
        .def_readonly("PAD", &PostfixNotationTokenDef::PAD)
        .def_readonly("END", &PostfixNotationTokenDef::END)
        .def_readonly("START", &PostfixNotationTokenDef::START)
        .def_readonly("BB", &PostfixNotationTokenDef::BB)
        .def_readonly("RXN", &PostfixNotationTokenDef::RXN);

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

    py::implicitly_convertible<std::shared_ptr<PostfixNotationFeaturizer>,
                               std::shared_ptr<Featurizer>>();
}
