#include "indexer.hpp"

#include <boost/python.hpp>

namespace py = boost::python;
using namespace synthesis_backend;

BOOST_PYTHON_MODULE(indexer) {
    py::def("get_suitable_reactant_indices_for_mol",
            (std::vector<size_t> (*)(const Reaction_sptr &, const Mol_sptr &))(
                &get_suitable_reactant_indices),
            (py::arg("reaction"), py::arg("mol")));

    py::def("get_suitable_reactant_indices_for_synthesis",
            (std::vector<size_t> (*)(const Reaction_sptr &,
                                     const Synthesis_sptr &))(
                &get_suitable_reactant_indices),
            (py::arg("reaction"), py::arg("synthesis")));
}
