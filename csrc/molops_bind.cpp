#include "molops.hpp"

#include <boost/python.hpp>

namespace py = boost::python;
using namespace synthesis_backend;

py::object murcko_scaffold_bind(const Mol_sptr &mol) {
    auto o = murcko_scaffold(mol);
    if (o) {
        return py::object(o.value());
    } else {
        return py::object(); // None
    }
}

py::list brics_fragments_bind(const Mol_sptr &mol) {
    auto frags = brics_fragments(mol);
    py::list out;
    for (const auto &frag : frags) {
        out.append(frag);
    }
    return out;
}

BOOST_PYTHON_MODULE(molops) {
    py::def("murcko_scaffold", &murcko_scaffold_bind, py::arg("mol"));
    py::def("brics_fragments", &brics_fragments_bind, py::arg("mol"));
}
