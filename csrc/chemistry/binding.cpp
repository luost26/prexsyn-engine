#include "binding.hpp"

#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>

#include "chemistry.hpp"

namespace py = pybind11;
using namespace prexsyn;

void def_module_chemistry(pybind11::module &m) {
    py::class_<Molecule, py::smart_holder>(m, "Molecule")
        .def_static("from_smiles", &Molecule::from_smiles)
        .def("smiles", &Molecule::smiles)
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
