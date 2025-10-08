
#include <boost/python/class.hpp>
#include <boost/python/import.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <memory>

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

#include "pydict_builder.hpp"

#include "product_pharmacophore_bind.hpp"
#include "product_property_bind.hpp"
#include "product_structure_bind.hpp"
#include "synthesis_bind.hpp"

namespace py = boost::python;
namespace np = boost::python::numpy;

using namespace prexsyn_engine;
BOOST_PYTHON_MODULE(featurizer) {
    np::initialize();

    py::class_<PyDictBuilder, std::shared_ptr<PyDictBuilder>>("PyDictBuilder",
                                                              py::init<>())
        .def("get", &PyDictBuilder::get)
        .def("erase_type", &PyDictBuilder::erase_type,
             py::return_internal_reference<>());

    std::ignore =
        py::class_<Featurizer, std::shared_ptr<Featurizer>, boost::noncopyable>(
            "Featurizer", py::no_init);

    auto featurizer_set_class =
        py::class_<FeaturizerSet, std::shared_ptr<FeaturizerSet>>(
            "FeaturizerSet", py::init<>())
            .def("__call__", &FeaturizerSet::operator(),
                 (py::arg("synthesis"), py::arg("builder")));
    py::register_ptr_to_python<std::shared_ptr<FeaturizerSet>>();
    auto def_featurizer_set_add_method = [&]<typename F>() {
        featurizer_set_class.def(
            "add",
            (FeaturizerSet & (FeaturizerSet::*)(std::shared_ptr<F>)) &
                FeaturizerSet::add,
            py::return_internal_reference<>(), (py::arg("featurizer")));
    };

    product_structure_bind();
    def_featurizer_set_add_method
        .template operator()<ProductStructureFeaturizer>();

    product_property_bind();
    def_featurizer_set_add_method
        .template operator()<ProductRDKitPropertyFeaturizer>();

    synthesis_bind();
    def_featurizer_set_add_method
        .template operator()<PostfixNotationFeaturizer>();

    product_pharmacophore_bind();
    def_featurizer_set_add_method
        .template operator()<ProductPharmacophoreFeaturizer>();
}
