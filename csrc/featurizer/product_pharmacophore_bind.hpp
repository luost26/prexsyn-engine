#pragma once

#include "product_pharmacophore.hpp"

#include <boost/python.hpp>
#include <boost/python/list.hpp>
#include <boost/python/object_fwd.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/with_custodian_and_ward.hpp>
#include <filesystem>
#include <omp.h>

namespace py = boost::python;
using namespace prexsyn_engine;

static void product_pharmacophore_featurizer_call_bind(
    ProductPharmacophoreFeaturizer &featurizer, py::object obj,
    FeatureBuilder &builder) {
    if (py::extract<Synthesis>(obj).check()) {
        featurizer(py::extract<Synthesis>(obj)(), builder);
    } else if (py::extract<PharmacophoreGraph>(obj).check()) {
        featurizer(py::extract<PharmacophoreGraph>(obj)(), builder);
    } else {
        throw std::runtime_error("Unsupported object type for featurization.");
    }
}

static PharmacophoreGraph
get_graph_bind(ProductPharmacophoreFeaturizer &featurizer, const Mol_sptr &mol,
               py::object features_obj) {
    if (!features_obj.is_none()) {
        if (py::extract<py::list>(features_obj).check()) {
            std::vector<RDKit::FeatSPtr> features;
            auto features_list = py::extract<py::list>(features_obj)();
            for (long i = 0; i < py::len(features_list); ++i) {
                features.push_back(
                    py::extract<RDKit::FeatSPtr>(features_list[i])());
            }
            return featurizer.get_graph(mol, features);
        } else if (py::extract<std::vector<RDKit::FeatSPtr>>(features_obj)
                       .check()) {
            auto features =
                py::extract<std::vector<RDKit::FeatSPtr>>(features_obj)();
            return featurizer.get_graph(mol, features);
        } else {
            throw std::runtime_error("Invalid features type provided.");
        }
    } else {
        return featurizer.get_graph(mol);
    }
}

static py::list get_graphs_bind(ProductPharmacophoreFeaturizer &featurizer,
                                const py::list &mols) {
    auto n = py::len(mols);
    std::vector<Mol_sptr> mols_vec;
    for (ssize_t i = 0; i < n; ++i) {
        mols_vec.push_back(py::extract<Mol_sptr>(mols[i])());
    }

    std::vector<PharmacophoreGraph> graphs(n);
#pragma omp parallel for
    for (ssize_t i = 0; i < n; ++i) {
        graphs[i] = featurizer.get_graph(mols_vec[i]);
    }

    py::list result;
    for (const auto &graph : graphs) {
        result.append(graph);
    }
    return result;
}

inline void product_pharmacophore_bind() {
    py::class_<ProductPharmacophoreFeaturizerOption>(
        "ProductPharmacophoreFeaturizerOption", py::init<>())
        .def_readwrite("name", &ProductPharmacophoreFeaturizerOption::name)
        .def_readwrite("feature_def",
                       &ProductPharmacophoreFeaturizerOption::feature_def)
        .def_readwrite("type_index_offset",
                       &ProductPharmacophoreFeaturizerOption::type_index_offset)
        .def_readwrite(
            "random_num_nodes_max",
            &ProductPharmacophoreFeaturizerOption::random_num_nodes_max)
        .def_readwrite(
            "random_num_nodes_min",
            &ProductPharmacophoreFeaturizerOption::random_num_nodes_min)
        .def_readwrite(
            "random_num_edges_max",
            &ProductPharmacophoreFeaturizerOption::random_num_edges_max)
        .def_readwrite("bond_weights",
                       &ProductPharmacophoreFeaturizerOption::bond_weights)
        .def_readwrite(
            "default_bond_weight",
            &ProductPharmacophoreFeaturizerOption::default_bond_weight);

    py::class_<std::vector<RDKit::FeatSPtr>>("FeatSPtrVector")
        .def(py::vector_indexing_suite<std::vector<RDKit::FeatSPtr>, true>());

    py::class_<ProductPharmacophoreFeaturizer>(
        "ProductPharmacophoreFeaturizer",
        py::init<const ProductPharmacophoreFeaturizerOption &>(
            (py::arg("option") = ProductPharmacophoreFeaturizerOption())))
        .def("__call__", &product_pharmacophore_featurizer_call_bind,
             (py::arg("obj"), py::arg("builder")))
        .def("get_features", &ProductPharmacophoreFeaturizer::get_features,
             (py::arg("mol")),
             py::return_value_policy<
                 py::return_by_value,
                 py::with_custodian_and_ward_postcall<0, 1>>())
        .def("get_graph", &get_graph_bind,
             (py::arg("mol"), py::arg("features") = py::object()),
             py::return_value_policy<py::return_by_value>())
        .def("get_graphs", &get_graphs_bind, (py::arg("mols")),
             py::return_value_policy<py::return_by_value>());

    py::register_ptr_to_python<
        std::shared_ptr<ProductPharmacophoreFeaturizer>>();
}
