#include "pharmacophore.hpp"

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost/python/numpy/ndarray.hpp>
#include <boost/python/object_fwd.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <omp.h>

#include "utils/pickling.hpp"

namespace py = boost::python;
namespace np = boost::python::numpy;
using namespace synthesis_backend;

py::list get_nodes(const PharmacophoreGraph &graph) {
    py::list nodes;
    for (const auto &node : graph.nodes) {
        nodes.append(node);
    }
    return nodes;
}

py::list get_edges(const PharmacophoreGraph &graph) {
    py::list edges;
    for (const auto &edge : graph.edge_weights) {
        edges.append(
            py::make_tuple(edge.first.first, edge.first.second, edge.second));
    }
    return edges;
}

PharmacophoreGraph
create_pharmacophore_graph_bind(const Mol_sptr &mol, const py::object &features,
                                const BondWeights &bond_weights,
                                float default_bond_weight) {
    if (py::extract<RDKit::MolChemicalFeatureFactory>(features).check()) {
        RDKit::MolChemicalFeatureFactory feature_factory =
            py::extract<RDKit::MolChemicalFeatureFactory>(features);
        return create_pharmacophore_graph(mol, feature_factory, bond_weights,
                                          default_bond_weight);
    }

    RDKit::FeatSPtrList feature_list;
    if (py::extract<py::list>(features).check()) {
        py::list feature_list_py = py::extract<py::list>(features);
        for (int i = 0; i < py::len(feature_list_py); ++i) {
            if (!py::extract<RDKit::FeatSPtr>(feature_list_py[i]).check()) {
                throw std::runtime_error("Invalid feature type in list.");
            }
            feature_list.push_back(
                py::extract<RDKit::FeatSPtr>(feature_list_py[i]));
        }
    } else if (py::extract<py::tuple>(features).check()) {
        py::tuple feature_tuple = py::extract<py::tuple>(features);
        for (int i = 0; i < py::len(feature_tuple); ++i) {
            if (!py::extract<RDKit::FeatSPtr>(feature_tuple[i]).check()) {
                throw std::runtime_error("Invalid feature type in tuple.");
            }
            feature_list.push_back(
                py::extract<RDKit::FeatSPtr>(feature_tuple[i]));
        }
    } else if (py::extract<RDKit::FeatSPtrList>(features).check()) {
        feature_list = py::extract<RDKit::FeatSPtrList>(features);
    } else {
        throw std::runtime_error("Invalid features type provided.");
    }

    return create_pharmacophore_graph(mol, feature_list, bond_weights,
                                      default_bond_weight);
}

py::tuple subgraph_similarity_bind(const PharmacophoreGraph &g1,
                                   const PharmacophoreGraph &g2,
                                   const py::list &subgraph_nodes_1,
                                   const py::list &subgraph_nodes_2,
                                   float lambda) {

    std::vector<size_t> nodes_1, nodes_2;
    for (int i = 0; i < py::len(subgraph_nodes_1); ++i) {
        nodes_1.push_back(py::extract<size_t>(subgraph_nodes_1[i]));
    }
    for (int i = 0; i < py::len(subgraph_nodes_2); ++i) {
        nodes_2.push_back(py::extract<size_t>(subgraph_nodes_2[i]));
    }

    auto result = subgraph_similarity(g1, g2, nodes_1, nodes_2, lambda);
    return py::make_tuple(result.first, result.second);
}

py::tuple pharmacophore_similarity_bind(const PharmacophoreGraph &g1,
                                        const PharmacophoreGraph &g2,
                                        float lambda) {
    auto result = pharmacophore_similarity(g1, g2, lambda);
    return py::make_tuple(result.first, result.second);
}

py::tuple pairwise_pharmacophore_similarity(const py::list &graphs1,
                                            const py::list &graphs2,
                                            float lambda) {
    auto n1 = py::len(graphs1), n2 = py::len(graphs2);
    np::ndarray result1 =
        np::zeros(py::make_tuple(n1, n2), np::dtype::get_builtin<float>());
    np::ndarray result2 =
        np::zeros(py::make_tuple(n1, n2), np::dtype::get_builtin<float>());

    std::vector<PharmacophoreGraph> graphs1_vec, graphs2_vec;
    for (int i = 0; i < n1; ++i) {
        if (!py::extract<PharmacophoreGraph>(graphs1[i]).check()) {
            throw std::runtime_error("Invalid graph type in graphs1 at index " +
                                     std::to_string(i));
        }
        graphs1_vec.push_back(py::extract<PharmacophoreGraph>(graphs1[i])());
    }
    for (int j = 0; j < n2; ++j) {
        if (!py::extract<PharmacophoreGraph>(graphs2[j]).check()) {
            throw std::runtime_error("Invalid graph type in graphs2 at index " +
                                     std::to_string(j));
        }
        graphs2_vec.push_back(py::extract<PharmacophoreGraph>(graphs2[j])());
    }

#pragma omp parallel for
    for (auto ij = 0; ij < n1 * n2; ++ij) {
        auto i = ij / n2;
        auto j = ij % n2;

        auto sim_pair =
            pharmacophore_similarity(graphs1_vec[i], graphs2_vec[j], lambda);
        auto ptr1 = reinterpret_cast<float *>(result1.get_data()) + ij;
        auto ptr2 = reinterpret_cast<float *>(result2.get_data()) + ij;
        *ptr1 = sim_pair.first;
        *ptr2 = sim_pair.second;
    }

    return py::make_tuple(result1, result2);
}

BOOST_PYTHON_MODULE(pharmacophore) {
    py::class_<PharmacophoreNode>(
        "PharmacophoreNode",
        py::init<const std::string &, const std::string &>())
        .def(py::init<const RDKit::MolChemicalFeature &>())
        .def_readonly("family", &PharmacophoreNode::family)
        .def_readonly("type", &PharmacophoreNode::type);

    py::class_<PharmacophoreGraph>("PharmacophoreGraph", py::init<>())
        .def("set_edge", &PharmacophoreGraph::set_edge,
             (py::arg("u"), py::arg("v"), py::arg("weight") = 1.0f))
        .def("unset_edge", &PharmacophoreGraph::unset_edge,
             (py::arg("u"), py::arg("v")))
        .def("has_edge", &PharmacophoreGraph::has_edge,
             (py::arg("u"), py::arg("v")))
        .def("get_nodes", &get_nodes,
             py::return_value_policy<py::return_by_value>())
        .def("get_edges", &get_edges,
             py::return_value_policy<py::return_by_value>())
        .def("get_min_spanning_tree",
             &PharmacophoreGraph::get_min_spanning_tree,
             py::return_value_policy<py::return_by_value>())
        .def("get_subgraph_with_n_edges",
             &PharmacophoreGraph::get_subgraph_with_n_edges,
             (py::arg("n_edges")),
             py::return_value_policy<py::return_by_value>())
        .def_pickle(pickle_suite<PharmacophoreGraph>());

    py::class_<BondWeights>("BondWeights")
        .def(py::map_indexing_suite<BondWeights>());

    py::scope().attr("bond_weights_EMPTY") = bond_weights_EMPTY;
    py::scope().attr("bond_weights_RELATIVE_BOND_LENGTH") =
        bond_weights_RELATIVE_BOND_LENGTH;

    py::def("create_pharmacophore_graph", &create_pharmacophore_graph_bind,
            (py::arg("mol"), py::arg("features"),
             py::arg("bond_weights") = bond_weights_EMPTY,
             py::arg("default_bond_weight") = 1.0f),
            py::return_value_policy<py::return_by_value>());

    py::def("subgraph_similarity", &subgraph_similarity_bind,
            (py::arg("g1"), py::arg("g2"), py::arg("subgraph_nodes_1"),
             py::arg("subgraph_nodes_2"), py::arg("lambda") = 1.0f),
            py::return_value_policy<py::return_by_value>());

    py::def("pharmacophore_similarity", &pharmacophore_similarity_bind,
            (py::arg("g1"), py::arg("g2"), py::arg("lambda") = 1.0f),
            py::return_value_policy<py::return_by_value>());

    py::def("pairwise_pharmacophore_similarity",
            &pairwise_pharmacophore_similarity,
            (py::arg("graphs1"), py::arg("graphs2"), py::arg("lambda") = 1.0f),
            py::return_value_policy<py::return_by_value>());
}
