#include "pharmacophore.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

#include <GraphMol/MolOps.h>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include "types.hpp"
#include "utils/assert.hpp"

namespace prexsyn_engine {

PharmacophoreNode::PharmacophoreNode(const std::string &family,
                                     const std::string &type)
    : family(family), type(type) {}

PharmacophoreNode::PharmacophoreNode(const RDKit::MolChemicalFeature &feature)
    : family(feature.getFamily()), type(feature.getType()) {}

PharmacophoreGraph::PharmacophoreGraph() {}

PharmacophoreGraph::PharmacophoreGraph(
    const std::vector<PharmacophoreNode> &nodes)
    : nodes(nodes), adj_list(nodes.size()) {}

PharmacophoreGraph::PharmacophoreGraph(const RDKit::FeatSPtrList &features) {
    for (const auto &feature : features) {
        Ensures(feature != nullptr);
        nodes.emplace_back(*feature);
    }
    adj_list.resize(nodes.size());
}

PharmacophoreGraph::PharmacophoreGraph(
    const std::vector<RDKit::FeatSPtr> &features)
    : adj_list(features.size()) {
    for (const auto &feature : features) {
        Ensures(feature != nullptr);
        nodes.emplace_back(*feature);
    }
}

std::string PharmacophoreGraph::pickle() const {
    std::ostringstream buffer;
    boost::archive::binary_oarchive oa(buffer);

    // nodes
    oa << nodes.size();
    for (const auto &node : nodes) {
        oa << node.family << node.type;
    }

    // adjacency list
    oa << adj_list.size();
    for (const auto &neighbors : adj_list) {
        oa << neighbors.size();
        for (const auto &neighbor : neighbors) {
            oa << neighbor;
        }
    }

    // edge weights
    oa << edge_weights.size();
    for (const auto &edge : edge_weights) {
        oa << edge.first.first << edge.first.second << edge.second;
    }

    return buffer.str();
}

PharmacophoreGraph *PharmacophoreGraph::unpickle(const std::string &data) {
    std::istringstream buffer(data);
    boost::archive::binary_iarchive ia(buffer);
    PharmacophoreGraph *graph = new PharmacophoreGraph();

    // nodes
    size_t num_nodes;
    ia >> num_nodes;
    for (size_t i = 0; i < num_nodes; ++i) {
        std::string family, type;
        ia >> family >> type;
        graph->nodes.emplace_back(family, type);
    }

    // adjacency list
    size_t num_adj_lists;
    ia >> num_adj_lists;
    graph->adj_list.resize(num_adj_lists);
    for (size_t i = 0; i < num_adj_lists; ++i) {
        size_t num_neighbors;
        ia >> num_neighbors;
        for (size_t j = 0; j < num_neighbors; ++j) {
            size_t neighbor;
            ia >> neighbor;
            graph->adj_list[i].insert(neighbor);
        }
    }

    // edge weights
    size_t num_edges;
    ia >> num_edges;
    for (size_t i = 0; i < num_edges; ++i) {
        size_t u, v;
        float weight;
        ia >> u >> v >> weight;
        graph->edge_weights[{u, v}] = weight;
    }

    return graph;
}

void PharmacophoreGraph::set_edge(size_t u, size_t v, float weight) {
    Ensures(u < nodes.size() && v < nodes.size());
    edge_weights[{u, v}] = weight;
    edge_weights[{v, u}] = weight;
    if (adj_list.size() <= u) {
        adj_list.resize(u + 1);
    }
    if (adj_list.size() <= v) {
        adj_list.resize(v + 1);
    }
    adj_list[u].insert(v);
    adj_list[v].insert(u);
}

void PharmacophoreGraph::unset_edge(size_t u, size_t v) {
    Ensures(u < nodes.size() && v < nodes.size());
    edge_weights.erase({u, v});
    edge_weights.erase({v, u});
    if (u < adj_list.size()) {
        adj_list[u].erase(v);
    }
    if (v < adj_list.size()) {
        adj_list[v].erase(u);
    }
}

bool PharmacophoreGraph::has_edge(size_t u, size_t v) const {
    if (u >= adj_list.size()) {
        return false;
    }
    const auto &neighbors = adj_list[u];
    return neighbors.find(v) != neighbors.end();
}

PharmacophoreGraph PharmacophoreGraph::get_min_spanning_tree() const {
    PharmacophoreGraph mst(nodes);

    std::set<std::pair<float, std::pair<size_t, size_t>>> edges;
    for (const auto &edge : edge_weights) {
        edges.insert({edge.second, edge.first});
    }

    std::vector<size_t> disjoint_set(nodes.size(), -1);
    for (size_t i = 0; i < nodes.size(); ++i) {
        disjoint_set[i] = i;
    }

    for (const auto &edge : edges) {
        float weight = edge.first;
        size_t u = edge.second.first;
        size_t v = edge.second.second;

        size_t root_u = u;
        while (disjoint_set[root_u] != root_u) {
            root_u = disjoint_set[root_u];
        }

        size_t root_v = v;
        while (disjoint_set[root_v] != root_v) {
            root_v = disjoint_set[root_v];
        }

        if (root_u != root_v) {
            mst.set_edge(u, v, weight);
            disjoint_set[root_v] = root_u; // Union operation
        }
    }
    return mst;
}

PharmacophoreGraph
PharmacophoreGraph::get_subgraph_with_n_edges(size_t n_edges) const {

    PharmacophoreGraph subgraph(nodes);

    std::set<std::tuple<float, size_t, size_t>> mst_edges;
    PharmacophoreGraph mst = get_min_spanning_tree();
    for (const auto &edge : mst.edge_weights) {
        if (edge.first.first < edge.first.second) { // Dont include self-loops
            mst_edges.insert(
                {edge.second, edge.first.first, edge.first.second});
        }
    }

    std::set<std::tuple<float, size_t, size_t>> other_edges;
    for (const auto &edge : edge_weights) {
        if (edge.first.first < edge.first.second &&
            !mst.has_edge(edge.first.first, edge.first.second)) {
            other_edges.insert(
                {edge.second, edge.first.first, edge.first.second});
        }
    }

    std::vector<std::tuple<float, size_t, size_t>> sorted_edges(
        mst_edges.begin(), mst_edges.end());
    sorted_edges.insert(sorted_edges.end(), other_edges.begin(),
                        other_edges.end());

    for (size_t i = 0; i < n_edges && i < sorted_edges.size(); ++i) {
        const auto &edge = sorted_edges[i];
        size_t u = std::get<1>(edge);
        size_t v = std::get<2>(edge);
        float weight = std::get<0>(edge);
        subgraph.set_edge(u, v, weight);
    }

    return subgraph;
}

bool PharmacophoreGraph::is_subgraph_connected(
    const std::vector<size_t> &subgraph) const {
    if (subgraph.empty()) {
        return false;
    }

    std::set<size_t> visited;
    std::vector<size_t> stack;
    stack.push_back(subgraph[0]);
    visited.insert(subgraph[0]);
    std::set<size_t> subgraph_set(subgraph.begin(), subgraph.end());

    while (!stack.empty()) {
        size_t node = stack.back();
        stack.pop_back();

        for (const auto &neighbor : adj_list[node]) {
            if (subgraph_set.find(neighbor) != subgraph_set.end() &&
                visited.find(neighbor) == visited.end()) {
                visited.insert(neighbor);
                stack.push_back(neighbor);
            }
        }
    }

    return visited.size() == subgraph.size();
}

float average_atoms_to_atoms_distance(
    const Mol_sptr &mol,
    const RDKit::MolChemicalFeature::AtomPtrContainer &atoms_1,
    const RDKit::MolChemicalFeature::AtomPtrContainer &atoms_2,
    const std::map<RDKit::Bond::BondType, float> &bond_weights,
    float default_bond_weight) {
    Ensures(mol != nullptr);
    float total_distance = 0.0f;
    int count = 0;
    for (const auto &atom_1 : atoms_1) {
        for (const auto &atom_2 : atoms_2) {
            Ensures(atom_1 != nullptr && atom_2 != nullptr);
            if (atom_1->getIdx() == atom_2->getIdx()) {
                count++;
                continue;
            }
            auto path = RDKit::MolOps::getShortestPath(*mol, atom_1->getIdx(),
                                                       atom_2->getIdx());

            if (path.size() < 2) {
                continue;
            }

            auto it1 = path.begin();
            auto it2 = std::next(path.begin());
            float distance = 0.0f;

            while (it2 != path.end()) {
                auto bond = mol->getBondBetweenAtoms(*it1, *it2);
                if (bond) {
                    auto bond_type = bond->getBondType();
                    auto weight_it = bond_weights.find(bond_type);
                    if (weight_it != bond_weights.end()) {
                        distance += weight_it->second;
                    } else {
                        distance += default_bond_weight;
                    }
                }
                it1 = it2;
                it2 = std::next(it2);
            }

            total_distance += distance;
            count++;
        }
    }

    return count > 0 ? total_distance / count : -1.0f;
}

std::shared_ptr<RDKit::MolChemicalFeatureFactory>
get_preset_feature_factory(const std::string &preset_name) {
    std::shared_ptr<RDKit::MolChemicalFeatureFactory> factory;
    factory.reset(
        RDKit::buildFeatureFactory(preset_feature_defs.at(preset_name)));
    if (!factory) {
        throw std::runtime_error("Failed to create feature factory");
    }
    return factory;
}

PharmacophoreGraph create_pharmacophore_graph(
    const Mol_sptr &mol,
    const RDKit::MolChemicalFeatureFactory &feature_factory,
    const BondWeights &bond_weights, float default_bond_weight) {
    Ensures(mol != nullptr);
    RDKit::FeatSPtrList features = feature_factory.getFeaturesForMol(*mol);
    return create_pharmacophore_graph(mol, features, bond_weights,
                                      default_bond_weight);
}

template <typename FeatSeq>
PharmacophoreGraph create_pharmacophore_graph(
    const Mol_sptr &mol, const FeatSeq &features,
    const std::map<RDKit::Bond::BondType, float> &bond_weights,
    float default_bond_weight) {

    PharmacophoreGraph graph(features);
    std::vector<RDKit::FeatSPtr> feature_vec(features.begin(), features.end());

    for (size_t i = 0; i < graph.nodes.size(); ++i) {
        for (size_t j = i + 1; j < graph.nodes.size(); ++j) {
            const auto &feature_i = feature_vec[i];
            const auto &feature_j = feature_vec[j];

            float distance = average_atoms_to_atoms_distance(
                mol, feature_i->getAtoms(), feature_j->getAtoms(), bond_weights,
                default_bond_weight);

            if (distance >= 0.0f) {
                graph.set_edge(i, j, distance);
            }
        }
    }

    return graph;
}

template PharmacophoreGraph create_pharmacophore_graph(
    const Mol_sptr &mol, const RDKit::FeatSPtrList &features,
    const BondWeights &bond_weights, float default_bond_weight);

template PharmacophoreGraph create_pharmacophore_graph(
    const Mol_sptr &mol, const std::vector<RDKit::FeatSPtr> &features,
    const BondWeights &bond_weights, float default_bond_weight);

std::pair<float, float>
subgraph_similarity(const PharmacophoreGraph &g1, const PharmacophoreGraph &g2,
                    std::vector<size_t> subgraph_nodes_1,
                    std::vector<size_t> subgraph_nodes_2, float lambda) {
    Ensures(subgraph_nodes_1.size() == subgraph_nodes_2.size());

    std::multiset<std::string> families_1, families_2;
    for (const auto &node : subgraph_nodes_1) {
        Ensures(node < g1.nodes.size());
        families_1.insert(g1.nodes[node].family);
    }
    for (const auto &node : subgraph_nodes_2) {
        Ensures(node < g2.nodes.size());
        families_2.insert(g2.nodes[node].family);
    }
    if (families_1 != families_2) {
        return {0.0f, 0.0f}; // Families do not match
    }

    size_t subgraph_size = subgraph_nodes_1.size();
    int num_edges_1 = 0, num_edges_2 = 0;
    for (size_t i = 0; i < subgraph_size; ++i) {
        size_t u1 = subgraph_nodes_1[i];
        const auto &neighbors_1 = g1.adj_list[u1];
        for (const auto &v1 : subgraph_nodes_1) {
            if (neighbors_1.find(v1) != neighbors_1.end()) {
                num_edges_1++;
            }
        }

        size_t u2 = subgraph_nodes_2[i];
        const auto &neighbors_2 = g2.adj_list[u2];
        for (const auto &v2 : subgraph_nodes_2) {
            if (neighbors_2.find(v2) != neighbors_2.end()) {
                num_edges_2++;
            }
        }
    }

    float score = 0.0f;
    std::sort(subgraph_nodes_2.begin(), subgraph_nodes_2.end());
    do {
        float current_score = 0.0f;
        bool family_matched = true;
        for (size_t i = 0; i < subgraph_size; ++i) {
            if (g1.nodes[subgraph_nodes_1[i]].family !=
                g2.nodes[subgraph_nodes_2[i]].family) {
                family_matched = false;
                break;
            }
        }
        if (!family_matched) {
            continue;
        }

        for (size_t i = 0; i < subgraph_size; ++i) {
            for (size_t j = 0; j < subgraph_size; ++j) {
                size_t u1 = subgraph_nodes_1[i];
                size_t v1 = subgraph_nodes_1[j];
                size_t u2 = subgraph_nodes_2[i];
                size_t v2 = subgraph_nodes_2[j];

                if (g1.has_edge(u1, v1) && g2.has_edge(u2, v2)) {
                    auto w1 = g1.edge_weights.at({u1, v1});
                    auto w2 = g2.edge_weights.at({u2, v2});
                    auto e = std::exp(-std::abs(w1 - w2) * lambda);
                    current_score += e;
                }
            }
        }
        score = std::max(score, current_score);
    } while (std::next_permutation(subgraph_nodes_2.begin(),
                                   subgraph_nodes_2.end()));

    return {score / num_edges_1, score / num_edges_2};
}

typedef std::multiset<std::string> FamilyMultiset;
typedef std::vector<size_t> SubgraphNodes;

FamilyMultiset get_subgraph_family_multiset(const PharmacophoreGraph &g,
                                            const SubgraphNodes &subgraph) {
    FamilyMultiset family_set;
    for (const auto &node : subgraph) {
        Ensures(node < g.nodes.size());
        family_set.insert(g.nodes[node].family);
    }
    return family_set;
}

std::map<FamilyMultiset, std::vector<SubgraphNodes>>
get_subgraphs_by_family(const PharmacophoreGraph &g) {
    std::vector<SubgraphNodes> subgraphs;

    auto add_subgraph = [&](const SubgraphNodes &subgraph) {
        if (g.is_subgraph_connected(subgraph))
            subgraphs.push_back(subgraph);
    };

    for (size_t i = 0; i < g.nodes.size(); ++i) {
        for (size_t j = i + 1; j < g.nodes.size(); ++j) {
            add_subgraph({i, j});
            for (size_t k = j + 1; k < g.nodes.size(); ++k) {
                add_subgraph({i, j, k});
                for (size_t l = k + 1; l < g.nodes.size(); ++l) {
                    add_subgraph({i, j, k, l});
                }
            }
        }
    }

    std::map<FamilyMultiset, std::vector<SubgraphNodes>> subgraphs_map;
    for (const auto &subgraph : subgraphs) {
        auto family_set = get_subgraph_family_multiset(g, subgraph);
        if (subgraphs_map.find(family_set) == subgraphs_map.end()) {
            subgraphs_map[family_set] = {};
        }
        subgraphs_map[family_set].push_back(subgraph);
    }

    return subgraphs_map;
}

std::pair<float, float> pharmacophore_similarity(const PharmacophoreGraph &g1,
                                                 const PharmacophoreGraph &g2,
                                                 float lambda) {
    std::map<FamilyMultiset, std::vector<SubgraphNodes>>
        subgraph_map_1 = get_subgraphs_by_family(g1),
        subgraph_map_2 = get_subgraphs_by_family(g2);

    int n_subgraphs_1 = 0, n_subgraphs_2 = 0;
    for (const auto &pair : subgraph_map_1) {
        n_subgraphs_1 += pair.second.size();
    }
    for (const auto &pair : subgraph_map_2) {
        n_subgraphs_2 += pair.second.size();
    }

    std::set<FamilyMultiset> all_families;
    for (const auto &pair : subgraph_map_1) {
        all_families.insert(pair.first);
    }
    for (const auto &pair : subgraph_map_2) {
        all_families.insert(pair.first);
    }

    float score_1 = 0.0f, score_2 = 0.0f;
    for (const auto &family : all_families) {
        const auto &it1 = subgraph_map_1.find(family);
        const auto &it2 = subgraph_map_2.find(family);
        if (it1 == subgraph_map_1.end() || it2 == subgraph_map_2.end())
            continue;

        std::vector<float> s1(it1->second.size(), 0.0f);
        std::vector<float> s2(it2->second.size(), 0.0f);
        for (size_t i = 0; i < it1->second.size(); ++i) {
            const auto &subgraph_1 = it1->second[i];
            for (size_t j = 0; j < it2->second.size(); ++j) {
                const auto &subgraph_2 = it2->second[j];
                auto [sim1, sim2] =
                    subgraph_similarity(g1, g2, subgraph_1, subgraph_2, lambda);
                s1[i] = std::max(s1[i], sim1);
                s2[j] = std::max(s2[j], sim2);
            }
        }
        score_1 += std::accumulate(s1.begin(), s1.end(), 0.0f);
        score_2 += std::accumulate(s2.begin(), s2.end(), 0.0f);
    }
    return {score_1 / n_subgraphs_1, score_2 / n_subgraphs_2};
}
} // namespace prexsyn_engine
