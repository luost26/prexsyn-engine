#include "product_pharmacophore.hpp"

#include <algorithm>
#include <random>
#include <tuple>

#include <GraphMol/MolChemicalFeatures/MolChemicalFeature.h>

#include "../utils/algorithm.hpp"
#include "../utils/logging.hpp"

namespace prexsyn_engine {

ProductPharmacophoreFeaturizer::ProductPharmacophoreFeaturizer(
    const ProductPharmacophoreFeaturizerOption &option)
    : name(option.name), type_index_offset(option.type_index_offset),
      random_num_nodes_max(option.random_num_nodes_max),
      random_num_nodes_min(option.random_num_nodes_min),
      random_num_edges_max(option.random_num_edges_max),
      bond_weights(option.bond_weights),
      default_bond_weight(option.default_bond_weight) {
    feature_factory.reset(RDKit::buildFeatureFactory(option.feature_def));

    for (auto feature_it = feature_factory->beginFeatureDefs();
         feature_it != feature_factory->endFeatureDefs(); ++feature_it) {

        logger()->info("Pharmacophoric family: {}, type: {}",
                       (*feature_it)->getFamily(), (*feature_it)->getType());
        std::string family = (*feature_it)->getFamily();
        std::string type = (*feature_it)->getType();
        family_type_to_index[std::make_tuple(family, type)] =
            family_type_to_index.size() + type_index_offset;
        if (family_to_index.find(family) == family_to_index.end()) {
            family_to_index[family] =
                family_to_index.size() + type_index_offset;
        }
    }
}

void ProductPharmacophoreFeaturizer::operator()(const Synthesis &syn,
                                                FeatureBuilder &dict) {
    std::mt19937 rng{std::random_device{}()};
    auto product = random_choice(syn.top(), rng);

    std::vector<Bool> node_exists(random_num_nodes_max, false);
    std::vector<Long> family_indices(random_num_nodes_max);
    std::vector<Long> type_indices(random_num_nodes_max);

    auto _features = feature_factory->getFeaturesForMol(*product);
    std::vector<RDKit::FeatSPtr> features(_features.begin(), _features.end());
    std::shuffle(features.begin(), features.end(), rng);
    size_t num_nodes = std::uniform_int_distribution<size_t>(
        random_num_nodes_min, random_num_nodes_max)(rng);
    if (features.size() > num_nodes) {
        features.resize(num_nodes);
    }

    for (size_t i = 0; i < features.size(); ++i) {
        node_exists[i] = true;
        auto feature_sptr = features[i];
        family_indices[i] = family_to_index.at(feature_sptr->getFamily());
        type_indices[i] = family_type_to_index.at(std::make_tuple(
            feature_sptr->getFamily(), feature_sptr->getType()));
    }

    PharmacophoreGraph graph = create_pharmacophore_graph(
        product, features, bond_weights, default_bond_weight);

    auto mst = graph.get_min_spanning_tree();
    std::vector<std::pair<size_t, size_t>> mst_edges{}, other_edges{};
    for (const auto &edge : graph.edge_weights) {
        if (edge.first.first > edge.first.second) {
            continue;
        }
        if (mst.has_edge(edge.first.first, edge.first.second)) {
            mst_edges.push_back(edge.first);
        } else {
            other_edges.push_back(edge.first);
        }
    }
    std::shuffle(other_edges.begin(), other_edges.end(), rng);

    std::vector<std::pair<size_t, size_t>> selected_edges;
    selected_edges.insert(selected_edges.end(), mst_edges.begin(),
                          mst_edges.end());
    selected_edges.insert(selected_edges.end(), other_edges.begin(),
                          other_edges.end());

    std::vector<Bool> edge_exists(random_num_edges_max, false);
    std::vector<Long> edge_indices_u(random_num_edges_max);
    std::vector<Long> edge_indices_v(random_num_edges_max);
    std::vector<Float> edge_weights(random_num_edges_max);

    for (size_t i = 0; i < random_num_edges_max && i < selected_edges.size();
         ++i) {
        auto edge = selected_edges[i];
        edge_exists[i] = true;
        edge_indices_u[i] = edge.first;
        edge_indices_v[i] = edge.second;
        edge_weights[i] = graph.edge_weights.at(edge);
    }

    // dict.add(name + "pharmacophore_types", type_indices);
    dict.add(name + ".node_features", family_indices);
    dict.add(name + ".node_exists", node_exists);
    dict.add(name + ".edge_features", edge_weights);
    dict.add(name + ".edge_u", edge_indices_u);
    dict.add(name + ".edge_v", edge_indices_v);
    dict.add(name + ".edge_exists", edge_exists);
}

void ProductPharmacophoreFeaturizer::operator()(const PharmacophoreGraph &graph,
                                                FeatureBuilder &dict) {
    std::vector<Bool> node_exists;
    std::vector<Long> family_indices;
    std::vector<Long> type_indices;
    for (const auto &node : graph.nodes) {
        node_exists.push_back(true);
        family_indices.push_back(family_to_index.at(node.family));
        type_indices.push_back(
            family_type_to_index.at(std::make_tuple(node.family, node.type)));
    }

    std::vector<Bool> edge_exists;
    std::vector<Long> edge_indices_u;
    std::vector<Long> edge_indices_v;
    std::vector<Float> edge_weights;
    for (const auto &edge : graph.edge_weights) {
        if (edge.first.first > edge.first.second) {
            continue;
        }
        edge_exists.push_back(true);
        edge_indices_u.push_back(edge.first.first);
        edge_indices_v.push_back(edge.first.second);
        edge_weights.push_back(edge.second);
    }

    dict.add(name + ".node_features", family_indices);
    dict.add(name + ".node_exists", node_exists);
    dict.add(name + ".edge_features", edge_weights);
    dict.add(name + ".edge_u", edge_indices_u);
    dict.add(name + ".edge_v", edge_indices_v);
    dict.add(name + ".edge_exists", edge_exists);
}

const std::map<std::string, size_t> &
ProductPharmacophoreFeaturizer::get_family_to_index() const {
    return family_to_index;
}

const std::map<ProductPharmacophoreFeaturizer::FamilyTypeTuple, size_t> &
ProductPharmacophoreFeaturizer::get_family_type_to_index() const {
    return family_type_to_index;
}

std::vector<RDKit::FeatSPtr>
ProductPharmacophoreFeaturizer::get_features(const Mol_sptr &mol) const {
    Ensures(mol != nullptr);
    auto features = feature_factory->getFeaturesForMol(*mol);
    return std::vector<RDKit::FeatSPtr>(features.begin(), features.end());
}

PharmacophoreGraph ProductPharmacophoreFeaturizer::get_graph(
    const Mol_sptr &mol,
    const std::optional<std::vector<RDKit::FeatSPtr>> &features) const {
    Ensures(mol != nullptr);
    auto features_v =
        features.has_value() ? features.value() : get_features(mol);
    return create_pharmacophore_graph(mol, features_v, bond_weights,
                                      default_bond_weight);
}

} // namespace prexsyn_engine
