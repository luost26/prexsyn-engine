#pragma once

#include <cstddef>
#include <filesystem>
#include <map>
#include <memory>
#include <optional>
#include <tuple>

#include <GraphMol/MolChemicalFeatures/MolChemicalFeatureFactory.h>

#include "../pharmacophore.hpp"
#include "builder.hpp"
#include "featurizer.hpp"

namespace prexsyn_engine {
struct ProductPharmacophoreFeaturizerOption {
    std::string name = "product_pharmacophore";
    std::string feature_def = fdef_BASE;

    size_t type_index_offset = 1;
    size_t random_num_nodes_max = 8;
    size_t random_num_nodes_min = 3;
    size_t random_num_edges_max = 16;

    BondWeights bond_weights = bond_weights_RELATIVE_BOND_LENGTH;
    float default_bond_weight = 1.0f;
};

class ProductPharmacophoreFeaturizer : public Featurizer {
  private:
    std::shared_ptr<RDKit::MolChemicalFeatureFactory> feature_factory;

    std::string name;

    typedef std::tuple<std::string, std::string> FamilyTypeTuple;
    const size_t type_index_offset;
    std::map<std::string, size_t> family_to_index;
    std::map<FamilyTypeTuple, size_t> family_type_to_index;

    size_t random_num_nodes_max;
    size_t random_num_nodes_min;
    size_t random_num_edges_max;

    BondWeights bond_weights;
    float default_bond_weight;

  public:
    ProductPharmacophoreFeaturizer(
        const ProductPharmacophoreFeaturizerOption &option);
    void operator()(const Synthesis &syn, FeatureBuilder &dict) override;
    void operator()(const PharmacophoreGraph &graph, FeatureBuilder &dict);
    const std::map<std::string, size_t> &get_family_to_index() const;
    const std::map<FamilyTypeTuple, size_t> &get_family_type_to_index() const;
    std::vector<RDKit::FeatSPtr> get_features(const Mol_sptr &mol) const;
    PharmacophoreGraph
    get_graph(const Mol_sptr &,
              const std::optional<std::vector<RDKit::FeatSPtr>> & =
                  std::nullopt) const;
};
} // namespace prexsyn_engine
