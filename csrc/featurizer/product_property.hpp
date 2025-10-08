#pragma once

#include <map>

#include <GraphMol/Descriptors/Property.h>

#include "builder.hpp"
#include "featurizer.hpp"

namespace prexsyn_engine {

static const std::vector<std::string> supported_rdkit_properties{
    "amw",
    "lipinskiHBA",
    "lipinskiHBD",
    "NumRotatableBonds",
    "NumHBD",
    "NumHBA",
    "NumHeavyAtoms",
    "NumAtoms",
    "NumHeteroatoms",
    "NumAmideBonds",
    "FractionCSP3",
    "NumRings",
    "NumAromaticRings",
    "NumAliphaticRings",
    "NumSaturatedRings",
    "NumHeterocycles",
    "NumAromaticHeterocycles",
    "NumSaturatedHeterocycles",
    "NumAliphaticHeterocycles",
    "NumSpiroAtoms",
    "NumBridgeheadAtoms",
    "NumAtomStereoCenters",
    "NumUnspecifiedAtomStereoCenters",
    "labuteASA",
    "tpsa",
    "CrippenClogP",
    "CrippenMR",

    "exactmw",
    "chi0v",
    "chi1v",
    "chi2v",
    "chi3v",
    "chi4v",
    "chi0n",
    "chi1n",
    "chi2n",
    "chi3n",
    "chi4n",
    "hallKierAlpha",
    "kappa1",
    "kappa2",
    "kappa3",
    "Phi",
};

struct ProductRDKitPropertyFeaturizerOption {
    std::string name = "product_rdkit_properties";
    unsigned int num_evaluated_properties = 4;

    unsigned int rdkit_property_index_offset = 1;
    std::vector<std::string> rdkit_properties = supported_rdkit_properties;
};

class ProductRDKitPropertyFeaturizer : public Featurizer {
    ProductRDKitPropertyFeaturizerOption option;
    RDKit::Descriptors::Properties rdkit_properties;
    std::map<std::string, size_t> supported_property_name_to_index;

  public:
    ProductRDKitPropertyFeaturizer(
        const ProductRDKitPropertyFeaturizerOption &option = {});
    size_t max_property_index() const;
    void operator()(const Synthesis &syn, FeatureBuilder &dict) override;
};
} // namespace prexsyn_engine
