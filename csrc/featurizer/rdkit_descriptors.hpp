#pragma once

#include <map>

#include <GraphMol/Descriptors/Property.h>

#include "../feature/builder.hpp"
#include "base.hpp"

namespace prexsyn_engine {
static const std::vector<std::string> SUPPORTED_RDKIT_DESCRIPTORS = {
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

class RDKitDescriptorsFeaturizer : public Featurizer {
    RDKit::Descriptors::Properties rdkit_properties_object;
    std::map<std::string, size_t> descriptor_name_to_index;
    std::vector<std::string> descriptor_names;

  public:
    std::string name;
    unsigned int num_evaluated_descriptors;

    RDKitDescriptorsFeaturizer(
        const std::string &name, unsigned int num_evaluated_descriptors = 4,
        const std::vector<std::string> &descriptor_names =
            SUPPORTED_RDKIT_DESCRIPTORS);
    size_t max_descriptor_index() const;
    size_t get_descriptor_index(const std::string &name) const;
    std::vector<std::string> get_descriptor_names() const;
    void operator()(const Synthesis &syn, FeatureBuilder &dict) override;
};
} // namespace prexsyn_engine
