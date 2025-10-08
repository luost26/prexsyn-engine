#pragma once

#include <GraphMol/Descriptors/Property.h>

#include "builder.hpp"
#include "featurizer.hpp"

namespace prexsyn_engine {

struct ProductStructureFeaturizerOption {
    std::vector<std::string> fp_types = {"ecfp4", "fcfp4"};
    std::string embedder_name_template = "product_{}_fingerprint";

    bool scaffold = true;
    std::string scaffold_fp_type = "ecfp4";
    std::string scaffold_embedder_name = "product_scaffold_fingerprint";

    unsigned int num_fragments = 8;
    std::string fragment_fp_type = "ecfp4";
    std::string fragment_embedder_name = "product_fragment_fingerprints";
};

class ProductStructureFeaturizer : public Featurizer {
    ProductStructureFeaturizerOption option;

  public:
    ProductStructureFeaturizer(
        const ProductStructureFeaturizerOption &option = {});
    void operator()(const Synthesis &syn, FeatureBuilder &dict) override;
};

} // namespace prexsyn_engine
