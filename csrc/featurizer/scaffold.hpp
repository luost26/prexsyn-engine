#pragma once

#include "../feature/builder.hpp"
#include "../fingerprints.hpp"
#include "base.hpp"

namespace prexsyn_engine {
class MurckoScaffoldFeaturizer : public Featurizer {
    FpFunc<float> fp_func;

  public:
    std::string name;

    MurckoScaffoldFeaturizer(const std::string &name,
                             const std::string &fp_type);
    void operator()(const Synthesis &syn, FeatureBuilder &dict) override;
};
} // namespace prexsyn_engine
