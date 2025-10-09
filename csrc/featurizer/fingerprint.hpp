#pragma once

#include "../feature/builder.hpp"
#include "../fingerprints.hpp"
#include "featurizer.hpp"

namespace prexsyn_engine {

class FingerprintFeaturizer : public Featurizer {
    FpFunc<float> fp_func;

  public:
    std::string name;

    FingerprintFeaturizer(const std::string &name, const std::string &fp_type);
    void operator()(const Synthesis &syn, FeatureBuilder &dict) override;
};
} // namespace prexsyn_engine
