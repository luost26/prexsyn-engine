#pragma once

#include "../feature/builder.hpp"
#include "../fingerprints.hpp"
#include "base.hpp"

namespace prexsyn_engine {
class BRICSFragmentsFeaturizer : public Featurizer {

    std::string name;
    unsigned int max_num_fragments;
    FpFunc<float> fp_func;

  public:
    BRICSFragmentsFeaturizer(const std::string &name,
                             const std::string &fp_type,
                             unsigned int max_num_fragments = 8);
    void operator()(const Synthesis &syn, FeatureBuilder &dict) override;
};
} // namespace prexsyn_engine
