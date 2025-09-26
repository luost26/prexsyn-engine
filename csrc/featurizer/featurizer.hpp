#pragma once

#include "../synthesis.hpp"
#include "builder.hpp"

namespace synthesis_backend {
class Featurizer {
  public:
    virtual ~Featurizer() = default;
    virtual void operator()(const Synthesis &, FeatureBuilder &) = 0;
};

class FeaturizerSet : public Featurizer {
    std::vector<std::shared_ptr<Featurizer>> featurizers_;

  public:
    template <typename F> FeaturizerSet &add(std::shared_ptr<F> featurizer) {
        featurizers_.push_back(
            std::dynamic_pointer_cast<Featurizer>(featurizer));
        return *this;
    }
    void operator()(const Synthesis &synthesis,
                    FeatureBuilder &builder) override {
        for (auto &featurizer : featurizers_) {
            (*featurizer)(synthesis, builder);
        }
    }
};
} // namespace synthesis_backend
