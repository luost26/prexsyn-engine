#pragma once

#include "../feature/builder.hpp"
#include "../synthesis.hpp"

namespace prexsyn_engine {
class Featurizer {
  public:
    virtual ~Featurizer() = default;
    virtual void operator()(const Synthesis &, FeatureBuilder &) = 0;
};

class FeaturizerSet : public Featurizer {
    std::vector<std::shared_ptr<Featurizer>> featurizers_;

  public:
    FeaturizerSet &add(std::shared_ptr<Featurizer> featurizer);
    void operator()(const Synthesis &synthesis,
                    FeatureBuilder &builder) override;
};
} // namespace prexsyn_engine
