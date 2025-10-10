#include "base.hpp"

namespace prexsyn_engine {
FeaturizerSet &FeaturizerSet::add(std::shared_ptr<Featurizer> featurizer) {
    featurizers_.push_back(featurizer);
    return *this;
}

void FeaturizerSet::operator()(const Synthesis &synthesis,
                               FeatureBuilder &builder) {
    for (auto &featurizer : featurizers_) {
        (*featurizer)(synthesis, builder);
    }
}
} // namespace prexsyn_engine
