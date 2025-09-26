#pragma once

#include "../synthesis.hpp"
#include "dtype.hpp"
#include "featurizer.hpp"

namespace synthesis_backend {
struct PostfixNotationTokenDef {
    int PAD = 0;
    int END = 1;
    int START = 2;
    int BB = 3;
    int RXN = 4;
};

struct PostfixNotationFeaturizerOption {
    unsigned int length = 16;
    PostfixNotationTokenDef token_def{};
};

using TypeToken = Long;
using BuildingBlockToken = Long;
using ReactionToken = Long;

class PostfixNotationFeaturizer : public Featurizer {
  public:
    PostfixNotationFeaturizerOption option;
    PostfixNotationFeaturizer(
        const PostfixNotationFeaturizerOption &option = {});
    void operator()(const Synthesis &synthesis,
                    FeatureBuilder &builder) override;
};
} // namespace synthesis_backend
