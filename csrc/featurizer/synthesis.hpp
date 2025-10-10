#pragma once

#include "../feature/dtype.hpp"
#include "../synthesis.hpp"
#include "base.hpp"

namespace prexsyn_engine {
struct PostfixNotationTokenDef {
    int PAD = 0;
    int END = 1;
    int START = 2;
    int BB = 3;
    int RXN = 4;
};

using TypeToken = Long;
using BuildingBlockToken = Long;
using ReactionToken = Long;

class PostfixNotationFeaturizer : public Featurizer {
  public:
    unsigned int max_length;
    PostfixNotationTokenDef token_def;

    PostfixNotationFeaturizer(unsigned int max_length = 16,
                              const PostfixNotationTokenDef &token_def = {});
    void operator()(const Synthesis &synthesis,
                    FeatureBuilder &builder) override;
};
} // namespace prexsyn_engine
