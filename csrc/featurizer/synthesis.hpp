#pragma once

#include "../feature/dtype.hpp"
#include "../synthesis.hpp"
#include "base.hpp"

namespace prexsyn_engine {

const int DEFAULT_PAD = 0;
const int DEFAULT_END = 1;
const int DEFAULT_START = 2;
const int DEFAULT_BB = 3;
const int DEFAULT_RXN = 4;

struct PostfixNotationTokenDef {
    int PAD;
    int END;
    int START;
    int BB;
    int RXN;
    PostfixNotationTokenDef(int pad = DEFAULT_PAD, int end = DEFAULT_END,
                            int start = DEFAULT_START, int bb = DEFAULT_BB,
                            int rxn = DEFAULT_RXN);
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
