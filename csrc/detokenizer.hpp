#pragma once

#include <span>

#include "container/building_block_list.hpp"
#include "container/reaction_list.hpp"
#include "featurizer/synthesis.hpp"
#include "synthesis.hpp"

namespace synthesis_backend {
class Detokenizer {
  private:
    std::shared_ptr<BuildingBlockList> building_blocks = nullptr;
    std::shared_ptr<ReactionList> reactions = nullptr;
    PostfixNotationTokenDef token_def;

  public:
    Detokenizer(
        const std::shared_ptr<BuildingBlockList> &building_blocks,
        const std::shared_ptr<ReactionList> &reactions,
        const PostfixNotationTokenDef &token_def = PostfixNotationTokenDef());
    void detokenize_one(const std::span<TypeToken> &token_types,
                        const std::span<BuildingBlockToken> &bb_indices,
                        const std::span<ReactionToken> &rxn_indices,
                        Synthesis_sptr &synthesis) const;
    Synthesis_sptr
    detokenize_one(const std::span<TypeToken> &token_types,
                   const std::span<BuildingBlockToken> &bb_indices,
                   const std::span<ReactionToken> &rxn_indices) const;
    SynthesisVector
    detokenize_many(const std::span<TypeToken> &token_types,
                    const std::span<BuildingBlockToken> &bb_indices,
                    const std::span<ReactionToken> &rxn_indices,
                    size_t stride) const;
};
} // namespace synthesis_backend
