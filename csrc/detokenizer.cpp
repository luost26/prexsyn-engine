#include "detokenizer.hpp"

#include <omp.h>

namespace prexsyn_engine {
Detokenizer::Detokenizer(
    const std::shared_ptr<BuildingBlockList> &building_blocks,
    const std::shared_ptr<ReactionList> &reactions,
    const PostfixNotationTokenDef &token_def)
    : building_blocks(building_blocks), reactions(reactions),
      token_def(token_def) {
    Ensures(building_blocks != nullptr);
    Ensures(reactions != nullptr);
}

void Detokenizer::detokenize_one(
    const std::span<TypeToken> &token_types,
    const std::span<BuildingBlockToken> &bb_indices,
    const std::span<ReactionToken> &rxn_indices,
    Synthesis_sptr &synthesis) const {
    Ensures(token_types.size() == bb_indices.size());
    Ensures(token_types.size() == rxn_indices.size());
    Ensures(synthesis != nullptr);

    for (size_t i = 0; i < token_types.size(); ++i) {
        auto token_type = token_types[i];
        if (token_type == token_def.BB) {
            synthesis->push(building_blocks->get(bb_indices[i]));
        } else if (token_type == token_def.RXN) {
            try {
                synthesis->push(reactions->get(rxn_indices[i]));
            } catch (const push_reaction_exception &e) {
                continue;
            }
        } else if (token_type == token_def.END) {
            break;
        } else if (token_type == token_def.START ||
                   token_type == token_def.PAD) {
            continue;
        } else {
            throw std::runtime_error("Unknown token type " +
                                     std::to_string(token_type) + " at index " +
                                     std::to_string(i));
        }
    }
}

Synthesis_sptr
Detokenizer::detokenize_one(const std::span<TypeToken> &token_types,
                            const std::span<BuildingBlockToken> &bb_indices,
                            const std::span<ReactionToken> &rxn_indices) const {
    Synthesis_sptr synthesis = std::make_shared<Synthesis>();
    detokenize_one(token_types, bb_indices, rxn_indices, synthesis);
    return synthesis;
}

SynthesisVector
Detokenizer::detokenize_many(const std::span<TypeToken> &token_types,
                             const std::span<BuildingBlockToken> &bb_indices,
                             const std::span<ReactionToken> &rxn_indices,
                             size_t stride) const {
    Ensures(token_types.size() == bb_indices.size());
    Ensures(token_types.size() == rxn_indices.size());
    Ensures(token_types.size() % stride == 0);

    size_t count = token_types.size() / stride;
    SynthesisVector results;
    results.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        results.emplace_back(new Synthesis());
    }
#pragma omp parallel for
    for (size_t i = 0; i < count; ++i) {
        detokenize_one(token_types.subspan(i * stride, stride),
                       bb_indices.subspan(i * stride, stride),
                       rxn_indices.subspan(i * stride, stride), results[i]);
    }
    return results;
}
} // namespace prexsyn_engine
