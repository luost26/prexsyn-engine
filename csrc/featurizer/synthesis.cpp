#include "synthesis.hpp"

#include <optional>
#include <random>
#include <variant>
#include <vector>

#include "../utils/assert.hpp"

namespace prexsyn_engine {

PostfixNotationTokenDef::PostfixNotationTokenDef(int pad, int end, int start,
                                                 int bb, int rxn)
    : PAD(pad), END(end), START(start), BB(bb), RXN(rxn) {}

PostfixNotationFeaturizer::PostfixNotationFeaturizer(
    unsigned int max_length, const PostfixNotationTokenDef &token_def)
    : max_length(max_length), token_def(token_def) {}

void PostfixNotationFeaturizer::operator()(const Synthesis &synthesis,
                                           FeatureBuilder &builder) {
    std::vector<TypeToken> token_types;
    std::vector<BuildingBlockToken> bb_indices;
    std::vector<ReactionToken> rxn_indices;

    auto add_token =
        [&](TypeToken token_type,
            std::optional<BuildingBlockToken> bb_index = std::nullopt,
            std::optional<ReactionToken> rxn_index = std::nullopt) {
            if (token_types.size() >= max_length) {
                return;
            }
            token_types.push_back(token_type);
            if (token_type == token_def.BB) {
                Ensures(bb_index.has_value());
                bb_indices.push_back(bb_index.value());
            } else {
                bb_indices.push_back(0);
            }
            if (token_type == token_def.RXN) {
                Ensures(rxn_index.has_value());
                rxn_indices.push_back(rxn_index.value());
            } else {
                rxn_indices.push_back(0);
            }
        };

    add_token(token_def.START);
    const auto &pfn = synthesis.get_postfix_notation();
    for (size_t i = 0; i < pfn.size(); ++i) {
        const auto &item = pfn[i];
        if (std::holds_alternative<Mol_sptr>(item)) {
            auto mol = std::get<Mol_sptr>(item);
            Ensures(mol->hasProp("building_block_index"));
            add_token(token_def.BB, mol->getProp<int>("building_block_index"));
        } else {
            auto rxn = std::get<Reaction_sptr>(item);
            Ensures(rxn->hasProp("reaction_index"));
            add_token(token_def.RXN, std::nullopt,
                      rxn->getProp<int>("reaction_index"));
        }
    }
    add_token(token_def.END);

    if (token_types.size() > max_length) {
        std::mt19937 rng{std::random_device{}()};
        size_t offset = std::uniform_int_distribution<size_t>(
            0, token_types.size() - max_length)(rng);
        token_types.erase(token_types.begin(), token_types.begin() + offset);
        token_types.resize(max_length);
        bb_indices.erase(bb_indices.begin(), bb_indices.begin() + offset);
        bb_indices.resize(max_length);
        rxn_indices.erase(rxn_indices.begin(), rxn_indices.begin() + offset);
        rxn_indices.resize(max_length);
    } else {
        for (size_t i = token_types.size(); i < max_length; ++i) {
            add_token(token_def.PAD);
        }
    }

    builder.add("synthesis.token_types", token_types);
    builder.add("synthesis.bb_indices", bb_indices);
    builder.add("synthesis.rxn_indices", rxn_indices);
}
} // namespace prexsyn_engine
