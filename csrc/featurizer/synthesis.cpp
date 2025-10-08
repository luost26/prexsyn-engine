#include "synthesis.hpp"

#include <optional>
#include <random>
#include <variant>
#include <vector>

#include "../utils/assert.hpp"

namespace prexsyn_engine {

PostfixNotationFeaturizer::PostfixNotationFeaturizer(
    const PostfixNotationFeaturizerOption &option)
    : option(option) {}

void PostfixNotationFeaturizer::operator()(const Synthesis &synthesis,
                                           FeatureBuilder &builder) {
    std::vector<TypeToken> token_types;
    std::vector<BuildingBlockToken> bb_indices;
    std::vector<ReactionToken> rxn_indices;

    auto add_token =
        [&](TypeToken token_type,
            std::optional<BuildingBlockToken> bb_index = std::nullopt,
            std::optional<ReactionToken> rxn_index = std::nullopt) {
            if (token_types.size() >= option.length) {
                return;
            }
            token_types.push_back(token_type);
            if (token_type == option.token_def.BB) {
                Ensures(bb_index.has_value());
                bb_indices.push_back(bb_index.value());
            } else {
                bb_indices.push_back(0);
            }
            if (token_type == option.token_def.RXN) {
                Ensures(rxn_index.has_value());
                rxn_indices.push_back(rxn_index.value());
            } else {
                rxn_indices.push_back(0);
            }
        };

    add_token(option.token_def.START);
    const auto &pfn = synthesis.get_postfix_notation();
    for (size_t i = 0; i < pfn.size(); ++i) {
        const auto &item = pfn[i];
        if (std::holds_alternative<Mol_sptr>(item)) {
            auto mol = std::get<Mol_sptr>(item);
            Ensures(mol->hasProp("building_block_index"));
            add_token(option.token_def.BB,
                      mol->getProp<int>("building_block_index"));
        } else {
            auto rxn = std::get<Reaction_sptr>(item);
            Ensures(rxn->hasProp("reaction_index"));
            add_token(option.token_def.RXN, std::nullopt,
                      rxn->getProp<int>("reaction_index"));
        }
    }
    add_token(option.token_def.END);

    if (token_types.size() > option.length) {
        std::mt19937 rng{std::random_device{}()};
        size_t offset = std::uniform_int_distribution<size_t>(
            0, token_types.size() - option.length)(rng);
        token_types.erase(token_types.begin(), token_types.begin() + offset);
        token_types.resize(option.length);
        bb_indices.erase(bb_indices.begin(), bb_indices.begin() + offset);
        bb_indices.resize(option.length);
        rxn_indices.erase(rxn_indices.begin(), rxn_indices.begin() + offset);
        rxn_indices.resize(option.length);
    } else {
        for (size_t i = token_types.size(); i < option.length; ++i) {
            add_token(option.token_def.PAD);
        }
    }

    builder.add("synthesis.token_types", token_types);
    builder.add("synthesis.bb_indices", bb_indices);
    builder.add("synthesis.rxn_indices", rxn_indices);
}
} // namespace prexsyn_engine
