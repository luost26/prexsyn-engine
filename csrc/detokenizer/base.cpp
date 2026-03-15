#include "base.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <stdexcept>

#include "../chemspace/chemspace.hpp"
#include "../descriptor/descriptor.hpp"

namespace prexsyn::detokenizer {

std::unique_ptr<chemspace::Synthesis>
detokenize(const std::span<const std::int64_t> &tokens,
           const std::shared_ptr<chemspace::ChemicalSpace> &cs,
           const descriptor::TokenDef &token_def) {
    auto length = tokens.size() / 3;
    if (tokens.size() != length * 3) {
        throw std::invalid_argument("Token size must be a multiple of 3");
    }

    auto syn = cs->new_synthesis();
    for (size_t i = 0; i < length; ++i) {
        auto token_type = tokens[i * 3];
        auto bb_idx = tokens[(i * 3) + 1];
        auto rxn_idx = tokens[(i * 3) + 2];
        if (token_type == token_def.bb) {
            syn->add_building_block(bb_idx);
        } else if (token_type == token_def.rxn) {
            syn->add_reaction(rxn_idx);
        } else if (token_type == token_def.end) {
            break;
        }
    }

    return syn;
}

} // namespace prexsyn::detokenizer
