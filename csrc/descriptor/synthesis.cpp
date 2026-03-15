#include "synthesis.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <span>

#include "../chemspace/chemspace.hpp"

namespace prexsyn::descriptor {

void SynthesisPostfixNotation::operator()(const chemspace::Synthesis &syn,
                                          std::span<std::byte> &out_raw) const {
    check_size(out_raw);
    auto out = cast<std::int64_t>(out_raw);

    const auto &tokens = syn.postfix_notation().tokens();

    std::fill(out.begin(), out.end(), token_def_.pad);
    out[0] = token_def_.start;

    size_t i = 1;
    for (; i <= std::min(tokens.size(), max_length_ - 1); ++i) {
        const auto &token = tokens[i - 1];
        if (token.type == chemspace::PostfixNotation::Token::Type::BuildingBlock) {
            out[i * 3] = token_def_.bb;
            out[(i * 3) + 1] = static_cast<std::int64_t>(token.index);
        } else if (token.type == chemspace::PostfixNotation::Token::Type::Reaction) {
            out[i * 3] = token_def_.rxn;
            out[(i * 3) + 2] = static_cast<std::int64_t>(token.index);
        }
    }

    if (i < max_length_) {
        out[i * 3] = token_def_.end;
    }
}

} // namespace prexsyn::descriptor
