#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>

#include "../chemspace/chemspace.hpp"
#include "base.hpp"

namespace prexsyn::descriptor {

struct TokenDef {
    std::int64_t pad = 0;
    std::int64_t end = 1;
    std::int64_t start = 2;
    std::int64_t bb = 3;
    std::int64_t rxn = 4;
};

class SynthesisPostfixNotation : public SynthesisDescriptor<std::int64_t> {
private:
    TokenDef token_def_;
    size_t max_length_;

public:
    SynthesisPostfixNotation(const TokenDef &token_def, size_t max_length)
        : token_def_(token_def), max_length_(max_length) {
        if (max_length_ < 4) {
            throw std::invalid_argument("max_length must be at least 4");
        }
    };

    size_t size() const override { return max_length_ * 3; /* type + bb_idx + rxn_idx */ }

    void operator()(const chemspace::Synthesis &syn, std::span<std::int64_t> &out) const override;
};

} // namespace prexsyn::descriptor
