#pragma once

#include <cstddef>
#include <vector>

namespace prexsyn::chemspace {

class PostfixNotation {
public:
    struct Token {
        size_t index;
        enum Type { BuildingBlock = 1, Reaction = 2 } type; // NOLINT(performance-enum-size)

        template <typename Archive> void serialize(Archive &ar, const unsigned int /* version */) {
            ar & index;
            ar & type;
        }
    };

private:
    std::vector<Token> tokens_;

public:
    template <typename Archive> void serialize(Archive &ar, const unsigned int /* version */) {
        ar & tokens_;
    }

    const std::vector<Token> &tokens() const { return tokens_; }

    void append(size_t index, Token::Type type) {
        tokens_.push_back(Token{.index = index, .type = type});
    }

    void extend(const std::vector<Token> &new_tokens) {
        tokens_.insert(tokens_.end(), new_tokens.begin(), new_tokens.end());
    }
};

} // namespace prexsyn::chemspace
