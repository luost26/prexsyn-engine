#pragma once

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace prexsyn::chemspace {

class PostfixNotation {
public:
    struct Token {
        enum Type { BuildingBlock = 1, Reaction = 2 }; // NOLINT(performance-enum-size)
        size_t index;
        Type type;

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

    void pop_back() {
        if (tokens_.empty()) {
            throw std::out_of_range("Cannot pop_back from an empty PostfixNotation");
        }
        tokens_.pop_back();
    }
};

} // namespace prexsyn::chemspace
