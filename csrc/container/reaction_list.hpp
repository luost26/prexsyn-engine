#pragma once

#include <filesystem>
#include <optional>

#include "../types.hpp"

namespace synthesis_backend {

using ReactionIndex = size_t;

class ReactionList {
    ReactionVector reactions;

    std::optional<ReactionIndex> add(const Reaction_sptr &reaction);

  public:
    ReactionList(const ReactionVector &reactions);
    static ReactionList *from_txt(const std::filesystem::path &path);
    const ReactionVector &get_reactions() const;
    void save(const std::filesystem::path &path) const;
    static ReactionList *load(const std::filesystem::path &path);
    Reaction_sptr get(size_t index) const;
    ReactionIndex size() const;
};
} // namespace synthesis_backend
