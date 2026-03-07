#pragma once

#include <cstddef>
#include <map>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "../chemistry/chemistry.hpp"

namespace prexsyn::chemspace {

struct ReactionEntry {
    std::shared_ptr<Reaction> reaction;
    std::string name;
};

struct ReactionItem : ReactionEntry {
    size_t index{};
};

class ReactionLibrary {
public:
    using Index = size_t;

private:
    std::vector<ReactionItem> reactions_;
    std::map<std::string, Index> name_to_index_;

    ReactionLibrary() = default;

public:
    size_t size() const { return reactions_.size(); }
    const ReactionItem &get(Index) const;
    const ReactionItem &get(const std::string &) const;
    Index add(const ReactionEntry &);

    auto begin() const noexcept { return reactions_.begin(); }
    auto end() const noexcept { return reactions_.end(); }

    struct Match {
        Index reaction_index;
        std::string_view reaction_name;
        Reaction::ReactantIndex reactant_index;
        std::string_view reactant_name;
        size_t count;
    };
    std::vector<Match> match_reactants(const Molecule &molecule) const;
};

} // namespace prexsyn::chemspace
