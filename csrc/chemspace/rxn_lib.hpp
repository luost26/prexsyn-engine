#pragma once

#include <cstddef>
#include <istream>
#include <map>
#include <memory>
#include <ostream>
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

    template <typename Archive> void serialize(Archive &ar, const unsigned int /* version */) {
        if constexpr (Archive::is_saving::value) {
            ar << reaction->serialize();
        } else {
            std::string rxn_data;
            ar >> rxn_data;
            reaction = Reaction::deserialize(rxn_data);
        }
        ar & name;
        ar & index;
    }
};

class ReactionLibrary {
public:
    using Index = size_t;

private:
    std::vector<ReactionItem> reactions_;
    std::map<std::string, Index> name_to_index_;

public:
    ReactionLibrary() = default;

    static std::unique_ptr<ReactionLibrary> deserialize(std::istream &);
    void serialize(std::ostream &) const;

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
