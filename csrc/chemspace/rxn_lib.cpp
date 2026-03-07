#include "rxn_lib.hpp"

#include <cstddef>
#include <istream>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../chemistry/chemistry.hpp"
#include "../utility/serialization.hpp"

namespace prexsyn::chemspace {

std::unique_ptr<ReactionLibrary> ReactionLibrary::deserialize(std::istream &data) {
    boost::archive::binary_iarchive ia(data);
    size_t num_items = 0;
    ia >> num_items;
    auto rxn_lib = std::make_unique<ReactionLibrary>();
    rxn_lib->reactions_.reserve(num_items);
    for (size_t i = 0; i < num_items; ++i) {
        ReactionItem item;
        ia >> item;
        rxn_lib->reactions_.push_back(std::move(item));
    }
    ia >> rxn_lib->name_to_index_;
    return rxn_lib;
}

void ReactionLibrary::serialize(std::ostream &stream) const {
    boost::archive::binary_oarchive oa(stream);
    oa << reactions_.size();
    for (const auto &item : reactions_) {
        oa << item;
    }
    oa << name_to_index_;
}

const ReactionItem &ReactionLibrary::get(Index index) const {
    if (index >= reactions_.size()) {
        throw std::out_of_range("Reaction index out of range");
    }
    return reactions_[index];
}

const ReactionItem &ReactionLibrary::get(const std::string &name) const {
    auto it = name_to_index_.find(name);
    if (it == name_to_index_.end()) {
        throw std::out_of_range("Reaction name not found: " + name);
    }
    return reactions_[it->second];
}

ReactionLibrary::Index ReactionLibrary::add(const ReactionEntry &entry) {
    if (name_to_index_.contains(entry.name)) {
        throw std::invalid_argument("Reaction with the same name already exists: " + entry.name);
    }
    auto new_index = reactions_.size();
    reactions_.push_back(ReactionItem{entry, new_index});
    name_to_index_[entry.name] = new_index;
    return new_index;
}

std::vector<ReactionLibrary::Match>
ReactionLibrary::match_reactants(const Molecule &molecule) const {
    std::vector<ReactionLibrary::Match> matches;
    for (size_t i = 0; i < reactions_.size(); ++i) {
        const auto &rxn = reactions_[i];
        auto rxn_matches = rxn.reaction->match_reactants(molecule);
        for (const auto &match : rxn_matches) {
            matches.push_back({
                .reaction_index = i,
                .reaction_name = rxn.name,
                .reactant_index = match.index,
                .reactant_name = match.name,
                .count = match.count,
            });
        }
    }
    return matches;
}

} // namespace prexsyn::chemspace
