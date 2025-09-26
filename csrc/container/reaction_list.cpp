#include "reaction_list.hpp"

#include <memory>

#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/ChemReactions/ReactionPickler.h>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include "../types.hpp"
#include "../utils/logging.hpp"

namespace synthesis_backend {

std::optional<ReactionIndex> ReactionList::add(const Reaction_sptr &reaction) {
    if (!reaction) {
        return std::nullopt;
    }
    Reaction_sptr reaction_copy =
        std::make_shared<Reaction_sptr::element_type>(*reaction);
    reaction_copy->initReactantMatchers();
    if (!reaction_copy->isInitialized()) {
        return std::nullopt;
    }
    auto next_index = reactions.size();
    reaction_copy->setProp<int>("reaction_index", (int)next_index);
    reactions.push_back(reaction_copy);
    return next_index;
}

ReactionList::ReactionList(const ReactionVector &reactions) {
    int original_index = 0;
    for (const auto &reaction : reactions) {
        reaction->setProp<int>("original_index", original_index);
        add(reaction);
    }
}

ReactionList *ReactionList::from_txt(const std::filesystem::path &path) {
    logger()->info("ReactionList: Loading reactions from text file {}",
                   path.string());
    std::ifstream ifs(path);
    if (!ifs) {
        throw std::runtime_error("Failed to open file for reading");
    }
    std::string line;
    std::vector<Reaction_sptr> reactions;
    while (std::getline(ifs, line)) {
        if (line.empty()) {
            continue;
        }
        reactions.emplace_back(RDKit::RxnSmartsToChemicalReaction(line));
    }
    ifs.close();
    logger()->info("ReactionList: Loaded {} reactions", reactions.size());
    return new ReactionList(reactions);
}

const ReactionVector &ReactionList::get_reactions() const { return reactions; }

void ReactionList::save(const std::filesystem::path &path) const {
    std::ofstream ofs(path, std::ios::binary);
    if (!ofs) {
        throw std::runtime_error("Failed to open file for writing");
    }
    boost::archive::binary_oarchive oa(ofs);
    size_t num_reactions = reactions.size();
    oa << num_reactions;
    for (const auto &reaction : reactions) {
        std::string pickle;
        RDKit::ReactionPickler::pickleReaction(*reaction, pickle,
                                               RDKit::PicklerOps::AllProps);
        oa << pickle;
    }
    ofs.close();
}

ReactionList *ReactionList::load(const std::filesystem::path &path) {
    logger()->info("ReactionList: Loading reactions from cache {}",
                   path.string());
    std::ifstream ifs(path, std::ios::binary);
    if (!ifs) {
        throw std::runtime_error("Failed to open file for reading");
    }
    boost::archive::binary_iarchive ia(ifs);
    size_t num_reactions;
    ia >> num_reactions;
    auto object = new ReactionList({});
    object->reactions.reserve(num_reactions);
    for (size_t i = 0; i < num_reactions; ++i) {
        std::string pickle;
        ia >> pickle;
        Reaction_sptr reaction =
            std::make_shared<Reaction_sptr::element_type>();
        RDKit::ReactionPickler::reactionFromPickle(pickle, reaction.get());
        if (!reaction) {
            throw std::runtime_error("Failed to unpickle reaction");
        }
        auto result = object->add(reaction);
        if (!result) {
            throw std::runtime_error("Failed to add unpickled reaction");
        }
    }
    ifs.close();
    logger()->info("ReactionList: Loaded {} reactions",
                   object->reactions.size());
    return object;
}

Reaction_sptr ReactionList::get(size_t index) const {
    if (index >= reactions.size()) {
        throw std::out_of_range("Index out of range in ReactionList::get " +
                                std::to_string(index) +
                                " >= " + std::to_string(reactions.size()));
    }
    return reactions[index];
}

ReactionIndex ReactionList::size() const { return reactions.size(); }
} // namespace synthesis_backend
