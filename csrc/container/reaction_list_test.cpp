#include "reaction_list.hpp"

#include <gtest/gtest.h>

#include "../utils/logging.hpp"

using namespace synthesis_backend;

TEST(reaction_list, load_txt) {
    auto reaction_list = *ReactionList::from_txt(
        std::filesystem::path("../data/reactions/hartenfeller_button.txt"));
    int cnt = 0;
    for (const auto &reaction : reaction_list.get_reactions()) {
        logger()->info("{}: {}", cnt++, reaction);
    }
    logger()->info("Number of reactions: {}",
                   reaction_list.get_reactions().size());
}
