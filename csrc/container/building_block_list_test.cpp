#include "building_block_list.hpp"

#include <gtest/gtest.h>
#include <spdlog/spdlog.h>

#include "../utils/logging.hpp"

using namespace synthesis_backend;

TEST(building_block_list, load_sdf) {
    spdlog::set_level(spdlog::level::debug);
    auto bb_list = *BuildingBlockList::from_sdf(
        std::filesystem::path("../data/building_blocks/mcule_subset.sdf"));
    for (int i = 0; i < 10; i++) {
        logger()->info("Building block #{}: {}", i, bb_list.get(i));
    }
}
