#include "chemspace.hpp"

#include <filesystem>

#include <gtest/gtest.h>
#include <spdlog/spdlog.h>

#include "container/building_block_list.hpp"

using namespace synthesis_backend;

TEST(chemspace, DISABLED_debug1) {
    auto sdf_path =
        std::filesystem::path("../data/building_blocks/"
                              "Enamine_Rush-Delivery_Building_Blocks-"
                              "US_223244cmpd_20231001.sdf");
    auto cache_path =
        std::filesystem::path("../data/building_blocks/"
                              "Enamine_Rush-Delivery_Building_Blocks-"
                              "US_223244cmpd_20231001.cache");
    if (!std::filesystem::exists(cache_path)) {
        auto bb_list = BuildingBlockList::from_sdf(sdf_path);
        bb_list->save(cache_path);
    }

    auto builder = ChemicalSpaceDefinitionBuilder();
    builder.building_blocks_from_cache(cache_path)
        .reactions_from_txt("../data/reactions/hartenfeller_button.txt")
        .secondary_building_blocks_from_single_reaction()
        .build_secondary_index();
}
