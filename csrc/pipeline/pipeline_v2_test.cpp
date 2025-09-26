#include "pipeline_v2.hpp"

#include <gtest/gtest.h>
#include <spdlog/common.h>

#include "../chemspace.hpp"
#include "../featurizer/product_property.hpp"
#include "../featurizer/product_structure.hpp"
#include "../featurizer/synthesis.hpp"
#include "../utils/logging.hpp"

using namespace synthesis_backend;

TEST(pipeline_v2, case1) {
    logger()->set_level(spdlog::level::info);
    auto csd =
        ChemicalSpaceDefinitionBuilder()
            .building_blocks_from_sdf(
                "../data/building_blocks/mcule_subset.sdf")
            .reactions_from_txt("../data/reactions/hartenfeller_button.txt")
            .secondary_building_blocks_from_single_reaction()
            .build_primary_index()
            .build_secondary_index()
            .build();
    auto featurizer = std::make_shared<FeaturizerSet>();
    featurizer->add(std::make_shared<PostfixNotationFeaturizer>())
        .add(std::make_shared<ProductStructureFeaturizer>())
        .add(std::make_shared<ProductRDKitPropertyFeaturizer>());

    auto ppl =
        DataPipelineV2<32768>(32, csd, SynthesisGeneratorOption(), featurizer);
    ppl.start();
    size_t count = 0;
    size_t n = 128;
    size_t step = 0;
    for (auto i = 0; i < 1000000; i++) {
        ppl.get(
            [&](const std::vector<
                typename DataBuffer<32768>::ReadTransaction::ReadEntry> &) {
                count += n;
                step++;
                logger()->info("Step={}, Entries={}", step, count);
            },
            128);
    }
}
