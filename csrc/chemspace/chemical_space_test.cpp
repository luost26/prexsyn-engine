#include <filesystem>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

#include <gtest/gtest.h>

#include "chemspace.hpp"

namespace {

using prexsyn::chemspace::ChemicalSpace;
using prexsyn::chemspace::IntermediateLibrary;

std::filesystem::path find_project_root() {
    auto current = std::filesystem::absolute(__FILE__);
    while (current.has_parent_path()) {
        current = current.parent_path();
        if (std::filesystem::exists(current / "resources/test/chemspace_small_1/rxn.txt") &&
            std::filesystem::exists(current / "resources/test/chemspace_small_1/bb.sdf")) {
            return current;
        }
    }
    throw std::runtime_error("Could not locate project root from __FILE__");
}

std::unique_ptr<ChemicalSpace> make_test_chemical_space() {
    const auto root = find_project_root();
    const auto rxn_path = root / "resources/test/chemspace_small_1/rxn.txt";
    const auto bb_path = root / "resources/test/chemspace_small_1/bb.sdf";

    auto rxn_lib = prexsyn::chemspace::rxn_lib_from_plain_text(rxn_path);
    auto bb_lib = prexsyn::chemspace::bb_lib_from_sdf(bb_path);
    return std::make_unique<ChemicalSpace>(std::move(bb_lib), std::move(rxn_lib),
                                           std::make_unique<IntermediateLibrary>());
}

} // namespace

TEST(ChemicalSpaceTest, EndToEndWorkflowMatchesMainExample) {
    auto chemspace = make_test_chemical_space();

    chemspace->build_reactant_lists_for_building_blocks();
    EXPECT_GT(chemspace->building_block_reactant_lists().num_matches(), 0U);

    chemspace->generate_intermediates();
    EXPECT_GT(chemspace->int_lib().size(), 0U);

    chemspace->build_reactant_lists_for_intermediates();
    EXPECT_GT(chemspace->intermediate_reactant_lists().num_matches(), 0U);

    std::stringstream ss;
    chemspace->serialize(ss);

    ss.seekg(0);
    const auto stats = ChemicalSpace::peek(ss);
    EXPECT_EQ(stats.num_building_blocks, chemspace->bb_lib().size());
    EXPECT_EQ(stats.num_reactions, chemspace->rxn_lib().size());
    EXPECT_EQ(stats.num_intermediates, chemspace->int_lib().size());

    ss.seekg(0);
    auto deserialized_chemspace = ChemicalSpace::deserialize(ss);
    ASSERT_NE(deserialized_chemspace, nullptr);
    EXPECT_EQ(deserialized_chemspace->bb_lib().size(), chemspace->bb_lib().size());
    EXPECT_EQ(deserialized_chemspace->rxn_lib().size(), chemspace->rxn_lib().size());
    EXPECT_EQ(deserialized_chemspace->int_lib().size(), chemspace->int_lib().size());

    std::ostringstream os;
    deserialized_chemspace->print_reactant_lists(os);
    EXPECT_NE(os.str().find("ReactionA"), std::string::npos);

    auto syn = deserialized_chemspace->new_synthesis();
    const auto bb_result_1 = syn->add_building_block("EN300-250786");
    ASSERT_TRUE(bb_result_1) << bb_result_1.message;
    const auto bb_result_2 = syn->add_building_block("EN300-101318");
    ASSERT_TRUE(bb_result_2) << bb_result_2.message;

    const auto rxn_result = syn->add_reaction("ReactionA");
    ASSERT_TRUE(rxn_result) << rxn_result.message;
    EXPECT_EQ(syn->synthesis().stack_size(), 1U);

    const auto top = syn->synthesis().stack_top();
    ASSERT_GT(top->size(), 0U);
    EXPECT_FALSE(top->at(0)->smiles().empty());

    const auto products = syn->products();
    ASSERT_FALSE(products.empty());
    EXPECT_FALSE(products.at(0)->smiles().empty());

    prexsyn::chemspace::ChemicalSpaceSynthesis syn2 = *syn;
    EXPECT_EQ(syn2.synthesis().stack_size(), syn->synthesis().stack_size());
    EXPECT_EQ(syn2.count_building_blocks(), syn->count_building_blocks());
    EXPECT_EQ(syn2.count_reactions(), syn->count_reactions());
}
