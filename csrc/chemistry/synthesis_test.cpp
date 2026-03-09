#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "chemistry.hpp"

namespace {

using prexsyn::Molecule;
using prexsyn::Reaction;
using prexsyn::Synthesis;
using prexsyn::SynthesisError;

constexpr const char *kReactionSmarts = "[NH2:1][c:2][c:3][C:4](=[O:5])[O][C].[C:6][NH2:7]>>[NH:1]"
                                        "1[c:2][c:3][C:4](=[O:5])[N:7]([C:6])C1=O";
constexpr const char *kExpectedProductSmiles = "CN(C)CCCn1c(=O)[nH]c2csnc2c1=O";

std::shared_ptr<Reaction> make_test_reaction() {
    return Reaction::from_smarts(kReactionSmarts, {"A", "B"});
}

std::shared_ptr<Molecule> make_reactant_a() { return Molecule::from_smiles("COC(=O)c1nscc1N"); }

std::shared_ptr<Molecule> make_reactant_b() { return Molecule::from_smiles("CN(C)CCCN"); }

std::shared_ptr<Molecule> make_non_matching_reactant() { return Molecule::from_smiles("CC"); }

} // namespace

TEST(SynthesisTest, PushReactionBuildsExpectedTopNodeAndPrecursors) {
    Synthesis synthesis;
    const auto reactant_a = make_reactant_a();
    const auto reactant_b = make_reactant_b();
    const auto reaction = make_test_reaction();

    synthesis.push(reactant_a);
    synthesis.push(reactant_b);
    synthesis.push(reaction);

    ASSERT_EQ(synthesis.stack_size(), 1);
    const auto &top = synthesis.stack_top();
    ASSERT_EQ(top->size(), 1);
    EXPECT_EQ(top->at(0)->smiles(), kExpectedProductSmiles);

    const auto precursors = top->precursors(0);
    ASSERT_EQ(precursors.size(), 2);

    EXPECT_EQ(precursors.at(0).precursor_index, 0);
    EXPECT_EQ(precursors.at(0).reactant_name, "B");
    EXPECT_EQ(precursors.at(0).item_index, 0);
    EXPECT_EQ(precursors.at(0).molecule->smiles(), reactant_b->smiles());

    EXPECT_EQ(precursors.at(1).precursor_index, 1);
    EXPECT_EQ(precursors.at(1).reactant_name, "A");
    EXPECT_EQ(precursors.at(1).item_index, 0);
    EXPECT_EQ(precursors.at(1).molecule->smiles(), reactant_a->smiles());
}

TEST(SynthesisTest, PushReactionThrowsWhenStackHasTooFewReactants) {
    Synthesis synthesis;
    synthesis.push(make_reactant_a());

    EXPECT_THROW({ synthesis.push(make_test_reaction()); }, SynthesisError);
}

TEST(SynthesisTest, PushReactionThrowsWhenReactionDoesNotProduceAnyProducts) {
    Synthesis synthesis;
    synthesis.push(make_non_matching_reactant());
    synthesis.push(make_non_matching_reactant());

    EXPECT_THROW({ synthesis.push(make_test_reaction()); }, SynthesisError);
}

TEST(SynthesisTest, UndoRestoresPrecursorNodesInOriginalStackOrder) {
    Synthesis synthesis;
    const auto reactant_a = make_reactant_a();
    const auto reactant_b = make_reactant_b();
    const auto reaction = make_test_reaction();

    synthesis.push(reactant_a);
    synthesis.push(reactant_b);
    synthesis.push(reaction);

    ASSERT_EQ(synthesis.nodes().size(), 3);
    ASSERT_EQ(synthesis.stack_size(), 1);

    synthesis.undo();

    ASSERT_EQ(synthesis.nodes().size(), 2);
    ASSERT_EQ(synthesis.stack_size(), 2);
    EXPECT_EQ(synthesis.stack_top(0)->at(0)->smiles(), reactant_b->smiles());
    EXPECT_EQ(synthesis.stack_top(1)->at(0)->smiles(), reactant_a->smiles());
}

TEST(SynthesisTest, PrecursorMoleculesThrowsOnInvalidItemIndex) {
    Synthesis synthesis;
    synthesis.push(make_reactant_a());
    synthesis.push(make_reactant_b());
    synthesis.push(make_test_reaction());

    const auto &top = synthesis.stack_top();
    ASSERT_EQ(top->size(), 1);
    EXPECT_THROW(
        {
            const auto unused = top->precursors(1);
            (void)unused;
        },
        std::out_of_range);
}
