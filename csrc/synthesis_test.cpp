#include "synthesis.hpp"

#include <gtest/gtest.h>
#include <memory>
#include <set>

#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include "types.hpp"

using namespace synthesis_backend;

Synthesis get_sample_1() {
    Mol_sptr mol1(RDKit::SmilesToMol("COC(=O)c1nscc1N"));
    Mol_sptr mol2(RDKit::SmilesToMol("CN(C)CCCN"));
    Reaction_sptr rxn(RDKit::RxnSmartsToChemicalReaction(
        "[NH2:1][c:2][c:3][C:4](=[O:5])[O][C].[C:6][NH2:7]>>[NH:1]1[c:2][c:3]["
        "C:4](=[O:5])[N:7]([C:6])C1=O"));
    rxn->initReactantMatchers();

    Synthesis synthesis;
    synthesis.push(mol1);
    synthesis.push(mol2);
    synthesis.push(rxn);
    return synthesis;
}

TEST(synthesis_test, synthesis_1) {
    Synthesis synthesis(get_sample_1());

    for (const auto &product : synthesis.top()) {
        std::string smiles = RDKit::MolToSmiles(*product);
        EXPECT_EQ(smiles, "CN(C)CCCn1c(=O)[nH]c2csnc2c1=O");
        std::cout << smiles << std::endl;
    }
}

TEST(synthesis_test, pfn_pickle) {
    Mol_sptr mol1(RDKit::SmilesToMol("COC(=O)c1nscc1N"));
    Mol_sptr mol2(RDKit::SmilesToMol("CN(C)CCCN"));
    Reaction_sptr rxn(RDKit::RxnSmartsToChemicalReaction(
        "[NH2:1][c:2][c:3][C:4](=[O:5])[O][C].[C:6][NH2:7]>>[NH:1]1[c:2][c:3]["
        "C:4](=[O:5])[N:7]([C:6])C1=O"));
    rxn->initReactantMatchers();
    PostfixNotation pfn;
    pfn.append(mol1);
    pfn.append(mol2);
    pfn.append(rxn);

    std::string pickle = pfn.pickle();
    auto pfn2 =
        std::unique_ptr<PostfixNotation>(PostfixNotation::unpickle(pickle));
    EXPECT_NE(pfn2, nullptr);

    EXPECT_EQ(pfn2->type(0), PostfixNotation::ItemType::Molecule);
    EXPECT_EQ(pfn2->type(1), PostfixNotation::ItemType::Molecule);
    EXPECT_EQ(pfn2->type(2), PostfixNotation::ItemType::Reaction);

    EXPECT_EQ(RDKit::MolToSmiles(*std::get<Mol_sptr>((*pfn2)[0])),
              RDKit::MolToSmiles(*mol1));
    EXPECT_EQ(RDKit::MolToSmiles(*std::get<Mol_sptr>((*pfn2)[1])),
              RDKit::MolToSmiles(*mol2));
    EXPECT_EQ(RDKit::ChemicalReactionToRxnSmarts(
                  *std::get<Reaction_sptr>((*pfn2)[2])),
              RDKit::ChemicalReactionToRxnSmarts(*rxn));
}

TEST(synthesis_test, synthesis_pickle) {
    Synthesis synthesis(get_sample_1());
    std::string pickle = synthesis.pickle();
    auto synthesis2 = std::unique_ptr<Synthesis>(Synthesis::unpickle(pickle));
    EXPECT_NE(synthesis2, nullptr);

    EXPECT_EQ(synthesis.stack_size(), synthesis2->stack_size());

    std::set<std::string> smiles_set1;
    for (const auto &product : synthesis.top()) {
        std::string smiles = RDKit::MolToSmiles(*product);
        smiles_set1.insert(smiles);
    }
    std::set<std::string> smiles_set2;
    for (const auto &product : synthesis2->top()) {
        std::string smiles = RDKit::MolToSmiles(*product);
        smiles_set2.insert(smiles);
    }
    EXPECT_EQ(smiles_set1, smiles_set2);
}
