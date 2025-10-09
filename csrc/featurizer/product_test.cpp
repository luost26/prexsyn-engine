#include "rdkit_descriptors.hpp"

#include <gtest/gtest.h>

#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include "../feature/builder.hpp"
#include "../synthesis.hpp"
#include "../types.hpp"

using namespace prexsyn_engine;

static Synthesis get_sample_1() {
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

class DummyBuilder : public FeatureBuilder {
  public:
#define DEFINE_ADD_METHOD(T)                                                   \
    void add(const std::string &key, const T &) {                              \
        std::cout << "Adding " #T << key << std::endl;                         \
    }
    DEFINE_ADD_METHOD(long)
    DEFINE_ADD_METHOD(float)
    DEFINE_ADD_METHOD(bool)
    DEFINE_ADD_METHOD(std::vector<long>)
    DEFINE_ADD_METHOD(std::vector<float>)
    DEFINE_ADD_METHOD(std::vector<bool>)
    DEFINE_ADD_METHOD(std::vector<std::vector<long>>)
    DEFINE_ADD_METHOD(std::vector<std::vector<float>>)
    DEFINE_ADD_METHOD(std::vector<std::vector<bool>>)
#undef DEFINE_ADD_METHOD
};

TEST(featurizer__product, rdkit_properties) {
    Synthesis synthesis(get_sample_1());
    RDKitDescriptorsFeaturizer featurizer("rdkit_descriptors");
    DummyBuilder builder;
    featurizer(synthesis, builder);
}
