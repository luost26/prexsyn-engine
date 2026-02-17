#pragma once

#include <GraphMol/ChemReactions/Reaction.h>

namespace prexsyn {

class Reaction {
private:
    std::shared_ptr<RDKit::ChemicalReaction> rdkit_rxn_;

public:
    Reaction(std::shared_ptr<RDKit::ChemicalReaction> rdkit_rxn)
        : rdkit_rxn_(std::move(rdkit_rxn)) {
        if (!rdkit_rxn_) {
            throw std::runtime_error("RDKit reaction pointer is null");
        }
    }
    static std::unique_ptr<Reaction> from_smarts(const std::string &smarts);

    const RDKit::ChemicalReaction &rdkit_rxn() const { return *rdkit_rxn_; }
    RDKit::ChemicalReaction &rdkit_rxn() { return *rdkit_rxn_; }
    std::shared_ptr<RDKit::ChemicalReaction> rdkit_rxn_ptr() const { return rdkit_rxn_; }
};

} // namespace prexsyn
