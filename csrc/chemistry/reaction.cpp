#include "reaction.hpp"

#include <GraphMol/ChemReactions/ReactionParser.h>

namespace prexsyn {

std::unique_ptr<Reaction> Reaction::from_smarts(const std::string &smarts) {
    std::shared_ptr<RDKit::ChemicalReaction> rdkit_rxn(RDKit::RxnSmartsToChemicalReaction(smarts));
    if (!rdkit_rxn) {
        throw std::runtime_error("Failed to parse SMARTS: " + smarts);
    }
    return std::make_unique<Reaction>(std::move(rdkit_rxn));
}

} // namespace prexsyn
