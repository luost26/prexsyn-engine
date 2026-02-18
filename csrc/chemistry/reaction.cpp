#include "reaction.hpp"

#include <algorithm>
#include <cstddef>
#include <memory>
#include <numeric>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <GraphMol/ChemReactions/ReactionParser.h>

#include "molecule.hpp"

namespace prexsyn {

std::shared_ptr<Molecule> ReactionOutcome::main_product() const {
    auto max_it = std::ranges::max_element(products, [](const auto &a, const auto &b) {
        return a->num_heavy_atoms() < b->num_heavy_atoms();
    });
    if (max_it == products.end()) {
        throw ReactionError("Reaction outcome has no products");
    }
    return *max_it;
}

std::unique_ptr<Reaction> Reaction::from_smarts(const std::string &smarts) {
    std::shared_ptr<RDKit::ChemicalReaction> rdkit_rxn(RDKit::RxnSmartsToChemicalReaction(smarts));
    if (!rdkit_rxn) {
        throw ReactionError("Failed to parse SMARTS: " + smarts);
    }
    return std::make_unique<Reaction>(std::move(rdkit_rxn));
}

std::vector<ReactionOutcome>
Reaction::apply(const std::vector<std::shared_ptr<Molecule>> &reactants, bool ignore_errors,
                const std::optional<std::vector<size_t>> &reactant_indices) const {
    std::vector<ReactionOutcome> outcomes;

    std::vector<RDKit::ROMOL_SPTR> rdk_reactants;
    rdk_reactants.resize(reactants.size());
    if (reactant_indices.has_value()) {
        for (size_t i = 0; i < reactants.size(); ++i) {
            auto reactant_index = reactant_indices->at(i);
            rdk_reactants.at(reactant_index) = reactants.at(i)->rdkit_mol_ptr();
        }
    } else {
        for (const auto &m : reactants) {
            rdk_reactants.push_back(m->rdkit_mol_ptr());
        }
    }

    auto rdk_outcomes = rdkit_rxn_->runReactants(rdk_reactants);
    for (const auto &rdk_outcome : rdk_outcomes) {
        ReactionOutcome outcome;
        if (reactant_indices.has_value()) {
            outcome.reactant_indices = *reactant_indices;
        } else {
            outcome.reactant_indices.resize(reactants.size());
            std::iota(outcome.reactant_indices.begin(), outcome.reactant_indices.end(), 0);
        }

        for (const auto &rdk_prod : rdk_outcome) {
            try {
                auto prod = Molecule::from_unsanitized_rdkit(rdk_prod);
                outcome.products.push_back(std::move(prod));
            } catch (const MoleculeError &e) {
                if (ignore_errors) {
                    continue;
                } else {
                    throw ReactionError("Failed to sanitize product molecule: " +
                                        std::string(e.what()));
                }
            }
        }

        outcomes.push_back(std::move(outcome));
    }

    return outcomes;
}

} // namespace prexsyn
