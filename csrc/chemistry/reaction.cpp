#include "reaction.hpp"

#include <algorithm>
#include <cstddef>
#include <map>
#include <memory>
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

std::unique_ptr<Reaction> Reaction::from_smarts(const std::string &smarts,
                                                const std::vector<std::string> &reactant_names) {
    std::shared_ptr<RDKit::ChemicalReaction> rdkit_rxn(RDKit::RxnSmartsToChemicalReaction(smarts));
    if (!rdkit_rxn) {
        throw ReactionError("Failed to parse SMARTS: " + smarts);
    }
    rdkit_rxn->initReactantMatchers();
    if (!rdkit_rxn->isInitialized()) {
        throw ReactionError("RDKit reaction is not initialized after parsing SMARTS: " + smarts);
    }

    return std::make_unique<Reaction>(std::move(rdkit_rxn), reactant_names);
}

std::vector<Reaction::ReactantMatch> Reaction::match_reactants(const Molecule &molecule) const {
    std::vector<Reaction::ReactantMatch> matches;
    for (size_t i = 0; i < num_reactants(); ++i) {
        const auto &tmpl = rdkit_rxn_->getReactants()[i];
        auto res = RDKit::SubstructMatch(molecule.rdkit_mol(), *tmpl);
        if (!res.empty()) {
            matches.push_back({
                .index = i,
                .name = reactant_names_.at(i),
                .count = res.size(),
            });
        }
    }
    return matches;
}

std::vector<ReactionOutcome>
Reaction::apply(const std::map<std::string, std::shared_ptr<Molecule>> &reactants,
                bool ignore_errors) const {
    if (reactants.size() != reactant_names_.size()) {
        throw ReactionError(
            "Number of reactants provided does not match number of reactant templates");
    }

    std::vector<ReactionOutcome> outcomes;

    std::vector<RDKit::ROMOL_SPTR> rdk_reactants;
    rdk_reactants.resize(reactants.size());
    for (const auto &[name, mol] : reactants) {
        if (!reactant_name_to_index_.contains(name)) {
            throw ReactionError("Unknown reactant name: " + name);
        }
        size_t index = reactant_name_to_index_.at(name);
        rdk_reactants[index] = mol->rdkit_mol_ptr();
    }

    auto rdk_outcomes = rdkit_rxn_->runReactants(rdk_reactants);
    for (const auto &rdk_outcome : rdk_outcomes) {
        ReactionOutcome outcome;
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

std::vector<ReactionOutcomeWithReactantAssignment>
Reaction::apply(const std::vector<std::shared_ptr<Molecule>> &reactants, bool ignore_errors) const {
    if (reactants.size() != reactant_names_.size()) {
        throw ReactionError(
            "Number of reactants provided does not match number of reactant templates");
    }

    std::vector<std::string> name_perm = reactant_names_;
    std::sort(name_perm.begin(), name_perm.end());

    std::vector<ReactionOutcomeWithReactantAssignment> results;
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-do-while)
    do {
        std::map<std::string, std::shared_ptr<Molecule>> reactant_map;
        for (size_t i = 0; i < name_perm.size(); ++i) {
            reactant_map[name_perm.at(i)] = reactants.at(i);
        }
        auto outcomes = apply(reactant_map, ignore_errors);
        for (auto &outcome : outcomes) {
            ReactionOutcomeWithReactantAssignment outcome_with_assignment;
            outcome_with_assignment.products = std::move(outcome.products);
            outcome_with_assignment.reactant_names = name_perm;
            results.push_back(std::move(outcome_with_assignment));
        }
    } while (std::next_permutation(name_perm.begin(), name_perm.end()));

    return results;
}

} // namespace prexsyn
