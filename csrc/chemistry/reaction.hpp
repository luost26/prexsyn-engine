#pragma once

#include <cstddef>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <GraphMol/ChemReactions/Reaction.h>

#include "molecule.hpp"

namespace prexsyn {

class ReactionError : public std::runtime_error {
public:
    explicit ReactionError(const std::string &message) : std::runtime_error(message) {}
};

struct ReactionOutcome {
    std::vector<std::shared_ptr<Molecule>> products;
    std::vector<size_t> reactant_indices;

    bool empty() const { return products.empty(); }
    size_t num_products() const { return products.size(); }
    std::shared_ptr<Molecule> main_product() const;
};

class Reaction {
private:
    std::shared_ptr<RDKit::ChemicalReaction> rdkit_rxn_;

public:
    Reaction(std::shared_ptr<RDKit::ChemicalReaction> rdkit_rxn)
        : rdkit_rxn_(std::move(rdkit_rxn)) {
        if (!rdkit_rxn_) {
            throw ReactionError("RDKit reaction pointer is null");
        }
    }
    static std::unique_ptr<Reaction> from_smarts(const std::string &smarts);

    const RDKit::ChemicalReaction &rdkit_rxn() const { return *rdkit_rxn_; }
    RDKit::ChemicalReaction &rdkit_rxn() { return *rdkit_rxn_; }
    std::shared_ptr<RDKit::ChemicalReaction> rdkit_rxn_ptr() const { return rdkit_rxn_; }

    std::vector<ReactionOutcome>
    apply(const std::vector<std::shared_ptr<Molecule>> &reactants, bool ignore_errors = false,
          const std::optional<std::vector<size_t>> &reactant_indices = std::nullopt) const;
};

} // namespace prexsyn
