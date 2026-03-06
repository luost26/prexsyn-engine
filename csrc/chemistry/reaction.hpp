#pragma once

#include <cstddef>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
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

    bool empty() const { return products.empty(); }
    size_t num_products() const { return products.size(); }
    std::shared_ptr<Molecule> main_product() const;
};

struct ReactionOutcomeWithReactantAssignment : public ReactionOutcome {
    std::vector<std::string> reactant_names;
};

class Reaction {
public:
    using ReactantIndex = size_t;

private:
    std::shared_ptr<RDKit::ChemicalReaction> rdkit_rxn_;
    std::vector<std::string> reactant_names_;
    std::map<std::string, ReactantIndex> reactant_name_to_index_;

public:
    Reaction(std::shared_ptr<RDKit::ChemicalReaction> rdkit_rxn,
             const std::vector<std::string> &reactant_names)
        : rdkit_rxn_(std::move(rdkit_rxn)), reactant_names_(reactant_names) {

        if (!rdkit_rxn_) {
            throw ReactionError("RDKit reaction pointer is null");
        }

        if (rdkit_rxn_->getNumReactantTemplates() != reactant_names.size()) {
            throw ReactionError(
                "Number of reactant names does not match number of reactant templates");
        }

        for (const auto &name : reactant_names) {
            if (reactant_name_to_index_.contains(name)) {
                throw ReactionError("Duplicate reactant name: " + name);
            }
            size_t index = reactant_name_to_index_.size();
            reactant_name_to_index_[name] = index;
        }
    }
    static std::unique_ptr<Reaction> from_smarts(const std::string &smarts,
                                                 const std::vector<std::string> &reactant_names);

    const RDKit::ChemicalReaction &rdkit_rxn() const { return *rdkit_rxn_; }
    RDKit::ChemicalReaction &rdkit_rxn() { return *rdkit_rxn_; }
    std::shared_ptr<RDKit::ChemicalReaction> rdkit_rxn_ptr() const { return rdkit_rxn_; }

    size_t num_reactants() const { return reactant_names_.size(); }
    const std::map<std::string, ReactantIndex> &reactant_name_to_index() const {
        return reactant_name_to_index_;
    }

    struct ReactantMatch {
        ReactantIndex index;
        std::string_view name;
        size_t count;
    };
    std::vector<ReactantMatch> match_reactants(const Molecule &) const;

    std::vector<ReactionOutcome>
    apply(const std::map<std::string, std::shared_ptr<Molecule>> &reactants,
          bool ignore_errors = false) const;

    std::vector<ReactionOutcomeWithReactantAssignment>
    apply(const std::vector<std::shared_ptr<Molecule>> &reactants,
          bool ignore_errors = false) const;
};

} // namespace prexsyn
