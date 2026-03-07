#pragma once

#include <cstddef>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
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
             const std::vector<std::string> &reactant_names);
    static std::unique_ptr<Reaction> from_smarts(const std::string &smarts,
                                                 const std::vector<std::string> &reactant_names);
    static std::unique_ptr<Reaction> from_smarts(const std::string &smarts);
    static std::unique_ptr<Reaction> from_rdkit_pickle(const std::string &,
                                                       const std::vector<std::string> &);

    static std::unique_ptr<Reaction> deserialize(const std::string &);
    std::string serialize() const;

    const RDKit::ChemicalReaction &rdkit_rxn() const { return *rdkit_rxn_; }
    RDKit::ChemicalReaction &rdkit_rxn() { return *rdkit_rxn_; }
    std::shared_ptr<RDKit::ChemicalReaction> rdkit_rxn_ptr() const { return rdkit_rxn_; }
    std::string rdkit_pickle() const;

    size_t num_reactants() const { return reactant_names_.size(); }
    const std::vector<std::string> &reactant_names() const { return reactant_names_; }
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
