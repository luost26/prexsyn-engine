#pragma once

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "molecule.hpp"
#include "reaction.hpp"

namespace prexsyn {

class SynthesisError : public std::runtime_error {
public:
    explicit SynthesisError(const std::string &message) : std::runtime_error(message) {}
};

class SynthesisNode {
private:
    struct Item {
        std::shared_ptr<Molecule> molecule;
        std::vector<std::string> reactant_names;
        std::vector<size_t> precursor_item_indices;
    };

    std::vector<Item> items_;
    std::shared_ptr<Reaction> reaction_;
    std::vector<std::shared_ptr<SynthesisNode>> precursor_nodes_;

    SynthesisNode() = default;

public:
    static std::unique_ptr<SynthesisNode> from_molecule(const std::shared_ptr<Molecule> &);
    static std::unique_ptr<SynthesisNode>
    from_reaction(const std::shared_ptr<Reaction> &,
                  const std::vector<std::shared_ptr<SynthesisNode>> &);

    void add_reaction_outcome(const ReactionOutcomeWithReactantAssignment &outcome,
                              const std::vector<size_t> &precursor_item_indices);

    size_t size() const { return items_.size(); }
    const auto &precursor_nodes() const { return precursor_nodes_; }
    std::shared_ptr<Molecule> at(size_t i) const { return items_.at(i).molecule; }

    struct PrecursorMolecule {
        size_t precursor_index;
        std::string reactant_name;
        std::shared_ptr<SynthesisNode> precursor_node;
        size_t item_index;
        std::shared_ptr<Molecule> molecule;
    };
    std::vector<PrecursorMolecule> precursors(size_t index) const;
};

class Synthesis {
private:
    std::vector<std::shared_ptr<SynthesisNode>> nodes_;
    std::vector<std::shared_ptr<SynthesisNode>> stack_;

public:
    const std::vector<std::shared_ptr<SynthesisNode>> &nodes() const { return nodes_; }
    size_t stack_size() const { return stack_.size(); }
    std::shared_ptr<SynthesisNode> stack_top(size_t i = 0) const;

    void push(const std::shared_ptr<Molecule> &);
    void push(const std::shared_ptr<Reaction> &);
    void undo();
};

} // namespace prexsyn
