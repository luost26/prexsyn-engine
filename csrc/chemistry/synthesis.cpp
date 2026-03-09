#include "synthesis.hpp"

#include <cstddef>
#include <functional>
#include <memory>
#include <ranges>
#include <string>
#include <utility>
#include <vector>

#include "molecule.hpp"
#include "reaction.hpp"

namespace prexsyn {

std::unique_ptr<SynthesisNode> SynthesisNode::from_molecule(const std::shared_ptr<Molecule> &mol) {
    std::unique_ptr<SynthesisNode> node{new SynthesisNode()};
    node->items_.push_back({mol, {}, {}});
    return node;
}

std::unique_ptr<SynthesisNode>
SynthesisNode::from_reaction(const std::shared_ptr<Reaction> &rxn,
                             const std::vector<std::shared_ptr<SynthesisNode>> &prec) {
    std::unique_ptr<SynthesisNode> node{new SynthesisNode()};
    node->reaction_ = rxn;
    node->precursor_nodes_ = prec;
    return node;
}

void SynthesisNode::add_reaction_outcome(const ReactionOutcomeWithReactantAssignment &outcome,
                                         const std::vector<size_t> &precursor_item_indices) {
    items_.push_back({outcome.main_product(), outcome.reactant_names, precursor_item_indices});
}

std::vector<SynthesisNode::PrecursorMolecule> SynthesisNode::precursors(size_t index) const {
    const auto &item = items_.at(index);
    std::vector<PrecursorMolecule> result;

    for (size_t i = 0; i < precursor_nodes_.size(); ++i) {
        const auto &reactant_name = item.reactant_names.at(i);
        const auto &pre_item_index = item.precursor_item_indices.at(i);
        const auto &node = precursor_nodes_.at(i);
        result.push_back({.precursor_index = i,
                          .reactant_name = reactant_name,
                          .precursor_node = node,
                          .item_index = pre_item_index,
                          .molecule = node->at(pre_item_index)});
    }

    return result;
}

std::shared_ptr<SynthesisNode> Synthesis::stack_top(size_t i) const {
    if (i >= stack_.size()) {
        throw SynthesisError("Stack index out of range");
    }
    return stack_.at(stack_.size() - 1 - i);
}

void Synthesis::push(const std::shared_ptr<Molecule> &molecule) {
    nodes_.push_back(SynthesisNode::from_molecule(molecule));
    stack_.push_back(nodes_.back());
}

static void cartesian_product(const std::vector<size_t> &sizes,
                              const std::function<void(const std::vector<size_t> &)> &callback) {
    std::vector<size_t> indices(sizes.size(), 0);
    while (true) {
        callback(indices);
        size_t i = 0;
        while (i < sizes.size()) {
            indices[i]++;
            if (indices[i] < sizes[i]) {
                break;
            }
            indices[i] = 0;
            i++;
        }
        if (i == sizes.size()) {
            break;
        }
    }
};

void Synthesis::push(const std::shared_ptr<Reaction> &reaction) {
    if (stack_.size() < reaction->num_reactants()) {
        throw SynthesisError("Not enough reactants on the stack for the reaction, got " +
                             std::to_string(stack_.size()) + " but need " +
                             std::to_string(reaction->num_reactants()));
    }

    std::vector<std::shared_ptr<SynthesisNode>> precursor_nodes;
    precursor_nodes.reserve(reaction->num_reactants());
    std::vector<size_t> sizes;
    auto stack_iter = stack_.end();
    for (size_t i = 0; i < reaction->num_reactants(); ++i) {
        const auto &node = *(--stack_iter);
        precursor_nodes.push_back(node);
        sizes.push_back(node->size());
    }

    auto new_node = SynthesisNode::from_reaction(reaction, precursor_nodes);
    cartesian_product(sizes, [&](const std::vector<size_t> &item_indices) -> void {
        std::vector<std::shared_ptr<Molecule>> reactants;
        reactants.reserve(precursor_nodes.size());
        for (size_t i = 0; i < precursor_nodes.size(); ++i) {
            reactants.push_back(precursor_nodes.at(i)->at(item_indices.at(i)));
        }
        auto outcomes = reaction->apply(reactants, /*ignore_errors=*/true);
        for (const auto &outcome : outcomes) {
            new_node->add_reaction_outcome(outcome, item_indices);
        }
    });

    if (new_node->size() == 0) {
        throw SynthesisError("The reaction did not produce any products.");
    }

    nodes_.push_back(std::move(new_node));
    stack_.erase(stack_iter, stack_.end());
    stack_.push_back(nodes_.back());
}

void Synthesis::undo() {
    if (stack_.empty()) {
        throw SynthesisError("Cannot undo because the stack is empty");
    }
    auto top_node = stack_.back();
    stack_.pop_back();
    nodes_.pop_back();
    for (const auto &it : std::ranges::reverse_view(top_node->precursor_nodes())) {
        stack_.push_back(it);
    }
}

} // namespace prexsyn
