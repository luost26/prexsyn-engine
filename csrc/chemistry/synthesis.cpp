#include "synthesis.hpp"

#include <algorithm>
#include <memory>
#include <utility>

#include "molecule.hpp"
#include "reaction.hpp"

namespace prexsyn::synthesis {

std::shared_ptr<Molecule> Outcome::main_product() const {
    auto max_it = std::ranges::max_element(products, [](const auto &a, const auto &b) {
        return a->num_heavy_atoms() < b->num_heavy_atoms();
    });
    if (max_it == products.end()) {
        throw SynthesisError("Outcome has no products");
    }
    return *max_it;
}

void Synthesis::add(const std::shared_ptr<Molecule> &mol) {
    Node node{.outcomes = {Outcome{.products = {mol}, .precursor_outcomes = {}}},
              .precursor_nodes = {}};
    nodes_.push_back(std::move(node));
}

} // namespace prexsyn::synthesis
