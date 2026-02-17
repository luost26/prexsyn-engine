#pragma once

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "molecule.hpp"
#include "reaction.hpp"

namespace prexsyn::synthesis {

class SynthesisError : public std::runtime_error {
public:
    explicit SynthesisError(const std::string &message) : std::runtime_error(message) {}
};

struct Outcome {
    std::vector<std::shared_ptr<Molecule>> products;
    std::vector<size_t> precursor_outcomes;

    bool empty() const { return products.empty(); }
    size_t num_products() const { return products.size(); }
    std::shared_ptr<Molecule> main_product() const;
};

struct Node {
    std::vector<Outcome> outcomes;
    std::vector<size_t> precursor_nodes;
};

class Synthesis {
    std::vector<Node> nodes_;
    std::vector<size_t> stack_;

public:
    void add(const std::shared_ptr<Molecule> &);
    void add(const std::shared_ptr<Reaction> &);
};

} // namespace prexsyn::synthesis
