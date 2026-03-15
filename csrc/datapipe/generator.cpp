#include "generator.hpp"

#include <cstddef>
#include <cstdlib>
#include <memory>
#include <optional>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

#include "../chemistry/chemistry.hpp"
#include "../chemspace/chemspace.hpp"
#include "../utility/random.hpp"

namespace prexsyn::datapipe {

Generator::Generator(std::shared_ptr<chemspace::ChemicalSpace> cs, const Config &config,
                     size_t random_seed)
    : cs_(std::move(cs)), config_(config), synthesis_(nullptr), rng_(random_seed) {
    if (cs_ == nullptr) {
        throw std::invalid_argument("null pointer for chemical space");
    }
    if (cs_->bb_lib().size() == 0) {
        throw std::invalid_argument("empty building block library");
    }
    if (cs_->rxn_lib().size() == 0) {
        throw std::invalid_argument("empty reaction library");
    }
}

bool Generator::is_limit_exceeded() const {
    if (synthesis_ == nullptr) {
        return false;
    }
    auto products = synthesis_->products();
    if (products.empty()) {
        return false;
    }
    const auto &product = products.front();
    return product->num_heavy_atoms() > config_.heavy_atom_limit ||
           synthesis_->count_building_blocks() > config_.max_building_blocks;
}

void Generator::clear_synthesis() { synthesis_.reset(); }

void Generator::init_synthesis() {
    synthesis_ = std::make_unique<chemspace::Synthesis>(*cs_);

    auto num_bb = cs_->bb_lib().size();
    std::uniform_int_distribution<size_t> dist(0, num_bb - 1);
    size_t index = dist(rng_);
    synthesis_->add_building_block(index);
}

void Generator::grow_synthesis() {
    auto products = synthesis_->products();
    if (products.empty()) {
        clear_synthesis();
        return;
    }
    auto product = random_choice(products, rng_);

    auto matches = cs_->rxn_lib().match_reactants(*product);
    if (matches.empty()) {
        clear_synthesis();
        return;
    }
    auto match = random_choice(matches, rng_);

    auto rxn = cs_->rxn_lib().get(match.reaction_index);
    chemspace::Synthesis::Result result;
    for (size_t i = 0; i < rxn.reaction->num_reactants(); ++i) {
        if (i == match.reactant_index) {
            continue;
        }
        const auto &rlist_bb = cs_->building_block_reactant_lists().get(match.reaction_index, i);
        const auto &rlist_int = cs_->intermediate_reactant_lists().get(match.reaction_index, i);
        const auto &[choice, index] = random_choice(rlist_bb, rlist_int, rng_);
        if (choice == which_vector::first) {
            result = synthesis_->add_building_block(index);
        } else {
            result = synthesis_->add_intermediate(index);
        }
        if (!result) {
            clear_synthesis();
            return;
        }
    }

    result = synthesis_->add_reaction(match.reaction_index);
    if (!result) {
        clear_synthesis();
        return;
    }
}

std::optional<std::shared_ptr<chemspace::Synthesis>> Generator::try_next() {
    if (synthesis_ == nullptr || is_limit_exceeded()) {
        init_synthesis();
    } else {
        grow_synthesis();
    }

    if (synthesis_ == nullptr) {
        return std::nullopt;
    } else {
        return std::make_shared<chemspace::Synthesis>(*synthesis_);
    }
}

std::shared_ptr<chemspace::Synthesis> Generator::next() {
    constexpr const int max_attempts = 100;
    for (int attempt = 0; attempt < max_attempts; ++attempt) {
        auto result = try_next();
        if (result.has_value()) {
            return result.value();
        }
    }
    throw std::runtime_error("too many failed attempts to find next synthesis");
}

std::pair<std::shared_ptr<chemspace::Synthesis>, std::shared_ptr<Molecule>>
Generator::next_with_product() {
    auto syn = next();
    auto product = random_choice(syn->products(), rng_);
    return {syn, product};
}

} // namespace prexsyn::datapipe
