#include "chemical_space.hpp"

#include <algorithm>
#include <cstddef>
#include <exception>
#include <istream>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <omp.h>

#include "../chemistry/chemistry.hpp"
#include "../utility/logging.hpp"
#include "../utility/serialization.hpp"
#include "bb_lib.hpp"
#include "int_lib.hpp"
#include "postfix_notation.hpp"
#include "rxn_lib.hpp"
#include "synthesis.hpp"

namespace prexsyn::chemspace {

void ReactantLists::init(const ReactionLibrary &rxn_lib) {
    r2b_.clear();
    r2b_.resize(rxn_lib.size());
    for (size_t rxn_idx = 0; rxn_idx < rxn_lib.size(); ++rxn_idx) {
        const auto &rxn_item = rxn_lib.get(rxn_idx);
        size_t num_reactants = rxn_item.reaction->num_reactants();
        r2b_[rxn_idx].resize(num_reactants);
    }
}

void ReactantLists::set(MolIndex bb, ReactionLibrary::Index rxn, Reaction::ReactantIndex rnt) {
    if (rxn >= r2b_.size()) {
        throw std::out_of_range("Reaction index out of range. Did you forget to call init?");
    }
    if (rnt >= r2b_[rxn].size()) {
        throw std::out_of_range("Reactant index out of range. Did you forget to call init?");
    }
    r2b_[rxn][rnt].push_back(bb);
    num_matches_++;
}

std::unique_ptr<ChemicalSpace> ChemicalSpace::deserialize(std::istream &is) {
    {
        boost::archive::binary_iarchive ia(is);
        size_t bb_lib_size = 0, rxn_lib_size = 0, int_lib_size = 0;
        ia >> bb_lib_size >> rxn_lib_size >> int_lib_size;
    }
    auto bb_lib = BuildingBlockLibrary::deserialize(is);
    auto rxn_lib = ReactionLibrary::deserialize(is);
    auto int_lib = IntermediateLibrary::deserialize(is);
    {
        boost::archive::binary_iarchive ia(is);
        ReactantMatchingConfig matching_config;
        ia >> matching_config;
        auto chemspace = std::make_unique<ChemicalSpace>(std::move(bb_lib), std::move(rxn_lib),
                                                         std::move(int_lib), matching_config);
        ia >> chemspace->rnt_bb_mapping_;
        ia >> chemspace->rnt_int_mapping_;
        return chemspace;
    }
}

ChemicalSpace::PeekStats ChemicalSpace::peek(std::istream &is) {
    boost::archive::binary_iarchive ia(is);
    PeekStats stats;
    ia >> stats.num_building_blocks >> stats.num_reactions >> stats.num_intermediates;
    return stats;
}

void ChemicalSpace::serialize(std::ostream &os) const {
    {
        // For peeking
        boost::archive::binary_oarchive oa(os);
        oa << bb_lib_->size() << rxn_lib_->size() << int_lib_->size();
    }
    bb_lib_->serialize(os);
    rxn_lib_->serialize(os);
    int_lib_->serialize(os);
    {
        boost::archive::binary_oarchive oa(os);
        oa << reactant_matching_config_;
        oa << rnt_bb_mapping_;
        oa << rnt_int_mapping_;
    }
}

void ChemicalSpace::generate_intermediates() {
    if (rnt_bb_mapping_.num_matches() == 0) {
        logger()->warn(
            "No building block reactant matches found. Please build reactant lists first.");
        return;
    }

    logger()->info("Starting to generate intermediates...");

    int_lib_->clear();
    std::vector<std::pair<BuildingBlockLibrary::Index, ReactionLibrary::Index>> bb_rxn_pairs;

    for (size_t rxn_idx = 0; rxn_idx < rnt_bb_mapping_.r2b_.size(); ++rxn_idx) {
        const auto &reactant_lists = rnt_bb_mapping_.r2b_[rxn_idx];
        if (reactant_lists.size() != 1) {
            continue;
        }
        const auto &bb_indices = reactant_lists[0];
        for (const auto &bb_idx : bb_indices) {
            bb_rxn_pairs.emplace_back(bb_idx, rxn_idx);
        }
    }

#pragma omp parallel for schedule(dynamic)
    for (const auto &[bb_idx, rxn_idx] : bb_rxn_pairs) {
        const auto &bb_item = bb_lib_->get(bb_idx);
        const auto &rxn_item = rxn_lib_->get(rxn_idx);
        PostfixNotation pfn{};
        pfn.append(bb_idx, PostfixNotation::Token::BuildingBlock);
        pfn.append(rxn_idx, PostfixNotation::Token::Reaction);

        try {
            auto outcomes = rxn_item.reaction->apply(std::vector{bb_item.molecule}, true);
            for (size_t i = 0; i < outcomes.size(); ++i) {
                std::string identifier =
                    bb_item.identifier + "@" + rxn_item.name + ":" + std::to_string(i);
#pragma omp critical
                {
                    int_lib_->add({
                        .postfix_notation = pfn,
                        .molecule = outcomes.at(i).main_product(),
                        .identifier = identifier,
                    });
                    if (int_lib_->size() % 10000 == 0) {
                        logger()->info("Generated {} intermediates...", int_lib_->size());
                    }
                }
            }
        } catch (const std::exception &e) {
            logger()->warn(
                "Error generating intermediate for building block {} and reaction {}: {}",
                bb_item.identifier, rxn_item.name, e.what());
            continue;
        }
    }

    logger()->info("Done. Intermediates: {}", int_lib_->size());
}

void ChemicalSpace::build_reactant_lists_for_building_blocks() {
    logger()->info("Starting to build reactant-building block lists...");
    rnt_bb_mapping_.init(*rxn_lib_);
#pragma omp parallel for schedule(dynamic)
    for (const auto &bb : *bb_lib_) {
        auto matches = rxn_lib_->match_reactants(*bb.molecule);
        for (const auto &match : matches) {
            if (!reactant_matching_config_.check(match)) {
                continue;
            }
#pragma omp critical
            {
                rnt_bb_mapping_.set(bb.index, match.reaction_index, match.reactant_index);
            }
        }
    }
    logger()->info("Done. Reactant-building block matches: {}", rnt_bb_mapping_.num_matches());
}

void ChemicalSpace::build_reactant_lists_for_intermediates() {
    logger()->info("Starting to build reactant-intermediate lists...");
    rnt_int_mapping_.init(*rxn_lib_);
#pragma omp parallel for schedule(dynamic)
    for (const auto &intm : *int_lib_) {
        auto matches = rxn_lib_->match_reactants(*intm.molecule);
        for (const auto &match : matches) {
            if (!reactant_matching_config_.check(match)) {
                continue;
            }
#pragma omp critical
            {
                rnt_int_mapping_.set(intm.index, match.reaction_index, match.reactant_index);
            }
        }
    }
    logger()->info("Done. Reactant-intermediate matches: {}", rnt_int_mapping_.num_matches());
}

void ChemicalSpace::print_reactant_lists(std::ostream &os) const {
    for (size_t rxn_idx = 0; rxn_idx < rnt_bb_mapping_.r2b_.size(); ++rxn_idx) {
        const auto &rxn_item = rxn_lib_->get(rxn_idx);
        os << "- " << rxn_item.name << " (index=" << rxn_idx << "):\n";
        for (size_t rnt_idx = 0; rnt_idx < rnt_bb_mapping_.r2b_[rxn_idx].size(); ++rnt_idx) {
            const auto &bb_indices = rnt_bb_mapping_.r2b_[rxn_idx][rnt_idx];
            const auto &int_indices = rnt_int_mapping_.r2b_[rxn_idx][rnt_idx];
            os << "    " << rxn_item.reaction->reactant_names()[rnt_idx] << ": "
               << bb_indices.size() << " building blocks, " << int_indices.size()
               << " intermediates ";

            os << " # [";
            for (size_t i = 0; i < std::min(bb_indices.size(), size_t(5)); ++i) {
                const auto &bb_item = bb_lib_->get(bb_indices[i]);
                os << bb_item.identifier << ", ";
            }
            if (bb_indices.size() > 5) {
                os << "...";
            }
            os << "]\n";
        }
    }
}

std::unique_ptr<ChemicalSpaceSynthesis> ChemicalSpace::new_synthesis() const {
    return std::make_unique<ChemicalSpaceSynthesis>(*this);
}

} // namespace prexsyn::chemspace
