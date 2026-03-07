#include "chemical_space.hpp"

#include <algorithm>
#include <boost/archive/binary_oarchive.hpp>
#include <cstddef>
#include <istream>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <utility>

#include <omp.h>

#include "../chemistry/chemistry.hpp"
#include "../utility/logging.hpp"
#include "../utility/serialization.hpp"
#include "bb_lib.hpp"
#include "rxn_lib.hpp"

namespace prexsyn::chemspace {

void ReactantBuildingBlockMapping::init(const ReactionLibrary &rxn_lib) {
    r2b_.clear();
    r2b_.resize(rxn_lib.size());
    for (size_t rxn_idx = 0; rxn_idx < rxn_lib.size(); ++rxn_idx) {
        const auto &rxn_item = rxn_lib.get(rxn_idx);
        size_t num_reactants = rxn_item.reaction->num_reactants();
        r2b_[rxn_idx].resize(num_reactants);
    }
}

void ReactantBuildingBlockMapping::set(BuildingBlockLibrary::Index bb, ReactionLibrary::Index rxn,
                                       Reaction::ReactantIndex rnt) {
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
        size_t bb_lib_size = 0, rxn_lib_size = 0;
        ia >> bb_lib_size >> rxn_lib_size;
    }
    auto bb_lib = BuildingBlockLibrary::deserialize(is);
    auto rxn_lib = ReactionLibrary::deserialize(is);
    {
        boost::archive::binary_iarchive ia(is);
        ReactantMatchingConfig matching_config;
        ia >> matching_config;
        auto chemspace =
            std::make_unique<ChemicalSpace>(std::move(bb_lib), std::move(rxn_lib), matching_config);
        ia >> chemspace->rnt_bb_mapping_;
        return chemspace;
    }
}

void ChemicalSpace::serialize(std::ostream &os) const {
    {
        // For peeking
        boost::archive::binary_oarchive oa(os);
        oa << bb_lib_->size() << rxn_lib_->size();
    }
    bb_lib_->serialize(os);
    rxn_lib_->serialize(os);
    {
        boost::archive::binary_oarchive oa(os);
        oa << reactant_matching_config_;
        oa << rnt_bb_mapping_;
    }
}

void ChemicalSpace::build_reactant_building_block_mapping() {
    logger()->info("Starting to build reactant-building block mapping...");
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
    logger()->info("Done. Found {} matches.", rnt_bb_mapping_.num_matches());
}

void ChemicalSpace::print_reactant_building_block_mapping(std::ostream &os) const {
    for (size_t rxn_idx = 0; rxn_idx < rnt_bb_mapping_.r2b_.size(); ++rxn_idx) {
        const auto &rxn_item = rxn_lib_->get(rxn_idx);
        os << "- " << rxn_item.name << " (#" << rxn_idx << "):\n";
        for (size_t rnt_idx = 0; rnt_idx < rnt_bb_mapping_.r2b_[rxn_idx].size(); ++rnt_idx) {
            const auto &bb_indices = rnt_bb_mapping_.r2b_[rxn_idx][rnt_idx];
            os << "  - " << rxn_item.reaction->reactant_names()[rnt_idx] << ": "
               << bb_indices.size() << " ";

            os << "[";
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

} // namespace prexsyn::chemspace
