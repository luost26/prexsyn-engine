#pragma once

#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

#include "../chemistry/chemistry.hpp"
#include "building_block.hpp"
#include "reaction.hpp"

namespace prexsyn::chemspace {

struct ReactantMatchingConfig {
    size_t selectivity_cutoff = 1;

    bool check(const ReactionLibrary::Match &match) const {
        bool selectivity_ok = match.count <= selectivity_cutoff;
        return selectivity_ok;
    }
};

class ReactantBuildingBlockMapping {
private:
    // reaction -> reactant -> [building block]
    std::vector<std::vector<std::vector<BuildingBlockLibrary::Index>>> r2b_;
    friend class ChemicalSpace;

public:
    const auto &get_building_blocks(ReactionLibrary::Index i, Reaction::ReactantIndex j) const {
        return r2b_.at(i).at(j);
    }

    void clear() { r2b_.clear(); }

    void set(BuildingBlockLibrary::Index, ReactionLibrary::Index, Reaction::ReactantIndex);
};

class ChemicalSpace {
private:
    std::unique_ptr<BuildingBlockLibrary> bb_lib_;
    std::unique_ptr<ReactionLibrary> rxn_lib_;

    ReactantMatchingConfig reactant_matching_config_;
    ReactantBuildingBlockMapping rnt_bb_mapping_;

public:
    ChemicalSpace(std::unique_ptr<BuildingBlockLibrary> bb_lib,
                  std::unique_ptr<ReactionLibrary> rxn_lib)
        : bb_lib_(std::move(bb_lib)), rxn_lib_(std::move(rxn_lib)) {}

    const BuildingBlockLibrary &bb_lib() const { return *bb_lib_; }
    const ReactionLibrary &rxn_lib() const { return *rxn_lib_; }
    const ReactantMatchingConfig &reactant_matching_config() const {
        return reactant_matching_config_;
    }

    void set_reactant_matching_config(const ReactantMatchingConfig &config) {
        reactant_matching_config_ = config;
        rnt_bb_mapping_.clear();
    }

    void build_reactant_building_block_mapping();
};

} // namespace prexsyn::chemspace
