#pragma once

#include <cstddef>
#include <istream>
#include <memory>
#include <ostream>
#include <utility>
#include <vector>

#include "../chemistry/chemistry.hpp"
#include "bb_lib.hpp"
#include "rxn_lib.hpp"

namespace prexsyn::chemspace {

struct ReactantMatchingConfig {
    size_t selectivity_cutoff = 1;

    template <typename Archive> void serialize(Archive &ar, const unsigned int /* version */) {
        ar & selectivity_cutoff;
    }

    bool check(const ReactionLibrary::Match &match) const {
        bool selectivity_ok = match.count <= selectivity_cutoff;
        return selectivity_ok;
    }
};

class ReactantBuildingBlockMapping {
private:
    // reaction -> reactant -> [building block]
    std::vector<std::vector<std::vector<BuildingBlockLibrary::Index>>> r2b_;
    size_t num_matches_ = 0;
    friend class ChemicalSpace;

public:
    template <typename Archive> void serialize(Archive &ar, const unsigned int /* version */) {
        ar & r2b_;
        ar & num_matches_;
    }

    const auto &get_building_blocks(ReactionLibrary::Index i, Reaction::ReactantIndex j) const {
        return r2b_.at(i).at(j);
    }

    size_t num_matches() const { return num_matches_; }

    void init(const ReactionLibrary &);
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
                  std::unique_ptr<ReactionLibrary> rxn_lib,
                  const ReactantMatchingConfig &matching_config = {})
        : bb_lib_(std::move(bb_lib)), rxn_lib_(std::move(rxn_lib)),
          reactant_matching_config_(matching_config) {}

    static std::unique_ptr<ChemicalSpace> deserialize(std::istream &);
    void serialize(std::ostream &) const;

    const BuildingBlockLibrary &bb_lib() const { return *bb_lib_; }
    const ReactionLibrary &rxn_lib() const { return *rxn_lib_; }
    const ReactantMatchingConfig &reactant_matching_config() const {
        return reactant_matching_config_;
    }

    void build_reactant_building_block_mapping();
    void print_reactant_building_block_mapping(std::ostream &) const;
};

} // namespace prexsyn::chemspace
