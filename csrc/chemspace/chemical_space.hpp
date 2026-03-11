#pragma once

#include <cstddef>
#include <istream>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <utility>
#include <vector>

#include "../chemistry/chemistry.hpp"
#include "bb_lib.hpp"
#include "int_lib.hpp"
#include "rxn_lib.hpp"
#include "synthesis.hpp"

namespace prexsyn::chemspace {

struct ReactantMatchingConfig {
    size_t selectivity_cutoff = 2;

    template <typename Archive> void serialize(Archive &ar, const unsigned int /* version */) {
        ar & selectivity_cutoff;
    }

    bool check(const ReactionLibrary::Match &match) const {
        bool selectivity_ok = match.count <= selectivity_cutoff;
        return selectivity_ok;
    }
};

class ReactantLists {
public:
    using MolIndex = size_t;

private:
    // reaction -> reactant -> [building block]
    std::vector<std::vector<std::vector<MolIndex>>> r2b_;
    size_t num_matches_ = 0;
    friend class ChemicalSpace;

public:
    template <typename Archive> void serialize(Archive &ar, const unsigned int /* version */) {
        ar & r2b_;
        ar & num_matches_;
    }

    const auto &get(ReactionLibrary::Index i, Reaction::ReactantIndex j) const {
        return r2b_.at(i).at(j);
    }

    size_t num_matches() const { return num_matches_; }

    void init(const ReactionLibrary &);
    void set(MolIndex, ReactionLibrary::Index, Reaction::ReactantIndex);
};

class ChemicalSpace {
private:
    std::unique_ptr<BuildingBlockLibrary> bb_lib_;
    std::unique_ptr<ReactionLibrary> rxn_lib_;
    std::unique_ptr<IntermediateLibrary> int_lib_;

    ReactantMatchingConfig reactant_matching_config_;
    ReactantLists rnt_bb_mapping_, rnt_int_mapping_;

public:
    ChemicalSpace(std::unique_ptr<BuildingBlockLibrary> bb_lib,
                  std::unique_ptr<ReactionLibrary> rxn_lib,
                  std::unique_ptr<IntermediateLibrary> int_lib,
                  const ReactantMatchingConfig &matching_config = {})
        : bb_lib_(std::move(bb_lib)), rxn_lib_(std::move(rxn_lib)), int_lib_(std::move(int_lib)),
          reactant_matching_config_(matching_config) {
        if (bb_lib_ == nullptr || rxn_lib_ == nullptr || int_lib_ == nullptr) {
            throw std::invalid_argument("null pointer not allowed");
        }
        rnt_bb_mapping_.init(*rxn_lib_);
        rnt_int_mapping_.init(*rxn_lib_);
    }

    static std::unique_ptr<ChemicalSpace> deserialize(std::istream &);
    struct PeekStats {
        size_t num_reactions = 0;
        size_t num_building_blocks = 0;
        size_t num_intermediates = 0;
    };
    static PeekStats peek(std::istream &);
    void serialize(std::ostream &) const;

    const BuildingBlockLibrary &bb_lib() const { return *bb_lib_; }
    const ReactionLibrary &rxn_lib() const { return *rxn_lib_; }
    const IntermediateLibrary &int_lib() const { return *int_lib_; }
    const ReactantMatchingConfig &reactant_matching_config() const {
        return reactant_matching_config_;
    }
    const ReactantLists &building_block_reactant_lists() const { return rnt_bb_mapping_; }
    const ReactantLists &intermediate_reactant_lists() const { return rnt_int_mapping_; }

    void generate_intermediates();

    void build_reactant_lists_for_building_blocks();
    void build_reactant_lists_for_intermediates();
    void print_reactant_lists(std::ostream &) const;

    std::unique_ptr<ChemicalSpaceSynthesis> new_synthesis() const;
};

} // namespace prexsyn::chemspace
