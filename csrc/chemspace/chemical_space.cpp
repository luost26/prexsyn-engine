#include "chemical_space.hpp"

#include <omp.h>

#include "../chemistry/chemistry.hpp"
#include "bb_lib.hpp"
#include "rxn_lib.hpp"

namespace prexsyn::chemspace {

void ReactantBuildingBlockMapping::set(BuildingBlockLibrary::Index bb, ReactionLibrary::Index rxn,
                                       Reaction::ReactantIndex rnt) {
    if (rxn >= r2b_.size()) {
        r2b_.resize(rxn + 1);
    }
    if (rnt >= r2b_[rxn].size()) {
        r2b_[rxn].resize(rnt + 1);
    }
    r2b_[rxn][rnt].push_back(bb);
}

void ChemicalSpace::build_reactant_building_block_mapping() {
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
}

} // namespace prexsyn::chemspace
