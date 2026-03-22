#pragma once

namespace prexsyn::enumerator {

struct EnumeratorConfig {
    unsigned int max_building_blocks = 5;
    unsigned int heavy_atom_limit = 50;
    unsigned int selectivity_cutoff = 2;
    unsigned int max_outcomes_per_reaction = 8;
};

constexpr const static EnumeratorConfig kDefaultEnumeratorConfig{};

} // namespace prexsyn::enumerator
