#pragma once

namespace prexsyn::enumerator {

struct EnumeratorConfig {
    unsigned int max_building_blocks = 5;
    unsigned int heavy_atom_limit = 50;
};

constexpr const static EnumeratorConfig default_config = {
    .max_building_blocks = 5,
    .heavy_atom_limit = 50,
};

} // namespace prexsyn::enumerator
