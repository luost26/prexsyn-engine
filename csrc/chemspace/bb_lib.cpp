#include "bb_lib.hpp"

#include <stdexcept>
#include <string>
#include <utility>

namespace prexsyn::chemspace {

const BuildingBlockItem &BuildingBlockLibrary::get(Index index) const {
    if (index >= building_blocks_.size()) {
        throw std::out_of_range("Building block index out of range");
    }
    return building_blocks_[index];
}

const BuildingBlockItem &BuildingBlockLibrary::get(const std::string &identifier) const {
    auto it = identifier_to_index_.find(identifier);
    if (it == identifier_to_index_.end()) {
        throw std::out_of_range("Building block identifier not found: " + identifier);
    }
    return building_blocks_[it->second];
}

BuildingBlockLibrary::Index BuildingBlockLibrary::add(const BuildingBlockEntry &entry) {
    if (identifier_to_index_.contains(entry.identifier)) {
        throw std::invalid_argument("Building block with the same identifier already exists: " +
                                    entry.identifier);
    }
    auto new_index = building_blocks_.size();
    building_blocks_.push_back(BuildingBlockItem{entry, new_index});
    identifier_to_index_[entry.identifier] = new_index;
    return new_index;
}

} // namespace prexsyn::chemspace
