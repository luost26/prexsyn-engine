#include "bb_lib.hpp"

#include <cstddef>
#include <istream>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>

#include "../utility/serialization.hpp"

namespace prexsyn::chemspace {

std::unique_ptr<BuildingBlockLibrary> BuildingBlockLibrary::deserialize(std::istream &data) {
    boost::archive::binary_iarchive ia(data);
    size_t num_items = 0;
    ia >> num_items;
    auto bb_lib = std::make_unique<BuildingBlockLibrary>();
    bb_lib->building_blocks_.reserve(num_items);
    for (size_t i = 0; i < num_items; ++i) {
        BuildingBlockItem item;
        ia >> item;
        bb_lib->building_blocks_.push_back(std::move(item));
    }
    ia >> bb_lib->identifier_to_index_;
    return bb_lib;
}

void BuildingBlockLibrary::serialize(std::ostream &stream) const {
    boost::archive::binary_oarchive oa(stream);
    oa << building_blocks_.size();
    for (const auto &item : building_blocks_) {
        oa << item;
    }
    oa << identifier_to_index_;
}

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
        throw BuildingBlockLibraryError("duplicate identifier: " + entry.identifier);
    }
    auto new_index = building_blocks_.size();
    building_blocks_.push_back(BuildingBlockItem{entry, new_index});
    identifier_to_index_[entry.identifier] = new_index;
    return new_index;
}

} // namespace prexsyn::chemspace
