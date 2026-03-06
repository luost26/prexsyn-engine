#pragma once

#include <cstddef>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "../chemistry/chemistry.hpp"

namespace prexsyn::chemspace {

struct BuildingBlockEntry {
    std::shared_ptr<Molecule> molecule;
    std::string identifier;
    std::set<std::string> classifications;
};

struct BuildingBlockItem : BuildingBlockEntry {
    size_t index{};
};

class BuildingBlockLibrary {
public:
    using Index = size_t;

private:
    std::vector<BuildingBlockItem> building_blocks_;
    std::map<std::string, Index> identifier_to_index_;

    BuildingBlockLibrary() = default;

public:
    size_t size() const { return building_blocks_.size(); }
    const BuildingBlockItem &get(Index) const;
    const BuildingBlockItem &get(const std::string &) const;
    Index add(const BuildingBlockEntry &);

    auto begin() const noexcept { return building_blocks_.begin(); }
    auto end() const noexcept { return building_blocks_.end(); }
};

} // namespace prexsyn::chemspace
