#pragma once

#include <cstddef>
#include <istream>
#include <map>
#include <memory>
#include <ostream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "../chemistry/chemistry.hpp"

namespace prexsyn::chemspace {

class BuildingBlockLibraryError : public std::runtime_error {
public:
    explicit BuildingBlockLibraryError(const std::string &message) : std::runtime_error(message) {}
};

struct BuildingBlockEntry {
    std::shared_ptr<Molecule> molecule;
    std::string identifier;
    std::set<std::string> labels;
};

struct BuildingBlockItem : BuildingBlockEntry {
    size_t index{};

    template <typename Archive> void serialize(Archive &ar, const unsigned int /* version */) {
        if constexpr (Archive::is_saving::value) {
            ar << molecule->serialize();
        } else {
            std::string mol_data;
            ar >> mol_data;
            molecule = Molecule::deserialize(mol_data);
        }
        ar & identifier;
        ar & labels;
        ar & index;
    }
};

class BuildingBlockLibrary {
public:
    using Index = size_t;

private:
    std::vector<BuildingBlockItem> building_blocks_;
    std::map<std::string, Index> identifier_to_index_;

public:
    BuildingBlockLibrary() = default;

    static std::unique_ptr<BuildingBlockLibrary> deserialize(std::istream &);
    void serialize(std::ostream &) const;

    size_t size() const { return building_blocks_.size(); }
    const BuildingBlockItem &get(Index) const;
    const BuildingBlockItem &get(const std::string &) const;
    Index add(const BuildingBlockEntry &);

    auto begin() const noexcept { return building_blocks_.begin(); }
    auto end() const noexcept { return building_blocks_.end(); }
};

} // namespace prexsyn::chemspace
