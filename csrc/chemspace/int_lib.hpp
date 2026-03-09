#pragma once

#include <cstddef>
#include <istream>
#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "../chemistry/chemistry.hpp"
#include "postfix_notation.hpp"

namespace prexsyn::chemspace {

struct IntermediateEntry {
    PostfixNotation postfix_notation;
    std::shared_ptr<Molecule> molecule;
    std::string identifier;
};

struct IntermediateItem : IntermediateEntry {
    size_t index{};

    template <typename Archive> void serialize(Archive &ar, const unsigned int /* version */) {
        if constexpr (Archive::is_saving::value) {
            ar << molecule->serialize();
        } else {
            std::string mol_data;
            ar >> mol_data;
            molecule = Molecule::deserialize(mol_data);
        }
        ar & postfix_notation;
        ar & identifier;
        ar & index;
    }
};

class IntermediateLibrary {
public:
    using Index = size_t;

private:
    std::vector<IntermediateItem> intermediates_;
    std::map<std::string, Index> identifier_to_index_;

public:
    IntermediateLibrary() = default;

    static std::unique_ptr<IntermediateLibrary> deserialize(std::istream &);
    void serialize(std::ostream &) const;

    size_t size() const { return intermediates_.size(); }
    const IntermediateItem &get(Index) const;
    const IntermediateItem &get(const std::string &) const;
    Index add(const IntermediateEntry &);
    void clear();

    auto begin() const noexcept { return intermediates_.begin(); }
    auto end() const noexcept { return intermediates_.end(); }
};

} // namespace prexsyn::chemspace
