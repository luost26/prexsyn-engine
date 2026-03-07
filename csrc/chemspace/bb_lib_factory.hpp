#pragma once

#include <filesystem>
#include <memory>

#include "../chemistry/chemistry.hpp"
#include "bb_lib.hpp"

namespace prexsyn::chemspace {

struct BuildingBlockPreprocessor {
    bool largest_fragment_only = true;

    std::shared_ptr<Molecule> operator()(const std::shared_ptr<Molecule> &mol) const {
        auto m = mol;
        m = largest_fragment_only ? mol->largest_fragment() : m;
        return m;
    }
};

std::unique_ptr<BuildingBlockLibrary> bb_lib_from_sdf(const std::filesystem::path &,
                                                      const BuildingBlockPreprocessor & = {});

} // namespace prexsyn::chemspace
