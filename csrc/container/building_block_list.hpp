#pragma once

#include <GraphMol/FileParsers/MolSupplier.h>
#include <filesystem>
#include <optional>

#include "../types.hpp"

namespace prexsyn_engine {
struct BuildingBlockPreprocessingOption {
    bool remove_Hs = true;
    bool largest_fragment_only = true;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int version);
};

using BuildingBlockIndex = size_t;

class BuildingBlockList {
    MolVector building_blocks;
    BuildingBlockPreprocessingOption option;

    std::optional<Mol_sptr> preprocess(const Mol_sptr &) const;
    std::optional<BuildingBlockIndex> preprocess_and_add(const Mol_sptr &,
                                                         int original_index);

  public:
    BuildingBlockList(const MolVector &,
                      const BuildingBlockPreprocessingOption &);
    BuildingBlockList(RDKit::MolSupplier &,
                      const BuildingBlockPreprocessingOption &);
    static BuildingBlockList *
    from_sdf(const std::filesystem::path &,
             const BuildingBlockPreprocessingOption & = {});
    const MolVector &get_building_blocks() const;
    void save(const std::filesystem::path &) const;
    static BuildingBlockList *load(const std::filesystem::path &);
    static size_t peek_size(const std::filesystem::path &);
    Mol_sptr get(size_t index) const;
    BuildingBlockIndex size() const;
};
} // namespace prexsyn_engine
