#pragma once

#include <memory>
#include <random>
#include <variant>

#include "../container/building_block_list.hpp"
#include "../container/reaction_list.hpp"
#include "../synthesis.hpp"
#include "../types.hpp"
#include "indexer.hpp"

namespace prexsyn_engine {
class no_available_building_blocks : public std::runtime_error {
  public:
    explicit no_available_building_blocks(const std::string &message)
        : std::runtime_error(message) {}
};

class no_available_reactions : public std::runtime_error {
  public:
    explicit no_available_reactions(const std::string &message)
        : std::runtime_error(message) {}
};

class ChemicalSpaceDefinition {
    std::shared_ptr<BuildingBlockList> primary_building_blocks = nullptr;
    std::shared_ptr<SynthesisVector> secondary_building_blocks = nullptr;
    std::shared_ptr<ReactionList> reactions = nullptr;
    std::shared_ptr<ReactionToMolecular<Mol_sptr>> primary_index = nullptr;
    std::shared_ptr<ReactionToMolecular<Synthesis_sptr>> secondary_index =
        nullptr;
    friend class ChemicalSpaceDefinitionBuilder;
    friend class SynthesisGenerator;

  public:
    const std::shared_ptr<BuildingBlockList>
    get_primary_building_blocks() const;
    const std::shared_ptr<SynthesisVector>
    get_secondary_building_blocks() const;
    const std::shared_ptr<ReactionList> get_reactions() const;
    void save(const std::filesystem::path &) const;
    std::variant<Mol_sptr, Synthesis_sptr>
    random_building_block(std::mt19937 &rng) const;
    std::variant<Mol_sptr, Synthesis_sptr>
    random_building_block(const ReactionReactantIndexTuple &,
                          std::mt19937 &rng) const;
    std::vector<ReactionReactantIndexTuple>
    available_reactions(const Mol_sptr &) const;
};

class ChemicalSpaceDefinitionBuilder {
    std::unique_ptr<BuildingBlockList> primary_building_blocks = nullptr;
    std::unique_ptr<SynthesisVector> secondary_building_blocks = nullptr;
    std::unique_ptr<ReactionList> reactions = nullptr;
    std::unique_ptr<ReactionToMolecular<Mol_sptr>> primary_index = nullptr;
    std::unique_ptr<ReactionToMolecular<Synthesis_sptr>> secondary_index =
        nullptr;

  public:
    ChemicalSpaceDefinitionBuilder() = default;
    ChemicalSpaceDefinitionBuilder(const ChemicalSpaceDefinitionBuilder &) =
        delete;

    ChemicalSpaceDefinitionBuilder &
    building_blocks_from_sdf(const std::filesystem::path &,
                             const BuildingBlockPreprocessingOption & = {});
    ChemicalSpaceDefinitionBuilder &
    building_blocks_from_cache(const std::filesystem::path &);

    ChemicalSpaceDefinitionBuilder &
    reactions_from_txt(const std::filesystem::path &);
    ChemicalSpaceDefinitionBuilder &
    reactions_from_cache(const std::filesystem::path &);

    ChemicalSpaceDefinitionBuilder &
    secondary_building_blocks_from_single_reaction();
    ChemicalSpaceDefinitionBuilder &
    secondary_building_blocks_from_cache(const std::filesystem::path &);
    ChemicalSpaceDefinitionBuilder &no_secondary_building_blocks();

    ChemicalSpaceDefinitionBuilder &build_primary_index();
    ChemicalSpaceDefinitionBuilder &build_secondary_index();

    ChemicalSpaceDefinitionBuilder &
    all_from_cache(const std::filesystem::path &);

    std::shared_ptr<ChemicalSpaceDefinition> build();

    SynthesisVector::size_type count_secondary_building_blocks() const;
    Synthesis_sptr
        get_secondary_building_block(SynthesisVector::size_type) const;
};

struct SynthesisGeneratorOption {
    unsigned int num_reactions_cutoff = 5;
    unsigned int num_product_atoms_cutoff = 50;
};

class SynthesisGenerator {
    std::shared_ptr<ChemicalSpaceDefinition> csd;
    std::unique_ptr<Synthesis> synthesis;
    std::mt19937 rng;

    void reset_synthesis();
    void init_synthesis();
    bool needs_reinit() const;

  public:
    SynthesisGeneratorOption option;
    SynthesisGenerator(const std::shared_ptr<ChemicalSpaceDefinition> &,
                       const SynthesisGeneratorOption & = {},
                       std::mt19937::result_type seed = std::random_device{}());
    SynthesisGenerator(const SynthesisGenerator &other);
    Synthesis_sptr next();
};
} // namespace prexsyn_engine
