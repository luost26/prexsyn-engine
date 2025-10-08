#pragma once

#include <filesystem>

#include "synthesis.hpp"
#include "types.hpp"

namespace prexsyn_engine {
typedef size_t ReactionIndex;
typedef size_t ReactantIndex;
typedef size_t MolecularIndex;

struct ReactionReactantIndexTuple {
    ReactionIndex reaction;
    ReactantIndex reactant;
};

std::vector<ReactantIndex>
get_suitable_reactant_indices(const Reaction_sptr &reaction,
                              const Mol_sptr &mol);
std::vector<ReactantIndex>
get_suitable_reactant_indices(const Reaction_sptr &reaction,
                              const Synthesis_sptr &synthesis);

template <typename Molecular_sptr> class ReactionToMolecular {
    std::vector<std::vector<std::vector<MolecularIndex>>> index;
    ReactionToMolecular() = default;

  public:
    ReactionToMolecular(const std::vector<Molecular_sptr> &,
                        const ReactionVector &);
    ReactionToMolecular(const ReactionToMolecular &) = default;

    void save(const std::filesystem::path &) const;
    static ReactionToMolecular<Molecular_sptr> *
    load(const std::filesystem::path &);

    size_t num_reactions() const;
    size_t num_reactants(ReactionIndex) const;
    const std::vector<MolecularIndex> &
        get_molecular_indices(ReactionIndex, ReactantIndex) const;
};
} // namespace prexsyn_engine
