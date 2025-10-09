#include "indexer.hpp"

#include <filesystem>
#include <omp.h>
#include <typeinfo>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>

#include "../synthesis.hpp"
#include "../types.hpp"
#include "../utils/assert.hpp"
#include "../utils/logging.hpp"

namespace prexsyn_engine {

std::vector<size_t> get_suitable_reactant_indices(const Reaction_sptr &reaction,
                                                  const Mol_sptr &mol) {
    std::vector<size_t> indices;
    for (size_t i = 0; i < reaction->getNumReactantTemplates(); ++i) {
        const auto &tmpl = *(reaction->getReactants()[i]);
        RDKit::MatchVectType res;
        RDKit::SubstructMatch(*mol, tmpl, res);
        if (!res.empty()) {
            indices.push_back(i);
        }
    }
    return indices;
}

std::vector<size_t>
get_suitable_reactant_indices(const Reaction_sptr &reaction,
                              const Synthesis_sptr &synthesis) {
    std::vector<size_t> indices;
    for (size_t i = 0; i < reaction->getNumReactantTemplates(); ++i) {
        const auto &tmpl = *(reaction->getReactants()[i]);
        for (const auto &product : synthesis->top()) {
            RDKit::MatchVectType res;
            bool matched = RDKit::SubstructMatch(*product, tmpl, res, false);
            if (matched) {
                indices.push_back(i);
                break;
            }
        }
    }
    return indices;
}

template <typename Molecular_sptr>
ReactionToMolecular<Molecular_sptr>::ReactionToMolecular(
    const std::vector<Molecular_sptr> &molecules,
    const ReactionVector &reactions) {
    logger()->info("Start creating index: Reaction[{}] -> {}[{}]",
                   reactions.size(), typeid(Molecular_sptr).name(),
                   molecules.size());
    index.resize(reactions.size());
    for (size_t i = 0; i < reactions.size(); ++i) {
        index[i].resize(reactions[i]->getNumReactantTemplates());
    }

    size_t count_links = 0;
#pragma omp parallel for
    for (size_t mol_idx = 0; mol_idx < molecules.size(); ++mol_idx) {
        const auto &mol = molecules[mol_idx];
        for (size_t rxn_idx = 0; rxn_idx < reactions.size(); ++rxn_idx) {
            std::vector<size_t> reactant_indices =
                get_suitable_reactant_indices(reactions[rxn_idx], mol);
#pragma omp critical
            {
                for (const auto &reactant_index : reactant_indices) {
                    index[rxn_idx][reactant_index].push_back(mol_idx);
                }
                count_links += reactant_indices.size();
            }
        }
    }
    logger()->info("Index created: {} links", count_links);
}

template <typename Molecular_sptr>
void ReactionToMolecular<Molecular_sptr>::save(
    const std::filesystem::path &filename) const {
    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        throw std::runtime_error("Could not open file for writing: " +
                                 filename.string());
    }

    boost::archive::binary_oarchive oa(ofs);
    oa << index;
    ofs.close();
}

template <typename Molecular_sptr>
ReactionToMolecular<Molecular_sptr> *ReactionToMolecular<Molecular_sptr>::load(
    const std::filesystem::path &filename) {
    std::ifstream ifs(filename);
    if (!ifs.is_open()) {
        throw std::runtime_error("Could not open file for reading: " +
                                 filename.string());
    }

    boost::archive::binary_iarchive ia(ifs);
    auto object = new ReactionToMolecular<Molecular_sptr>();
    ia >> object->index;
    ifs.close();
    return object;
}

template <typename Molecular_sptr>
size_t ReactionToMolecular<Molecular_sptr>::num_reactions() const {
    return index.size();
}

template <typename Molecular_sptr>
size_t ReactionToMolecular<Molecular_sptr>::num_reactants(
    ReactionIndex reaction_index) const {
    Ensures(reaction_index < index.size());
    return index[reaction_index].size();
}

template <typename Molecular_sptr>
const std::vector<MolecularIndex> &
ReactionToMolecular<Molecular_sptr>::get_molecular_indices(
    ReactionIndex reaction_index, ReactantIndex reactant_index) const {
    Ensures(reaction_index < index.size());
    Ensures(reactant_index < index[reaction_index].size());
    return index[reaction_index][reactant_index];
}

template class ReactionToMolecular<Mol_sptr>;
template class ReactionToMolecular<Synthesis_sptr>;
} // namespace prexsyn_engine
