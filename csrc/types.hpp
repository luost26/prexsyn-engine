#pragma once

#include <memory>
#include <set>

#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/GraphMol.h>

namespace prexsyn_engine {
typedef RDKit::ROMOL_SPTR Mol_sptr;
typedef std::shared_ptr<RDKit::ChemicalReaction> Reaction_sptr;

typedef std::set<Mol_sptr> MolSet;
typedef std::vector<Mol_sptr> MolVector;
typedef std::vector<Reaction_sptr> ReactionVector;
} // namespace prexsyn_engine
