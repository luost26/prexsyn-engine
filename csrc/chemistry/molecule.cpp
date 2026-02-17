#include "molecule.hpp"

#include <GraphMol/SmilesParse/SmilesParse.h>

namespace prexsyn {

std::unique_ptr<Molecule> Molecule::from_smiles(const std::string &smiles) {
    RDKit::ROMOL_SPTR rdkit_mol(RDKit::SmilesToMol(smiles));
    if (!rdkit_mol) {
        throw std::runtime_error("Failed to parse SMILES: " + smiles);
    }
    return std::make_unique<Molecule>(std::move(rdkit_mol));
}

} // namespace prexsyn
