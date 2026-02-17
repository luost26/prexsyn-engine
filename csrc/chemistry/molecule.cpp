#include "molecule.hpp"

#include <memory>
#include <string>
#include <utility>

#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

namespace prexsyn {

std::unique_ptr<Molecule> Molecule::from_smiles(const std::string &smiles) {
    RDKit::ROMOL_SPTR rdkit_mol(RDKit::SmilesToMol(smiles));
    if (!rdkit_mol) {
        throw MoleculeError("Failed to parse SMILES: " + smiles);
    }
    return std::make_unique<Molecule>(std::move(rdkit_mol));
}

std::unique_ptr<Molecule> Molecule::from_unsanitized_rdkit(const RDKit::ROMOL_SPTR &rdkit_mol) {
    if (!rdkit_mol) {
        throw MoleculeError("RDKit molecule pointer is null");
    }
    auto sanitized = std::make_unique<RDKit::RWMol>(*rdkit_mol);
    try {
        RDKit::MolOps::sanitizeMol(*sanitized);
    } catch (const RDKit::MolSanitizeException &e) {
        throw MoleculeError("Failed to sanitize RDKit molecule: " + std::string(e.what()));
    }
    return std::make_unique<Molecule>(std::move(sanitized));
}

} // namespace prexsyn
