#pragma once

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include <GraphMol/GraphMol.h>

namespace prexsyn {

class MoleculeError : public std::runtime_error {
public:
    explicit MoleculeError(const std::string &message) : std::runtime_error(message) {}
};

class Molecule {
private:
    RDKit::ROMOL_SPTR rdkit_mol_;

public:
    Molecule(RDKit::ROMOL_SPTR rdkit_mol) : rdkit_mol_(std::move(rdkit_mol)) {
        if (!rdkit_mol_) {
            throw MoleculeError("RDKit molecule pointer is null");
        }
    }
    static std::unique_ptr<Molecule> from_smiles(const std::string &smiles);
    static std::unique_ptr<Molecule> from_unsanitized_rdkit(const RDKit::ROMOL_SPTR &rdkit_mol);

    const RDKit::ROMol &rdkit_mol() const { return *rdkit_mol_; }
    RDKit::ROMol &rdkit_mol() { return *rdkit_mol_; }
    RDKit::ROMOL_SPTR rdkit_mol_ptr() const { return rdkit_mol_; }

    unsigned int num_heavy_atoms() const { return rdkit_mol_->getNumHeavyAtoms(); }
};

} // namespace prexsyn
