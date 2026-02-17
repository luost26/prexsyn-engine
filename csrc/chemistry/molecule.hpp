#pragma once

#include <memory>

#include <GraphMol/GraphMol.h>

namespace prexsyn {

class Molecule {
private:
    RDKit::ROMOL_SPTR rdkit_mol_;

public:
    Molecule(RDKit::ROMOL_SPTR rdkit_mol) : rdkit_mol_(std::move(rdkit_mol)) {
        if (!rdkit_mol_) {
            throw std::runtime_error("RDKit molecule pointer is null");
        }
    }
    static std::unique_ptr<Molecule> from_smiles(const std::string &smiles);

    const RDKit::ROMol &rdkit_mol() const { return *rdkit_mol_; }
    RDKit::ROMol &rdkit_mol() { return *rdkit_mol_; }
    RDKit::ROMOL_SPTR rdkit_mol_ptr() const { return rdkit_mol_; }
};

} // namespace prexsyn
