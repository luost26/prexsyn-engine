#pragma once

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

namespace prexsyn {

class MoleculeError : public std::runtime_error {
public:
    explicit MoleculeError(const std::string &message) : std::runtime_error(message) {}
};

class Molecule {
private:
    RDKit::ROMOL_SPTR rdkit_mol_;
    std::string smiles_;

public:
    Molecule(RDKit::ROMOL_SPTR rdkit_mol) : rdkit_mol_(std::move(rdkit_mol)) {
        if (!rdkit_mol_) {
            throw MoleculeError("RDKit molecule pointer is null");
        }
        smiles_ = RDKit::MolToSmiles(*rdkit_mol_);
    }
    static std::unique_ptr<Molecule> from_smiles(const std::string &smiles);
    static std::unique_ptr<Molecule> from_unsanitized_rdkit(const RDKit::ROMOL_SPTR &rdkit_mol);
    static std::unique_ptr<Molecule> from_rdkit_pickle(const std::string &);

    static std::unique_ptr<Molecule> deserialize(const std::string &data) {
        return from_rdkit_pickle(data);
    }
    std::string serialize() const { return rdkit_pickle(); }

    const RDKit::ROMol &rdkit_mol() const { return *rdkit_mol_; }
    RDKit::ROMol &rdkit_mol() { return *rdkit_mol_; }
    RDKit::ROMOL_SPTR rdkit_mol_ptr() const { return rdkit_mol_; }
    std::string rdkit_pickle() const;

    unsigned int num_heavy_atoms() const { return rdkit_mol_->getNumHeavyAtoms(); }
    const std::string &smiles() const { return smiles_; }

    std::unique_ptr<Molecule> largest_fragment() const;
};

} // namespace prexsyn
