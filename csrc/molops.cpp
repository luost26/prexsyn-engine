#include "molops.hpp"

#include <GraphMol/ChemTransforms/ChemTransforms.h>

#include "types.hpp"
#include "utils/assert.hpp"

namespace prexsyn_engine {

std::optional<Mol_sptr> sanitize(const Mol_sptr &mol) {
    Ensures(mol != nullptr);
    auto sanitized = std::make_unique<RDKit::RWMol>(*mol);
    try {
        RDKit::MolOps::sanitizeMol(*sanitized);
    } catch (const RDKit::MolSanitizeException &e) {
        return std::nullopt;
    }
    return Mol_sptr(new RDKit::ROMol(*sanitized));
}

std::optional<Mol_sptr> murcko_scaffold(const Mol_sptr &mol) {
    Ensures(mol != nullptr);
    // MurckoDecompose would actually mutate the mol's prop dict, which is not
    // thread-safe and causes double-free errors
    Mol_sptr mol_copy{new RDKit::ROMol(*mol)};
    Mol_sptr scaffold{RDKit::MurckoDecompose(*mol_copy)};
    if (!scaffold) {
        return std::nullopt;
    }
    return sanitize(scaffold);
}

MolVector brics_fragments(const Mol_sptr &mol) {
    Ensures(mol != nullptr);
    Mol_sptr fragmented_mol{RDKit::MolFragmenter::fragmentOnBRICSBonds(*mol)};
    auto frags = RDKit::MolOps::getMolFrags(*fragmented_mol);
    MolVector out;
    for (const auto &frag : frags) {
        auto o = sanitize(frag);
        if (o) {
            out.push_back(o.value());
        }
    }
    return out;
}
} // namespace prexsyn_engine
