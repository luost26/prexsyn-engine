#include "morgan.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <utility>

#include <DataStructs/ExplicitBitVect.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>
#include <GraphMol/Fingerprints/RDKitFPGenerator.h>

#include "../chemistry/chemistry.hpp"

namespace prexsyn::descriptor {

std::unique_ptr<MorganFingerprint> MorganFingerprint::ecfp4() {
    std::unique_ptr<Generator> generator{
        RDKit::MorganFingerprint::getMorganGenerator<std::uint64_t>(2, false, false, true, false,
                                                                    true, nullptr, nullptr)};
    return std::make_unique<MorganFingerprint>(std::move(generator));
}

std::unique_ptr<MorganFingerprint> MorganFingerprint::fcfp4() {
    auto feature_inv_gen = RDKit::MorganFingerprint::MorganFeatureAtomInvGenerator();
    std::unique_ptr<Generator> generator{
        RDKit::MorganFingerprint::getMorganGenerator<std::uint64_t>(
            2, false, false, true, false, true, &feature_inv_gen, nullptr)};
    return std::make_unique<MorganFingerprint>(std::move(generator));
}

void MorganFingerprint::operator()(const Molecule &mol, std::span<std::byte> &out) const {
    check_size(out);
    auto out_t = cast<bool>(out);

    // getFingerprint returns a raw pointer with released ownership
    std::unique_ptr<ExplicitBitVect> fp{generator_->getFingerprint(mol.rdkit_mol())};
    for (size_t i = 0; i < out.size(); ++i) {
        out_t[i] = (*fp)[i];
    }
}

} // namespace prexsyn::descriptor
