#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <utility>
#include <vector>

#include <GraphMol/Fingerprints/MorganGenerator.h>

#include "../chemistry/chemistry.hpp"
#include "base.hpp"

namespace prexsyn::descriptor {

class MorganFingerprint : public MoleculeDescriptor {
private:
    using Generator = RDKit::FingerprintGenerator<std::uint64_t>;
    std::unique_ptr<Generator> generator_;

    std::unique_ptr<RDKit::AtomInvariantsGenerator> atom_invariants_generator_ = nullptr;
    std::unique_ptr<RDKit::BondInvariantsGenerator> bond_invariants_generator_ = nullptr;

public:
    MorganFingerprint(std::unique_ptr<Generator> generator) : generator_(std::move(generator)) {};

    static std::unique_ptr<MorganFingerprint> ecfp4();
    static std::unique_ptr<MorganFingerprint> fcfp4();

    std::vector<size_t> size() const override { return {generator_->getOptions()->d_fpSize}; }
    DataType::T dtype() const override { return DataType::bool8; }

    void operator()(const Molecule &mol, std::span<std::byte> &out) const override;
};

} // namespace prexsyn::descriptor
