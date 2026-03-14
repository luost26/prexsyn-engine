#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <utility>

#include <GraphMol/Fingerprints/MorganGenerator.h>

#include "../chemistry/chemistry.hpp"
#include "base.hpp"

namespace prexsyn::descriptor {

class MorganFingerprint : public MoleculeDescriptor {
private:
    using Generator = RDKit::FingerprintGenerator<std::uint64_t>;
    std::unique_ptr<Generator> generator_;

public:
    MorganFingerprint(std::unique_ptr<Generator> generator) : generator_(std::move(generator)) {};

    static std::unique_ptr<MorganFingerprint> ecfp4();
    static std::unique_ptr<MorganFingerprint> fcfp4();

    size_t size() const override { return generator_->getOptions()->d_fpSize; }
    DataType::T dtype() const override { return DataType::bool8; }

    void operator()(const Molecule &, std::span<std::byte> &) const override;
};

} // namespace prexsyn::descriptor
