#pragma once

#include <cstddef>
#include <span>
#include <vector>

#include "../chemistry/chemistry.hpp"
#include "../chemspace/chemspace.hpp"

namespace prexsyn::descriptor {

template <typename D>
    requires(std::is_floating_point_v<D> || std::is_integral_v<D>)
class MoleculeDescriptor {
public:
    MoleculeDescriptor() = default;
    MoleculeDescriptor(const MoleculeDescriptor &) = default;
    MoleculeDescriptor(MoleculeDescriptor &&) = default;
    MoleculeDescriptor &operator=(const MoleculeDescriptor &) = default;
    MoleculeDescriptor &operator=(MoleculeDescriptor &&) = default;

    virtual ~MoleculeDescriptor() = default;
    virtual size_t size() const = 0;
    virtual void operator()(const Molecule &, std::span<D> &out) const = 0;

    std::vector<D> operator()(const Molecule &mol) const {
        std::vector<D> out(size());
        (*this)(mol, out);
        return out;
    }
};

template <typename D>
    requires(std::is_floating_point_v<D> || std::is_integral_v<D>)
class SynthesisDescriptor {
public:
    SynthesisDescriptor() = default;
    SynthesisDescriptor(const SynthesisDescriptor &) = default;
    SynthesisDescriptor(SynthesisDescriptor &&) = default;
    SynthesisDescriptor &operator=(const SynthesisDescriptor &) = default;
    SynthesisDescriptor &operator=(SynthesisDescriptor &&) = default;

    virtual ~SynthesisDescriptor() = default;
    virtual size_t size() const = 0;
    virtual void operator()(const chemspace::Synthesis &, std::span<D> &out) const = 0;

    std::vector<D> operator()(const chemspace::Synthesis &syn) const {
        std::vector<D> out(size());
        (*this)(syn, out);
        return out;
    }
};

} // namespace prexsyn::descriptor
