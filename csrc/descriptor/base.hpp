#pragma once

#include <cstddef>
#include <functional>
#include <numeric>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

// IWYU pragma: begin_exports
#include "../utility/data_type.hpp"
// IWYU pragma: end_exports
#include "../chemistry/chemistry.hpp"
#include "../chemspace/chemspace.hpp"

namespace prexsyn::descriptor {

template <typename T>
    requires std::is_same_v<T, Molecule> || std::is_same_v<T, chemspace::Synthesis>
class Descriptor {
public:
    Descriptor() = default;
    Descriptor(const Descriptor &) = default;
    Descriptor(Descriptor &&) = default;
    Descriptor &operator=(const Descriptor &) = default;
    Descriptor &operator=(Descriptor &&) = default;

    virtual ~Descriptor() = default;
    virtual DataType::T dtype() const = 0;
    virtual std::vector<size_t> size() const = 0;
    size_t num_elements() const {
        auto s = size();
        return std::accumulate(s.begin(), s.end(), 1, std::multiplies<size_t>());
    }
    size_t size_in_bytes() const { return num_elements() * DataType::get_size(dtype()); }

    template <SupportedDataType U> std::span<U> cast(std::span<std::byte> &out) const {
        if (DataType::get_dtype<U>() != dtype()) {
            throw std::runtime_error("Data type mismatch");
        }
        return std::span<U>(reinterpret_cast<U *>(out.data()),
                            out.size() / DataType::get_size(dtype()));
    }

    virtual void operator()(const T &, std::span<std::byte> &) const = 0;

    void check_size(std::span<std::byte> &out) const {
        if (out.size() != size_in_bytes()) {
            throw std::runtime_error("descriptor size mismatch, expected " +
                                     std::to_string(size_in_bytes()) + " bytes but got " +
                                     std::to_string(out.size()) + " bytes");
        }
    }
};

template class Descriptor<Molecule>;
using MoleculeDescriptor = Descriptor<Molecule>;

template class Descriptor<chemspace::Synthesis>;
using SynthesisDescriptor = Descriptor<chemspace::Synthesis>;

} // namespace prexsyn::descriptor
