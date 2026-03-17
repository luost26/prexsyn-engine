#pragma once

#include <cstdint>
#include <stdexcept>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "data_type.hpp"

namespace prexsyn {

inline auto data_type_to_numpy_dtype(const DataType::T &d) {
    namespace py = pybind11;

    switch (d) {
    case DataType::float32:
        return py::dtype::of<float>();
    case DataType::int64:
        return py::dtype::of<std::int64_t>();
    case DataType::bool8:
        return py::dtype::of<bool>();
    default:
        throw std::runtime_error("Unsupported data type");
    }
}

}; // namespace prexsyn
