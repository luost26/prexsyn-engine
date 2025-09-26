#pragma once

#include <concepts>
#include <stdexcept>
#include <vector>

#include "../utils/assert.hpp"

namespace synthesis_backend {
template <typename T>
    requires std::integral<T> || std::floating_point<T>
static std::vector<long> shape_of(const std::vector<std::vector<T>> &vec) {
    Ensures(!vec.empty());
    std::vector<long> shape{(long)vec.size()};
    size_t second_dim_size = vec[0].size();
    for (const auto &v : vec) {
        Ensures(v.size() == second_dim_size);
    }
    shape.push_back((long)second_dim_size);
    return shape;
}

template <typename T>
    requires std::integral<T> || std::floating_point<T>
static std::vector<long> shape_of(const std::vector<T> &vec) {
    return {(long)vec.size()};
}

template <typename T>
    requires std::integral<T> || std::floating_point<T>
static std::vector<long> shape_of(const T &) {
    return {};
}

template <typename T>
concept supported_dtype = std::is_same_v<T, float> || std::is_same_v<T, long> ||
                          std::is_same_v<T, bool>;

using Long = long;
using Float = float;
using Bool = bool;

struct DType {
    enum Type {
        Float = 0,
        Long = 1,
        Bool = 2,
    };

    template <supported_dtype T> static Type get_type() {
        if constexpr (std::is_same_v<T, float>) {
            return Float;
        } else if constexpr (std::is_same_v<T, long>) {
            return Long;
        } else if constexpr (std::is_same_v<T, bool>) {
            return Bool;
        } else {
            throw std::runtime_error("Unsupported type");
        }
    }

    static size_t get_size(Type type) {
        switch (type) {
        case Float:
            return sizeof(float);
        case Long:
            return sizeof(long);
        case Bool:
            return sizeof(bool);
        default:
            throw std::runtime_error("Unsupported type");
        }
    }
};

} // namespace synthesis_backend
