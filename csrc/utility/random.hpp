#pragma once

#include <cstddef>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

namespace prexsyn {

template <typename T, typename RNG> const T &random_choice(const std::vector<T> &vec, RNG &rng) {
    if (vec.empty()) {
        throw std::out_of_range("Cannot choose from an empty vector");
    }
    std::uniform_int_distribution<size_t> dist(0, vec.size() - 1);
    return vec[dist(rng)];
}

enum class which_vector : std::uint8_t { first, second };

template <typename T, typename RNG>
std::pair<which_vector, const T &> random_choice(const std::vector<T> &v1, const std::vector<T> &v2,
                                                 RNG &rng) {
    if (v1.empty() && v2.empty()) {
        throw std::out_of_range("Cannot choose from two empty vectors");
    }

    auto total = v1.size() + v2.size();
    std::uniform_int_distribution<size_t> dist(0, total - 1);
    auto idx = dist(rng);
    if (idx < v1.size()) {
        return {which_vector::first, v1[idx]};
    } else {
        return {which_vector::second, v2[idx - v1.size()]};
    }
}

} // namespace prexsyn
