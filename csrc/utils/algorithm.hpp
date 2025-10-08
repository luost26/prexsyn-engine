#pragma once

#include <exception>
#include <functional>
#include <random>
#include <set>
#include <stdexcept>
#include <variant>
#include <vector>

namespace prexsyn_engine {

class stop_iteration : public std::exception {};

template <typename Container>
void generate_combinations(
    const std::vector<Container> &vect_of_containers,
    std::function<void(std::vector<typename Container::value_type>)> callback) {

    using ValueType = typename Container::value_type;

    std::vector<ValueType> current_combination;
    for (const Container &container : vect_of_containers) {
        if (container.empty()) {
            // If any container is empty, we cannot form a combination
            return;
        }
        current_combination.push_back(*container.begin());
    }

    std::function<void(size_t)> backtrack = [&](size_t depth) {
        if (depth == vect_of_containers.size()) {
            callback(current_combination);
            return;
        }
        for (const ValueType &value : vect_of_containers[depth]) {
            current_combination[depth] = value;
            backtrack(depth + 1);
        }
    };
    backtrack(0);
}

template <typename T>
T random_choice(const std::vector<T> &vec, std::mt19937 &rng) {
    if (vec.empty()) {
        throw std::out_of_range("Cannot choose from an empty vector");
    }
    std::uniform_int_distribution<> dist(0, vec.size() - 1);

    return vec[dist(rng)];
}

template <typename T>
T random_choice(const std::set<T> &set, std::mt19937 &rng) {
    if (set.empty()) {
        throw std::out_of_range("Cannot choose from an empty set");
    }

    std::uniform_int_distribution<> dist(0, set.size() - 1);

    auto it = set.begin();
    std::advance(it, dist(rng));
    return *it;
}

template <typename T1, typename T2>
std::variant<T1, T2>
random_choice(const std::vector<T1> &vec1, const std::vector<T2> &vec2,
              std::mt19937 &rng, bool *selected_2 = nullptr) {
    if (vec1.empty() && vec2.empty()) {
        throw std::out_of_range("Cannot choose from empty vectors");
    }
    std::uniform_int_distribution<> dist(0, vec1.size() + vec2.size() - 1);
    size_t index = dist(rng);
    if (index < vec1.size()) {
        if (selected_2) {
            *selected_2 = false;
        }
        return vec1[index];
    } else {
        if (selected_2) {
            *selected_2 = true;
        }
        index -= vec1.size();
        return vec2[index];
    }
}

template <typename T>
T random_choice(const std::vector<T> &vec1, const std::vector<T> &vec2,
                std::mt19937 &rng, bool *selected_2 = nullptr) {
    if (vec1.empty() && vec2.empty()) {
        throw std::out_of_range("Cannot choose from empty vectors");
    }
    std::uniform_int_distribution<> dist(0, vec1.size() + vec2.size() - 1);
    size_t index = dist(rng);
    if (index < vec1.size()) {
        if (selected_2) {
            *selected_2 = false;
        }
        return vec1[index];
    } else {
        if (selected_2) {
            *selected_2 = true;
        }
        index -= vec1.size();
        return vec2[index];
    }
}
} // namespace prexsyn_engine
