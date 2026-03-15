#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <type_traits>

namespace prexsyn {

template <typename T>
concept SupportedDataType =
    std::is_same_v<T, float> || std::is_same_v<T, std::int64_t> || std::is_same_v<T, bool>;

struct DataType {
    enum T : std::uint8_t {
        float32,
        int64,
        bool8,
    };

    template <typename U> static constexpr T get_dtype() {
        if constexpr (std::is_same_v<U, float>) {
            return T::float32;
        } else if constexpr (std::is_same_v<U, std::int64_t>) {
            return T::int64;
        } else if constexpr (std::is_same_v<U, bool>) {
            return T::bool8;
        } else {
            static_assert(!std::is_same_v<U, U>, "Unsupported data type");
        }
    }

    static size_t get_size(const T &t) {
        switch (t) {
        case T::float32:
            return 4;
        case T::int64:
            return 8;
        case T::bool8:
            return 1;
        }
        return 0;
    }

    static std::string to_string(const T &t) {
        switch (t) {
        case T::float32:
            return "float32";
        case T::int64:
            return "int64";
        case T::bool8:
            return "bool8";
        }
        return "unknown";
    }
};

} // namespace prexsyn
