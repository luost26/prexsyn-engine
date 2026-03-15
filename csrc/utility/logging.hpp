#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include <fmt/format.h>
#include <spdlog/logger.h>

namespace prexsyn {

using Logger = spdlog::logger;

std::shared_ptr<Logger> logger();

std::shared_ptr<Logger> create_logger(const std::string &module);

} // namespace prexsyn

template <> struct fmt::formatter<std::vector<size_t>> : fmt::formatter<std::string> {
    auto format(const std::vector<size_t> &vec, format_context &ctx) const -> decltype(ctx.out()) {
        std::string result = "[";
        for (size_t i = 0; i < vec.size(); ++i) {
            result += std::to_string(vec[i]);
            if (i < vec.size() - 1) {
                result += ", ";
            }
        }
        result += "]";
        return fmt::format_to(ctx.out(), "{}", result);
    }
};
