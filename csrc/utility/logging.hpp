#pragma once

#include <memory>

#include <spdlog/logger.h>

namespace prexsyn {

std::shared_ptr<spdlog::logger> logger();

} // namespace prexsyn
