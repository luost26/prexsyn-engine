#pragma once

#include <string>

#include <spdlog/logger.h>

namespace prexsyn {

using Logger = spdlog::logger;

std::shared_ptr<Logger> logger();

std::shared_ptr<Logger> create_logger(const std::string &module);

} // namespace prexsyn
