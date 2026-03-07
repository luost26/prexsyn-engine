#pragma once

#include <filesystem>
#include <memory>

#include "rxn_lib.hpp"

namespace prexsyn::chemspace {

std::unique_ptr<ReactionLibrary> rxn_lib_from_plain_text(const std::filesystem::path &,
                                                         bool ignore_errors = false);

}
