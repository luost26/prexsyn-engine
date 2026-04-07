#pragma once

#include <filesystem>
#include <memory>
#include <string>

#include "rxn_lib.hpp"

namespace prexsyn::chemspace {

std::unique_ptr<ReactionLibrary> rxn_lib_from_plain_text(const std::filesystem::path &,
                                                         bool ignore_errors = false);

struct ReactionCSVConfig {
    std::string name_column = "id";
    std::string smarts_column = "smarts";
    std::string reactant_name_column = "reactant_names";
    std::string reactant_name_delimiter = ";";
};

std::unique_ptr<ReactionLibrary> rxn_lib_from_csv(const std::filesystem::path &,
                                                  const ReactionCSVConfig &config = {},
                                                  bool ignore_errors = false);

std::unique_ptr<ReactionLibrary> rxn_lib_from_json(const std::filesystem::path &,
                                                   bool ignore_errors = false);

} // namespace prexsyn::chemspace
