#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>

#include "../chemspace/chemspace.hpp"
#include "../descriptor/descriptor.hpp"

namespace prexsyn::detokenizer {

std::unique_ptr<chemspace::Synthesis> detokenize(const std::span<const std::int64_t> &,
                                                 const std::shared_ptr<chemspace::ChemicalSpace> &,
                                                 const descriptor::TokenDef &,
                                                 std::optional<size_t> max_outcomes_per_reaction);

}
