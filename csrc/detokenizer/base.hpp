#pragma once

#include <cstdint>
#include <memory>
#include <span>

#include "../chemspace/chemspace.hpp"
#include "../descriptor/descriptor.hpp"

namespace prexsyn::detokenizer {

std::unique_ptr<chemspace::Synthesis> detokenize(const std::span<const std::int64_t> &,
                                                 const std::shared_ptr<chemspace::ChemicalSpace> &,
                                                 const descriptor::TokenDef &);

}
