#pragma once

#include <optional>

#include "types.hpp"
namespace prexsyn_engine {

std::optional<Mol_sptr> sanitize(const Mol_sptr &mol);
std::optional<Mol_sptr> murcko_scaffold(const Mol_sptr &mol);
MolVector brics_fragments(const Mol_sptr &mol);

} // namespace prexsyn_engine
