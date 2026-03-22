#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <vector>

#include "../chemspace/chemspace.hpp"
#include "../descriptor/descriptor.hpp"

namespace prexsyn::detokenizer {

class MultiThreadedDetokenizer {
private:
    std::shared_ptr<chemspace::ChemicalSpace> cs_;
    descriptor::TokenDef token_def_;
    std::optional<size_t> max_outcomes_per_reaction_;

public:
    MultiThreadedDetokenizer(const std::shared_ptr<chemspace::ChemicalSpace> &cs,
                             const descriptor::TokenDef &token_def,
                             std::optional<size_t> max_outcomes_per_reaction = std::nullopt)
        : cs_(cs), token_def_(token_def), max_outcomes_per_reaction_(max_outcomes_per_reaction) {}

    std::vector<std::unique_ptr<chemspace::Synthesis>>
    operator()(size_t batch_size, const std::span<const std::int64_t> &) const;
};

} // namespace prexsyn::detokenizer
