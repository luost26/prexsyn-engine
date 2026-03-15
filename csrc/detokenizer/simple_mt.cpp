#include "simple_mt.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <stdexcept>
#include <vector>

#include <omp.h>

#include "../chemspace/chemspace.hpp"
#include "base.hpp"

namespace prexsyn::detokenizer {

std::vector<std::unique_ptr<chemspace::Synthesis>>
MultiThreadedDetokenizer::operator()(size_t batch_size,
                                     const std::span<const std::int64_t> &tokens) const {

    auto seqlen = (tokens.size() / batch_size) / 3;
    if (tokens.size() != batch_size * seqlen * 3) {
        throw std::invalid_argument("size must be batch_size * seqlen * 3");
    }

    std::vector<std::unique_ptr<chemspace::Synthesis>> out;
    out.resize(batch_size);

#pragma omp parallel for
    for (size_t i = 0; i < batch_size; ++i) {
        out[i] = detokenize(tokens.subspan(i * seqlen * 3, seqlen * 3), cs_, token_def_);
    }

    return out;
}

} // namespace prexsyn::detokenizer
