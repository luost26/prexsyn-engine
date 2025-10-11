#include "scaffold.hpp"

#include "../molops.hpp"
#include "../utils/algorithm.hpp"

namespace prexsyn_engine {
MurckoScaffoldFeaturizer::MurckoScaffoldFeaturizer(const std::string &name,
                                                   const std::string &fp_type)
    : fp_func(get_fp_func<float>(fp_type)), name(name) {}

void MurckoScaffoldFeaturizer::operator()(const Synthesis &syn,
                                          FeatureBuilder &dict) {
    std::mt19937 rng{std::random_device{}()};
    auto product = random_choice(syn.top(), rng);
    auto scaffold = murcko_scaffold(product);
    // Fallback to product if scaffold is empty
    auto fp = fp_func(scaffold.has_value() ? scaffold : product);
    dict.add(name + ".fingerprint", fp);
}
} // namespace prexsyn_engine
