#include "fingerprint.hpp"

#include "../utils/algorithm.hpp"

namespace prexsyn_engine {
FingerprintFeaturizer::FingerprintFeaturizer(const std::string &name,
                                             const std::string &fp_type)
    : fp_func(get_fp_func<float>(fp_type)), name(name) {}

void FingerprintFeaturizer::operator()(const Synthesis &syn,
                                       FeatureBuilder &dict) {
    std::mt19937 rng{std::random_device{}()};
    auto product = random_choice(syn.top(), rng);
    auto fp = fp_func(product);
    dict.add(name + ".fingerprint", fp);
}
} // namespace prexsyn_engine
