#include "substructures.hpp"

#include "../fingerprints.hpp"
#include "../molops.hpp"
#include "../utils/algorithm.hpp"

namespace prexsyn_engine {
BRICSFragmentsFeaturizer::BRICSFragmentsFeaturizer(
    const std::string &name, const std::string &fp_type,
    unsigned int max_num_fragments)
    : name(name), max_num_fragments(max_num_fragments),
      fp_func(get_fp_func<float>(fp_type)) {}

void BRICSFragmentsFeaturizer::operator()(const Synthesis &syn,
                                          FeatureBuilder &dict) {
    std::mt19937 rng{std::random_device{}()};
    auto product = random_choice(syn.top(), rng);
    auto fragments = brics_fragments(product);
    std::shuffle(fragments.begin(), fragments.end(),
                 std::mt19937{std::random_device{}()});
    std::vector<std::vector<float>> frag_fps;
    std::vector<bool> frag_exists;
    for (size_t i = 0; i < max_num_fragments; ++i) {
        if (i < fragments.size()) {
            frag_fps.push_back(fp_func(fragments[i]));
            frag_exists.push_back(true);
        } else {
            frag_fps.push_back(fp_func(std::nullopt));
            frag_exists.push_back(false);
        }
    }
    dict.add(name + ".fingerprints", frag_fps);
    dict.add(name + ".fingerprint_exists", frag_exists);
}
} // namespace prexsyn_engine
