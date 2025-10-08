#include "product_structure.hpp"

#include <algorithm>
#include <random>

#include "../fingerprints.hpp"
#include "../molops.hpp"
#include "../utils/algorithm.hpp"

namespace prexsyn_engine {

ProductStructureFeaturizer::ProductStructureFeaturizer(
    const ProductStructureFeaturizerOption &option)
    : option(option) {
    if (option.embedder_name_template.find("{}") == std::string::npos) {
        throw std::invalid_argument(
            "Embedder name template must contain placeholder '{}'");
    }
}

void ProductStructureFeaturizer::operator()(const Synthesis &syn,
                                            FeatureBuilder &dict) {
    std::mt19937 rng{std::random_device{}()};
    auto product = random_choice(syn.top(), rng);
    for (const auto &fp_type : option.fp_types) {
        auto fp = fp_func<float>(fp_type)(product);
        std::string name = option.embedder_name_template;
        name.replace(name.find("{}"), 2, fp_type);
        dict.add(name + ".fingerprint", fp);
    }

    if (option.scaffold) {
        auto scaffold = murcko_scaffold(product);
        auto fp = fp_func<float>(option.scaffold_fp_type)(
            scaffold.has_value() ? scaffold
                                 : product); // Use product if no scaffold
        dict.add(option.scaffold_embedder_name + ".fingerprint", fp);
    }

    if (option.num_fragments > 0) {
        auto fragments = brics_fragments(product);
        std::shuffle(fragments.begin(), fragments.end(),
                     std::mt19937{std::random_device{}()});
        std::vector<std::vector<float>> frag_fps;
        std::vector<bool> frag_exists;
        auto fp_f = fp_func<float>(option.fragment_fp_type);
        for (size_t i = 0; i < option.num_fragments; ++i) {
            if (i < fragments.size()) {
                frag_fps.push_back(fp_f(fragments[i]));
                frag_exists.push_back(true);
            } else {
                frag_fps.push_back(fp_f(std::nullopt));
                frag_exists.push_back(false);
            }
        }
        dict.add(option.fragment_embedder_name + ".fingerprints", frag_fps);
        dict.add(option.fragment_embedder_name + ".fingerprint_exists",
                 frag_exists);
    }
}

} // namespace prexsyn_engine
