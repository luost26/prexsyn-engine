#include "rdkit_descriptors.hpp"

#include "../utils/algorithm.hpp"

namespace prexsyn_engine {
RDKitDescriptorsFeaturizer::RDKitDescriptorsFeaturizer(
    const std::string &name, unsigned int num_evaluated_descriptors,
    const std::vector<std::string> &descriptor_names)
    : name(name), num_evaluated_descriptors(num_evaluated_descriptors),
      descriptor_names(descriptor_names) {
    for (size_t i = 1; i <= SUPPORTED_RDKIT_DESCRIPTORS.size(); i++) {
        auto name = SUPPORTED_RDKIT_DESCRIPTORS[i - 1];
        descriptor_name_to_index[name] = i;
    }

    for (auto desc_name : descriptor_names) {
        Ensures(descriptor_name_to_index.find(desc_name) !=
                descriptor_name_to_index.end());
    }
}

size_t RDKitDescriptorsFeaturizer::max_property_index() const {
    return 1 + SUPPORTED_RDKIT_DESCRIPTORS.size();
}

static float calc_property(const std::string &name, const RDKit::ROMol &mol) {
    return (*RDKit::Descriptors::Properties::getProperty(name))(mol);
}

void RDKitDescriptorsFeaturizer::operator()(const Synthesis &syn,
                                            FeatureBuilder &dict) {
    std::mt19937 rng{std::random_device{}()};
    auto product = random_choice(syn.top(), rng);

    std::vector<Long> prop_indices(num_evaluated_descriptors);
    std::vector<Float> prop_values(num_evaluated_descriptors);
    auto prop_indices_it = prop_indices.begin();
    auto prop_values_it = prop_values.begin();

    auto add_prop = [&](size_t index, float value) {
        if (prop_indices_it == prop_indices.end() ||
            prop_values_it == prop_values.end()) {
            throw std::runtime_error(
                "ProductRDKitPropertyFeaturizer: "
                "property indices or values vector is full");
        }
        *prop_indices_it++ = index;
        *prop_values_it++ = value;
    };

    std::vector<std::string> descs_to_eval;
    std::sample(descriptor_names.begin(), descriptor_names.end(),
                std::back_inserter(descs_to_eval), num_evaluated_descriptors,
                std::mt19937{std::random_device{}()});
    for (const auto &desc_name : descs_to_eval) {
        add_prop(descriptor_name_to_index[desc_name],
                 calc_property(desc_name, *product));
    }

    dict.add(name + ".types", prop_indices);
    dict.add(name + ".values", prop_values);
}

} // namespace prexsyn_engine
