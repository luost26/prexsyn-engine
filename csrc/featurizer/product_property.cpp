#include "product_property.hpp"

#include <algorithm>
#include <random>

#include "../utils/algorithm.hpp"
#include "../utils/assert.hpp"

namespace synthesis_backend {

ProductRDKitPropertyFeaturizer::ProductRDKitPropertyFeaturizer(
    const ProductRDKitPropertyFeaturizerOption &option)
    : option(option) {
    for (size_t i = 0; i < supported_rdkit_properties.size(); i++) {
        auto name = supported_rdkit_properties[i];
        supported_property_name_to_index[name] =
            i + option.rdkit_property_index_offset;
    }

    for (auto prop_name : option.rdkit_properties) {
        Ensures(supported_property_name_to_index.find(prop_name) !=
                supported_property_name_to_index.end());
    }
}

size_t ProductRDKitPropertyFeaturizer::max_property_index() const {
    return option.rdkit_property_index_offset +
           supported_rdkit_properties.size();
}

static float calc_property(const std::string &name, const RDKit::ROMol &mol) {
    return (*RDKit::Descriptors::Properties::getProperty(name))(mol);
}

void ProductRDKitPropertyFeaturizer::operator()(const Synthesis &syn,
                                                FeatureBuilder &dict) {
    std::mt19937 rng{std::random_device{}()};
    auto product = random_choice(syn.top(), rng);

    std::vector<Long> prop_indices(option.num_evaluated_properties);
    std::vector<Float> prop_values(option.num_evaluated_properties);
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

    std::vector<std::string> props_to_eval;
    std::sample(option.rdkit_properties.begin(), option.rdkit_properties.end(),
                std::back_inserter(props_to_eval),
                option.num_evaluated_properties,
                std::mt19937{std::random_device{}()});
    for (const auto &prop_name : props_to_eval) {
        add_prop(supported_property_name_to_index[prop_name],
                 calc_property(prop_name, *product));
    }

    dict.add(option.name + ".types", prop_indices);
    dict.add(option.name + ".values", prop_values);
}

} // namespace synthesis_backend
