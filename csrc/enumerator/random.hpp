#pragma once

#include <cstddef>
#include <memory>
#include <optional>
#include <random>
#include <utility>

#include "../chemistry/chemistry.hpp"
#include "../chemspace/chemspace.hpp"
#include "base.hpp"

namespace prexsyn::enumerator {

class RandomEnumerator {
public:
    using Config = EnumeratorConfig;

private:
    std::shared_ptr<chemspace::ChemicalSpace> cs_;
    Config config_;

    std::unique_ptr<chemspace::Synthesis> synthesis_;
    std::mt19937 rng_;

    bool not_growable() const;

    void clear_synthesis();
    void init_synthesis();
    void grow_synthesis();
    std::optional<std::shared_ptr<chemspace::Synthesis>> try_next();

public:
    RandomEnumerator(std::shared_ptr<chemspace::ChemicalSpace> cs,
                     const Config &config = default_config,
                     std::optional<size_t> random_seed = std::nullopt);

    std::shared_ptr<chemspace::Synthesis> next();
    std::pair<std::shared_ptr<chemspace::Synthesis>, std::shared_ptr<Molecule>> next_with_product();
};

} // namespace prexsyn::enumerator
